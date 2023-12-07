// system
#include <unordered_map>
#include <set>
#include <chrono>
#include <regex>
#include <filesystem>
#include <random>

namespace fs = std::filesystem;

// PCL
#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl/io/ply_io.h>
#include <pcl/common/io.h>
#include <pcl/registration/transformation_estimation_svd.h>
#include <pcl/search/kdtree.h>
#include <pcl/octree/octree.h>
#include <pcl/filters/uniform_sampling.h>
#include <pcl/common/distances.h>
#include <pcl/common/common.h>
#include <pcl/common/intersections.h>

#include <Eigen/Dense>

#include <spdlog/spdlog.h>

// local
#include "Common.h"
#include "Registration.h"
#include "util.h"
#include "PLANE_Extraction.h"

namespace GPSCO
{
	Params::Params()
	{
		// Parameter initialization
		min_support_points = 1000;
		SmoothnessThreshold = 2.0;
		CurvatureThreshold = 1.0;
		parallel_thresh = 5.0;
		coplanar_thresh = 2.0;
		dist_thresh = 0.5;
	}

	Params::~Params()
	{

	}

	MatchGroup::MatchGroup()
	{
		group1_index = -1;
		group2_index = -1;
		match_num = -1;
		score = -1.0;
	}

	MatchGroup::~MatchGroup()
	{
	}

	Group_Base::Group_Base()
	{
		group1_index = -1;
		group2_index = -1;
		group_angle = -1.0;
	}

	Group_Base::~Group_Base()
	{
	}

	Group_Three::Group_Three()
	{
		group1_index = -1;
		group2_index = -1;
		group3_index = -1;
		group_angle12 = -1.0;
		group_angle13 = -1.0;
		group_angle23 = -1.0;
	}

	Group_Three::~Group_Three()
	{
	}

	Group_Three_Pair::Group_Three_Pair()
	{
		group_src_index1 = -1;
		group_src_index2 = -1;
		group_src_index3 = -1;
		group_tgt_index1 = -1;
		group_tgt_index2 = -1;
		group_tgt_index3 = -1;

		Intersection_points_src.reset(new GPSCO::cloud);
		Intersection_points_tgt.reset(new GPSCO::cloud);

		match_num = 0;
		score = 0.0f;
	}

	Group_Three_Pair::~Group_Three_Pair()
	{
	}

	RT_Info::RT_Info()
	{
		Intersection_points_src.reset(new GPSCO::cloud);
		Intersection_points_tgt.reset(new GPSCO::cloud);
		plane_match_num = 0;
		confidence = 0.0f;
		rms = 0.0f;
	}

	RT_Info::~RT_Info()
	{
	}

	void RT_Info::Compute_RMS()
	{
		rms = 0.0f;
		cloudptr trans_cloud(new cloud);
		Intersection_points_src->width = Intersection_points_src->points.size();
		Intersection_points_src->height = 1;
		pcl::transformPointCloud(*Intersection_points_src, *trans_cloud, rt);
		for (int i = 0; i < trans_cloud->points.size(); i++)
		{
			rms += sqrt(pcl::squaredEuclideanDistance(trans_cloud->points[i], Intersection_points_tgt->points[i]));
		}
		rms /= trans_cloud->size();
	}

	namespace Registration
	{
		bool Regis(
			const GPSCO::cloudptr cloud_src,
			const GPSCO::cloudptr cloud_tgt,
			GPSCO::Params& params,
			Eigen::Matrix4f& RT)
		{
			if (cloud_src->empty())
			{
				spdlog::error("Source point cloud is empty.");
				return false;
			}
			if (cloud_tgt->empty())
			{
				spdlog::error("Target point cloud is empty.");
				return false;
			}

			// Planar extraction
			std::vector<GPSCO::PLANE> planes_src;
			std::vector<GPSCO::PLANE> planes_tgt;
			if (!GPSCO::PLANE_Extraction::PLANE_Tetect_RegionGrow(cloud_src, params.min_support_points,
				params.SmoothnessThreshold, params.CurvatureThreshold, planes_src))
			{
				spdlog::error("Source point cloud plane extraction failed.");
				return false;
			}
			if (!GPSCO::PLANE_Extraction::PLANE_Tetect_RegionGrow(cloud_tgt, params.min_support_points,
				params.SmoothnessThreshold, params.CurvatureThreshold, planes_tgt))
			{
				spdlog::error("Target point cloud plane extraction failed.");
				return false;
			}

			// Planar clustering
			std::vector<std::vector<std::vector<GPSCO::PLANE>>> PlaneGroups_src;
			std::vector<std::vector<std::vector<GPSCO::PLANE>>> PlaneGroups_tgt;
			if (!Plane_Cluster(planes_src, params.parallel_thresh, params.coplanar_thresh, PlaneGroups_src))
			{
				spdlog::error("Source point cloud plane clustering failed.");
				return false;
			}
			if (!Plane_Cluster(planes_tgt, params.parallel_thresh, params.coplanar_thresh, PlaneGroups_tgt))
			{
				spdlog::error("Target point cloud plane clustering failed.");
				return false;
			}

			// Export Plane Groups And Pland Clusters
//			Export_groups(PlaneGroups_src, "D:\\Code\\CLion\\GPSCO\\results\\");
//			Export_cluster(PlaneGroups_src, "D:\\Code\\CLion\\GPSCO\\results\\");

			return true;
		}

		bool Plane_Cluster(
			const std::vector<GPSCO::PLANE>& Planes,
			float parallel_thresh,
			float coplanar_thresh,
			std::vector<std::vector<std::vector<GPSCO::PLANE>>>& PlaneGroups)
		{
			if (Planes.empty())
			{
				spdlog::info("The number of planes to be clustered is empty.");
				return false;
			}

			else
			{
				// Cluster parallel planes
				PlaneGroups.clear();
				std::vector<std::vector<GPSCO::PLANE>> PlaneGroups_parallel;
				for (int i = 0; i < Planes.size(); ++i)
				{
					bool is_parallel = false;
					for (int j = 0; j < PlaneGroups_parallel.size(); ++j)
					{
						float plane_angle =
							pcl::getAngle3D(Planes[i].normal, PlaneGroups_parallel[j][0].normal, true);

						if ((plane_angle < parallel_thresh) || (std::fabs(180 - plane_angle) < parallel_thresh))
						{
							is_parallel = true;
							PlaneGroups_parallel[j].push_back(Planes[i]);
						}
					}

					// non-parallel
					if (!is_parallel)
					{
						PlaneGroups_parallel.push_back(std::vector<GPSCO::PLANE>{ Planes[i] });
					}
				}

				spdlog::info("The num of groups is {}.", PlaneGroups_parallel.size());

				// Cluster coplanar planes
				for (int i = 0; i < PlaneGroups_parallel.size(); ++i)
				{
					int planenum = PlaneGroups_parallel[i].size();
					std::vector<GPSCO::PLANE> planearray(planenum); // store the sequential index number of the plane
					// ----------------------------
					// Count the number of lines connecting the centroid of the plane
					// and the normal vector in the same direction to determine the order of the planes.
					// Take the first plane normal vector as the positive direction
					for (int j = 0; j < planenum; ++j)
					{
						int num = 0; // the Number in the same direction
						for (int k = 0; k < planenum; ++k)
						{
							if (j != k)
							{
								Eigen::Vector3f direc = PlaneGroups_parallel[i][j].centroid
									- PlaneGroups_parallel[i][k].centroid;
								if ((PlaneGroups_parallel[i][0].normal).dot(direc) > 0)
								{ num++; }
							}
						}
						planearray[num] = PlaneGroups_parallel[i][j];
					}

					std::vector<std::vector<GPSCO::PLANE>> PlaneGroups_coplanar;
					PlaneGroups_coplanar.push_back({{ planearray[0] }});
					for (int j = 1; j < planearray.size(); ++j)
					{
						bool is_coplanar = true;
						for (int k = 0; k < PlaneGroups_coplanar.back().size(); ++k)
						{
							// Is the direction of the line connecting two centroids
							// approximately perpendicular to the normal vector of the plane
							Eigen::Vector3f centroidAB =
								planearray[j].centroid - PlaneGroups_coplanar.back()[k].centroid;
							double angle1 = pcl::getAngle3D(centroidAB, planearray[j].normal, true);
							double angle2 = pcl::getAngle3D(centroidAB, PlaneGroups_coplanar.back()[k].normal, true);
							if (!(std::fabs(90 - angle1) < coplanar_thresh && std::fabs(90 - angle2) < coplanar_thresh))
							{
								is_coplanar = false;
								break;
							}
						}

						if (!is_coplanar)
						{
							PlaneGroups_coplanar.push_back({{ planearray[j] }});
						}
						else
						{
							PlaneGroups_coplanar.back().push_back(planearray[j]);
						}
					}

					PlaneGroups.push_back(PlaneGroups_coplanar);
				}

				int num_cluster = 0;
				for (int i = 0; i < PlaneGroups.size(); ++i)
				{
					num_cluster += PlaneGroups[i].size();
				}
				spdlog::info("The number of plane clusters is {}.", num_cluster);

				return true;
			}
		}

		bool Export_groups(
			std::vector<std::vector<std::vector<GPSCO::PLANE>>>& PlaneGroups,
			std::string outpath)
		{
			if (PlaneGroups.empty())
			{
				spdlog::error("The Plane Groups is empty!");
				return false;
			}
			else
			{
				std::mt19937 random;
				std::ofstream outfile; // Creating an Output File Stream Object
				for (int i = 0; i < PlaneGroups.size(); ++i)
				{
					auto rgb = random();
					auto red_ = (rgb >> 16) & 0xff;
					auto green_ = (rgb >> 8) & 0xff;
					auto blue_ = rgb & 0xff;

					std::string file = outpath + "g" + std::to_string(i + 1) + ".txt";
					outfile.open(file);
					for (const auto& cluster : PlaneGroups[i])
					{
						for (auto& plane : cluster)
						{
							for (const auto& pt : plane.points->points)
							{
								outfile << pt.x << " " << pt.y << " " << pt.z << " " <<
										red_ << " " << green_ << " " << blue_ << std::endl;
							}
						}
					}
					outfile.close();
				}
				spdlog::info("Successful export to {}", outpath);
				return true;
			}
		}

		bool Export_cluster(
			std::vector<std::vector<std::vector<GPSCO::PLANE>>>& PlaneGroups,
			std::string outpath)
		{
			if (PlaneGroups.empty())
			{
				spdlog::error("The Plane Clusters is empty!");
				return false;
			}
			else
			{
				std::mt19937 random;
				std::ofstream outfile; // Creating an Output File Stream Object
				int cluster_idx = 0;
				for (int i = 0; i < PlaneGroups.size(); ++i)
				{
					for (const auto& cluster : PlaneGroups[i])
					{
						auto rgb = random();
						auto red_ = (rgb >> 16) & 0xff;
						auto green_ = (rgb >> 8) & 0xff;
						auto blue_ = rgb & 0xff;

						cluster_idx++;
						std::string file = outpath + "cluster" + std::to_string(cluster_idx) + ".txt";
						outfile.open(file);

						for (auto& plane : cluster)
						{
							for (const auto& pt : plane.points->points)
							{
								outfile << pt.x << " " << pt.y << " " << pt.z << " " <<
										red_ << " " << green_ << " " << blue_ << std::endl;
							}
						}

						outfile.close();
					}
				}
				spdlog::info("Successful export to {}", outpath);
				return true;
			}
		}
	}
}