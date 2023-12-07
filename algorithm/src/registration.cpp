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
			std::vector<std::vector<std::vector<int>>> PlaneGroups_src;
			std::vector<std::vector<std::vector<int>>> PlaneGroups_tgt;
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

			// Initial match based moving alignment
			std::vector<std::vector<std::vector<MatchGroup>>> group_table; // The relationship table between plane groups
			Compute_Group_Table(planes_src, planes_tgt, PlaneGroups_src, PlaneGroups_tgt, params.dist_thresh, group_table);

			return true;
		}

		bool Plane_Cluster(
			const std::vector<GPSCO::PLANE>& Planes,
			float parallel_thresh,
			float coplanar_thresh,
			std::vector<std::vector<std::vector<int>>>& PlaneGroups)
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
				std::vector<std::vector<int>> PlaneGroups_parallel;
				for (int i = 0; i < Planes.size(); ++i)
				{
					bool is_parallel = false;
					for (int j = 0; j < PlaneGroups_parallel.size(); ++j)
					{
						float plane_angle =
							pcl::getAngle3D(Planes[i].normal, Planes[PlaneGroups_parallel[j][0]].normal, true);

						if ((plane_angle < parallel_thresh) || (std::fabs(180 - plane_angle) < parallel_thresh))
						{
							is_parallel = true;
							PlaneGroups_parallel[j].push_back(i);
						}
					}

					// non-parallel
					if (!is_parallel)
					{
						PlaneGroups_parallel.push_back(std::vector<int>{ i });
					}
				}

				spdlog::info("The num of groups is {}.", PlaneGroups_parallel.size());

				// Cluster coplanar planes
				for (int i = 0; i < PlaneGroups_parallel.size(); ++i)
				{
					int planenum = PlaneGroups_parallel[i].size();
					std::vector<int> planearray(planenum, -1); // store the sequential index number of the plane
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
								Eigen::Vector3f direc = Planes[PlaneGroups_parallel[i][j]].centroid
									- Planes[PlaneGroups_parallel[i][k]].centroid;
								if ((Planes[PlaneGroups_parallel[i][0]].normal).dot(direc) > 0)
								{ num++; }
							}
						}
						planearray[num] = PlaneGroups_parallel[i][j];
					}

					std::vector<std::vector<int>> PlaneGroups_coplanar;
					PlaneGroups_coplanar.push_back({{ planearray[0] }});
					for (int j = 1; j < planearray.size(); ++j)
					{
						bool is_coplanar = true;
						for (int k = 0; k < PlaneGroups_coplanar.back().size(); ++k)
						{
							// Is the direction of the line connecting two centroids
							// approximately perpendicular to the normal vector of the plane
							Eigen::Vector3f centroidAB =
								Planes[planearray[j]].centroid - Planes[PlaneGroups_coplanar.back()[k]].centroid;
							double angle1 = pcl::getAngle3D(centroidAB, Planes[planearray[j]].normal, true);
							double angle2 = pcl::getAngle3D(centroidAB, Planes[PlaneGroups_coplanar.back()[k]].normal, true);
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
			const std::vector<GPSCO::PLANE>& Planes,
			std::vector<std::vector<std::vector<int>>>& PlaneGroups,
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
						for (auto& planeidx : cluster)
						{
							for (const auto& pt : Planes[planeidx].points->points)
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
			const std::vector<GPSCO::PLANE>& Planes,
			std::vector<std::vector<std::vector<int>>>& PlaneGroups,
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

						for (auto& planeidx : cluster)
						{
							for (const auto& pt : Planes[planeidx].points->points)
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

		bool Compute_Group_Table(
			const std::vector<GPSCO::PLANE>& Planes_src,
			const std::vector<GPSCO::PLANE>& Planes_tgt,
			const std::vector<std::vector<std::vector<int>>>& PlaneGroups_src,
			const std::vector<std::vector<std::vector<int>>>& PlaneGroups_tgt,
			float dist_thresh,
			std::vector<std::vector<std::vector<MatchGroup>>>& group_table)
		{
			if (PlaneGroups_src.size() < 3 || PlaneGroups_tgt.size() < 3)
			{
				if (PlaneGroups_src.empty() || PlaneGroups_tgt.empty())
				{
					spdlog::error("Plane Groups is empty.");
					return false;
				}
				else
				{
					spdlog::error("Insufficient number of plane groups.");
					return false;
				}
			}
			else
			{
				// initialization for group_table
				group_table.clear();
				// Resize the outermost dimension to have PlaneGroups_src.size() elements
				group_table.resize(PlaneGroups_src.size());
				// Resize each inner vector to have PlaneGroups_tgt.size() elements
				for (auto& row : group_table)
				{
					row.resize(PlaneGroups_tgt.size());
				}

				// Use the moving alignment method to estimate the matching scores between plane groups
				// and record the correspondence between the planes
				for (int i = 0; i < PlaneGroups_src.size(); ++i)
				{
					for (int j = 0; j < PlaneGroups_tgt.size(); ++j)
					{
						std::vector<MatchGroup> match_vec;
						if (Moving_alignment(i, j, Planes_src, Planes_tgt,
							PlaneGroups_src[i], PlaneGroups_tgt[j], dist_thresh, match_vec))
						{
							group_table[i][j] = match_vec;
						}
						else
						{
							spdlog::error("Moving_alignment failure!");
							return false;
						}
					}
				}

				spdlog::info("The relationship table between plane groups is successfully obtained.");
				return true;
			}
		}

		bool Moving_alignment(
			int group1_index,
			int group2_index,
			const std::vector<GPSCO::PLANE>& Planes_src,
			const std::vector<GPSCO::PLANE>& Planes_tgt,
			const std::vector<std::vector<int>>& planegroup_src,
			const std::vector<std::vector<int>>& planegroup_tgt,
			float dist_thresh,
			std::vector<MatchGroup>& match_vec)
		{
			if (Planes_src.empty() || Planes_tgt.empty()
				|| planegroup_src.empty() || planegroup_tgt.empty() || dist_thresh <= 0)
			{
				spdlog::error("There was a problem entering data.");
				return false;
			}

			int num_src = planegroup_src.size();
			int num_tgt = planegroup_tgt.size();

			if (num_src == 1)
			{
				for (int i = 0; i < num_tgt; ++i)
				{
					MatchGroup matchpair;
					matchpair.group1_index = group1_index;
					matchpair.group2_index = group2_index;
					matchpair.planepairs.push_back(std::make_pair(0, i));
					matchpair.match_num = 1;
					match_vec.push_back(matchpair);
				}
				return true;
			}
			if (num_tgt == 1)
			{
				for (int i = 0; i < num_src; ++i)
				{
					MatchGroup matchpair;
					matchpair.group1_index = group1_index;
					matchpair.group2_index = group2_index;
					matchpair.planepairs.push_back(std::make_pair(i, 0));
					matchpair.match_num = 1;
					match_vec.push_back(matchpair);
				}
				return true;
			}

			// Calculate the average normal vector of the plane group
			// to avoid individual plane normal vectors affecting the distance calculation
			Eigen::Vector3f Normal_avg_src(0.0, 0.0, 0.0), Normal_avg_tgt(0.0, 0.0, 0.0);
			int n_planesrc, n_planetgt;
			n_planesrc = n_planetgt = 0;
			for (const auto& cluster : planegroup_src)
			{
				for (const auto& planeidx : cluster)
				{
					n_planesrc++;
					if (Normal_avg_src.dot(Planes_src[planeidx].normal) >= 0)
						Normal_avg_src += Planes_src[planeidx].normal;
					else
						Normal_avg_src -= Planes_src[planeidx].normal;
				}
			}
			Normal_avg_src /= n_planesrc;
			Normal_avg_src.normalized();

			for (const auto& cluster : planegroup_tgt)
			{
				for (const auto& planeidx : cluster)
				{
					n_planetgt++;
					if (Normal_avg_tgt.dot(Planes_tgt[planeidx].normal) >= 0)
						Normal_avg_tgt += Planes_tgt[planeidx].normal;
					else
						Normal_avg_tgt -= Planes_tgt[planeidx].normal;
				}
			}
			Normal_avg_tgt /= n_planetgt;
			Normal_avg_tgt.normalized();

			// Apply for a dynamic array as a plane distance table
			std::vector<std::vector<float>> d_src(num_src, std::vector<float>(num_src, 0));
			std::vector<std::vector<float>> d_tgt(num_tgt, std::vector<float>(num_tgt, 0));

			for (int i = 0; i < num_src; i++)
			{
				for (int j = 0; j < num_src; j++)
				{
					if (i < j)
					{
						d_src[i][j] = d_src[j][i]
							= GPSCO::Registration::dist_TwoClsuter(
								Planes_src,
								Normal_avg_src,
								planegroup_src[i],
								planegroup_src[j]);
					}
				}
			}
			for (int i = 0; i < num_tgt; i++)
			{
				for (int j = 0; j < num_tgt; j++)
				{
					if (i < j)
					{
						d_tgt[i][j] = d_tgt[j][i]
							= GPSCO::Registration::dist_TwoClsuter(
								Planes_tgt,
								Normal_avg_tgt,
								planegroup_tgt[i],
								planegroup_tgt[j]);
					}
				}
			}

			// distance arrays
			std::vector<std::vector<float>> disarrays_src;
			std::vector<std::vector<float>> disarrays_tgt;

			for (int i = 0; i < num_src; i++)
			{
				std::vector<float> disarray_src;
				for (int j = 0; j < num_src; j++)
				{
					if (i < j)
					{
						disarray_src.push_back(d_src[i][j]);
					}
					else if (i > j)
					{
						disarray_src.push_back(-d_src[i][j]);
					}
				}
				disarrays_src.push_back(disarray_src);
			}
			for (int i = 0; i < num_tgt; i++)
			{
				std::vector<float> disarray_tgt;
				for (int j = 0; j < num_tgt; j++)
				{
					if (i < j)
					{
						disarray_tgt.push_back(d_tgt[i][j]);
					}
					else if (i > j)
					{
						disarray_tgt.push_back(-d_tgt[i][j]);
					}
				}
				disarrays_tgt.push_back(disarray_tgt);
			}

			// moving alignment
			for (int i = 0; i < num_src; i++)
			{
				for (int j = 0; j < num_tgt; j++)
				{
					// src_i aligns with tgt_j
					int src_index = -1;
					int tgt_index = -1;
					MatchGroup matchpair;
					matchpair.planepairs.push_back(std::make_pair(i, j));
					for (int m = 0; m < disarrays_src[i].size(); m++)
					{
						for (int n = 0; n < disarrays_tgt[j].size(); n++)
						{
							if (std::abs(disarrays_src[i][m] - disarrays_tgt[j][n]) <= dist_thresh)
							{
								if (m < i)
								{ src_index = m; }
								else
								{ src_index = m + 1; }
								if (n < j)
								{ tgt_index = n; }
								else
								{ tgt_index = n + 1; }

								matchpair.planepairs.push_back(std::make_pair(src_index, tgt_index));
							}
						}
					}
					match_vec.push_back(matchpair);
				}
			}

			// Do it again in reverse
			disarrays_tgt.clear();
			for (int i = 0; i < num_tgt; i++)
			{
				std::vector<float> disarray_tgt;
				for (int j = 0; j < num_tgt; j++)
				{
					if (i < j)
					{
						disarray_tgt.push_back(-d_tgt[i][j]);
					}
					else if (i > j)
					{
						disarray_tgt.push_back(d_tgt[i][j]);
					}
				}
				disarrays_tgt.push_back(disarray_tgt);
			}
			for (int i = 0; i < num_src; i++)
			{
				for (int j = 0; j < num_tgt; j++)
				{
					// src_i aligns with tgt_j
					int src_index = -1;
					int tgt_index = -1;
					MatchGroup matchpair;
					matchpair.planepairs.push_back(std::make_pair(i, j));
					for (int m = 0; m < disarrays_src[i].size(); m++)
					{
						for (int n = 0; n < disarrays_tgt[j].size(); n++)
						{
							if (std::abs(disarrays_src[i][m] - disarrays_tgt[j][n]) <= dist_thresh)
							{
								if (m < i)
								{ src_index = m; }
								else
								{ src_index = m + 1; }
								if (n < j)
								{ tgt_index = n; }
								else
								{ tgt_index = n + 1; }

								matchpair.planepairs.push_back(std::make_pair(src_index, tgt_index));
							}
						}
					}
					match_vec.push_back(matchpair);
				}
			}

			// Decompose the results that do not meet the conditions to avoid interfering with the score judgment
			std::vector<MatchGroup> match_vec_temp;
			for (const auto& match_temp : match_vec)
			{
				if (match_temp.planepairs.size() > 1)
				{
					std::vector<std::vector<std::pair<int, int>>> vec_;
					vec_.push_back(match_temp.planepairs);
					// Move planepairs[0] to the correct sorting position
					auto pair0 = vec_[0][0];
					vec_[0].erase(vec_[0].begin());
					if (pair0.first >= vec_[0].back().first)
						vec_[0].push_back(pair0);
					else
					{
						for (int loc = 0; loc < vec_[0].size(); loc++)
						{
							if (pair0.first <= vec_[0][loc].first)
							{
								vec_[0].insert(vec_[0].begin() + loc, pair0);
								break;
							}
						}
					}

					bool IsSplit = true;
					while (IsSplit)
					{
						IsSplit = false;
						int vec_num = vec_.size();
						for (int vec_idx = 0; vec_idx < vec_num; vec_idx++)
						{
							for (int pair_idx = 0; pair_idx < vec_[vec_idx].size() - 1; pair_idx++)
							{
								// Split
								if (d_src[vec_[vec_idx][pair_idx].first][vec_[vec_idx][pair_idx + 1].first] < 0.3 ||
									d_tgt[vec_[vec_idx][pair_idx].second][vec_[vec_idx][pair_idx + 1].second] < 0.3)
								{
									IsSplit = true;
									vec_.push_back(vec_[vec_idx]);
									vec_[vec_idx].erase(vec_[vec_idx].begin() + pair_idx);
									vec_.back().erase(vec_.back().begin() + pair_idx + 1);

									// Jump out of double loop
									pair_idx = vec_[vec_idx].size();
									vec_idx = vec_num;
								}
							}
						}
					}

					for (const auto& vec : vec_)
					{
						MatchGroup temp;
						temp.group1_index = match_temp.group1_index;
						temp.group2_index = match_temp.group2_index;
						temp.planepairs = vec;
						temp.match_num = vec.size();
						match_vec_temp.push_back(temp);
					}
				}
				else
				{
					match_vec_temp.push_back(match_temp);
				}
			}

			// Remove duplicates, such as (1-2, 2-4) and (2-4, 1-2)
			std::set<std::vector<std::pair<int, int>>> unique_combinations;
			for (auto& matchgroup_vec_temp_ : match_vec_temp)
			{
				std::sort(matchgroup_vec_temp_.planepairs.begin(), matchgroup_vec_temp_.planepairs.end());
				unique_combinations.insert(matchgroup_vec_temp_.planepairs);
			}

			match_vec.clear();
			for (const auto& combination : unique_combinations)
			{
				MatchGroup matchpair;
				matchpair.group1_index = group1_index;
				matchpair.group2_index = group2_index;
				matchpair.planepairs = combination;
				matchpair.match_num = combination.size();
				match_vec.push_back(matchpair);
			}

			std::sort(match_vec.begin(), match_vec.end());

			return true;
		}

		float dist_TwoClsuter(
			const std::vector<GPSCO::PLANE>& Planes,
			const Eigen::Vector3f Normal_avg,
			const std::vector<int>& cluster1,
			const std::vector<int>& cluster2)
		{
			float dist = 0.0;
			for (int i = 0; i < cluster1.size(); i++)
			{
				for (int j = 0; j < cluster2.size(); j++)
				{
					float dist1 = GPSCO::distancePointToPlane(
						Planes[cluster1[i]].centroid,
						Normal_avg,
						-Normal_avg.dot(Planes[cluster2[j]].centroid));
					float dist2 = GPSCO::distancePointToPlane(
						Planes[cluster2[j]].centroid,
						Normal_avg,
						-Normal_avg.dot(Planes[cluster1[i]].centroid));
					dist += (dist1 + dist2);
				}
			}
			dist /= (2 * cluster1.size() * cluster2.size());

			return dist;
		}

		bool Get_Group_Three(
			const std::vector<GPSCO::PLANE>& Planes,
			const std::vector<std::vector<std::vector<int>>>& PlaneGroups,
			float angle_min,
			float angle_max,
			std::vector<GPSCO::Group_Three>& group_three_vector)
		{
			if (Planes.empty())
			{
				return false;
			}

			group_three_vector.clear();
			std::vector<GPSCO::Group_Base> group_base;
			// Get plane group pairs that meet the threshold
			for (int i = 0; i < PlaneGroups.size(); i++)
			{
				for (int j = 0; j < PlaneGroups.size(); j++)
				{
					if (i < j)
					{
						float group_angle = pcl::getAngle3D(Planes[PlaneGroups[i][0][0]].normal,
							Planes[PlaneGroups[j][0][0]].normal, true);
						if (angle_min < group_angle && group_angle < angle_max)
						{
							GPSCO::Group_Base group_base_temp;
							group_base_temp.group1_index = i;
							group_base_temp.group2_index = j;
							// make the angle acute
							if (group_angle < 90)
								group_base_temp.group_angle = group_angle;
							else
								group_base_temp.group_angle = 180 - group_angle;
							group_base.push_back(group_base_temp);
						}
					}
				}
			}
			if (group_base.empty())
			{
				spdlog::error("The number of plane group pairs that meet the threshold is 0.");
				return false;
			}

			// Get three plane group that meet the threshold
			// Hash checks whether it already exists to avoid duplication
			std::unordered_map<Array3I, int, HashArray3I> hash;
			for (const auto& group_base_temp : group_base)
			{
				for (int i = 0; i < PlaneGroups.size(); i++)
				{
					if (group_base_temp.group1_index != i && group_base_temp.group2_index != i)
					{
						Array3I idx;
						idx[0] = group_base_temp.group1_index;
						idx[1] = group_base_temp.group2_index;
						idx[2] = i;
						std::sort(idx.begin(), idx.end());
						// Never appeared
						if (hash.count(idx) == 0)
						{
							hash.insert(std::make_pair(idx, 1));
							auto plane_index1 = PlaneGroups[group_base_temp.group1_index][0][0];
							auto plane_index2 = PlaneGroups[group_base_temp.group2_index][0][0];
							float angle1 = pcl::getAngle3D(Planes[plane_index1].normal,
								Planes[PlaneGroups[i][0][0]].normal, true);
							float angle2 = pcl::getAngle3D(Planes[plane_index2].normal,
								Planes[PlaneGroups[i][0][0]].normal, true);
							if ((angle_min < angle1 && angle1 < angle_max) &&
								(angle_min < angle2 && angle2 < angle_max))
							{
								Group_Three group_three_temp;
								group_three_temp.group1_index = group_base_temp.group1_index;
								group_three_temp.group2_index = group_base_temp.group2_index;
								group_three_temp.group3_index = i;
								group_three_temp.group_angle12 = group_base_temp.group_angle;
								if (angle1 < 90)
									group_three_temp.group_angle13 = angle1;
								else
									group_three_temp.group_angle13 = 180 - angle1;
								if (angle2 < 90)
									group_three_temp.group_angle23 = angle2;
								else
									group_three_temp.group_angle23 = 180 - angle2;

								group_three_vector.push_back(group_three_temp);
							}
						}
					}
				}
			}

			if (group_three_vector.empty())
			{
				spdlog::error("The number of three plane groups that meet the threshold is 0.");
				return false;
			}

			return true;
		}

		bool Get_transformation_matrix(
			const std::vector<GPSCO::PLANE>& Planes_src,
			const std::vector<GPSCO::PLANE>& Planes_tgt,
			const std::vector<std::vector<std::vector<int>>>& PlaneGroups_src,
			const std::vector<std::vector<std::vector<int>>>& PlaneGroups_tgt,
			std::vector<std::vector<std::vector<MatchGroup>>>& group_table,
			const std::vector<GPSCO::Group_Three>& group_three_vector_src,
			const std::vector<GPSCO::Group_Three>& group_three_vector_tgt,
			Eigen::Matrix4f& final_rt)
		{
			float angle_thresh = 5.0; // Conditions for matching between two candidates
			std::vector<GPSCO::Group_Three_Pair> group_three_pair_vec; // Three pairs of non-parallel plane groups
			for (const auto& group_three_src : group_three_vector_src)
			{
				for (const auto& group_three_tgt : group_three_vector_tgt)
				{
					// Six possible matches
					// type 1: 1-1,2-2,3-3
					if (std::fabs(group_three_src.group_angle12 - group_three_tgt.group_angle12) < angle_thresh &&
						std::fabs(group_three_src.group_angle13 - group_three_tgt.group_angle13) < angle_thresh &&
						std::fabs(group_three_src.group_angle23 - group_three_tgt.group_angle23) < angle_thresh)
					{
						GPSCO::Group_Three_Pair group_three_pair_temp;
						group_three_pair_temp.group_src_index1 = group_three_src.group1_index;
						group_three_pair_temp.group_src_index2 = group_three_src.group2_index;
						group_three_pair_temp.group_src_index3 = group_three_src.group3_index;
						group_three_pair_temp.group_tgt_index1 = group_three_tgt.group1_index;
						group_three_pair_temp.group_tgt_index2 = group_three_tgt.group2_index;
						group_three_pair_temp.group_tgt_index3 = group_three_tgt.group3_index;
						auto score1 = group_table[group_three_pair_temp.group_src_index1][group_three_pair_temp.group_tgt_index1][0].match_num;
						auto score2 = group_table[group_three_pair_temp.group_src_index2][group_three_pair_temp.group_tgt_index2][0].match_num;
						auto score3 = group_table[group_three_pair_temp.group_src_index3][group_three_pair_temp.group_tgt_index3][0].match_num;
						group_three_pair_temp.match_score = score1 + score2 + score3;
						group_three_pair_vec.push_back(group_three_pair_temp);
					}
					// type 2: 1-1,2-3,3-2
					if (std::fabs(group_three_src.group_angle12 - group_three_tgt.group_angle13) < angle_thresh &&
						std::fabs(group_three_src.group_angle13 - group_three_tgt.group_angle12) < angle_thresh &&
						std::fabs(group_three_src.group_angle23 - group_three_tgt.group_angle23) < angle_thresh)
					{
						GPSCO::Group_Three_Pair group_three_pair_temp;
						group_three_pair_temp.group_src_index1 = group_three_src.group1_index;
						group_three_pair_temp.group_src_index2 = group_three_src.group2_index;
						group_three_pair_temp.group_src_index3 = group_three_src.group3_index;
						group_three_pair_temp.group_tgt_index1 = group_three_tgt.group1_index;
						group_three_pair_temp.group_tgt_index2 = group_three_tgt.group3_index;
						group_three_pair_temp.group_tgt_index3 = group_three_tgt.group2_index;
						auto score1 = group_table[group_three_pair_temp.group_src_index1][group_three_pair_temp.group_tgt_index1][0].match_num;
						auto score2 = group_table[group_three_pair_temp.group_src_index2][group_three_pair_temp.group_tgt_index2][0].match_num;
						auto score3 = group_table[group_three_pair_temp.group_src_index3][group_three_pair_temp.group_tgt_index3][0].match_num;
						group_three_pair_temp.match_score = score1 + score2 + score3;
						group_three_pair_vec.push_back(group_three_pair_temp);
					}
					// type 3: 1-2,2-1,3-3
					if (std::fabs(group_three_src.group_angle12 - group_three_tgt.group_angle12) < angle_thresh &&
						std::fabs(group_three_src.group_angle13 - group_three_tgt.group_angle23) < angle_thresh &&
						std::fabs(group_three_src.group_angle23 - group_three_tgt.group_angle13) < angle_thresh)
					{
						GPSCO::Group_Three_Pair group_three_pair_temp;
						group_three_pair_temp.group_src_index1 = group_three_src.group1_index;
						group_three_pair_temp.group_src_index2 = group_three_src.group2_index;
						group_three_pair_temp.group_src_index3 = group_three_src.group3_index;
						group_three_pair_temp.group_tgt_index1 = group_three_tgt.group3_index;
						group_three_pair_temp.group_tgt_index2 = group_three_tgt.group1_index;
						group_three_pair_temp.group_tgt_index3 = group_three_tgt.group2_index;
						auto score1 = group_table[group_three_pair_temp.group_src_index1][group_three_pair_temp.group_tgt_index1][0].match_num;
						auto score2 = group_table[group_three_pair_temp.group_src_index2][group_three_pair_temp.group_tgt_index2][0].match_num;
						auto score3 = group_table[group_three_pair_temp.group_src_index3][group_three_pair_temp.group_tgt_index3][0].match_num;
						group_three_pair_temp.match_score = score1 + score2 + score3;
						group_three_pair_vec.push_back(group_three_pair_temp);
					}
					// type 4: 1-2,2-3,3-1
					if (std::fabs(group_three_src.group_angle12 - group_three_tgt.group_angle23) < angle_thresh &&
						std::fabs(group_three_src.group_angle13 - group_three_tgt.group_angle12) < angle_thresh &&
						std::fabs(group_three_src.group_angle23 - group_three_tgt.group_angle13) < angle_thresh)
					{
						GPSCO::Group_Three_Pair group_three_pair_temp;
						group_three_pair_temp.group_src_index1 = group_three_src.group1_index;
						group_three_pair_temp.group_src_index2 = group_three_src.group2_index;
						group_three_pair_temp.group_src_index3 = group_three_src.group3_index;
						group_three_pair_temp.group_tgt_index1 = group_three_tgt.group2_index;
						group_three_pair_temp.group_tgt_index2 = group_three_tgt.group3_index;
						group_three_pair_temp.group_tgt_index3 = group_three_tgt.group1_index;
						auto score1 = group_table[group_three_pair_temp.group_src_index1][group_three_pair_temp.group_tgt_index1][0].match_num;
						auto score2 = group_table[group_three_pair_temp.group_src_index2][group_three_pair_temp.group_tgt_index2][0].match_num;
						auto score3 = group_table[group_three_pair_temp.group_src_index3][group_three_pair_temp.group_tgt_index3][0].match_num;
						group_three_pair_temp.match_score = score1 + score2 + score3;
						group_three_pair_vec.push_back(group_three_pair_temp);
					}
					// type 5: 1-3,2-1,3-2
					if (std::fabs(group_three_src.group_angle12 - group_three_tgt.group_angle13) < angle_thresh &&
						std::fabs(group_three_src.group_angle13 - group_three_tgt.group_angle23) < angle_thresh &&
						std::fabs(group_three_src.group_angle23 - group_three_tgt.group_angle12) < angle_thresh)
					{
						GPSCO::Group_Three_Pair group_three_pair_temp;
						group_three_pair_temp.group_src_index1 = group_three_src.group1_index;
						group_three_pair_temp.group_src_index2 = group_three_src.group2_index;
						group_three_pair_temp.group_src_index3 = group_three_src.group3_index;
						group_three_pair_temp.group_tgt_index1 = group_three_tgt.group3_index;
						group_three_pair_temp.group_tgt_index2 = group_three_tgt.group1_index;
						group_three_pair_temp.group_tgt_index3 = group_three_tgt.group2_index;
						auto score1 = group_table[group_three_pair_temp.group_src_index1][group_three_pair_temp.group_tgt_index1][0].match_num;
						auto score2 = group_table[group_three_pair_temp.group_src_index2][group_three_pair_temp.group_tgt_index2][0].match_num;
						auto score3 = group_table[group_three_pair_temp.group_src_index3][group_three_pair_temp.group_tgt_index3][0].match_num;
						group_three_pair_temp.match_score = score1 + score2 + score3;
						group_three_pair_vec.push_back(group_three_pair_temp);
					}
					// type 6: 1-3,2-2,3-1
					if (std::fabs(group_three_src.group_angle12 - group_three_tgt.group_angle23) < angle_thresh &&
						std::fabs(group_three_src.group_angle13 - group_three_tgt.group_angle13) < angle_thresh &&
						std::fabs(group_three_src.group_angle23 - group_three_tgt.group_angle12) < angle_thresh)
					{
						GPSCO::Group_Three_Pair group_three_pair_temp;
						group_three_pair_temp.group_src_index1 = group_three_src.group1_index;
						group_three_pair_temp.group_src_index2 = group_three_src.group2_index;
						group_three_pair_temp.group_src_index3 = group_three_src.group3_index;
						group_three_pair_temp.group_tgt_index1 = group_three_tgt.group3_index;
						group_three_pair_temp.group_tgt_index2 = group_three_tgt.group2_index;
						group_three_pair_temp.group_tgt_index3 = group_three_tgt.group1_index;
						auto score1 = group_table[group_three_pair_temp.group_src_index1][group_three_pair_temp.group_tgt_index1][0].match_num;
						auto score2 = group_table[group_three_pair_temp.group_src_index2][group_three_pair_temp.group_tgt_index2][0].match_num;
						auto score3 = group_table[group_three_pair_temp.group_src_index3][group_three_pair_temp.group_tgt_index3][0].match_num;
						group_three_pair_temp.match_score = score1 + score2 + score3;
						group_three_pair_vec.push_back(group_three_pair_temp);
					}
				}
			}

			// Sorted by match score, highest to lowest
			std::sort(group_three_pair_vec.begin(), group_three_pair_vec.end());

			return true;
		}
	}
}