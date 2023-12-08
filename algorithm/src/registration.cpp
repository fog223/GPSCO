// system
#include <unordered_map>
#include <set>
#include <chrono>
#include <random>
// PCL
#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl/io/ply_io.h>
#include <pcl/common/io.h>
#include <pcl/common/common.h>
#include <pcl/common/distances.h>
#include <pcl/common/intersections.h>
#include <pcl/registration/transformation_estimation_svd.h>
#include <pcl/search/kdtree.h>
#include <pcl/octree/octree.h>
#include <pcl/filters/uniform_sampling.h>
// Eigen
#include <Eigen/Dense>
// spdlog
#include <spdlog/spdlog.h>
// local
#include "Common.h"
#include "PLANE_Extraction.h"
#include "Registration.h"

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

		match_score = 0.0f;
	}

	Group_Three_Pair::~Group_Three_Pair()
	{
	}

	RT_Info::RT_Info()
	{
		Intersection_points_src.reset(new GPSCO::cloud);
		Intersection_points_tgt.reset(new GPSCO::cloud);
		confidence = 0.0f;
		Dmean = 0.0f;
	}

	RT_Info::~RT_Info()
	{
	}

	void RT_Info::Compute_Dmean()
	{
		Dmean = 0.0f;
		cloudptr trans_cloud(new cloud);
		Intersection_points_src->width = Intersection_points_src->points.size();
		Intersection_points_src->height = 1;
		pcl::transformPointCloud(*Intersection_points_src, *trans_cloud, rt);
		for (int i = 0; i < trans_cloud->points.size(); i++)
		{
			Dmean += sqrt(pcl::squaredEuclideanDistance(trans_cloud->points[i], Intersection_points_tgt->points[i]));
		}
		Dmean /= trans_cloud->size();
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

			// Get candidate
			float angle_min = 30;
			float angle_max = 150;
			std::vector<GPSCO::Group_Three> group_three_vector_src;
			std::vector<GPSCO::Group_Three> group_three_vector_tgt;
			Get_Group_Three(planes_src, PlaneGroups_src, angle_min, angle_max, group_three_vector_src);
			Get_Group_Three(planes_tgt, PlaneGroups_tgt, angle_min, angle_max, group_three_vector_tgt);

			// GPSCO determines the transformation matrix
			if (Get_transformation_matrix(planes_src, planes_tgt, PlaneGroups_src, PlaneGroups_tgt,
				group_table, group_three_vector_src, group_three_vector_tgt, RT))
			{
				spdlog::info("Transformation matrix obtained successfully.");
				return true;
			}
			else
			{
				spdlog::info("Transformation Matrix obtained failed.");
				return false;
			}
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

			// Global Planar Structural Constraint Optimal(highest score)
			std::vector<GPSCO::RT_Info> RT_vector;
			for (int i = 0; i < group_three_pair_vec.size(); i++)
			{
				if (group_three_pair_vec[i].match_score == group_three_pair_vec[0].match_score)
				{
					int index1_src = group_three_pair_vec[i].group_src_index1;
					int index1_tgt = group_three_pair_vec[i].group_tgt_index1;
					for (int j = 0; j < group_table[index1_src][index1_tgt].size(); j++)
					{
						if (group_table[index1_src][index1_tgt][j].match_num == group_table[index1_src][index1_tgt][0].match_num)
						{
							int index2_src = group_three_pair_vec[i].group_src_index2;
							int index2_tgt = group_three_pair_vec[i].group_tgt_index2;
							for (int m = 0; m < group_table[index2_src][index2_tgt].size(); m++)
							{
								if (group_table[index2_src][index2_tgt][m].match_num == group_table[index2_src][index2_tgt][0].match_num)
								{
									int index3_src = group_three_pair_vec[i].group_src_index3;
									int index3_tgt = group_three_pair_vec[i].group_tgt_index3;
									for (int n = 0; n < group_table[index3_src][index3_tgt].size(); n++)
									{
										if (group_table[index3_src][index3_tgt][n].match_num == group_table[index3_src][index3_tgt][0].match_num)
										{
											RT_Info rt_info;
											rt_info.plane_pairs1 = group_table[index1_src][index1_tgt][j];
											rt_info.plane_pairs2 = group_table[index2_src][index2_tgt][m];
											rt_info.plane_pairs3 = group_table[index3_src][index3_tgt][n];
											RT_vector.push_back(rt_info);
										}
										else
										{ break; }
									}
								}
								else
								{ break; }
							}
						}
						else
						{ break; }
					}
				}
				else
				{ break; }
			}

			if (RT_vector.empty())
			{
				spdlog::error("RT_vector is empty.");
				return false;
			}

			// Verification
			std::vector<Eigen::Matrix4f> rt_bucket; // Clustering for transformation matrices to avoid repeated validations
			for (auto& rt_temp : RT_vector)
			{
				// Planar correspondences to point-pair correspondences
				for (const auto& pair1 : rt_temp.plane_pairs1.planepairs)
				{
					for (const auto& pair2 : rt_temp.plane_pairs2.planepairs)
					{
						for (const auto& pair3 : rt_temp.plane_pairs3.planepairs)
						{
							Eigen::Vector3f inter_pt_src(0.0, 0.0, 0.0), inter_pt_tgt(0.0, 0.0, 0.0);
							Eigen::Vector4f plane_a_src, plane_b_src, plane_c_src;
							Eigen::Vector4f plane_a_tgt, plane_b_tgt, plane_c_tgt;
							std::vector<Eigen::Vector3f> intersection_points_src; // Used for averaging
							std::vector<Eigen::Vector3f> intersection_points_tgt; // Used for averaging
							for (const auto& plane_a_index : PlaneGroups_src[rt_temp.plane_pairs1.group1_index][pair1.first])
							{
								plane_a_src << Planes_src[plane_a_index].coefficients[0],
									Planes_src[plane_a_index].coefficients[1],
									Planes_src[plane_a_index].coefficients[2],
									Planes_src[plane_a_index].coefficients[3];
								for (const auto& plane_b_index : PlaneGroups_src[rt_temp.plane_pairs2.group1_index][pair2.first])
								{
									plane_b_src << Planes_src[plane_b_index].coefficients[0],
										Planes_src[plane_b_index].coefficients[1],
										Planes_src[plane_b_index].coefficients[2],
										Planes_src[plane_b_index].coefficients[3];
									for (const auto& plane_c_index : PlaneGroups_src[rt_temp.plane_pairs3.group1_index][pair3.first])
									{
										plane_c_src << Planes_src[plane_c_index].coefficients[0],
											Planes_src[plane_c_index].coefficients[1],
											Planes_src[plane_c_index].coefficients[2],
											Planes_src[plane_c_index].coefficients[3];

										if (pcl::threePlanesIntersection(plane_a_src, plane_b_src, plane_c_src, inter_pt_src))
										{
											intersection_points_src.push_back(inter_pt_src);
										}
									}
								}
							}
							for (const auto& plane_a_index : PlaneGroups_tgt[rt_temp.plane_pairs1.group2_index][pair1.second])
							{
								plane_a_tgt << Planes_tgt[plane_a_index].coefficients[0],
									Planes_tgt[plane_a_index].coefficients[1],
									Planes_tgt[plane_a_index].coefficients[2],
									Planes_tgt[plane_a_index].coefficients[3];
								for (const auto& plane_b_index : PlaneGroups_tgt[rt_temp.plane_pairs2.group2_index][pair2.second])
								{
									plane_b_tgt << Planes_tgt[plane_b_index].coefficients[0],
										Planes_tgt[plane_b_index].coefficients[1],
										Planes_tgt[plane_b_index].coefficients[2],
										Planes_tgt[plane_b_index].coefficients[3];
									for (const auto& plane_c_index : PlaneGroups_tgt[rt_temp.plane_pairs3.group2_index][pair3.second])
									{
										plane_c_tgt << Planes_tgt[plane_c_index].coefficients[0],
											Planes_tgt[plane_c_index].coefficients[1],
											Planes_tgt[plane_c_index].coefficients[2],
											Planes_tgt[plane_c_index].coefficients[3];

										if (pcl::threePlanesIntersection(plane_a_tgt, plane_b_tgt, plane_c_tgt, inter_pt_tgt))
										{
											intersection_points_tgt.push_back(inter_pt_tgt);
										}
									}
								}
							}
							if (!intersection_points_src.empty() && !intersection_points_tgt.empty())
							{
								// Average
								inter_pt_src.setZero();
								for (int num = 0; num < intersection_points_src.size(); num++)
								{
									inter_pt_src += intersection_points_src[num];
								}
								inter_pt_src /= intersection_points_src.size();
								rt_temp.Intersection_points_src->push_back(
									pcl::PointXYZ(inter_pt_src[0], inter_pt_src[1], inter_pt_src[2]));

								inter_pt_tgt.setZero();
								for (int num = 0; num < intersection_points_tgt.size(); num++)
								{
									inter_pt_tgt += intersection_points_tgt[num];
								}
								inter_pt_tgt /= intersection_points_tgt.size();
								rt_temp.Intersection_points_tgt->push_back(
									pcl::PointXYZ(inter_pt_tgt[0], inter_pt_tgt[1], inter_pt_tgt[2]));
							}
						}
					}
				}

				if (rt_temp.Intersection_points_src->size() != rt_temp.Intersection_points_tgt->size())
				{
					spdlog::error("Inconsistent number of match points.");
					return false;
				}
				else if (rt_temp.Intersection_points_src->size() < 3)
				{
					// Insufficient number of point pairs,
					// only matching planes can be applied to compute transformation matrices
					// Two pairs of corresponding planes can determine the rotation matrix
					// Three pairs of corresponding planes can uniquely determine a transformation matrix
					spdlog::warn("Insufficient number of matching point pairs.");
					return false;
					// ...(Add more when needed, laughing)
				}
				else
				{
					// SVD
					pcl::registration::TransformationEstimationSVD<pcl::PointXYZ, pcl::PointXYZ> SVD;
					SVD.estimateRigidTransformation(
						*rt_temp.Intersection_points_src,
						*rt_temp.Intersection_points_tgt,
						rt_temp.rt);

					// Point-to-point average distance for constraints.
					rt_temp.Compute_Dmean();
					if (rt_temp.Dmean > 0.1)
						continue;

					// Avoiding repetitive transformation matrix calculations
					bool IsCluster = false;
					for (const auto& rt : rt_bucket)
					{
						Eigen::Matrix3f R = rt.block<3, 3>(0, 0);
						Eigen::Vector3f T = rt.block<3, 1>(0, 3);

						Eigen::AngleAxisf aaDiff(R.transpose() * rt_temp.rt.block<3, 3>(0, 0));
						float angleDiff = aaDiff.angle();

						// Convert radians to angles
						float angleDiffDegrees = angleDiff * 180.0 / M_PI;

						if (angleDiffDegrees < 5.0)
						{
							// Calculate the distance difference
							if ((T - rt_temp.rt.block<3, 1>(0, 3)).norm() < 0.2)
							{
								IsCluster = true;
								break;
							}
						}
					}
					if (IsCluster)
						continue;
					else
						rt_bucket.push_back(rt_temp.rt);

					// Evaluate, the proportional sum of the overlapping planes
					Evaluate(Planes_src, Planes_tgt, rt_temp);
				}
			}

			// The highest scored transformation matrix is the final result
			std::sort(RT_vector.begin(), RT_vector.end());
			final_rt = RT_vector[0].rt;

			return true;
		}

		void Evaluate(
			const std::vector<GPSCO::PLANE>& planes_src,
			const std::vector<GPSCO::PLANE>& planes_tgt,
			RT_Info& rt_temp)
		{
			float dist_plane = 0.05; // The distance at which two planes may overlap
			float dist_pointpair = 0.3; // The maximum distance between point pairs that overlap
			std::vector<int> pointIdxNKNSearch(1);
			std::vector<float> pointNKNSquaredDistance(1);
			for (const auto& plane_src : planes_src)
			{
				for (const auto& plane_tgt : planes_tgt)
				{
					auto normal_trans = rt_temp.rt.block<3, 3>(0, 0) * plane_src.normal;
					auto centroid_trans = rt_temp.rt.block<3, 3>(0, 0) * plane_src.centroid + rt_temp.rt.block<3, 1>(0, 3);
					// Determine if there is an overlap
					auto angle = GPSCO::angleBetweenVectors(normal_trans, plane_tgt.normal);
					if (angle < 2.0)
					{
						auto dist1_temp =
							GPSCO::distancePointToPlane(centroid_trans, plane_tgt.normal, plane_tgt.coefficients[3]);
						auto dist2_temp =
							GPSCO::distancePointToPlane(plane_tgt.centroid, normal_trans, -normal_trans.dot(centroid_trans));
						if ((dist1_temp + dist2_temp) / 2.0 < dist_plane)
						{
							GPSCO::cloudptr transformed_cloud(new GPSCO::cloud);
							pcl::transformPointCloud(*plane_src.patches, *transformed_cloud, rt_temp.rt);
							pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
							kdtree.setInputCloud(transformed_cloud);

							// Minimum distance between two plane point clouds
							plane_tgt.kdtree.nearestKSearch(transformed_cloud->points[0], 1, pointIdxNKNSearch, pointNKNSquaredDistance);
							kdtree.nearestKSearch(plane_tgt.patches->points[pointIdxNKNSearch[0]], 1, pointIdxNKNSearch, pointNKNSquaredDistance);
							if (pointNKNSquaredDistance[0] < dist_pointpair)
							{
								float overlap_ptnum = 0;
								for (const auto& pt : transformed_cloud->points)
								{
									plane_tgt.kdtree.nearestKSearch(pt, 1, pointIdxNKNSearch, pointNKNSquaredDistance);
									if (pointNKNSquaredDistance[0] < dist_pointpair)
									{
										overlap_ptnum++;
									}
								}
								rt_temp.confidence += (2 * overlap_ptnum) / (plane_src.patches->points.size() + plane_tgt.patches->points.size());
							}
						}
					}
				}
			}
		}
	}
}