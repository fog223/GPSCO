// system
#include <unordered_map>
#include <set>
#include <random>
// PCL
#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl/io/ply_io.h>
#include <pcl/common/common.h>
#include <pcl/common/distances.h>
#include <pcl/common/intersections.h>
#include <pcl/registration/transformation_estimation_svd.h>
#include <pcl/octree/octree.h>
// Eigen
#include <Eigen/Dense>
// spdlog
#include <spdlog/spdlog.h>
// local
#include "common.h"
#include "plane_extraction.h"
#include "registration.h"
#include "costfunction.h"

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
		area = 0.0f;
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
		area_score = 0.0f;
	}

	Group_Three_Pair::~Group_Three_Pair()
	{
	}

	RT_Info::RT_Info()
	{
		Intersection_points_src.reset(new GPSCO::cloud);
		Intersection_points_tgt.reset(new GPSCO::cloud);
		ini_rt.setIdentity();
		rt.setIdentity();
		confidence = 0.0f;
	}

	RT_Info::~RT_Info()
	{
	}

	bool Registration::Regis()
	{
		auto time_start = clock();

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
		auto time_begin = clock();
		if (!GPSCO::PLANE_Extraction::PLANE_Tetect_RegionGrow(cloud_src, options.min_support_points, options.max_plane_num,
			options.SmoothnessThreshold, options.CurvatureThreshold, options.SegSize, Planes_src))
		{
			spdlog::error("Source point cloud plane extraction failed.");
			return false;
		}
		if (!GPSCO::PLANE_Extraction::PLANE_Tetect_RegionGrow(cloud_tgt, options.min_support_points, options.max_plane_num,
			options.SmoothnessThreshold, options.CurvatureThreshold, options.SegSize, Planes_tgt))
		{
			spdlog::error("Target point cloud plane extraction failed.");
			return false;
		}
		auto time_end = clock();
		time_plane_extra = (double)(time_end - time_begin) / CLOCKS_PER_SEC;

		// large to small
		std::sort(Planes_src.begin(), Planes_src.end());
		std::sort(Planes_tgt.begin(), Planes_tgt.end());

		// Planar clustering
		time_begin = clock();
		if (!Plane_Cluster(Planes_src, PlaneGroups_src))
		{
			spdlog::error("Source point cloud plane clustering failed.");
			return false;
		}
		if (!Plane_Cluster(Planes_tgt, PlaneGroups_tgt))
		{
			spdlog::error("Target point cloud plane clustering failed.");
			return false;
		}
		time_end = clock();
		time_plane_cluster = (double)(time_end - time_begin) / CLOCKS_PER_SEC;

		// Compute the mean normal vector of a planar group
		Compute_AvgPlaneGroupNorm(Planes_src, PlaneGroups_src, AvgNorm_Group_src);
		Compute_AvgPlaneGroupNorm(Planes_tgt, PlaneGroups_tgt, AvgNorm_Group_tgt);

		// Export Plane, Plane Groups or Plane Clusters
//		GPSCO::PLANE_Extraction::Export_plane(Planes_src, "D:\\Code\\CLion\\GPSCO\\results\\src_");
//		GPSCO::PLANE_Extraction::Export_plane(Planes_tgt, "D:\\Code\\CLion\\GPSCO\\results\\tgt_");
//		Export_Cluster("D:\\Code\\CLion\\GPSCO\\results");
//		Export_Groups("D:\\Code\\CLion\\GPSCO\\results");

		// Initial match based moving alignment
		time_begin = clock();
		if (!Compute_Group_Table())
		{
			spdlog::error("Initial match failed.");
			return false;
		}
		else
			spdlog::info("Initial match completed.");
		time_end = clock();
		time_Match = (double)(time_end - time_begin) / CLOCKS_PER_SEC;

		// Get candidate
		if (Get_Group_Three(Planes_src, PlaneGroups_src, group_three_vector_src)
			&& Get_Group_Three(Planes_tgt, PlaneGroups_tgt, group_three_vector_tgt))
		{
			spdlog::info("Get candidate item completed.");
			// GPSCO determines the transformation matrix
			time_begin = clock();
			if (Get_transformation_matrix())
			{
				IsSuccessful = true;
				spdlog::info("Transformation matrix obtained successfully.");

				time_end = clock();
				time_Verify = (double)(time_end - time_begin) / CLOCKS_PER_SEC;

				auto time_finsh = clock();
				time_Total = (double)(time_finsh - time_start) / CLOCKS_PER_SEC;

				IsSuccess = true;
				return true;
			}
		}

		auto time_finsh = clock();
		time_Total = (double)(time_finsh - time_start) / CLOCKS_PER_SEC;

		spdlog::error("Transformation Matrix obtained failed.");
		return false;
	}

	bool Registration::Plane_Cluster(const std::vector<GPSCO::PLANE>& Planes,
		std::vector<std::vector<std::vector<int>>>& PlaneGroups)
	{
		PlaneGroups.clear();

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

					if ((plane_angle < options.parallel_thresh) || (std::fabs(180 - plane_angle) < options.parallel_thresh))
					{
						is_parallel = true;
						PlaneGroups_parallel[j].push_back(i);
					}
				}

				// too much
				if (PlaneGroups_parallel.size() > 9)
					continue;

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
						if (!(std::fabs(90 - angle1) < options.coplanar_thresh && std::fabs(90 - angle2) < options.coplanar_thresh))
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

	bool Registration::Compute_AvgPlaneGroupNorm(
		const std::vector<GPSCO::PLANE>& Planes,
		const std::vector<std::vector<std::vector<int>>>& PlaneGroups,
		std::vector<Eigen::Vector3f>& AvgNorm_Group)
	{
		for (const auto& planegroup : PlaneGroups)
		{
			int num_group = 0;
			Eigen::Vector3f avgnorm(0, 0, 0);
			for (const auto& cluster_idx : planegroup)
			{
				for (const auto& plane_idx : cluster_idx)
				{
					num_group++;
					if (avgnorm.dot(Planes[plane_idx].normal) > 0)
						avgnorm += Planes[plane_idx].normal;
					else
						avgnorm -= Planes[plane_idx].normal;
				}
			}
			avgnorm /= num_group;
			avgnorm.normalize();

			AvgNorm_Group.push_back(avgnorm);
		}

		return true;
	}

	bool Registration::Compute_Group_Table()
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
					if (Moving_alignment(i, j, match_vec))
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

//				spdlog::info("The relationship table between plane groups is successfully obtained.");
			return true;
		}
	}

	bool Registration::Moving_alignment(int group1_index, int group2_index, std::vector<MatchGroup>& match_vec)
	{
		auto planegroup_src = PlaneGroups_src[group1_index];
		auto planegroup_tgt = PlaneGroups_tgt[group2_index];

		if (planegroup_src.empty() || planegroup_tgt.empty())
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

		Eigen::Vector3f Normal_avg_src = AvgNorm_Group_src[group1_index];
		Eigen::Vector3f Normal_avg_tgt = AvgNorm_Group_tgt[group2_index];

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
						= Dist_TwoClsuter(Planes_src, Normal_avg_src, planegroup_src[i], planegroup_src[j]);
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
						= Dist_TwoClsuter(Planes_tgt, Normal_avg_tgt, planegroup_tgt[i], planegroup_tgt[j]);
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
						if (std::abs(disarrays_src[i][m] - disarrays_tgt[j][n]) <= options.e_pl2pldist)
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
						if (std::abs(disarrays_src[i][m] - disarrays_tgt[j][n]) <= options.e_pl2pldist)
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
							if (d_src[vec_[vec_idx][pair_idx].first][vec_[vec_idx][pair_idx + 1].first] < options.min_pl2pldist ||
								d_tgt[vec_[vec_idx][pair_idx].second][vec_[vec_idx][pair_idx + 1].second] < options.min_pl2pldist)
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

					if (vec_.size() > 10000)
						break;
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
			matchpair.group1_index = group1_index;//int
			matchpair.group2_index = group2_index;
			matchpair.planepairs = combination;//std::vector<std::pair<int, int>>
			matchpair.match_num = combination.size();
			for (const auto& cluster : combination)
			{
				float area_src = 0.0f;
				float area_tgt = 0.0f;
				for (const auto planeidx_src : planegroup_src[cluster.first])
				{
					area_src += Planes_src[planeidx_src].points->size();
				}
				for (const auto planeidx_tgt : planegroup_tgt[cluster.second])
				{
					area_tgt += Planes_tgt[planeidx_tgt].points->size();
				}
				matchpair.area += std::min(area_src, area_tgt);
			}
			match_vec.push_back(matchpair);
		}

		std::sort(match_vec.begin(), match_vec.end());

		return true;
	}

	float Registration::Dist_TwoClsuter(
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

	bool Registration::Get_Group_Three(
		const std::vector<GPSCO::PLANE>& Planes,
		const std::vector<std::vector<std::vector<int>>>& PlaneGroups,
		std::vector<GPSCO::Group_Three>& group_three_vector)
	{
		bool IsTimeOut = false;
		bool Cancel = false;

		std::thread TimeWait([&]()
		{
		  for (int i = 0; i < options.wait_time; ++i)
		  {
			  if (Cancel)
			  {
				  return;
			  }
			  std::this_thread::sleep_for(std::chrono::milliseconds(500));
		  }
		  IsTimeOut = true;
		});
		TimeWait.detach();

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
				if (IsTimeOut)
				{
					spdlog::warn("Timeout_1!");
					return false;
				}

				if (i < j)
				{
					float group_angle = pcl::getAngle3D(Planes[PlaneGroups[i][0][0]].normal,
						Planes[PlaneGroups[j][0][0]].normal, true);
					if (options.angle_min < group_angle && group_angle < options.angle_max)
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
				if (IsTimeOut)
				{
					spdlog::warn("Timeout_2!");
					return false;
				}

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

						Eigen::Vector3f cross_normal = (Planes[plane_index1].normal).cross(Planes[plane_index2].normal);
						float angle = pcl::getAngle3D(cross_normal,
							Planes[PlaneGroups[i][0][0]].normal, true);

						if (angle < 90 - options.angle_min || angle > 90 + options.angle_min)
						{
							float angle1 = pcl::getAngle3D(Planes[plane_index1].normal,
								Planes[PlaneGroups[i][0][0]].normal, true);
							float angle2 = pcl::getAngle3D(Planes[plane_index2].normal,
								Planes[PlaneGroups[i][0][0]].normal, true);
							if ((options.angle_min < angle1 && angle1 < options.angle_max) &&
								(options.angle_min < angle2 && angle2 < options.angle_max))
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
		}

		if (group_three_vector.empty())
		{
			spdlog::error("The number of three plane groups that meet the threshold is 0.");
			return false;
		}

		Cancel = true;

		return true;
	}

	bool Registration::Get_transformation_matrix()
	{
		// Three pairs of non-parallel plane groups
		std::vector<GPSCO::Group_Three_Pair> group_three_pair_vec;
		for (const auto& group_three_src : group_three_vector_src)
		{
			for (const auto& group_three_tgt : group_three_vector_tgt)
			{
				// Six possible matches
				// type 1: 1-1,2-2,3-3
				if (std::fabs(group_three_src.group_angle12 - group_three_tgt.group_angle12) < options.e_angle &&
					std::fabs(group_three_src.group_angle13 - group_three_tgt.group_angle13) < options.e_angle &&
					std::fabs(group_three_src.group_angle23 - group_three_tgt.group_angle23) < options.e_angle)
				{
					GPSCO::Group_Three_Pair group_three_pair_temp;
					group_three_pair_temp.group_src_index1 = group_three_src.group1_index;
					group_three_pair_temp.group_src_index2 = group_three_src.group2_index;
					group_three_pair_temp.group_src_index3 = group_three_src.group3_index;
					group_three_pair_temp.group_tgt_index1 = group_three_tgt.group1_index;
					group_three_pair_temp.group_tgt_index2 = group_three_tgt.group2_index;
					group_three_pair_temp.group_tgt_index3 = group_three_tgt.group3_index;
					auto matchscore1 = group_table[group_three_pair_temp.group_src_index1][group_three_pair_temp.group_tgt_index1][0].match_num;
					auto matchscore2 = group_table[group_three_pair_temp.group_src_index2][group_three_pair_temp.group_tgt_index2][0].match_num;
					auto matchscore3 = group_table[group_three_pair_temp.group_src_index3][group_three_pair_temp.group_tgt_index3][0].match_num;
					group_three_pair_temp.match_score = matchscore1 + matchscore2 + matchscore3;
					group_three_pair_vec.push_back(group_three_pair_temp);
				}
				// type 2: 1-1,2-3,3-2
				if (std::fabs(group_three_src.group_angle12 - group_three_tgt.group_angle13) < options.e_angle &&
					std::fabs(group_three_src.group_angle13 - group_three_tgt.group_angle12) < options.e_angle &&
					std::fabs(group_three_src.group_angle23 - group_three_tgt.group_angle23) < options.e_angle)
				{
					GPSCO::Group_Three_Pair group_three_pair_temp;
					group_three_pair_temp.group_src_index1 = group_three_src.group1_index;
					group_three_pair_temp.group_src_index2 = group_three_src.group2_index;
					group_three_pair_temp.group_src_index3 = group_three_src.group3_index;
					group_three_pair_temp.group_tgt_index1 = group_three_tgt.group1_index;
					group_three_pair_temp.group_tgt_index2 = group_three_tgt.group3_index;
					group_three_pair_temp.group_tgt_index3 = group_three_tgt.group2_index;
					auto matchscore1 = group_table[group_three_pair_temp.group_src_index1][group_three_pair_temp.group_tgt_index1][0].match_num;
					auto matchscore2 = group_table[group_three_pair_temp.group_src_index2][group_three_pair_temp.group_tgt_index2][0].match_num;
					auto matchscore3 = group_table[group_three_pair_temp.group_src_index3][group_three_pair_temp.group_tgt_index3][0].match_num;
					group_three_pair_temp.match_score = matchscore1 + matchscore2 + matchscore3;
					group_three_pair_vec.push_back(group_three_pair_temp);
				}
				// type 3: 1-2,2-1,3-3
				if (std::fabs(group_three_src.group_angle12 - group_three_tgt.group_angle12) < options.e_angle &&
					std::fabs(group_three_src.group_angle13 - group_three_tgt.group_angle23) < options.e_angle &&
					std::fabs(group_three_src.group_angle23 - group_three_tgt.group_angle13) < options.e_angle)
				{
					GPSCO::Group_Three_Pair group_three_pair_temp;
					group_three_pair_temp.group_src_index1 = group_three_src.group1_index;
					group_three_pair_temp.group_src_index2 = group_three_src.group2_index;
					group_three_pair_temp.group_src_index3 = group_three_src.group3_index;
					group_three_pair_temp.group_tgt_index1 = group_three_tgt.group2_index;
					group_three_pair_temp.group_tgt_index2 = group_three_tgt.group1_index;
					group_three_pair_temp.group_tgt_index3 = group_three_tgt.group3_index;
					auto matchscore1 = group_table[group_three_pair_temp.group_src_index1][group_three_pair_temp.group_tgt_index1][0].match_num;
					auto matchscore2 = group_table[group_three_pair_temp.group_src_index2][group_three_pair_temp.group_tgt_index2][0].match_num;
					auto matchscore3 = group_table[group_three_pair_temp.group_src_index3][group_three_pair_temp.group_tgt_index3][0].match_num;
					group_three_pair_temp.match_score = matchscore1 + matchscore2 + matchscore3;
					group_three_pair_vec.push_back(group_three_pair_temp);
				}
				// type 4: 1-2,2-3,3-1
				if (std::fabs(group_three_src.group_angle12 - group_three_tgt.group_angle23) < options.e_angle &&
					std::fabs(group_three_src.group_angle13 - group_three_tgt.group_angle12) < options.e_angle &&
					std::fabs(group_three_src.group_angle23 - group_three_tgt.group_angle13) < options.e_angle)
				{
					GPSCO::Group_Three_Pair group_three_pair_temp;
					group_three_pair_temp.group_src_index1 = group_three_src.group1_index;
					group_three_pair_temp.group_src_index2 = group_three_src.group2_index;
					group_three_pair_temp.group_src_index3 = group_three_src.group3_index;
					group_three_pair_temp.group_tgt_index1 = group_three_tgt.group2_index;
					group_three_pair_temp.group_tgt_index2 = group_three_tgt.group3_index;
					group_three_pair_temp.group_tgt_index3 = group_three_tgt.group1_index;
					auto matchscore1 = group_table[group_three_pair_temp.group_src_index1][group_three_pair_temp.group_tgt_index1][0].match_num;
					auto matchscore2 = group_table[group_three_pair_temp.group_src_index2][group_three_pair_temp.group_tgt_index2][0].match_num;
					auto matchscore3 = group_table[group_three_pair_temp.group_src_index3][group_three_pair_temp.group_tgt_index3][0].match_num;
					group_three_pair_temp.match_score = matchscore1 + matchscore2 + matchscore3;
					group_three_pair_vec.push_back(group_three_pair_temp);
				}
				// type 5: 1-3,2-1,3-2
				if (std::fabs(group_three_src.group_angle12 - group_three_tgt.group_angle13) < options.e_angle &&
					std::fabs(group_three_src.group_angle13 - group_three_tgt.group_angle23) < options.e_angle &&
					std::fabs(group_three_src.group_angle23 - group_three_tgt.group_angle12) < options.e_angle)
				{
					GPSCO::Group_Three_Pair group_three_pair_temp;
					group_three_pair_temp.group_src_index1 = group_three_src.group1_index;
					group_three_pair_temp.group_src_index2 = group_three_src.group2_index;
					group_three_pair_temp.group_src_index3 = group_three_src.group3_index;
					group_three_pair_temp.group_tgt_index1 = group_three_tgt.group3_index;
					group_three_pair_temp.group_tgt_index2 = group_three_tgt.group1_index;
					group_three_pair_temp.group_tgt_index3 = group_three_tgt.group2_index;
					auto matchscore1 = group_table[group_three_pair_temp.group_src_index1][group_three_pair_temp.group_tgt_index1][0].match_num;
					auto matchscore2 = group_table[group_three_pair_temp.group_src_index2][group_three_pair_temp.group_tgt_index2][0].match_num;
					auto matchscore3 = group_table[group_three_pair_temp.group_src_index3][group_three_pair_temp.group_tgt_index3][0].match_num;
					group_three_pair_temp.match_score = matchscore1 + matchscore2 + matchscore3;
					group_three_pair_vec.push_back(group_three_pair_temp);
				}
				// type 6: 1-3,2-2,3-1
				if (std::fabs(group_three_src.group_angle12 - group_three_tgt.group_angle23) < options.e_angle &&
					std::fabs(group_three_src.group_angle13 - group_three_tgt.group_angle13) < options.e_angle &&
					std::fabs(group_three_src.group_angle23 - group_three_tgt.group_angle12) < options.e_angle)
				{
					GPSCO::Group_Three_Pair group_three_pair_temp;
					group_three_pair_temp.group_src_index1 = group_three_src.group1_index;
					group_three_pair_temp.group_src_index2 = group_three_src.group2_index;
					group_three_pair_temp.group_src_index3 = group_three_src.group3_index;
					group_three_pair_temp.group_tgt_index1 = group_three_tgt.group3_index;
					group_three_pair_temp.group_tgt_index2 = group_three_tgt.group2_index;
					group_three_pair_temp.group_tgt_index3 = group_three_tgt.group1_index;
					auto matchscore1 = group_table[group_three_pair_temp.group_src_index1][group_three_pair_temp.group_tgt_index1][0].match_num;
					auto matchscore2 = group_table[group_three_pair_temp.group_src_index2][group_three_pair_temp.group_tgt_index2][0].match_num;
					auto matchscore3 = group_table[group_three_pair_temp.group_src_index3][group_three_pair_temp.group_tgt_index3][0].match_num;
					group_three_pair_temp.match_score = matchscore1 + matchscore2 + matchscore3;
					group_three_pair_vec.push_back(group_three_pair_temp);
				}
			}
		}

		if (group_three_pair_vec.empty())
			return false;

		// Sorted by match score, highest to lowest
		std::sort(group_three_pair_vec.begin(), group_three_pair_vec.end());

		// Global Planar Structural Constraint Optimal(highest score)
		options.start = clock();
		int min_match_num = 3;
		if (!options.IsValidate)
		{
			min_match_num = group_three_pair_vec[0].match_score;
		}

		for (int i = 0; i < group_three_pair_vec.size(); i++)
		{
			if (group_three_pair_vec[i].match_score < min_match_num)
				break;
			Find_optimal_RT(group_three_pair_vec[i], min_match_num, RT_vector);

			clock_t end = clock();
			float seconds = static_cast<float>(end - options.start) / CLOCKS_PER_SEC;
			if (seconds > options.wait_time)
			{
				spdlog::warn("Timeout!");
				return false;
			}
		}

		if (RT_vector.empty())
		{
			spdlog::error("RT_vector is empty.");
			return false;
		}

		// Verification
		std::vector<Eigen::Matrix4f> rt_bucket;
		for (auto& rt_temp : RT_vector)
		{
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
			Evaluate(rt_temp);
		}

		// The highest scored transformation matrix is the final result
		std::sort(RT_vector.begin(), RT_vector.end());
		Final_RT = RT_vector[0].rt;

		return true;
	}

	bool Registration::Find_optimal_RT(
		const GPSCO::Group_Three_Pair& group_three_pair,
		int& min_match_num,
		std::vector<GPSCO::RT_Info>& RT_vector)
	{
		int idx_src_g1 = group_three_pair.group_src_index1;
		int idx_tgt_g1 = group_three_pair.group_tgt_index1;
		int idx_src_g2 = group_three_pair.group_src_index2;
		int idx_tgt_g2 = group_three_pair.group_tgt_index2;
		int idx_src_g3 = group_three_pair.group_src_index3;
		int idx_tgt_g3 = group_three_pair.group_tgt_index3;
		Eigen::Vector3f avgnorm_src_g1 = AvgNorm_Group_src[idx_src_g1];
		Eigen::Vector3f avgnorm_src_g2 = AvgNorm_Group_src[idx_src_g2];
		Eigen::Vector3f avgnorm_src_g3 = AvgNorm_Group_src[idx_src_g3];
		Eigen::Vector3f avgnorm_tgt_g1 = AvgNorm_Group_tgt[idx_tgt_g1];
		Eigen::Vector3f avgnorm_tgt_g2 = AvgNorm_Group_tgt[idx_tgt_g2];
		Eigen::Vector3f avgnorm_tgt_g3 = AvgNorm_Group_tgt[idx_tgt_g3];

		int max_gidx_src = idx_src_g1, mid_gidx_src = idx_src_g2, min_gidx_src = idx_src_g3;
		int max_gidx_tgt = idx_tgt_g1, mid_gidx_tgt = idx_tgt_g2, min_gidx_tgt = idx_tgt_g3;
		int num1 = group_table[idx_src_g1][idx_tgt_g1].size();
		int num2 = group_table[idx_src_g2][idx_tgt_g2].size();
		int num3 = group_table[idx_src_g3][idx_tgt_g3].size();
		if (num1 > num2)
		{
			// 3 -> 1 -> 2
			if (num3 > num1)
			{
				max_gidx_src = idx_src_g3;
				max_gidx_tgt = idx_tgt_g3;
				mid_gidx_src = idx_src_g1;
				mid_gidx_tgt = idx_tgt_g1;
				min_gidx_src = idx_src_g2;
				min_gidx_tgt = idx_tgt_g2;
			}
				// 1 -> 3 -> 2
			else if (num3 > num2)
			{
				max_gidx_src = idx_src_g1;
				max_gidx_tgt = idx_tgt_g1;
				mid_gidx_src = idx_src_g3;
				mid_gidx_tgt = idx_tgt_g3;
				min_gidx_src = idx_src_g2;
				min_gidx_tgt = idx_tgt_g2;
			}
		}
		else
		{
			// 3 -> 2 -> 1
			if (num3 > num2)
			{
				max_gidx_src = idx_src_g3;
				max_gidx_tgt = idx_tgt_g3;
				mid_gidx_src = idx_src_g2;
				mid_gidx_tgt = idx_tgt_g2;
				min_gidx_src = idx_src_g1;
				min_gidx_tgt = idx_tgt_g1;
			}
				// 2 -> 3 -> 1
			else if (num3 > num1)
			{
				max_gidx_src = idx_src_g2;
				max_gidx_tgt = idx_tgt_g2;
				mid_gidx_src = idx_src_g3;
				mid_gidx_tgt = idx_tgt_g3;
				min_gidx_src = idx_src_g1;
				min_gidx_tgt = idx_tgt_g1;
			}
				// 2 -> 1 -> 3
			else
			{
				max_gidx_src = idx_src_g2;
				max_gidx_tgt = idx_tgt_g2;
				mid_gidx_src = idx_src_g1;
				mid_gidx_tgt = idx_tgt_g1;
				min_gidx_src = idx_src_g3;
				min_gidx_tgt = idx_tgt_g3;
			}
		}

		auto gg1_vec = group_table[max_gidx_src][max_gidx_tgt];
		auto gg2_vec = group_table[mid_gidx_src][mid_gidx_tgt];
		auto gg3_vec = group_table[min_gidx_src][min_gidx_tgt];

		for (int i = 0; i < gg1_vec.size(); i++)
		{
			if ((gg1_vec[i].match_num + gg2_vec[0].match_num + gg3_vec[0].match_num) >= min_match_num)
			{
				auto p11_idx = PlaneGroups_src[max_gidx_src][gg1_vec[i].planepairs[0].first][0];
				auto p12_idx = PlaneGroups_tgt[max_gidx_tgt][gg1_vec[i].planepairs[0].second][0];
				Eigen::Vector4f plane_a_src, plane_a_tgt;
				plane_a_src << Planes_src[p11_idx].coefficients[0],
					Planes_src[p11_idx].coefficients[1],
					Planes_src[p11_idx].coefficients[2],
					Planes_src[p11_idx].coefficients[3];
				plane_a_tgt << Planes_tgt[p12_idx].coefficients[0],
					Planes_tgt[p12_idx].coefficients[1],
					Planes_tgt[p12_idx].coefficients[2],
					Planes_tgt[p12_idx].coefficients[3];

				bool lock1 = false;
				if (gg1_vec[i].match_num > 1)
				{
					lock1 = true;
					auto startp_src_idx = PlaneGroups_src[max_gidx_src][gg1_vec[i].planepairs[0].first][0];
					auto endp_src_idx = PlaneGroups_src[max_gidx_src][gg1_vec[i].planepairs.back().first][0];
					auto startp_tgt_idx = PlaneGroups_tgt[max_gidx_tgt][gg1_vec[i].planepairs[0].second][0];
					auto endp_tgt_idx = PlaneGroups_tgt[max_gidx_tgt][gg1_vec[i].planepairs.back().second][0];
					if ((Planes_src[endp_src_idx].centroid - Planes_src[startp_src_idx].centroid)
						.dot(avgnorm_src_g1) < 0)
						avgnorm_src_g1 = -avgnorm_src_g1;
					if ((Planes_tgt[endp_tgt_idx].centroid - Planes_tgt[startp_tgt_idx].centroid)
						.dot(avgnorm_tgt_g1) < 0)
						avgnorm_tgt_g1 = -avgnorm_tgt_g1;
				}

				for (int j = 0; j < gg2_vec.size(); j++)
				{
					if ((gg1_vec[i].match_num + gg2_vec[j].match_num + gg3_vec[0].match_num) >= min_match_num)
					{
						auto p21_idx = PlaneGroups_src[mid_gidx_src][gg2_vec[j].planepairs[0].first][0];
						auto p22_idx = PlaneGroups_tgt[mid_gidx_tgt][gg2_vec[j].planepairs[0].second][0];
						Eigen::Vector4f plane_b_src, plane_b_tgt;
						plane_b_src << Planes_src[p21_idx].coefficients[0],
							Planes_src[p21_idx].coefficients[1],
							Planes_src[p21_idx].coefficients[2],
							Planes_src[p21_idx].coefficients[3];
						plane_b_tgt << Planes_tgt[p22_idx].coefficients[0],
							Planes_tgt[p22_idx].coefficients[1],
							Planes_tgt[p22_idx].coefficients[2],
							Planes_tgt[p22_idx].coefficients[3];

						bool lock2 = false;
						if (gg2_vec[j].match_num > 1)
						{
							lock2 = true;
							auto startp_src_idx = PlaneGroups_src[mid_gidx_src][gg2_vec[j].planepairs[0].first][0];
							auto endp_src_idx = PlaneGroups_src[mid_gidx_src][gg2_vec[j].planepairs.back().first][0];
							auto startp_tgt_idx = PlaneGroups_tgt[mid_gidx_tgt][gg2_vec[j].planepairs[0].second][0];
							auto endp_tgt_idx = PlaneGroups_tgt[mid_gidx_tgt][gg2_vec[j].planepairs.back().second][0];
							if ((Planes_src[endp_src_idx].centroid - Planes_src[startp_src_idx].centroid)
								.dot(avgnorm_src_g2) < 0)
								avgnorm_src_g2 = -avgnorm_src_g2;
							if ((Planes_tgt[endp_tgt_idx].centroid - Planes_tgt[startp_tgt_idx].centroid)
								.dot(avgnorm_tgt_g2) < 0)
								avgnorm_tgt_g2 = -avgnorm_tgt_g2;
						}

						for (int k = 0; k < gg3_vec.size(); k++)
						{
							clock_t end = clock();
							float seconds = static_cast<float>(end - options.start) / CLOCKS_PER_SEC;
							if (seconds > options.wait_time)
							{
								return false;
							}

							if ((gg1_vec[i].match_num + gg2_vec[j].match_num + gg3_vec[k].match_num) >= min_match_num)
							{
								auto p31_idx = PlaneGroups_src[min_gidx_src][gg3_vec[k].planepairs[0].first][0];
								auto p32_idx = PlaneGroups_tgt[min_gidx_tgt][gg3_vec[k].planepairs[0].second][0];
								Eigen::Vector4f plane_c_src, plane_c_tgt;
								plane_c_src << Planes_src[p31_idx].coefficients[0],
									Planes_src[p31_idx].coefficients[1],
									Planes_src[p31_idx].coefficients[2],
									Planes_src[p31_idx].coefficients[3];
								plane_c_tgt << Planes_tgt[p32_idx].coefficients[0],
									Planes_tgt[p32_idx].coefficients[1],
									Planes_tgt[p32_idx].coefficients[2],
									Planes_tgt[p32_idx].coefficients[3];

								Eigen::Vector3f inter_pt_src(0.0, 0.0, 0.0), inter_pt_tgt(0.0, 0.0, 0.0);
								if (!(pcl::threePlanesIntersection(plane_a_src, plane_b_src, plane_c_src, inter_pt_src)
									&& pcl::threePlanesIntersection(plane_a_tgt, plane_b_tgt, plane_c_tgt, inter_pt_tgt)))
									break;

								bool lock3 = false;
								if (gg3_vec[k].match_num > 1)
								{
									lock3 = true;
									auto startp_src_idx = PlaneGroups_src[min_gidx_src][gg3_vec[k].planepairs[0].first][0];
									auto endp_src_idx = PlaneGroups_src[min_gidx_src][gg3_vec[k].planepairs.back().first][0];
									auto startp_tgt_idx = PlaneGroups_tgt[min_gidx_tgt][gg3_vec[k].planepairs[0].second][0];
									auto endp_tgt_idx = PlaneGroups_tgt[min_gidx_tgt][gg3_vec[k].planepairs.back().second][0];
									if ((Planes_src[endp_src_idx].centroid - Planes_src[startp_src_idx].centroid)
										.dot(avgnorm_src_g3) < 0)
										avgnorm_src_g3 = -avgnorm_src_g3;
									if ((Planes_tgt[endp_tgt_idx].centroid - Planes_tgt[startp_tgt_idx].centroid)
										.dot(avgnorm_tgt_g3) < 0)
										avgnorm_tgt_g3 = -avgnorm_tgt_g3;
								}

								RT_Info rt_info;
								rt_info.plane_pairs1 = gg1_vec[i];
								rt_info.plane_pairs2 = gg2_vec[j];
								rt_info.plane_pairs3 = gg3_vec[k];

								Get_Intersections(rt_info);
								if (rt_info.Intersection_points_src->points.empty())
									continue;

								if (!options.IsValidate && rt_info.Intersection_points_src->points.size() >= 3)
								{
									pcl::registration::TransformationEstimationSVD<pcl::PointXYZ, pcl::PointXYZ> SVD;
									SVD.estimateRigidTransformation(
										*rt_info.Intersection_points_src,
										*rt_info.Intersection_points_tgt,
										rt_info.rt);

									RT_vector.push_back(rt_info);
								}

								Eigen::Matrix3f rotMatrix1, rotMatrix2, rotMatrix3, rotMatrix4;
								bool possible1 = false, possible2 = false, possible3 = false, possible4 = false;
								if (lock1 && lock2)
								{
									GPSCO::Compute_rotMatrix(avgnorm_src_g1, avgnorm_src_g2, avgnorm_tgt_g1, avgnorm_tgt_g2, rotMatrix1);
									auto angle11 = pcl::getAngle3D(rotMatrix1 * avgnorm_src_g1, avgnorm_tgt_g1, true);
									auto angle12 = pcl::getAngle3D(rotMatrix1 * avgnorm_src_g2, avgnorm_tgt_g2, true);
									auto angle13 = pcl::getAngle3D(rotMatrix1 * avgnorm_src_g3, avgnorm_tgt_g3, true);
									if (lock3)
									{
										if (angle11 < 10 && angle12 < 10 && angle13 < 10)
											possible1 = true;
									}
									else
									{
										if (angle11 < 10 && angle12 < 10 && (angle13 < 10 || angle13 > 170))
											possible1 = true;
									}
								}
								else if (lock1 && lock3)
								{
									GPSCO::Compute_rotMatrix(avgnorm_src_g1, avgnorm_src_g3, avgnorm_tgt_g1, avgnorm_tgt_g3, rotMatrix1);
									auto angle11 = pcl::getAngle3D(rotMatrix1 * avgnorm_src_g1, avgnorm_tgt_g1, true);
									auto angle12 = pcl::getAngle3D(rotMatrix1 * avgnorm_src_g2, avgnorm_tgt_g2, true);
									auto angle13 = pcl::getAngle3D(rotMatrix1 * avgnorm_src_g3, avgnorm_tgt_g3, true);
									if (angle11 < 10 && (angle12 < 10 || angle12 > 170) && angle13 < 10)
									{
										possible1 = true;
									}
								}
								else if (lock2 && lock3)
								{
									GPSCO::Compute_rotMatrix(avgnorm_src_g2, avgnorm_src_g3, avgnorm_tgt_g2, avgnorm_tgt_g3, rotMatrix1);
									auto angle11 = pcl::getAngle3D(rotMatrix1 * avgnorm_src_g1, avgnorm_tgt_g1, true);
									auto angle12 = pcl::getAngle3D(rotMatrix1 * avgnorm_src_g2, avgnorm_tgt_g2, true);
									auto angle13 = pcl::getAngle3D(rotMatrix1 * avgnorm_src_g3, avgnorm_tgt_g3, true);
									if ((angle11 < 10 || angle11 > 170) && angle12 < 10 && angle13 < 10)
									{
										possible1 = true;
									}
								}
								else if (lock1)
								{
									GPSCO::Compute_rotMatrix(avgnorm_src_g1, avgnorm_src_g2, avgnorm_tgt_g1, avgnorm_tgt_g2, rotMatrix1);
									GPSCO::Compute_rotMatrix(avgnorm_src_g1, avgnorm_src_g2, avgnorm_tgt_g1, -avgnorm_tgt_g2, rotMatrix2);
									auto angle11 = pcl::getAngle3D(rotMatrix1 * avgnorm_src_g1, avgnorm_tgt_g1, true);
									auto angle12 = pcl::getAngle3D(rotMatrix1 * avgnorm_src_g2, avgnorm_tgt_g2, true);
									auto angle13 = pcl::getAngle3D(rotMatrix1 * avgnorm_src_g3, avgnorm_tgt_g3, true);
									auto angle21 = pcl::getAngle3D(rotMatrix2 * avgnorm_src_g1, avgnorm_tgt_g1, true);
									auto angle22 = pcl::getAngle3D(rotMatrix2 * avgnorm_src_g2, avgnorm_tgt_g2, true);
									auto angle23 = pcl::getAngle3D(rotMatrix2 * avgnorm_src_g3, avgnorm_tgt_g3, true);
									if (angle11 < 10 && (angle12 < 10 || angle12 > 170) && (angle13 < 10 || angle13 > 170))
									{
										possible1 = true;
									}
									if (angle21 < 10 && (angle22 < 10 || angle22 > 170) && (angle23 < 10 || angle23 > 170))
									{
										possible2 = true;
									}
								}
								else if (lock2)
								{
									GPSCO::Compute_rotMatrix(avgnorm_src_g2, avgnorm_src_g1, avgnorm_tgt_g2, avgnorm_tgt_g1, rotMatrix1);
									GPSCO::Compute_rotMatrix(avgnorm_src_g2, avgnorm_src_g1, avgnorm_tgt_g2, -avgnorm_tgt_g1, rotMatrix2);
									auto angle11 = pcl::getAngle3D(rotMatrix1 * avgnorm_src_g1, avgnorm_tgt_g1, true);
									auto angle12 = pcl::getAngle3D(rotMatrix1 * avgnorm_src_g2, avgnorm_tgt_g2, true);
									auto angle13 = pcl::getAngle3D(rotMatrix1 * avgnorm_src_g3, avgnorm_tgt_g3, true);
									auto angle21 = pcl::getAngle3D(rotMatrix2 * avgnorm_src_g1, avgnorm_tgt_g1, true);
									auto angle22 = pcl::getAngle3D(rotMatrix2 * avgnorm_src_g2, avgnorm_tgt_g2, true);
									auto angle23 = pcl::getAngle3D(rotMatrix2 * avgnorm_src_g3, avgnorm_tgt_g3, true);
									if ((angle11 < 10 || angle11 > 170) && angle12 < 10 && (angle13 < 10 || angle13 > 170))
									{
										possible1 = true;
									}
									if ((angle21 < 10 || angle21 > 170) && angle22 < 10 && (angle23 < 10 || angle23 > 170))
									{
										possible2 = true;
									}
								}
								else if (lock3)
								{
									GPSCO::Compute_rotMatrix(avgnorm_src_g3, avgnorm_src_g1, avgnorm_tgt_g3, avgnorm_tgt_g1, rotMatrix1);
									GPSCO::Compute_rotMatrix(avgnorm_src_g3, avgnorm_src_g1, avgnorm_tgt_g3, -avgnorm_tgt_g1, rotMatrix2);
									auto angle11 = pcl::getAngle3D(rotMatrix1 * avgnorm_src_g1, avgnorm_tgt_g1, true);
									auto angle12 = pcl::getAngle3D(rotMatrix1 * avgnorm_src_g2, avgnorm_tgt_g2, true);
									auto angle13 = pcl::getAngle3D(rotMatrix1 * avgnorm_src_g3, avgnorm_tgt_g3, true);
									auto angle21 = pcl::getAngle3D(rotMatrix2 * avgnorm_src_g1, avgnorm_tgt_g1, true);
									auto angle22 = pcl::getAngle3D(rotMatrix2 * avgnorm_src_g2, avgnorm_tgt_g2, true);
									auto angle23 = pcl::getAngle3D(rotMatrix2 * avgnorm_src_g3, avgnorm_tgt_g3, true);
									if ((angle11 < 10 || angle11 > 170) && (angle12 < 10 || angle12 > 170) && angle13 < 10)
									{
										possible1 = true;
									}
									if ((angle21 < 10 || angle21 > 170) && (angle22 < 10 || angle22 > 170) && angle23 < 10)
									{
										possible2 = true;
									}
								}
								else
								{
									GPSCO::Compute_rotMatrix(avgnorm_src_g1, avgnorm_src_g2, avgnorm_tgt_g1, avgnorm_tgt_g2, rotMatrix1);
									GPSCO::Compute_rotMatrix(avgnorm_src_g1, avgnorm_src_g2, avgnorm_tgt_g1, -avgnorm_tgt_g2, rotMatrix2);
									GPSCO::Compute_rotMatrix(avgnorm_src_g1, avgnorm_src_g2, -avgnorm_tgt_g1, avgnorm_tgt_g2, rotMatrix3);
									GPSCO::Compute_rotMatrix(avgnorm_src_g1, avgnorm_src_g2, -avgnorm_tgt_g1, -avgnorm_tgt_g2, rotMatrix4);
									auto angle11 = pcl::getAngle3D(rotMatrix1 * avgnorm_src_g1, avgnorm_tgt_g1, true);
									auto angle12 = pcl::getAngle3D(rotMatrix1 * avgnorm_src_g2, avgnorm_tgt_g2, true);
									auto angle13 = pcl::getAngle3D(rotMatrix1 * avgnorm_src_g3, avgnorm_tgt_g3, true);
									auto angle21 = pcl::getAngle3D(rotMatrix2 * avgnorm_src_g1, avgnorm_tgt_g1, true);
									auto angle22 = pcl::getAngle3D(rotMatrix2 * avgnorm_src_g2, avgnorm_tgt_g2, true);
									auto angle23 = pcl::getAngle3D(rotMatrix2 * avgnorm_src_g3, avgnorm_tgt_g3, true);
									auto angle31 = pcl::getAngle3D(rotMatrix3 * avgnorm_src_g1, avgnorm_tgt_g1, true);
									auto angle32 = pcl::getAngle3D(rotMatrix3 * avgnorm_src_g2, avgnorm_tgt_g2, true);
									auto angle33 = pcl::getAngle3D(rotMatrix3 * avgnorm_src_g3, avgnorm_tgt_g3, true);
									auto angle41 = pcl::getAngle3D(rotMatrix4 * avgnorm_src_g1, avgnorm_tgt_g1, true);
									auto angle42 = pcl::getAngle3D(rotMatrix4 * avgnorm_src_g2, avgnorm_tgt_g2, true);
									auto angle43 = pcl::getAngle3D(rotMatrix4 * avgnorm_src_g3, avgnorm_tgt_g3, true);
									if ((angle11 < 10 || angle11 > 170) && (angle12 < 10 || angle12 > 170) && (angle13 < 10 || angle13 > 170))
									{
										possible1 = true;
									}
									if ((angle21 < 10 || angle21 > 170) && (angle22 < 10 || angle22 > 170) && (angle23 < 10 || angle23 > 170))
									{
										possible2 = true;
									}
									if ((angle31 < 10 || angle31 > 170) && (angle32 < 10 || angle32 > 170) && (angle33 < 10 || angle33 > 170))
									{
										possible3 = true;
									}
									if ((angle41 < 10 || angle41 > 170) && (angle42 < 10 || angle42 > 170) && (angle43 < 10 || angle43 > 170))
									{
										possible4 = true;
									}
								}

								if (possible1)
								{
									rt_info.ini_rt(0, 0) = rotMatrix1(0, 0);
									rt_info.ini_rt(0, 1) = rotMatrix1(0, 1);
									rt_info.ini_rt(0, 2) = rotMatrix1(0, 2);
									rt_info.ini_rt(1, 0) = rotMatrix1(1, 0);
									rt_info.ini_rt(1, 1) = rotMatrix1(1, 1);
									rt_info.ini_rt(1, 2) = rotMatrix1(1, 2);
									rt_info.ini_rt(2, 0) = rotMatrix1(2, 0);
									rt_info.ini_rt(2, 1) = rotMatrix1(2, 1);
									rt_info.ini_rt(2, 2) = rotMatrix1(2, 2);
									auto t = inter_pt_tgt - rotMatrix1 * inter_pt_src;
									rt_info.ini_rt(0, 3) = t[0];
									rt_info.ini_rt(1, 3) = t[1];
									rt_info.ini_rt(2, 3) = t[2];
									rt_info.ini_rt(3, 0) = 0.0f;
									rt_info.ini_rt(3, 1) = 0.0f;
									rt_info.ini_rt(3, 2) = 0.0f;
									rt_info.ini_rt(3, 3) = 1.0f;

									// Optimisation of the transformation matrix
									Fine_RT(rt_info);

									// Matching plane are checked and calibrated
									int match_num = 0;
									Inspect_Structure(rt_info, match_num);

									if (min_match_num == match_num)
										RT_vector.push_back(rt_info);
									else if (min_match_num < match_num)
									{
										min_match_num = match_num;
										RT_vector.clear();
										RT_vector.push_back(rt_info);
									}
								}
								if (possible2)
								{
									rt_info.ini_rt(0, 0) = rotMatrix2(0, 0);
									rt_info.ini_rt(0, 1) = rotMatrix2(0, 1);
									rt_info.ini_rt(0, 2) = rotMatrix2(0, 2);
									rt_info.ini_rt(1, 0) = rotMatrix2(1, 0);
									rt_info.ini_rt(1, 1) = rotMatrix2(1, 1);
									rt_info.ini_rt(1, 2) = rotMatrix2(1, 2);
									rt_info.ini_rt(2, 0) = rotMatrix2(2, 0);
									rt_info.ini_rt(2, 1) = rotMatrix2(2, 1);
									rt_info.ini_rt(2, 2) = rotMatrix2(2, 2);
									auto t = inter_pt_tgt - rotMatrix2 * inter_pt_src;
									rt_info.ini_rt(0, 3) = t[0];
									rt_info.ini_rt(1, 3) = t[1];
									rt_info.ini_rt(2, 3) = t[2];
									rt_info.ini_rt(3, 0) = 0.0f;
									rt_info.ini_rt(3, 1) = 0.0f;
									rt_info.ini_rt(3, 2) = 0.0f;
									rt_info.ini_rt(3, 3) = 1.0f;

									// Optimisation of the transformation matrix
									Fine_RT(rt_info);

									// Matching plane are checked and calibrated
									int match_num = 0;
									Inspect_Structure(rt_info, match_num);

									if (min_match_num == match_num)
										RT_vector.push_back(rt_info);
									else if (min_match_num < match_num)
									{
										min_match_num = match_num;
										RT_vector.clear();
										RT_vector.push_back(rt_info);
									}
								}
								if (possible3)
								{
									rt_info.ini_rt(0, 0) = rotMatrix3(0, 0);
									rt_info.ini_rt(0, 1) = rotMatrix3(0, 1);
									rt_info.ini_rt(0, 2) = rotMatrix3(0, 2);
									rt_info.ini_rt(1, 0) = rotMatrix3(1, 0);
									rt_info.ini_rt(1, 1) = rotMatrix3(1, 1);
									rt_info.ini_rt(1, 2) = rotMatrix3(1, 2);
									rt_info.ini_rt(2, 0) = rotMatrix3(2, 0);
									rt_info.ini_rt(2, 1) = rotMatrix3(2, 1);
									rt_info.ini_rt(2, 2) = rotMatrix3(2, 2);
									auto t = inter_pt_tgt - rotMatrix3 * inter_pt_src;
									rt_info.ini_rt(0, 3) = t[0];
									rt_info.ini_rt(1, 3) = t[1];
									rt_info.ini_rt(2, 3) = t[2];
									rt_info.ini_rt(3, 0) = 0.0f;
									rt_info.ini_rt(3, 1) = 0.0f;
									rt_info.ini_rt(3, 2) = 0.0f;
									rt_info.ini_rt(3, 3) = 1.0f;

									// Optimisation of the transformation matrix
									Fine_RT(rt_info);

									// Matching plane are checked and calibrated
									int match_num = 0;
									Inspect_Structure(rt_info, match_num);

									if (min_match_num == match_num)
										RT_vector.push_back(rt_info);
									else if (min_match_num < match_num)
									{
										min_match_num = match_num;
										RT_vector.clear();
										RT_vector.push_back(rt_info);
									}
								}
								if (possible4)
								{
									rt_info.ini_rt(0, 0) = rotMatrix4(0, 0);
									rt_info.ini_rt(0, 1) = rotMatrix4(0, 1);
									rt_info.ini_rt(0, 2) = rotMatrix4(0, 2);
									rt_info.ini_rt(1, 0) = rotMatrix4(1, 0);
									rt_info.ini_rt(1, 1) = rotMatrix4(1, 1);
									rt_info.ini_rt(1, 2) = rotMatrix4(1, 2);
									rt_info.ini_rt(2, 0) = rotMatrix4(2, 0);
									rt_info.ini_rt(2, 1) = rotMatrix4(2, 1);
									rt_info.ini_rt(2, 2) = rotMatrix4(2, 2);
									auto t = inter_pt_tgt - rotMatrix4 * inter_pt_src;
									rt_info.ini_rt(0, 3) = t[0];
									rt_info.ini_rt(1, 3) = t[1];
									rt_info.ini_rt(2, 3) = t[2];
									rt_info.ini_rt(3, 0) = 0.0f;
									rt_info.ini_rt(3, 1) = 0.0f;
									rt_info.ini_rt(3, 2) = 0.0f;
									rt_info.ini_rt(3, 3) = 1.0f;

									// Optimisation of the transformation matrix
									Fine_RT(rt_info);

									// Matching plane are checked and calibrated
									int match_num = 0;
									Inspect_Structure(rt_info, match_num);

									if (min_match_num == match_num)
										RT_vector.push_back(rt_info);
									else if (min_match_num < match_num)
									{
										min_match_num = match_num;
										RT_vector.clear();
										RT_vector.push_back(rt_info);
									}
								}
							}
							else
								break;
						}
					}
					else
						break;
				}
			}
			else
				break;
		}

		return true;
	}

	bool Registration::Get_Intersections(GPSCO::RT_Info& rt_info)
	{
		// Planar correspondences to point-pair correspondences
		for (const auto& pair1 : rt_info.plane_pairs1.planepairs)
		{
			for (const auto& pair2 : rt_info.plane_pairs2.planepairs)
			{
				for (const auto& pair3 : rt_info.plane_pairs3.planepairs)
				{
					Eigen::Vector3f inter_pt_src(0.0, 0.0, 0.0), inter_pt_tgt(0.0, 0.0, 0.0);
					Eigen::Vector4f plane_a_src, plane_b_src, plane_c_src;
					Eigen::Vector4f plane_a_tgt, plane_b_tgt, plane_c_tgt;
					std::vector<Eigen::Vector3f> intersection_points_src; // Used for averaging
					std::vector<Eigen::Vector3f> intersection_points_tgt; // Used for averaging
					for (const auto& plane_a_index : PlaneGroups_src[rt_info.plane_pairs1.group1_index][pair1.first])
					{
						plane_a_src << Planes_src[plane_a_index].coefficients[0],
							Planes_src[plane_a_index].coefficients[1],
							Planes_src[plane_a_index].coefficients[2],
							Planes_src[plane_a_index].coefficients[3];
						for (const auto& plane_b_index : PlaneGroups_src[rt_info.plane_pairs2.group1_index][pair2.first])
						{
							plane_b_src << Planes_src[plane_b_index].coefficients[0],
								Planes_src[plane_b_index].coefficients[1],
								Planes_src[plane_b_index].coefficients[2],
								Planes_src[plane_b_index].coefficients[3];
							for (const auto& plane_c_index : PlaneGroups_src[rt_info.plane_pairs3.group1_index][pair3.first])
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
					for (const auto& plane_a_index : PlaneGroups_tgt[rt_info.plane_pairs1.group2_index][pair1.second])
					{
						plane_a_tgt << Planes_tgt[plane_a_index].coefficients[0],
							Planes_tgt[plane_a_index].coefficients[1],
							Planes_tgt[plane_a_index].coefficients[2],
							Planes_tgt[plane_a_index].coefficients[3];
						for (const auto& plane_b_index : PlaneGroups_tgt[rt_info.plane_pairs2.group2_index][pair2.second])
						{
							plane_b_tgt << Planes_tgt[plane_b_index].coefficients[0],
								Planes_tgt[plane_b_index].coefficients[1],
								Planes_tgt[plane_b_index].coefficients[2],
								Planes_tgt[plane_b_index].coefficients[3];
							for (const auto& plane_c_index : PlaneGroups_tgt[rt_info.plane_pairs3.group2_index][pair3.second])
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
						rt_info.Intersection_points_src->push_back(
							pcl::PointXYZ(inter_pt_src[0], inter_pt_src[1], inter_pt_src[2]));

						inter_pt_tgt.setZero();
						for (int num = 0; num < intersection_points_tgt.size(); num++)
						{
							inter_pt_tgt += intersection_points_tgt[num];
						}
						inter_pt_tgt /= intersection_points_tgt.size();
						rt_info.Intersection_points_tgt->push_back(
							pcl::PointXYZ(inter_pt_tgt[0], inter_pt_tgt[1], inter_pt_tgt[2]));
					}
				}
			}
		}

		return true;
	}

	bool Registration::Fine_RT(RT_Info& rt_info)
	{
//		// Before optimisation
//		std::cout << "ini_rt:" << std::endl;
//		std::cout << rt_info.ini_rt << std::endl;

		// Extract rotation part (quaternion)
		Eigen::Quaternionf quaternion(rt_info.ini_rt.block<3, 3>(0, 0));
		// Extract translation vector
		Eigen::Vector3f translation(rt_info.ini_rt.block<3, 1>(0, 3));

		// para_q(w, x, y, z), para_t
		double para_q[4] = { quaternion.w(), quaternion.x(), quaternion.y(), quaternion.z() };
		double para_t[3] = { translation[0], translation[1], translation[2] };

		// Ceres: Nonlinear least squares optimisation
		// using point-to-surface distance as residual equation
		ceres::LossFunction* loss_function_dist = new ceres::ScaledLoss(nullptr, 1.0, ceres::TAKE_OWNERSHIP);
		ceres::LossFunction* loss_function_angle = new ceres::ScaledLoss(nullptr, 1.0, ceres::TAKE_OWNERSHIP);
		ceres::Manifold* q_parameterization = new ceres::EigenQuaternionManifold();
		ceres::Problem::Options problem_options;
		ceres::Problem problem(problem_options);
		problem.AddParameterBlock(para_q, 4, q_parameterization);
		problem.AddParameterBlock(para_t, 3);

		auto g11_idx = rt_info.plane_pairs1.group1_index;
		auto g12_idx = rt_info.plane_pairs1.group2_index;
		for (const auto& plane_pair : rt_info.plane_pairs1.planepairs)
		{
			for (const auto& plene_src_idx : PlaneGroups_src[g11_idx][plane_pair.first])
			{
				Eigen::Vector3d pt_src = Planes_src[plene_src_idx].centroid.cast<double>();
				Eigen::Vector3d normal_src = Planes_src[plene_src_idx].normal.cast<double>();
				for (const auto& plene_tgt_idx : PlaneGroups_tgt[g12_idx][plane_pair.second])
				{
					Eigen::Vector3d normal_tgt = Planes_tgt[plene_tgt_idx].normal.cast<double>();
					double d = Planes_tgt[plene_tgt_idx].coefficients[3];
					ceres::CostFunction* cost_function_dist = GPSCO::Error_P2Plane::Create(pt_src, normal_tgt, d);
					problem.AddResidualBlock(cost_function_dist, loss_function_dist, para_q, para_t);
					ceres::CostFunction* cost_function_angle = GPSCO::Error_Normal2Normal::Create(normal_src, normal_tgt);
					problem.AddResidualBlock(cost_function_angle, loss_function_angle, para_q);
				}
			}
		}

		auto g21_idx = rt_info.plane_pairs2.group1_index;
		auto g22_idx = rt_info.plane_pairs2.group2_index;
		for (const auto& plane_pair : rt_info.plane_pairs2.planepairs)
		{
			for (const auto& plene_src_idx : PlaneGroups_src[g21_idx][plane_pair.first])
			{
				Eigen::Vector3d pt_src = Planes_src[plene_src_idx].centroid.cast<double>();
				Eigen::Vector3d normal_src = Planes_src[plene_src_idx].normal.cast<double>();
				for (const auto& plene_tgt_idx : PlaneGroups_tgt[g22_idx][plane_pair.second])
				{
					Eigen::Vector3d normal_tgt = Planes_tgt[plene_tgt_idx].normal.cast<double>();
					double d = Planes_tgt[plene_tgt_idx].coefficients[3];
					ceres::CostFunction* cost_function_dist = GPSCO::Error_P2Plane::Create(pt_src, normal_tgt, d);
					problem.AddResidualBlock(cost_function_dist, loss_function_dist, para_q, para_t);
					ceres::CostFunction* cost_function_angle = GPSCO::Error_Normal2Normal::Create(normal_src, normal_tgt);
					problem.AddResidualBlock(cost_function_angle, loss_function_angle, para_q);
				}
			}
		}

		auto g31_idx = rt_info.plane_pairs3.group1_index;
		auto g32_idx = rt_info.plane_pairs3.group2_index;
		for (const auto& plane_pair : rt_info.plane_pairs3.planepairs)
		{
			for (const auto& plene_src_idx : PlaneGroups_src[g31_idx][plane_pair.first])
			{
				Eigen::Vector3d pt_src = Planes_src[plene_src_idx].centroid.cast<double>();
				Eigen::Vector3d normal_src = Planes_src[plene_src_idx].normal.cast<double>();
				for (const auto& plene_tgt_idx : PlaneGroups_tgt[g32_idx][plane_pair.second])
				{
					Eigen::Vector3d normal_tgt = Planes_tgt[plene_tgt_idx].normal.cast<double>();
					double d = Planes_tgt[plene_tgt_idx].coefficients[3];
					ceres::CostFunction* cost_function_dist = GPSCO::Error_P2Plane::Create(pt_src, normal_tgt, d);
					problem.AddResidualBlock(cost_function_dist, loss_function_dist, para_q, para_t);
					ceres::CostFunction* cost_function_angle = GPSCO::Error_Normal2Normal::Create(normal_src, normal_tgt);
					problem.AddResidualBlock(cost_function_angle, loss_function_angle, para_q);
				}
			}
		}

		ceres::Solver::Options options;
		options.linear_solver_type = ceres::DENSE_QR;
		options.minimizer_progress_to_stdout = false;
		ceres::Solver::Summary summary;
		ceres::Solve(options, &problem, &summary);

		// A quaternion and a translation vector are combined into a transformation matrix
		Eigen::Quaternionf fine_q(para_q[0], para_q[1], para_q[2], para_q[3]);
		Eigen::Vector3f fine_t(para_t[0], para_t[1], para_t[2]);
		rt_info.rt.block<3, 3>(0, 0) = fine_q.toRotationMatrix();
		rt_info.rt.block<3, 1>(0, 3) = fine_t;

//		// --------------------------Test---------------------------------
//		// After optimisation
//		std::cout << "fine_rt:" << std::endl;
//		std::cout << rt_info.rt << std::endl;
//
//		// Calculate the rotation difference
//		Eigen::Matrix3f rotationDiff = rt_info.ini_rt.block<3, 3>(0, 0).transpose() * rt_info.rt.block<3, 3>(0, 0);
//		Eigen::AngleAxisf rotationAngleAxis(rotationDiff);
//		double rotationAngle = rotationAngleAxis.angle();
//
//		// Calculate the translation difference (Euclidean distance)
//		Eigen::Vector3f translationDiff = rt_info.ini_rt.block<3, 1>(0, 3) - rt_info.rt.block<3, 1>(0, 3);
//		double translationDistance = translationDiff.norm();
//
//		// Output the rotation angle and translation distance
//		std::cout << "Rotation Angle Difference: " << rotationAngle << " radians" << std::endl;
//		std::cout << "Translation Distance: " << translationDistance << std::endl;

		return true;
	}

	bool Registration::Inspect_Structure(const RT_Info& rt_temp, int& match_num)
	{
		match_num = 0;

		std::vector<int> pointIdxNKNSearch(1);
		std::vector<float> pointNKNSquaredDistance(1);

		auto g1_src_idx = rt_temp.plane_pairs1.group1_index;
		auto g1_tgt_idx = rt_temp.plane_pairs1.group2_index;
		for (const auto& plane_pair : rt_temp.plane_pairs1.planepairs)
		{
			bool IsMatch = false;
			for (const auto& p_src_idx : PlaneGroups_src[g1_src_idx][plane_pair.first])
			{
				Eigen::Vector3f pt_trans_src = rt_temp.rt.block<3, 3>(0, 0) * Planes_src[p_src_idx].centroid + rt_temp.rt.block<3, 1>(0, 3);

				for (const auto& p_tgt_idx : PlaneGroups_tgt[g1_tgt_idx][plane_pair.second])
				{
					// Minimum distance between two plane point clouds
					Planes_tgt[p_tgt_idx].kdtree.nearestKSearch(pcl::PointXYZ(pt_trans_src[0], pt_trans_src[1], pt_trans_src[2]),
						1, pointIdxNKNSearch, pointNKNSquaredDistance);
					Eigen::Vector3f pt_trans_tgt = rt_temp.rt.inverse().block<3, 3>(0, 0)
						* Planes_tgt[p_tgt_idx].patches->points[pointIdxNKNSearch[0]].getVector3fMap() + rt_temp.rt.inverse().block<3, 1>(0, 3);
					Planes_src[p_src_idx].kdtree.nearestKSearch(pcl::PointXYZ(pt_trans_tgt[0], pt_trans_tgt[1], pt_trans_tgt[2]),
						1, pointIdxNKNSearch, pointNKNSquaredDistance);
					if (sqrt(pointNKNSquaredDistance[0]) < options.max_dist_inspect)
					{
						IsMatch = true;
						break;
					}
				}

				if (IsMatch)
					break;
			}

			if (IsMatch)
				match_num++;
		}

		auto g2_src_idx = rt_temp.plane_pairs2.group1_index;
		auto g2_tgt_idx = rt_temp.plane_pairs2.group2_index;
		for (const auto& plane_pair : rt_temp.plane_pairs2.planepairs)
		{
			bool IsMatch = false;

			for (const auto& p_src_idx : PlaneGroups_src[g2_src_idx][plane_pair.first])
			{
				Eigen::Vector3f pt_trans_src = rt_temp.rt.block<3, 3>(0, 0) * Planes_src[p_src_idx].centroid + rt_temp.rt.block<3, 1>(0, 3);

				for (const auto& p_tgt_idx : PlaneGroups_tgt[g2_tgt_idx][plane_pair.second])
				{
					// Minimum distance between two plane point clouds
					Planes_tgt[p_tgt_idx].kdtree.nearestKSearch(pcl::PointXYZ(pt_trans_src[0], pt_trans_src[1], pt_trans_src[2]),
						1, pointIdxNKNSearch, pointNKNSquaredDistance);
					Eigen::Vector3f pt_trans_tgt = rt_temp.rt.inverse().block<3, 3>(0, 0)
						* Planes_tgt[p_tgt_idx].patches->points[pointIdxNKNSearch[0]].getVector3fMap() + rt_temp.rt.inverse().block<3, 1>(0, 3);
					Planes_src[p_src_idx].kdtree.nearestKSearch(pcl::PointXYZ(pt_trans_tgt[0], pt_trans_tgt[1], pt_trans_tgt[2]),
						1, pointIdxNKNSearch, pointNKNSquaredDistance);
					if (sqrt(pointNKNSquaredDistance[0]) < options.max_dist_inspect)
					{
						IsMatch = true;
						break;
					}
				}
			}

			if (IsMatch)
				match_num++;
		}

		auto g3_src_idx = rt_temp.plane_pairs3.group1_index;
		auto g3_tgt_idx = rt_temp.plane_pairs3.group2_index;
		for (const auto& plane_pair : rt_temp.plane_pairs3.planepairs)
		{
			bool IsMatch = false;
			for (const auto& p_src_idx : PlaneGroups_src[g3_src_idx][plane_pair.first])
			{
				Eigen::Vector3f pt_trans_src = rt_temp.rt.block<3, 3>(0, 0) * Planes_src[p_src_idx].centroid + rt_temp.rt.block<3, 1>(0, 3);

				for (const auto& p_tgt_idx : PlaneGroups_tgt[g3_tgt_idx][plane_pair.second])
				{
					// Minimum distance between two plane point clouds
					Planes_tgt[p_tgt_idx].kdtree.nearestKSearch(pcl::PointXYZ(pt_trans_src[0], pt_trans_src[1], pt_trans_src[2]),
						1, pointIdxNKNSearch, pointNKNSquaredDistance);
					Eigen::Vector3f pt_trans_tgt = rt_temp.rt.inverse().block<3, 3>(0, 0)
						* Planes_tgt[p_tgt_idx].patches->points[pointIdxNKNSearch[0]].getVector3fMap() + rt_temp.rt.inverse().block<3, 1>(0, 3);
					Planes_src[p_src_idx].kdtree.nearestKSearch(pcl::PointXYZ(pt_trans_tgt[0], pt_trans_tgt[1], pt_trans_tgt[2]),
						1, pointIdxNKNSearch, pointNKNSquaredDistance);
					if (sqrt(pointNKNSquaredDistance[0]) < options.max_dist_inspect)
					{
						IsMatch = true;
						break;
					}
				}
			}

			if (IsMatch)
				match_num++;
		}

		return true;
	}

	bool Registration::Evaluate(RT_Info& rt_temp)
	{
		std::vector<int> pointIdxNKNSearch(1);
		std::vector<float> pointNKNSquaredDistance(1);

		for (const auto& plane_src : Planes_src)
		{
			for (const auto& plane_tgt : Planes_tgt)
			{
				auto normal_trans = rt_temp.rt.block<3, 3>(0, 0) * plane_src.normal;
				auto centroid_trans = rt_temp.rt.block<3, 3>(0, 0) * plane_src.centroid + rt_temp.rt.block<3, 1>(0, 3);
				// Determine if there is an overlap
				auto angle = GPSCO::angleBetweenVectors(normal_trans, plane_tgt.normal);
				if (angle < options.overlap_angle)
				{
					auto dist1_temp =
						GPSCO::distancePointToPlane(centroid_trans, plane_tgt.normal, plane_tgt.coefficients[3]);
					auto dist2_temp =
						GPSCO::distancePointToPlane(plane_tgt.centroid, normal_trans, -normal_trans.dot(centroid_trans));
					if ((dist1_temp + dist2_temp) / 2.0 < options.overlap_dist)
					{
						GPSCO::cloudptr transformed_cloud(new GPSCO::cloud);
						pcl::transformPointCloud(*plane_src.patches, *transformed_cloud, rt_temp.rt);
						pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
						kdtree.setInputCloud(transformed_cloud);

						// Minimum distance between two plane point clouds
						plane_tgt.kdtree.nearestKSearch(transformed_cloud->points[0], 1, pointIdxNKNSearch, pointNKNSquaredDistance);
						kdtree.nearestKSearch(plane_tgt.patches->points[pointIdxNKNSearch[0]], 1, pointIdxNKNSearch, pointNKNSquaredDistance);
						if (pointNKNSquaredDistance[0] < options.max_dist_evaluate)
						{
							float overlap_ptnum = 0;
							for (const auto& pt : transformed_cloud->points)
							{
								plane_tgt.kdtree.nearestKSearch(pt, 1, pointIdxNKNSearch, pointNKNSquaredDistance);
								if (pointNKNSquaredDistance[0] < options.max_dist_evaluate)
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

		return true;
	}

	bool Registration::Export_Cluster(std::string outpath)
	{
		if (PlaneGroups_src.empty() || PlaneGroups_tgt.empty())
		{
			spdlog::error("The Plane Clusters is empty!");
			return false;
		}
		else
		{
			std::mt19937 random;
			std::ofstream outfile; // Creating an Output File Stream Object
			int cluster_idx = -1;
			for (int i = 0; i < PlaneGroups_src.size(); ++i)
			{
				for (const auto& cluster : PlaneGroups_src[i])
				{
					auto rgb = random();
					auto red_ = (rgb >> 16) & 0xff;
					auto green_ = (rgb >> 8) & 0xff;
					auto blue_ = rgb & 0xff;

					cluster_idx++;
					std::string file = outpath + "\\src_cluster_" + std::to_string(cluster_idx) + ".txt";
					outfile.open(file);

					for (auto& planeidx : cluster)
					{
						for (const auto& pt : Planes_src[planeidx].points->points)
						{
							outfile << pt.x << " " << pt.y << " " << pt.z << " " <<
									red_ << " " << green_ << " " << blue_ << std::endl;
						}
					}

					outfile.close();
				}
			}
			for (int i = 0; i < PlaneGroups_tgt.size(); ++i)
			{
				for (const auto& cluster : PlaneGroups_tgt[i])
				{
					auto rgb = random();
					auto red_ = (rgb >> 16) & 0xff;
					auto green_ = (rgb >> 8) & 0xff;
					auto blue_ = rgb & 0xff;

					cluster_idx++;
					std::string file = outpath + "\\tgt_cluster_" + std::to_string(cluster_idx) + ".txt";
					outfile.open(file);

					for (auto& planeidx : cluster)
					{
						for (const auto& pt : Planes_tgt[planeidx].points->points)
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

	bool Registration::Export_Groups(std::string outpath)
	{
		if (PlaneGroups_src.empty() || PlaneGroups_tgt.empty())
		{
			spdlog::error("The Plane Groups is empty!");
			return false;
		}
		else
		{
			std::mt19937 random;
			std::ofstream outfile; // Creating an Output File Stream Object
			for (int i = 0; i < PlaneGroups_src.size(); ++i)
			{
				auto rgb = random();
				auto red_ = (rgb >> 16) & 0xff;
				auto green_ = (rgb >> 8) & 0xff;
				auto blue_ = rgb & 0xff;

				std::string file = outpath + "\\src_g_" + std::to_string(i) + ".txt";
				outfile.open(file);
				for (const auto& cluster : PlaneGroups_src[i])
				{
					for (auto& planeidx : cluster)
					{
						for (const auto& pt : Planes_src[planeidx].points->points)
						{
							outfile << pt.x << " " << pt.y << " " << pt.z << " " <<
									red_ << " " << green_ << " " << blue_ << std::endl;
						}
					}
				}
				outfile.close();
			}
			for (int i = 0; i < PlaneGroups_tgt.size(); ++i)
			{
				auto rgb = random();
				auto red_ = (rgb >> 16) & 0xff;
				auto green_ = (rgb >> 8) & 0xff;
				auto blue_ = rgb & 0xff;

				std::string file = outpath + "\\tgt_g_" + std::to_string(i) + ".txt";
				outfile.open(file);
				for (const auto& cluster : PlaneGroups_tgt[i])
				{
					for (auto& planeidx : cluster)
					{
						for (const auto& pt : Planes_tgt[planeidx].points->points)
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
}