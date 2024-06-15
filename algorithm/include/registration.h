/**
 *  \file   registration.h
 *  \brief  GPSCO: Global Planar Structural Constraint Optimal
 *  \author fog
 *  \email  luochengwen22@mails.ucas.ac.cn
 *  \date   2023/12/6
 *  \note
 */
//

#pragma once
// system
#include <vector>
// Eigen
#include <Eigen/Dense>
// local
#include "common.h"

namespace GPSCO
{
	// Main parameters in GPSCO algorithm
	class Params
	{
	 public:
		Params();
		~Params();

	 public:
		// Plane extraction parameters
		int min_support_points;
		float SmoothnessThreshold;
		float CurvatureThreshold;
		// Cluster parameters
		float parallel_thresh;
		float coplanar_thresh;
		// Match parameters
		float dist_thresh;
	};

	// single matching plane group
	class MatchGroup
	{
	 public:
		MatchGroup();
		~MatchGroup();

		/// Sort, Descending
		/// \param b Another plane group
		/// \return
		bool
		operator<(const MatchGroup& b) const
		{
			return match_num > b.match_num;
		}

	 public:
		int group1_index;
		int group2_index;
		std::vector<std::pair<int, int>> planepairs; // planes_src and planes_tgt correspond one to one
		int match_num; // Matching score of two plane groups
		float area; // Sum of the areas of the matching planes (take the one with the smaller area)
	};

	// In the same point cloud, two plane groups that meet the requirements (30~150°)
	class Group_Base
	{
	 public:
		Group_Base();
		~Group_Base();

	 public:
		int group1_index;
		int group2_index;
		float group_angle;
	};

	// In the same point cloud, three plane groups that meet the requirements (30~150°)
	class Group_Three
	{
	 public:
		Group_Three();
		~Group_Three();

	 public:
		int group1_index;
		int group2_index;
		int group3_index;
		float group_angle12; // The angle between group1 and group2
		float group_angle13; // The angle between group1 and group1
		float group_angle23; // The angle between group2 and group3
	};

	// Three pairs of non-parallel planar groups in different point clouds that meet the requirements
	class Group_Three_Pair
	{
	 public:
		Group_Three_Pair();
		~Group_Three_Pair();

		// Sort, Descending
		bool operator<(const Group_Three_Pair& b) const
		{
			return match_score > b.match_score;
		}

	 public:
		int group_src_index1;
		int group_src_index2;
		int group_src_index3;
		int group_tgt_index1;
		int group_tgt_index2;
		int group_tgt_index3;

		float match_score; // Sum of matching scores for the three corresponding plane groups
		float area_score; // Sum of area scores for the three corresponding plane groups
	};

	// Information necessary to calculate the transformation matrix
	class RT_Info
	{
	 public:
		RT_Info();
		~RT_Info();

		/// Confidence in descending order, Dmean in ascending order while "confidence == b.confidence"
		/// \param b another transformation matrix
		/// \return
		bool
		operator<(const RT_Info& b) const
		{
			return confidence > b.confidence;
		}

	 public:
		MatchGroup plane_pairs1; // plane pairs in plane group 1
		MatchGroup plane_pairs2; // plane pairs in plane group 2
		MatchGroup plane_pairs3; // plane pairs in plane group 3

		GPSCO::cloudptr Intersection_points_src; // Intersection of three planes
		GPSCO::cloudptr Intersection_points_tgt; // Intersection of three planes
		std::vector<std::pair<int, int>> match_plane; // Matching planes that can overlap after final validation

		Eigen::Matrix4f ini_rt; // initial transformation matrix
		Eigen::Matrix4f rt; // transformation matrix
		float confidence; // possibility
	};

	class Registration
	{
	 public:
//		// typedef
//		typedef pcl::PointCloud<pcl::PointXYZ> cloud;
//		typedef pcl::PointCloud<pcl::PointXYZ>::Ptr cloudptr;

		struct Options
		{
			Options()
			{
			}

			// Whether or not to validate, Fewer parallel planes(true)
			bool IsValidate = false;
			// Plane extraction parameters
			float SmoothnessThreshold = 5.0f;
			float CurvatureThreshold = 2.0f;
			// constraint - Reduced computation and false detection.
			// adjusted depending on the number and density of points in the scanned point cloud and the scene.
			int min_support_points = 50;
			int max_plane_num = 30;
			// Average distance between neighbouring points in a point cloud
			float scale = 0.01;
			//
			float SegSize = 0.2;
			// Cluster parameters
			float parallel_thresh = 5.0f;
			float coplanar_thresh = 2.0f;
			// Match parameters
			float e_pl2pldist = 0.05f;
			//
			float min_pl2pldist = 0.3f;
			// angle
			float angle_min = 30;
			float angle_max = 150;
			// Maximum constraint distance between matching planes, Inspect_Structure
			float max_dist_inspect = 1.0;
			// Conditions for matching between two candidates
			float e_angle = 5.0;
			// The angle and distance at which two planes may overlap
			float overlap_angle = 2;
			float overlap_dist = 0.05;
			// The maximum distance between point pairs that overlap, Evaluate
			float max_dist_evaluate = 0.3;

			// wait time (Beyond that time, the registration is considered to have failed.)
			int wait_time = 10; // s
			clock_t start;
		};

		Registration(Options options_ = Options()) : options(options_)
		{
		}

		bool SetCloud(GPSCO::cloudptr cloud_src_, GPSCO::cloudptr cloud_tgt_)
		{
			if (cloud_src_->empty() || cloud_tgt_->empty())
			{
				spdlog::error("Input data is empty.");
				return false;
			}

			cloud_src = cloud_src_;
			cloud_tgt = cloud_tgt_;

			return true;
		}

		Eigen::Matrix4f GetRT()
		{
			return Final_RT;
		}

		/// GPSCO algorithm for point cloud registration
		/// \return Is the transformation matrix calculation successful
		bool Regis();

		/// Configuration items
		Options options;

		// running time
		float time_plane_extra = 0.0f;
		float time_plane_cluster = 0.0f;
		float time_Match = 0.0f;
		float time_Verify = 0.0f;
		float time_Total = 0.0f;

		// Is Success
		bool IsSuccess = false;

	 private:
		/// Planar clustering
		/// \param Planes plane to be clustered
		/// \param PlaneGroups Clustered plane group set
		/// \return Is clustering successful
		bool Plane_Cluster(
			const std::vector<GPSCO::PLANE>& Planes,
			std::vector<std::vector<std::vector<int>>>& PlaneGroups);

		/// Compute the mean normal vector of a planar group
		/// \param Planes
		/// \param PlaneGroups
		/// \param AvgNorm_Group
		/// \return
		bool Compute_AvgPlaneGroupNorm(
			const std::vector<GPSCO::PLANE>& Planes,
			const std::vector<std::vector<std::vector<int>>>& PlaneGroups,
			std::vector<Eigen::Vector3f>& AvgNorm_Group);

		/// Evaluate matching scores between plane groups based on moving alignment method
		bool Compute_Group_Table();

		/// Maximum number of aligned planes between two plane groups
		/// \param group1_index Subscript index value of group 1
		/// \param group2_index Subscript index value of group 2
		/// \param match_vec Storage plane correspondence
		bool Moving_alignment(
			int group1_index,
			int group2_index,
			std::vector<MatchGroup>& match_vec);

		/// The distance between two plane clusters
		/// \param Planes plane set
		/// \param Normal_avg Average normal vector of plane group
		/// \param cluster1 planar cluster 1
		/// \param cluster2 planar cluster 2
		/// \return distance between planar clusters
		float Dist_TwoClsuter(
			const std::vector<GPSCO::PLANE>& Planes,
			const Eigen::Vector3f Normal_avg,
			const std::vector<int>& cluster1,
			const std::vector<int>& cluster2);

		/// In the same point cloud, obtain three plane groups that meet the requirements.
		/// Angle(n1, n2) ~ (angle_min, angle_max). Generally set to 30~150
		/// Angle(n3, n1.corss(n2)) ~ (0, 90-angle_min) or (90+angle_min, 180)
		/// \param Planes plane set
		/// \param PlaneGroups plane group set
		/// \param group_three_vector Eligible candidates
		/// \return does it exist
		bool Get_Group_Three(
			const std::vector<GPSCO::PLANE>& Planes,
			const std::vector<std::vector<std::vector<int>>>& PlaneGroups,
			std::vector<GPSCO::Group_Three>& group_three_vector);

		/// The GPSCO algorithm is used to directly determine the optimal match.
		/// When there is more than one optimal situation, the transformation matrix is verified.
		bool Get_transformation_matrix();

		/// Finding the optimal transformation matrix
		bool Find_optimal_RT(
			const GPSCO::Group_Three_Pair& group_three_pair,
			int& min_match_num,
			std::vector<GPSCO::RT_Info>& RT_vector);

		bool Get_Intersections(GPSCO::RT_Info& rt_info);

		/// Optimisation transformation matrix
		bool Fine_RT(RT_Info& rt_info);

		/// reject the fake and preserve the genuine
		bool Inspect_Structure(
			const RT_Info& rt_temp,
			int& match_num);

		/// Evaluation of the transformation matrix, The score is the proportional sum of the overlapping planes
		/// \param rt_temp Transformation matrix to be evaluated
		bool Evaluate(RT_Info& rt_temp);

		/// Save the plane clusters separately, with different colour
		/// output: .txt
		bool Export_Cluster(std::string outpath);

		/// Save the plane groups separately, with different colour
		/// output: .txt
		bool Export_Groups(std::string outpath);

		// original point cloud
		GPSCO::cloudptr cloud_src;
		GPSCO::cloudptr cloud_tgt;

		std::vector<GPSCO::PLANE> Planes_src;
		std::vector<GPSCO::PLANE> Planes_tgt;

		std::vector<std::vector<std::vector<int>>> PlaneGroups_src;
		std::vector<std::vector<std::vector<int>>> PlaneGroups_tgt;

		std::vector<Eigen::Vector3f> AvgNorm_Group_src;
		std::vector<Eigen::Vector3f> AvgNorm_Group_tgt;

		std::vector<GPSCO::Group_Three> group_three_vector_src;
		std::vector<GPSCO::Group_Three> group_three_vector_tgt;

		std::vector<std::vector<std::vector<MatchGroup>>> group_table; // The relationship table between plane groups

		std::vector<GPSCO::RT_Info> RT_vector;

		Eigen::Matrix4f Final_RT = Eigen::Matrix4f::Identity();

		bool IsSuccessful = false;
	};
}
