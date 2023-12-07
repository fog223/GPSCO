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

		GPSCO::cloudptr Intersection_points_src; // 三平面交点
		GPSCO::cloudptr Intersection_points_tgt; // 三平面交点

		std::vector<Eigen::Matrix3f> Rs; // 候选的变换矩阵
		std::vector<Eigen::Vector3f> Ts; // 候选的平移向量

		float match_score; // Sum of matching scores for the three corresponding plane groups
	};

// 变换矩阵
	class RT_Info
	{
	 public:
		RT_Info();
		~RT_Info();

		//// 排序，从小到大
		// bool operator<(const RT_Info& b) const {
		//	return confidence < b.confidence;
		// }

//		/// \brief 排序，从大到小
//		/// \param b
//		/// \return
		bool
		operator<(const RT_Info& b) const
		{
			if (confidence == b.confidence)
				return plane_match_num > b.plane_match_num;
			else
				return confidence > b.confidence;
		}

		void Compute_RMS();

//		bool
//		operator<(const RT_Info& b) const
//		{
//			if (plane_match_num == b.plane_match_num)
//				return confidence > b.confidence;
//			else
//				return plane_match_num > b.plane_match_num;
//		}

	 public:
		MatchGroup plane_pairs1;
		MatchGroup plane_pairs2;
		MatchGroup plane_pairs3;

		GPSCO::cloudptr Intersection_points_src; // 三平面交点
		GPSCO::cloudptr Intersection_points_tgt; // 三平面交点

		Eigen::Matrix4f rt;     // 变换矩阵
		int plane_match_num; // 匹配的平面数目
		float confidence;     // 置信度
		float rms; // 平面交点间的rms
	};

	namespace Registration
	{
		/// GPSCO algorithm for point cloud registration
		/// \param cloud_src Source point cloud
		/// \param cloud_tgt target point cloud
		/// \param params main parameters in GPSCO algorithm
		/// \param RT transformation matrix
		/// \return Is the transformation matrix calculation successful
		bool Regis(
			const GPSCO::cloudptr cloud_src,
			const GPSCO::cloudptr cloud_tgt,
			GPSCO::Params& params,
			Eigen::Matrix4f& RT);

		/// Planar clustering
		/// \param Planes plane to be clustered
		/// \param parallel_thresh Parallel threshold for planar clustering
		/// \param coplanar_thresh Coplanar threshold for planar clustering
		/// \param PlaneGroups Clustered plane group set
		/// \return Is clustering successful
		bool Plane_Cluster(
			const std::vector<GPSCO::PLANE>& Planes,
			float parallel_thresh,
			float coplanar_thresh,
			std::vector<std::vector<std::vector<int>>>& PlaneGroups);

		/// Save the plane groups separately, with different colour
		/// \param Planes
		/// \param PlaneGroups
		/// \param outpath
		/// \return .txt
		bool Export_groups(
			const std::vector<GPSCO::PLANE>& Planes,
			std::vector<std::vector<std::vector<int>>>& PlaneGroups,
			std::string outpath);

		/// Save the plane clusters separately, with different colour
		/// \param Planes
		/// \param PlaneGroups
		/// \param outpath
		/// \return .txt
		bool Export_cluster(
			const std::vector<GPSCO::PLANE>& Planes,
			std::vector<std::vector<std::vector<int>>>& PlaneGroups,
			std::string outpath);

		/// Evaluate matching scores between plane groups based on moving alignment method
		/// \param Planes_src
		/// \param Planes_tgt
		/// \param PlaneGroups_src plane groups in source point cloud
		/// \param PlaneGroups_tgt plane groups in target point cloud
		/// \param dist_thresh matching threshold
		/// \param group_table A table of matching relationship between plane groups
		/// \return
		bool Compute_Group_Table(
			const std::vector<GPSCO::PLANE>& Planes_src,
			const std::vector<GPSCO::PLANE>& Planes_tgt,
			const std::vector<std::vector<std::vector<int>>>& PlaneGroups_src,
			const std::vector<std::vector<std::vector<int>>>& PlaneGroups_tgt,
			float dist_thresh,
			std::vector<std::vector<std::vector<MatchGroup>>>& group_table);

		/// Maximum number of aligned planes between two plane groups
		/// \param group1_index Subscript index value of group 1
		/// \param group2_index Subscript index value of group 2
		/// \param Planes_src
		/// \param Planes_tgt
		/// \param planegroup_src a plane group in source point cloud
		/// \param planegroup_tgt a plane group in target point cloud
		/// \param dist_thresh
		/// \param match_vec Storage plane correspondence
		/// \return
		bool Moving_alignment(
			int group1_index,
			int group2_index,
			const std::vector<GPSCO::PLANE>& Planes_src,
			const std::vector<GPSCO::PLANE>& Planes_tgt,
			const std::vector<std::vector<int>>& planegroup_src,
			const std::vector<std::vector<int>>& planegroup_tgt,
			float dist_thresh,
			std::vector<MatchGroup>& match_vec);

		// The distance between two plane clusters
		/// \param Planes plane set
		/// \param Normal_avg Average normal vector of plane group
		/// \param cluster1 planar cluster 1
		/// \param cluster2 planar cluster 2
		/// \return distance between planar clusters
		float dist_TwoClsuter(
			const std::vector<GPSCO::PLANE>& Planes,
			const Eigen::Vector3f Normal_avg,
			const std::vector<int>& cluster1,
			const std::vector<int>& cluster2);

		/// In the same point cloud, obtain three plane groups that meet the requirements.
		/// The angles between the plane groups are in the range of angle_min~angle_max. Generally set to 30~150
		/// \param Planes plane set
		/// \param PlaneGroups plane group set
		/// \param angle_min minimum angle
		/// \param angle_max maximum angle
		/// \param group_three_vector Eligible candidates
		/// \return does it exist
		bool Get_Group_Three(
			const std::vector<GPSCO::PLANE>& Planes,
			const std::vector<std::vector<std::vector<int>>>& PlaneGroups,
			float angle_min,
			float angle_max,
			std::vector<GPSCO::Group_Three>& group_three_vector);

		/// The GPSCO algorithm is used to directly determine the optimal match.
		/// When there is more than one optimal situation, the transformation matrix is verified.
		/// \param Planes_src
		/// \param Planes_tgt
		/// \param PlaneGroups_src
		/// \param PlaneGroups_tgt
		/// \param group_table
		/// \param group_three_vector_src candidates in source point cloud
		/// \param group_three_vector_tgt candidates in target point cloud
		/// \param final_rt final registration result
		/// \return
		bool Get_transformation_matrix(
			const std::vector<GPSCO::PLANE>& Planes_src,
			const std::vector<GPSCO::PLANE>& Planes_tgt,
			const std::vector<std::vector<std::vector<int>>>& PlaneGroups_src,
			const std::vector<std::vector<std::vector<int>>>& PlaneGroups_tgt,
			std::vector<std::vector<std::vector<MatchGroup>>>& group_table,
			const std::vector<GPSCO::Group_Three>& group_three_vector_src,
			const std::vector<GPSCO::Group_Three>& group_three_vector_tgt,
			Eigen::Matrix4f& final_rt);
	}
}
