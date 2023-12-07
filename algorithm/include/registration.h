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
			if (match_num == b.match_num)
				return score > b.score;
			else
				return match_num > b.match_num;
		}

	 public:
		int group1_index;
		int group2_index;
		std::vector<std::pair<int, int>> planepairs;
		int match_num;
		float score;
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

// 不同点云中, 符合配准阈值的平面组（三个平面）对
	class Group_Three_Pair
	{
	 public:
		Group_Three_Pair();
		~Group_Three_Pair();

		bool operator<(const Group_Three_Pair& b) const
		{
			if (match_num == b.match_num)
			{
				return score > b.score;
			}
			return match_num > b.match_num;
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

		int match_num; // 三组对应平面分数，出现的非1分个数，（例：2+2+2 > 5+5+1）
		float score; // 三个对应平面分数总和
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
			std::vector<std::vector<std::vector<GPSCO::PLANE>>>& PlaneGroups);

		/// Save the plane groups separately, with different colour
		/// \param PlaneGroups
		/// \param outpath
		/// \return .txt
		bool Export_groups(
			std::vector<std::vector<std::vector<GPSCO::PLANE>>>& PlaneGroups,
			std::string outpath);

		/// Save the plane clusters separately, with different colour
		/// \param PlaneGroups
		/// \param outpath
		/// \return .txt
		bool Export_cluster(
			std::vector<std::vector<std::vector<GPSCO::PLANE>>>& PlaneGroups,
			std::string outpath);
	}
}
