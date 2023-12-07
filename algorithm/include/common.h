/**
 *  \file   common.h
 *  \brief  Define the class of plane
 *  \author fog
 *  \email  luochengwen22@mails.ucas.ac.cn
 *  \date   2023/12/4
 *  \note
 */
//

#pragma once

// system
#include <array>
// pcl
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/kdtree/kdtree_flann.h>

namespace GPSCO
{
	// typedef
	typedef pcl::PointCloud<pcl::PointXYZ> cloud;
	typedef pcl::PointCloud<pcl::PointXYZ>::Ptr cloudptr;

	class PLANE
	{
	 public:
		PLANE();
		~PLANE();

		/// Sort, Descending
		/// \param b Another plane
		/// \return
		bool
		operator<(const PLANE& b) const
		{
			return points->size() > b.points->size();
		}

		/// Calculate plane normal vector, centroid, coefficients
		/// \return
		void
		ComputeProperties();

	 public:
		cloudptr points;
		Eigen::Vector3f normal;
		Eigen::Vector3f centroid;
		std::vector<float> coefficients; // Ax+By+Cz+D=0
	};

	// The distance from point to plane
	float distancePointToPlane(const Eigen::Vector3f& point, const Eigen::Vector3f& planeNormal, const float& d);

	// Hash
	using Array3I = std::array<int, 3>;
	struct HashArray3I
	{
		std::size_t operator()(const Array3I& rhs) const noexcept
		{
			return (std::hash<int>()(rhs[0])) ^ (std::hash<int>()(rhs[1])) ^ (std::hash<int>()(rhs[2]));
		}
	};
}