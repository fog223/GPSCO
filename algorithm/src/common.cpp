// system
#include <algorithm>
// EIGEN
#include <Eigen/Eigenvalues>
// pcl
#include <pcl/common/centroid.h>
#include <pcl/keypoints/uniform_sampling.h>
#include <pcl/registration/ia_fpcs.h>
// local
#include "common.h"

namespace GPSCO
{
	PLANE::PLANE()
	{
		points.reset(new GPSCO::cloud);
		patches.reset(new GPSCO::cloud);
		normal.setZero();
		centroid.setZero();
		coefficients.clear();
	}

	PLANE::~PLANE()
	{
	}

	void
	PLANE::ComputeProperties()
	{
//		// Approach 1
//		Eigen::Vector4f centroid_temp; // homogeneous coordinates
//		pcl::compute3DCentroid(*points, centroid_temp);
//		// get centroid
//		centroid = centroid_temp.head<3>();

		// Approach 2
		// To avoid noise interference, the median is used instead of the mean Centroid
		// Extract x, y, and z coordinates
		std::vector<float> x_values, y_values, z_values;
		for (const auto& point : points->points)
		{
			x_values.push_back(point.x);
			y_values.push_back(point.y);
			z_values.push_back(point.z);
		}

		// Sort coordinates
		std::sort(x_values.begin(), x_values.end());
		std::sort(y_values.begin(), y_values.end());
		std::sort(z_values.begin(), z_values.end());

		// Calculate medians
		float x_median, y_median, z_median;

		size_t num_points = points->points.size();
		if (num_points % 2 == 0)
		{
			size_t mid_idx = num_points / 2;
			x_median = (x_values[mid_idx - 1] + x_values[mid_idx]) / 2;
			y_median = (y_values[mid_idx - 1] + y_values[mid_idx]) / 2;
			z_median = (z_values[mid_idx - 1] + z_values[mid_idx]) / 2;
		}
		else
		{
			size_t mid_idx = num_points / 2;
			x_median = x_values[mid_idx];
			y_median = y_values[mid_idx];
			z_median = z_values[mid_idx];
		}

		// Get centroid
		centroid[0] = x_median;
		centroid[1] = y_median;
		centroid[2] = z_median;

		Eigen::Vector4f centroid_temp(centroid[0], centroid[1], centroid[2], 1.0); // homogeneous coordinates

		// Calculate the 3x3 covariance matrix
		Eigen::Matrix3f covariance_matrix;
		pcl::computeCovarianceMatrix(*points, centroid_temp, covariance_matrix);

		Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> es;
		es.compute(covariance_matrix);

		// eigen values (and vectors) are sorted in ascending order
		const auto& eVec = es.eigenvectors();

		// get normal
		normal = eVec.col(0).normalized(); // smallest eigenvalue

		// get coefficients
		coefficients.push_back(normal[0]);
		coefficients.push_back(normal[1]);
		coefficients.push_back(normal[2]);
		coefficients.push_back(-normal.dot(centroid));
	}

	void
	PLANE::Segment(float size)
	{
		pcl::UniformSampling<pcl::PointXYZ> us;
		us.setInputCloud(points);
		us.setRadiusSearch(size);
		us.filter(*patches);
	}

	void
	PLANE::build()
	{
		kdtree.setInputCloud(patches);
	}

	float
	distancePointToPlane(const Eigen::Vector3f& point, const Eigen::Vector3f& planeNormal, const float& d)
	{
		float distance = (point.dot(planeNormal) + d) / planeNormal.norm();
		return std::fabs(distance);
	}

	float
	angleBetweenVectors(const Eigen::Vector3f& v1, const Eigen::Vector3f& v2)
	{
		float dot = v1.dot(v2);
		float v1_norm = v1.norm();
		float v2_norm = v2.norm();
		double cos_angle = dot / (v1_norm * v2_norm);

		// Handling of results outside the range [-1, 1] due to numerical errors
		// std::max, std::min must be "double"
		cos_angle = std::max(-1.0, std::min(1.0, cos_angle));

		float radian_angle = std::acos(cos_angle);
		float degree_angle = radian_angle * 180.0 / M_PI;

		// Ensure that the angle is less than 90 degrees
		if (degree_angle > 90.0)
		{
			degree_angle = 180.0 - degree_angle;
		}

		return degree_angle;
	}

	float
	distancePointToPlane(const Eigen::Vector3f& point, const Eigen::Vector3f& planeNormal, const double& d)
	{
		float distance = (point.dot(planeNormal) + d) / planeNormal.norm();
		return std::fabs(distance);
	}

	bool Compute_density(cloudptr points, float& density)
	{
//		float dist_max = 100.0;
//		density = pcl::getMeanPointDensity<pcl::PointXYZ>(points, dist_max, 8);
		return true;
	}

	bool Compute_rotMatrix(Eigen::Vector3f normal1_src, Eigen::Vector3f normal2_src,
		Eigen::Vector3f normal1_tgt, Eigen::Vector3f normal2_tgt, Eigen::Matrix3f& rotMatrix)
	{
		Eigen::Matrix3f n_matrix;
		n_matrix.col(0) = normal1_src;
		n_matrix.col(1) = normal2_src;
		n_matrix.col(2) = normal1_src.cross(normal2_src);
		Eigen::Matrix3f n_inverse = n_matrix.inverse();

		Eigen::Matrix3f T;
		T.col(0) = normal1_tgt;
		T.col(1) = normal2_tgt;
		T.col(2) = normal1_tgt.cross(normal2_tgt);
		T *= n_inverse;

		// Orthogonalise T to obtain the rotation matrix R
		Eigen::JacobiSVD<Eigen::Matrix3f> svd(T, Eigen::ComputeFullU | Eigen::ComputeFullV);
		rotMatrix = svd.matrixU() * svd.matrixV().transpose();

		return true;
	}
}