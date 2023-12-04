// pcl
#include <pcl/common/impl/centroid.hpp>
// EIGEN
#include <Eigen/Eigenvalues>
// local
#include "Common.h"

namespace GPSCO
{
	PLANE::PLANE()
	{
		points.reset(new GPSCO::cloud);
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
		Eigen::Vector4f centroid_temp; // homogeneous coordinates
		pcl::compute3DCentroid(*points, centroid_temp);

		// get centroid
		centroid = centroid_temp.head<3>();

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
}