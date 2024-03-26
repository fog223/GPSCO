/**
 *  \file   costfunction.h
 *  \brief
 *  \author fog
 *  \email  luochengwen22@mails.ucas.ac.cn
 *  \date   2024/3/21
 *  \note   
 */
//

#ifndef GPSCO_ALGORITHM_INCLUDE_COSTFUNCTION_H_
#define GPSCO_ALGORITHM_INCLUDE_COSTFUNCTION_H_

#include <cmath>
#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include <eigen3/Eigen/Dense>

namespace GPSCO
{
// Residual between intersections of two planes
	struct Error_P2P
	{
		Error_P2P(Eigen::Vector3d pt_src_, Eigen::Vector3d pt_tgt_) : pt_src(pt_src_), pt_tgt(pt_tgt_)
		{
		}

		template<typename T>
		bool operator()(const T* q, const T* t, T* residual) const
		{
			Eigen::Quaternion<T> q_w{ q[3], q[0], q[1], q[2] }; // Note the real and imaginary parts of the quaternion
			Eigen::Matrix<T, 3, 1> t_w{ t[0], t[1], t[2] };
			Eigen::Matrix<T, 3, 1> pt_s{ T(pt_src.x()), T(pt_src.y()), T(pt_src.z()) };
			Eigen::Matrix<T, 3, 1> pt_w;
			pt_w = q_w * pt_s + t_w;

			Eigen::Matrix<T, 3, 1> pt_t(T(pt_tgt.x()), T(pt_tgt.y()), T(pt_tgt.z()));
			residual[0] = (pt_w - pt_t).norm();

			return true;
		}

		static ceres::CostFunction* Create(Eigen::Vector3d pt_src_, Eigen::Vector3d pt_tgt_)
		{
			//  the size of the residuals, 1
			//  the size of the first parameter block in the optimization problem, 4
			//  the size of the second parameter block in the optimization problem, 3
			return (new ceres::AutoDiffCostFunction<Error_P2P, 1, 4, 3>(
				new Error_P2P(pt_src_, pt_tgt_)));
		}

		Eigen::Vector3d pt_src, pt_tgt;
	};

// Residual of point-to-plane distance
	struct Error_P2Plane
	{
		Error_P2Plane(Eigen::Vector3d pt_src_, Eigen::Vector3d normal_tgt_, double d_)
			: pt_src(pt_src_), normal_tgt(normal_tgt_), d(d_)
		{
		}

		template<typename T>
		bool operator()(const T* q, const T* t, T* residual) const
		{
			Eigen::Quaternion<T> q_w{ q[0], q[1], q[2], q[3] }; // Note the real and imaginary parts of the quaternion
			Eigen::Matrix<T, 3, 1> t_w{ t[0], t[1], t[2] };
			Eigen::Matrix<T, 3, 1> pt_s{ T(pt_src.x()), T(pt_src.y()), T(pt_src.z()) };
			Eigen::Matrix<T, 3, 1> pt_w;
			pt_w = q_w * pt_s + t_w;

			Eigen::Matrix<T, 3, 1> norm(T(normal_tgt.x()), T(normal_tgt.y()), T(normal_tgt.z()));

			T dist = norm.dot(pt_w) + T(d);
			if (dist > T(0))
				residual[0] = dist;
			else
				residual[0] = -dist;

			return true;
		}

		static ceres::CostFunction* Create(Eigen::Vector3d pt_src_, Eigen::Vector3d normal_tgt_, double d_)
		{
			//  the size of the residuals, 1
			//  the size of the first parameter block in the optimization problem, 4
			//  the size of the second parameter block in the optimization problem, 3
			return (new ceres::AutoDiffCostFunction<Error_P2Plane, 1, 4, 3>(
				new Error_P2Plane(pt_src_, normal_tgt_, d_)));
		}

		Eigen::Vector3d pt_src, normal_tgt;
		double d;
	};

	// Angular difference between two normal vectors
	struct Error_Normal2Normal
	{
		Error_Normal2Normal(Eigen::Vector3d Nor_src_, Eigen::Vector3d Nor_tgt_)
			: Nor_src(Nor_src_), Nor_tgt(Nor_tgt_)
		{
		}

		template<typename T>
		bool operator()(const T* q, T* residual) const
		{
			Eigen::Quaternion<T> q_w{ q[0], q[1], q[2], q[3] }; // Note the real and imaginary parts of the quaternion
			Eigen::Matrix<T, 3, 1> nor_s{ T(Nor_src[0]), T(Nor_src[1]), T(Nor_src[2]) };
			Eigen::Matrix<T, 3, 1> nor_t{ T(Nor_tgt[0]), T(Nor_tgt[1]), T(Nor_tgt[2]) };
			Eigen::Matrix<T, 3, 1> nor_s_trans;
			nor_s_trans = q_w * nor_s;

			// Calculate the angle between nor_s_trans and nor_t
			if (nor_s_trans.dot(nor_t) > T(0))
			{
				T cos_angle = nor_s_trans.dot(nor_t) / (nor_s_trans.norm() * nor_t.norm());
				residual[0] = T(acos(cos_angle));
			}
			else
			{
				T cos_angle = nor_s_trans.dot(-nor_t) / (nor_s_trans.norm() * nor_t.norm());
				residual[0] = T(acos(cos_angle));
			}

			return true;
		}

		static ceres::CostFunction* Create(Eigen::Vector3d Nor_src_, Eigen::Vector3d Nor_tgt_)
		{
			//  the size of the residuals, 1
			//  the size of the first parameter block in the optimization problem, 4
			return (new ceres::AutoDiffCostFunction<Error_Normal2Normal, 1, 4>(
				new Error_Normal2Normal(Nor_src_, Nor_tgt_)));
		}

		Eigen::Vector3d Nor_src, Nor_tgt;
	};
}

#endif //GPSCO_ALGORITHM_INCLUDE_COSTFUNCTION_H_
