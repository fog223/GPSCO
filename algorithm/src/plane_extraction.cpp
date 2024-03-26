// system
#include <fstream>
#include <random>
#include <set>
// spdlog
#include <spdlog/spdlog.h>
// pcl
#include <pcl/io/ply_io.h>
#include <pcl/features/normal_3d_omp.h>
#include <pcl/segmentation/region_growing.h>
#include <pcl/kdtree/kdtree_flann.h>
// Eigen
#include <Eigen/Eigenvalues>
// local
#include "plane_extraction.h"

namespace GPSCO
{
	namespace PLANE_Extraction
	{
		bool PLANE_Tetect_RegionGrow(
			cloudptr cloud,
			int min_support_points,
			float SmoothnessThreshold,
			float CurvatureThreshold,
			std::vector<GPSCO::PLANE>& outPlanes)
		{
			if (cloud->empty())
			{
				spdlog::error("The cloud data is empty!");
				return false;
			}

			//------------------------------- Normal Estimation -------------------------------
			spdlog::info("Normal Estimation in progress...");
			pcl::search::Search<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>);
			pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>);

			// (normals) windows differ from ubuntu
#ifdef WIN32
			pcl::NormalEstimationOMP<pcl::PointXYZ, pcl::Normal> ne;
			ne.setSearchMethod(tree);
			ne.setInputCloud(cloud);
			ne.setKSearch(16);
			ne.compute(*normals);
#else
			// Compute normals and curvature for each point
			pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
			kdtree.setInputCloud(cloud);
			std::vector<int> pointIdxNKNSearch(16);
			std::vector<float> pointNKNSquaredDistance(16);
			pcl::Normal normal_pt;
			for (const auto& pt : cloud->points)
			{
				pcl::PointCloud<pcl::PointXYZ>::Ptr points(new pcl::PointCloud<pcl::PointXYZ>);
				if (kdtree.nearestKSearch(pt, 16, pointIdxNKNSearch, pointNKNSquaredDistance) > 0)
				{
					for (const auto& idx : pointIdxNKNSearch)
						points->push_back(cloud->points[idx]);
				}

				Eigen::Vector4f centroid_temp; // homogeneous coordinates
				pcl::compute3DCentroid(*points, centroid_temp);

				// Calculate the 3x3 covariance matrix
				Eigen::Matrix3f covariance_matrix;
				pcl::computeCovarianceMatrix(*points, centroid_temp, covariance_matrix);

				Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> es;
				es.compute(covariance_matrix);

				// eigen values (and vectors) are sorted in ascending order
				const auto& eVec = es.eigenvectors();

				// get normal
				auto normal = eVec.col(0).normalized(); // smallest eigenvalue
				normal_pt.normal_x = normal[0];
				normal_pt.normal_y = normal[1];
				normal_pt.normal_z = normal[2];

				if ((es.eigenvalues()[0] + es.eigenvalues()[1] + es.eigenvalues()[2]) == 0.0)
				{
					normal_pt.curvature = 0;
				}
				else
				{
					normal_pt.curvature = es.eigenvalues()[0] / (es.eigenvalues()[0] + es.eigenvalues()[1] + es.eigenvalues()[2]);
				}

				normals->push_back(normal_pt);
			}
#endif

			//------------------------------- Region growing -------------------------------
			spdlog::info("Region growing in progress...");
			pcl::RegionGrowing<pcl::PointXYZ, pcl::Normal> rg;
			std::vector<pcl::PointIndices> clusters;
			rg.setMinClusterSize(min_support_points);
			rg.setMaxClusterSize(1000000);
			rg.setSearchMethod(tree);
			rg.setNumberOfNeighbours(16);
			rg.setInputCloud(cloud);
			rg.setInputNormals(normals);
			rg.setSmoothnessThreshold(SmoothnessThreshold / 180.0 * M_PI);
			rg.setCurvatureThreshold(CurvatureThreshold);
			rg.extract(clusters);

			// Obtain planes for each of the different connectivity areas

			if (clusters.empty())
			{
				spdlog::warn("No planes extracted.");
				return false;
			}
			else
			{
				for (const auto& cluster : clusters)
				{
					GPSCO::PLANE plane;
					for (const auto& idx : cluster.indices)
					{
						plane.points->push_back(cloud->points[idx]);
					}
					plane.ComputeProperties();
					plane.Segment(0.2);
					plane.build();
					outPlanes.push_back(plane);
				}
				spdlog::info("The number of planes extracted is {}.", outPlanes.size());
				return true;
			}
		}

		bool Export_plane(
			std::vector<GPSCO::PLANE>& Planes,
			std::string outpath)
		{
			if (Planes.empty())
			{
				spdlog::error("The Planes is empty!");
				return false;
			}
			else
			{
				std::mt19937 random;
				std::ofstream outfile; // Creating an Output File Stream Object
				for (int i = 0; i < Planes.size(); ++i)
				{
					auto rgb = random();
					auto red_ = (rgb >> 16) & 0xff;
					auto green_ = (rgb >> 8) & 0xff;
					auto blue_ = rgb & 0xff;

					std::string file = outpath + std::to_string(i) + ".txt";
					outfile.open(file);
					for (const auto& pt : Planes[i].points->points)
					{
						outfile << pt.x << " " << pt.y << " " << pt.z << " " <<
								red_ << " " << green_ << " " << blue_ << std::endl;
					}
					outfile.close();
				}
				spdlog::info("Successful export to {}", outpath);
				return true;
			}
		}
	}
}
