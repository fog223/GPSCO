// system
#include <fstream>
#include <random>
#include <set>
// spdlog
#include <spdlog/spdlog.h>
// pcl
#include <pcl/features/normal_3d_omp.h>
#include <pcl/segmentation/region_growing.h>
#include <pcl/io/ply_io.h>
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
			pcl::NormalEstimationOMP<pcl::PointXYZ, pcl::Normal> ne;
			pcl::search::Search<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>);
			pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>);
			ne.setSearchMethod(tree);
			ne.setInputCloud(cloud);
			ne.setKSearch(16);
			ne.compute(*normals);

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

		bool Export_txt(
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

					std::string file = outpath + "p" + std::to_string(i + 1) + ".txt";
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
