/**
 *  \file   exam_overlap.cpp
 *  \brief  Calculate the overlap between point clouds
 *  \author fog
 *  \email  luochengwen22@mails.ucas.ac.cn
 *  \date   2023/12/8
 *  \note
 */
//

// system
#include <iostream>
#include <fstream>
#include <filesystem>
// spdlog
#include <spdlog/spdlog.h>
// pcl
#include <pcl/io/pcd_io.h>
#include <pcl/io/ply_io.h>
#include <pcl/kdtree/kdtree_flann.h>

namespace fs = std::filesystem;

int main()
{
	// Load point cloud
	std::vector<std::string> files_cloud;
	for (const auto& entry : fs::directory_iterator("D:\\Benchmark_HS\\HS_1\\2-AlignedPointCloud"))
	{
		const fs::path& filePath = entry.path();
		if (fs::is_regular_file(filePath) && filePath.extension() == ".ply")
		{
			files_cloud.push_back(filePath.string());
		}
	}
	std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> cloud_all;
	std::vector<pcl::KdTreeFLANN<pcl::PointXYZ>::Ptr> kdtree_all;
	for (const auto& path : files_cloud)
	{
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
		pcl::io::loadPLYFile<pcl::PointXYZ>(path, *cloud);
		cloud_all.push_back(cloud);
		// Build kdtree for each point cloud data
		pcl::KdTreeFLANN<pcl::PointXYZ>::Ptr kdtree(new pcl::KdTreeFLANN<pcl::PointXYZ>);
		kdtree->setInputCloud(cloud);
		kdtree_all.push_back(kdtree);
	}
	spdlog::info("Load pointcloud success!");

	// Determine which two stations need to be registered with each other.
	std::vector<std::string> files_transmatrix;
	for (const auto& entry : fs::directory_iterator("D:\\Benchmark_HS\\HS_1\\3-GroundTruth"))
	{
		const fs::path& filePath = entry.path();
		if (fs::is_regular_file(filePath) && filePath.extension() == ".txt")
		{
			files_transmatrix.push_back(filePath.string());
		}
	}

	std::vector<std::pair<int, int>> pair_regis_idx;
	for (const auto& path : files_transmatrix)
	{
		std::ifstream file(path);
		if (!file.is_open())
		{
			std::cout << "Can't Load File" << std::endl;
			continue;
		}

		// Determine the indexes of the two stations to be registered.
		size_t start = path.find("transformation") + 14; // Find the position after "transformation"
		size_t end1 = path.find("_", start); // find the position of the next "_"
		size_t end2 = path.find(".", end1); // find the position of the next "."
		std::string extractedNumber1 = path.substr(start, end1 - start);
		std::string extractedNumber2 = path.substr(end1 + 1, end2 - end1 - 1);
		int idx1 = std::stoi(extractedNumber1); // Convert the extracted string to an integer.
		int idx2 = std::stoi(extractedNumber2); // Convert the extracted string to an integer.
		pair_regis_idx.push_back(std::pair<int, int>(idx1 - 1, idx2 - 1));
	}

	std::sort(pair_regis_idx.begin(), pair_regis_idx.end());

	// compute overlap(IoU)
	spdlog::info("Start overlap calculation");

	std::ofstream outfile("D:\\Benchmark_HS\\HS_1\\overlap.txt");
	float dist_thresh = 0.02;
	std::vector<int> pointIdxNKNSearch(1);
	std::vector<float> pointNKNSquaredDistance(1);
	int process = 0;
	float overlap_avg = 0.0f;
	for (const auto& pair : pair_regis_idx)
	{
		process++;
		spdlog::info("{}/{}", process, pair_regis_idx.size());

		float num_overlap = 0;
		for (const auto& pt : cloud_all[pair.first]->points)
		{
			kdtree_all[pair.second]->nearestKSearch(pt, 1, pointIdxNKNSearch, pointNKNSquaredDistance);
			if (sqrt(pointNKNSquaredDistance[0]) < dist_thresh)
			{
				num_overlap++;
			}
		}
		auto overlap = (2 * num_overlap) / float(cloud_all[pair.first]->size() + cloud_all[pair.second]->size());
		outfile << pair.first + 1 << "_sampled.ply and " << pair.second + 1 << "_sampled.ply overlap: " <<
				std::fixed << std::setprecision(3) << overlap << std::endl;
		overlap_avg += overlap;
	}
	outfile.close();

	overlap_avg /= pair_regis_idx.size();
	spdlog::info("The average overlap between point clouds is: {}.", overlap_avg);
	spdlog::info("Overlap calculation complete!");

	return 0;
}