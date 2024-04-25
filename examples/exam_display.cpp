/**
 *  \file   exam_display.cpp
 *  \brief  Convert all successfully matched point clouds to the same coordinate system
 *  \author fog
 *  \email  luochengwen22@mails.ucas.ac.cn
 *  \date   2024/3/29
 *  \note
 */
//

// system
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <filesystem>
// spdlog
#include <spdlog/spdlog.h>
// pcl
#include <pcl/io/pcd_io.h>
#include <pcl/io/ply_io.h>
#include <pcl/common/transforms.h>

namespace fs = std::filesystem;

// Custom comparison function to sort based on numeric values in file names
bool numericStringCompare(const std::string& a, const std::string& b)
{
	// Extract file names from full paths
	std::string fileNameA = a.substr(a.find_last_of("/\\") + 1);
	std::string fileNameB = b.substr(b.find_last_of("/\\") + 1);

	// Extract numeric part from file names
	int numA = std::stoi(fileNameA);
	int numB = std::stoi(fileNameB);

	// Compare using numeric values
	return numA < numB;
}

int main(int argc, char** argv)
{
	if (argc != 3)
	{
		spdlog::error("Parameter input error!");
		return -1;
	}

	std::string root_dir = argv[1];

	std::vector<std::string> files_path;
	for (const auto& entry : fs::directory_iterator(root_dir + "/1-RawPointCloud"))
	{
		const fs::path& filePath = entry.path();
		if (fs::is_regular_file(filePath) && filePath.extension() == ".pcd")
		{
			files_path.push_back(filePath.string());
		}
	}

	// Sort by file name
	std::sort(files_path.begin(), files_path.end(), numericStringCompare);

	// All point cloud data
	std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> cloud_all;
	int process = 0;
	for (const auto& file : files_path)
	{
		process++;
		spdlog::info("Data processing: {}/{}", process, files_path.size());

		// Load point cloud data
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_(new pcl::PointCloud<pcl::PointXYZ>);
		if (pcl::io::loadPCDFile<pcl::PointXYZ>(file, *cloud_) == -1)
		{
			spdlog::error("{} load failure!", file);
			return 0;
		}
		else
		{
			cloud_all.push_back(cloud_);
		}
	}
	spdlog::info("Data load completed.");

	// load GPSCO_RT
	std::vector<Eigen::Matrix4f> Gpsco_RTs;
	std::vector<std::pair<int, int>> pair_regis_idx;
	std::ifstream file(argv[2]);
	if (!file.is_open())
	{
		spdlog::error("{} Can't Load File.", argv[2]);
		return 0;
	}
	else
	{
		int index1, index2;
		Eigen::Matrix4f rt_temp;
		while (!file.eof())
		{
			file >> index1 >> index2;
			for (int i = 0; i < 4; ++i)
			{
				for (int j = 0; j < 4; ++j)
				{
					file >> rt_temp(i, j);
				}
			}

			pair_regis_idx.push_back(std::pair<int, int>(index1, index2));
			Gpsco_RTs.push_back(rt_temp);
		}
		// Close the file
		file.close();
	}
	pair_regis_idx.pop_back();
	Gpsco_RTs.pop_back();

	spdlog::info("Gpsco_RTs load completed.");

//	// HS1
//	std::vector<std::vector<int>> trans_path(29);
//	trans_path[0] = std::vector < int > { 19, 31, 33, 35, 37, 38 };
//	trans_path[1] = std::vector < int > { 31, 33, 35, 37, 38 };
//	trans_path[2] = std::vector < int > { 33, 35, 37, 38 };
//	trans_path[3] = std::vector < int > { 35, 37, 38 };
//	trans_path[4] = std::vector < int > { 37, 38 };
//	trans_path[5] = std::vector < int > { 38 };
//	trans_path[6] = std::vector < int > { -39, 38 };
//	trans_path[7] = std::vector < int > { -41, -39, 38 };
//	trans_path[8] = std::vector < int > { -43, -41, -39, 38 };
//	trans_path[9] = std::vector < int > { -44, -43, -41, -39, 38 };
//	trans_path[10] = std::vector < int > { 1, 3, 5, 7 };
//	trans_path[11] = std::vector < int > { 3, 5, 7 };
//	trans_path[12] = std::vector < int > { 5, 7 };
//	trans_path[13] = std::vector < int > { 7 };
//	trans_path[14] = std::vector < int > {};
//	trans_path[15] = std::vector < int > { -9 };
//	trans_path[16] = std::vector < int > { -11, -9 };
//	trans_path[17] = std::vector < int > { -13, -11, -9 };
//	trans_path[18] = std::vector < int > { -15, -13, -11, -9 };
//	trans_path[19] = std::vector < int > { -17, -15, -13, -11, -9 };
//	trans_path[20] = std::vector < int > { -18, -15, -13, -11, -9 };
//	trans_path[21] = std::vector < int > { -16, -13, -11, -9 };
//	trans_path[22] = std::vector < int > { -14, -11, -9 };
//	trans_path[23] = std::vector < int > { -12, -9 };
//	trans_path[24] = std::vector < int > { -10 };
//	trans_path[25] = std::vector < int > { -8, 7 };
//	trans_path[26] = std::vector < int > { -6, 5, 7 };
//	trans_path[27] = std::vector < int > { -4, 3, 5, 7 };
//	trans_path[28] = std::vector < int > { -2, 1, 3, 5, 7 };

//	// HS2
//	std::vector<std::vector<int>> trans_path(15);
//	trans_path[0] = std::vector < int > { 6, 7, 10 };
//	trans_path[1] = std::vector < int > { 7, 10 };
//	trans_path[2] = std::vector < int > { 10 };
//	trans_path[3] = std::vector < int > { 12, 17 };
//	trans_path[4] = std::vector < int > { 13, 15, 17 };
//	trans_path[5] = std::vector < int > { 15, 17 };
//	trans_path[6] = std::vector < int > { 17 };
//	trans_path[7] = std::vector < int > {};
//	trans_path[8] = std::vector < int > { -19 };
//	trans_path[9] = std::vector < int > { -20, -19 };
//	trans_path[10] = std::vector < int > { 1, -21, -19 };
//	trans_path[11] = std::vector < int > { -21, -19 };
//	trans_path[12] = std::vector < int > { 3, -16, 17 };
//	trans_path[13] = std::vector < int > { -16, 17 };
//	trans_path[14] = std::vector < int > { -14, 15, 17 };

//	// Whu-Park
//	std::vector<std::vector<int>> trans_path(32);
//	trans_path[0] = std::vector < int > { 0, 1, 2, 3, 4, 5, 8, 9, 10, 11, 12 };
//	trans_path[1] = std::vector < int > { 1, 2, 3, 4, 5, 8, 9, 10, 11, 12 };
//	trans_path[2] = std::vector < int > { 2, 3, 4, 5, 8, 9, 10, 11, 12 };
//	trans_path[3] = std::vector < int > { 3, 4, 5, 8, 9, 10, 11, 12 };
//	trans_path[4] = std::vector < int > { 4, 5, 8, 9, 10, 11, 12 };
//	trans_path[5] = std::vector < int > { 5, 8, 9, 10, 11, 12 };
//	trans_path[6] = std::vector < int > { 6, 4, 5, 8, 9, 10, 11, 12 };
//	trans_path[7] = std::vector < int > { 7, 5, 8, 9, 10, 11, 12 };
//	trans_path[8] = std::vector < int > { 8, 9, 10, 11, 12 };
//	trans_path[9] = std::vector < int > { 9, 10, 11, 12 };
//	trans_path[10] = std::vector < int > { 10, 11, 12 };
//	trans_path[11] = std::vector < int > { 11, 12 };
//	trans_path[12] = std::vector < int > { 12 };
//	trans_path[13] = std::vector < int > {};
//	trans_path[14] = std::vector < int > { 13 };
//	trans_path[15] = std::vector < int > { 14, 13 };
//	trans_path[16] = std::vector < int > { 15, 14, 13 };
//	trans_path[17] = std::vector < int > { 16, 6, 4, 5, 8, 9, 10, 11, 12 };
//	trans_path[18] = std::vector < int > { 17, 6, 4, 5, 8, 9, 10, 11, 12 };
//	trans_path[19] = std::vector < int > { 18, 17, 6, 4, 5, 8, 9, 10, 11, 12 };
//	trans_path[20] = std::vector < int > { 19, 18, 17, 6, 4, 5, 8, 9, 10, 11, 12 };
//	trans_path[21] = std::vector < int > { 20, 15, 14, 13 };
//	trans_path[22] = std::vector < int > { 21, 23, 24, 26, 25, 28, 29, 30 };
//	trans_path[23] = std::vector < int > { 22, 23, 24, 26, 25, 28, 29, 30 };
//	trans_path[24] = std::vector < int > { 23, 24, 26, 25, 28, 29, 30 };
//	trans_path[25] = std::vector < int > { 24, 26, 25, 28, 29, 30 };
//	trans_path[26] = std::vector < int > { 25, 28, 29, 30 };
//	trans_path[27] = std::vector < int > { 26, 25, 28, 29, 30 };
//	trans_path[28] = std::vector < int > { 27, 26, 25, 28, 29, 30 };
//	trans_path[29] = std::vector < int > { 28, 29, 30 };
//	trans_path[30] = std::vector < int > { 29, 30 };
//	trans_path[31] = std::vector < int > { 30 };

//	// Whu-Campus
//	std::vector<std::vector<int>> trans_path(10);
//	trans_path[0] = std::vector<int>{ 0, 1 };
//	trans_path[1] = std::vector<int>{ 1 };
//	trans_path[2] = std::vector<int>{};
//	trans_path[3] = std::vector<int>{ 2 };
//	trans_path[4] = std::vector<int>{ 3, 2 };
//	trans_path[5] = std::vector<int>{ 4, 3, 2 };
//	trans_path[6] = std::vector<int>{ 5, 4, 3, 2 };
//	trans_path[7] = std::vector<int>{ 6 };
//	trans_path[8] = std::vector<int>{ 7, 6 };
//	trans_path[9] = std::vector<int>{ 8, 7, 6 };

//	// HS1,HS2,Whu-Park,Whu-Campus
//	std::mt19937 random;
//	for (int i = 0; i < cloud_all.size(); i++)
//	{
//		auto rgb = random();
//		auto red_ = (rgb >> 16) & 0xff;
//		auto green_ = (rgb >> 8) & 0xff;
//		auto blue_ = rgb & 0xff;
//
//		Eigen::Matrix4f trans_matrix = Eigen::Matrix4f::Identity();
//		for (int j = 0; j < trans_path[i].size(); j++)
//		{
//			if (trans_path[i][j] >= 0)
//				trans_matrix = Gpsco_RTs[trans_path[i][j]] * trans_matrix;
//			else
//				trans_matrix = Gpsco_RTs[-trans_path[i][j]].inverse() * trans_matrix;
//		}
//
//		pcl::PointCloud<pcl::PointXYZ>::Ptr trans_cloud(new pcl::PointCloud<pcl::PointXYZ>);
//		pcl::PointCloud<pcl::PointXYZRGB>::Ptr trans_cloud_rgb(new pcl::PointCloud<pcl::PointXYZRGB>);
//		pcl::transformPointCloud(*cloud_all[i], *trans_cloud, trans_matrix);
//		pcl::PointXYZRGB pt_rgb;
//		for (const auto& pt : trans_cloud->points)
//		{
//			pt_rgb.x = pt.x;
//			pt_rgb.y = pt.y;
//			pt_rgb.z = pt.z;
//			pt_rgb.r = red_;
//			pt_rgb.g = green_;
//			pt_rgb.b = blue_;
//
//			trans_cloud_rgb->push_back(pt_rgb);
//		}
//		trans_cloud_rgb->width = trans_cloud_rgb->size();
//		trans_cloud_rgb->height = 1;
//
//		std::string outpath = "D:/References/GPSCO/Results/Whu-Campus/" + std::to_string(i + 1) + "_trans.pcd";
//		pcl::io::savePCDFileASCII(outpath, *trans_cloud_rgb);
//	}

	// ETH-Apartment, ETH-Stairs
	std::mt19937 random;
	Eigen::Matrix4f trans_matrix = Eigen::Matrix4f::Identity();
	for (int i = 0; i < cloud_all.size(); i++)
	{
		auto rgb = random();
		auto red_ = (rgb >> 16) & 0xff;
		auto green_ = (rgb >> 8) & 0xff;
		auto blue_ = rgb & 0xff;

		if (i == 0)
		{
			pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_rgb(new pcl::PointCloud<pcl::PointXYZRGB>);
			pcl::PointXYZRGB pt_rgb;
			for (const auto& pt : cloud_all[0]->points)
			{
				pt_rgb.x = pt.x;
				pt_rgb.y = pt.y;
				pt_rgb.z = pt.z;
				pt_rgb.r = red_;
				pt_rgb.g = green_;
				pt_rgb.b = blue_;

				cloud_rgb->push_back(pt_rgb);
			}
			cloud_rgb->width = cloud_rgb->size();
			cloud_rgb->height = 1;

			std::string outpath = "D:/References/GPSCO/Results/ETH-Hauptgebaude/" + std::to_string(i) + "_trans.pcd";
			pcl::io::savePCDFileASCII(outpath, *cloud_rgb);
		}
		else
		{
			trans_matrix = trans_matrix * Gpsco_RTs[i - 1].inverse();
			pcl::PointCloud<pcl::PointXYZ>::Ptr trans_cloud(new pcl::PointCloud<pcl::PointXYZ>);
			pcl::PointCloud<pcl::PointXYZRGB>::Ptr trans_cloud_rgb(new pcl::PointCloud<pcl::PointXYZRGB>);
			pcl::transformPointCloud(*cloud_all[i], *trans_cloud, trans_matrix);
			pcl::PointXYZRGB pt_rgb;
			for (const auto& pt : trans_cloud->points)
			{
				pt_rgb.x = pt.x;
				pt_rgb.y = pt.y;
				pt_rgb.z = pt.z;
				pt_rgb.r = red_;
				pt_rgb.g = green_;
				pt_rgb.b = blue_;

				trans_cloud_rgb->push_back(pt_rgb);
			}
			trans_cloud_rgb->width = trans_cloud_rgb->size();
			trans_cloud_rgb->height = 1;

			std::string outpath = "D:/References/GPSCO/Results/ETH-Hauptgebaude/" + std::to_string(i) + "_trans.pcd";
			pcl::io::savePCDFileASCII(outpath, *trans_cloud_rgb);
		}
	}

	return 0;
}