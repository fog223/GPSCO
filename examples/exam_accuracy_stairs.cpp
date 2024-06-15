/**
 *  \file   exam_accuracy_stairs.cpp
 *  \brief  Evaluation of the accuracy of the algorithm GPSCO(ETH,ALS Datasets - Stairs)
 *  \author fog
 *  \email  luochengwen22@mails.ucas.ac.cn
 *  \date   2024/3/28
 *  \note
 */
//

// system
#include <iostream>
#include <vector>
#include <filesystem>
// spdlog
#include <spdlog/spdlog.h>
// pcl
#include <pcl/io/pcd_io.h>
#include <pcl/io/ply_io.h>
// local
#include "plane_extraction.h"
#include "registration.h"

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
	if (argc != 4)
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
	std::vector<GPSCO::cloudptr> cloud_all;
	int process = 0;
	for (const auto& file : files_path)
	{
		process++;
		spdlog::info("Data processing: {}/{}", process, files_path.size());

		// Load point cloud data
		GPSCO::cloudptr cloud_(new GPSCO::cloud);
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
	spdlog::info("Data processing completed.");

	// Get the pose of each frame
	std::vector<Eigen::Matrix4f> RTs_Pose;
	std::ifstream file(root_dir + "/pose_scanner_leica.csv");
	if (!file.is_open())
	{
		spdlog::error("{} Can't Load File.", argv[2]);
		return 0;
	}
	else
	{
		// Read the file line by line
		std::string line;
		Eigen::Matrix4f rt_temp;
		while (std::getline(file, line))
		{
			std::stringstream ss(line);
			std::vector<std::string> row;

			// Split each line by commas
			std::string cell;
			while (std::getline(ss, cell, ','))
			{
				row.push_back(cell);
			}

			if (row[0] != "poseId")
			{
				rt_temp(0, 0) = std::stof(row[2]);
				rt_temp(0, 1) = std::stof(row[3]);
				rt_temp(0, 2) = std::stof(row[4]);
				rt_temp(0, 3) = std::stof(row[5]);
				rt_temp(1, 0) = std::stof(row[6]);
				rt_temp(1, 1) = std::stof(row[7]);
				rt_temp(1, 2) = std::stof(row[8]);
				rt_temp(1, 3) = std::stof(row[9]);
				rt_temp(2, 0) = std::stof(row[10]);
				rt_temp(2, 1) = std::stof(row[11]);
				rt_temp(2, 2) = std::stof(row[12]);
				rt_temp(2, 3) = std::stof(row[13]);
				rt_temp(3, 0) = std::stof(row[14]);
				rt_temp(3, 1) = std::stof(row[15]);
				rt_temp(3, 2) = std::stof(row[16]);
				rt_temp(3, 3) = std::stof(row[17]);

				RTs_Pose.push_back(rt_temp);
			}
		}
		// Close the file
		file.close();
	}

	// Determine which two stations need to be aligned with each other.
	std::vector<Eigen::Matrix4f> RTs_GroundTruth;
	std::vector<std::pair<int, int>> pair_regis_idx;
	for (int i = 0; i < RTs_Pose.size() - 1; i++)
	{
		pair_regis_idx.push_back(std::pair<int, int>(i, i + 1));
		RTs_GroundTruth.push_back(RTs_Pose[i + 1].inverse() * RTs_Pose[i]);
	}

	spdlog::info("GroundTruth load completed.");

	// Save registration results to txt file
	std::ofstream outfile_rt(argv[2]);
	// Check if the file was opened successfully
	if (!outfile_rt.is_open())
	{
		std::cerr << "Error opening file: " << argv[2] << std::endl;
		return 1;
	}

	// Save running time and regis info to csv file
	std::ofstream outfile(argv[3]);
	// Check if the file was opened successfully
	if (!outfile.is_open())
	{
		std::cerr << "Error opening file: " << argv[3] << std::endl;
		return 1;
	}
	// Write the header row
	outfile << "Pair,plane_extra(s),plane_cluster(s),Match(s),Verify(s),total(s),e_R,e_T,Is_success\n";

	// SRR, R_rmse, T_rmse
	// Rotation threshold is 5°, translation threshold is 1m
	float R_thresh = 5;
	float T_thresh = 1.0;
	int num_rtmatch = 0; // correct match
	std::vector<float> R_error;
	float R_error_avg = 0.0f;
	std::vector<float> T_error;
	float T_error_avg = 0.0f;
	// running time
	float time_plane_extra = 0.0f;
	float time_plane_cluster = 0.0f;
	float time_Match = 0.0f;
	float time_Verify = 0.0f;
	float time_total = 0.0f;
	/// Registration
	int process_regis = 0;
	for (const auto& pair : pair_regis_idx)
	{
		process_regis++;
		spdlog::info("processing : {} and {}, {} / {}. ", pair.first, pair.second, process_regis, pair_regis_idx.size());

		GPSCO::Registration::Options options;
		options.min_support_points = 50;
		options.max_plane_num = 10;
		options.SmoothnessThreshold = 5.0;
		options.CurvatureThreshold = 2.0;
		options.parallel_thresh = 5.0;
		options.coplanar_thresh = 2.0;
		options.e_pl2pldist = 0.05;

		if ((pair.first == 15 && pair.second == 16)
			|| (pair.first == 16 && pair.second == 17)
			|| (pair.first == 26 && pair.second == 27)
			|| (pair.first == 27 && pair.second == 28))
			options.max_plane_num = 20;

		if ((pair.first == 28 && pair.second == 29)
			|| (pair.first == 29 && pair.second == 30))
			options.max_plane_num = 30;

		GPSCO::Registration regis_(options);

		regis_.SetCloud(cloud_all[pair.first], cloud_all[pair.second]);
		regis_.Regis();

		time_plane_extra += regis_.time_plane_extra;
		time_plane_cluster += regis_.time_plane_cluster;
		time_Match += regis_.time_Match;
		time_Verify += regis_.time_Verify;

		if(regis_.IsSuccess)
			time_total += regis_.time_Total;

		// Write the data rows
		outfile_rt << pair.first << " " << pair.second << std::endl;
		outfile_rt << regis_.GetRT() << std::endl;

		auto pair_str = std::to_string(pair.first) + " % " + std::to_string(pair.second);
		outfile << pair_str << "," << regis_.time_plane_extra << ","
				<< regis_.time_plane_cluster << "," << regis_.time_Match << ","
				<< regis_.time_Verify << "," << regis_.time_Total;

		Eigen::Matrix3f R1 = RTs_GroundTruth[process_regis - 1].block<3, 3>(0, 0);
		Eigen::Vector3f T1 = RTs_GroundTruth[process_regis - 1].block<3, 1>(0, 3);
		Eigen::Matrix3f R2 = regis_.GetRT().block<3, 3>(0, 0);
		Eigen::Vector3f T2 = regis_.GetRT().block<3, 1>(0, 3);

		Eigen::AngleAxisf aaDiff(R1.transpose() * R2);
		float angleDiff = aaDiff.angle();

		// Convert radians to angles
		float angleDiffDegrees = angleDiff * 180.0 / M_PI;
		// Calculate the distance difference
		float e_dist = (T1 - T2).norm();

		if (angleDiffDegrees < R_thresh && e_dist < T_thresh)
		{
			num_rtmatch++;
			R_error.push_back(angleDiffDegrees);
			R_error_avg += pow(angleDiffDegrees, 2);
			T_error.push_back(e_dist);
			T_error_avg += pow(e_dist, 2);

			outfile << "," << angleDiffDegrees << "," << e_dist << ",True\n";
		}
		else
		{
			spdlog::warn("{} and {} Registration failed!", pair.first, pair.second);
			spdlog::info("angleDiffDegrees: {}.", angleDiffDegrees);
			spdlog::info("e_dist: {}.", e_dist);

			outfile << "," << angleDiffDegrees << "," << e_dist << ",False\n";
		}
	}
	// Close the file
	outfile_rt.close();

	outfile << "\navg_plane_extra(s),avg_plane_cluster(s),avg_Match(s),avg_Verify(s),avg_total(s)\n";
	outfile << time_plane_extra / float(RTs_GroundTruth.size()) << ","
			<< time_plane_cluster / float(RTs_GroundTruth.size()) << ","
			<< time_Match / float(RTs_GroundTruth.size()) << ","
			<< time_Verify / float(RTs_GroundTruth.size()) << ","
			<< time_total / float(num_rtmatch) << "\n";

	spdlog::info("Registration completed.");

	// calculation accuracy
	spdlog::info("plane extrction Time_Avg: {}", time_plane_extra / float(files_path.size()));
	spdlog::info("plane cluster Time_Avg: {}", time_plane_cluster / float(files_path.size()));
	spdlog::info("Match Time_Avg: {}", time_Match / float(pair_regis_idx.size()));
	spdlog::info("Verify Time_Avg: {}", time_Verify / float(pair_regis_idx.size()));

//	float srr = float(num_rtmatch) / float(RTs_GroundTruth.size());
	spdlog::info("Successful Registration Rate(SRR): {} / {}.", num_rtmatch, RTs_GroundTruth.size());
	spdlog::info("Time: {}.", time_total / float(pair_regis_idx.size()));

	if (num_rtmatch > 0)
	{
		std::sort(R_error.begin(), R_error.end());
		std::sort(T_error.begin(), T_error.end());

		spdlog::info("Rotation error min: {}.", R_error.front());
		spdlog::info("Rotation error max: {}.", R_error.back());
		spdlog::info("Rotation error RMSE: {}.", sqrt(R_error_avg / R_error.size()));

		spdlog::info("Translation error min: {}.", T_error.front());
		spdlog::info("Translation error max: {}.", T_error.back());
		spdlog::info("Translation error RMSE: {}.", sqrt(T_error_avg / T_error.size()));

		outfile << "\nR_min,R_max,R_RMSE,T_min,T_max,T_RMSE,SSR\n";
		outfile << R_error.front() << "," << R_error.back() << "," << sqrt(R_error_avg / R_error.size()) << ","
				<< T_error.front() << "," << T_error.back() << "," << sqrt(T_error_avg / T_error.size()) << ","
				<< std::to_string(num_rtmatch) + " % " + std::to_string(RTs_GroundTruth.size()) << "\n";
	}
	else
	{
		spdlog::info("Rotation error min: NULL.");
		spdlog::info("Rotation error max: NULL.");
		spdlog::info("Rotation error RMSE: NULL.");

		spdlog::info("Translation error min: NULL.");
		spdlog::info("Translation error max: NULL.");
		spdlog::info("Translation error RMSE: NULL.");

		outfile << "\nR_min,R_max,R_RMSE,T_min,T_max,T_RMSE,SSR\n";
		outfile << "-,-,-,-,-,-,-\n";
	}
	// Close the file
	outfile.close();

	return 0;
}



