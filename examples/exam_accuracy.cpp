/**
 *  \file   exam_accuracy.cpp
 *  \brief  Evaluation of the accuracy of the algorithm GPSCO
 *  \author fog
 *  \email  luochengwen22@mails.ucas.ac.cn
 *  \date   2023/12/12
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
bool numericStringCompare(const std::string& a, const std::string& b) {
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
	if (argc != 2 && argc != 3)
	{
		spdlog::error("Parameter input error!");
		return -1;
	}

	std::string root_dir = argv[1];

	GPSCO::Params params;
	params.min_support_points = 1000;
	params.SmoothnessThreshold = 2.0;
	params.CurvatureThreshold = 1.0;
	params.parallel_thresh = 5.0;
	params.coplanar_thresh = 2.0;
	params.dist_thresh = 0.05;

	std::vector<std::string> files_path;
	for (const auto& entry : fs::directory_iterator(root_dir + "/1-RawPointCloud-Sampled"))
	{
		const fs::path& filePath = entry.path();
		if (fs::is_regular_file(filePath) && filePath.extension() == ".ply")
		{
			files_path.push_back(filePath.string());
		}
	}

	// Sort by file name
	std::sort(files_path.begin(), files_path.end(), numericStringCompare);

	// starting moment
	auto time_begin = clock();
	float time_plane_extra = 0.0f;
	float time_plane_cluster = 0.0f;
	float time_Match = 0.0f;
	float time_Verify = 0.0f;

	std::vector<GPSCO::cloudptr> cloud_all;                                  // All point cloud data
	std::vector<std::vector<GPSCO::PLANE>> planes_all;                       // Plane extracted from point cloud
	std::vector<std::vector<std::vector<std::vector<int>>>> PlaneGroups_all; // plane groups
	std::vector<std::vector<GPSCO::Group_Three>> group_three_vector_all;     // Candidates

	int process = 0;
	for (const auto& file : files_path)
	{
		process++;
		spdlog::info("Data processing: {}/{}", process, files_path.size());

		// Load point cloud data
		GPSCO::cloudptr cloud_(new GPSCO::cloud);
		if (pcl::io::loadPLYFile<pcl::PointXYZ>(file, *cloud_) == -1)
		{
			spdlog::error("{} load failure!", file);
			return 0;
		}
		else
		{
			cloud_all.push_back(cloud_);

			// planar extraction
			auto start = clock();
			std::vector<GPSCO::PLANE> outPlanes;
			if (GPSCO::PLANE_Extraction::PLANE_Tetect_RegionGrow(
				cloud_, params.min_support_points, params.SmoothnessThreshold, params.CurvatureThreshold, outPlanes))
			{
				auto end = clock();
				auto Time = (double)(end - start) / CLOCKS_PER_SEC;
				time_plane_extra += Time;

				// large to small
				std::sort(outPlanes.begin(), outPlanes.end());
				planes_all.push_back(outPlanes);

				// planar clustering
				start = clock();
				std::vector<std::vector<std::vector<int>>> PlaneGroups;
				if (GPSCO::Registration::Plane_Cluster(
					outPlanes, params.parallel_thresh, params.coplanar_thresh, PlaneGroups))
				{
					end = clock();
					Time = (double)(end - start) / CLOCKS_PER_SEC;
					time_plane_cluster += Time;

					PlaneGroups_all.push_back(PlaneGroups);

					// Find Candidate
					float threshold_min = 30.0;
					float threshold_max = 150.0;
					std::vector<GPSCO::Group_Three> group_three_vector;
					if (GPSCO::Registration::Get_Group_Three(
						outPlanes, PlaneGroups, threshold_min, threshold_max, group_three_vector))
					{
						group_three_vector_all.push_back(group_three_vector);
					}
					else
					{
						spdlog::error("Get_Group_Three failure: {}", file);
						return 0;
					}
				}
				else
				{
					spdlog::error("{} plane clustering failure!", file);
					return 0;
				}
			}
			else
			{
				spdlog::error("{} plane extraction failure!", file);
				return 0;
			}
		}
	}

	spdlog::info("Data processing completed.");

	// load GroundTruth
	std::vector<std::string> files_transmatrix;
	for (const auto& entry : fs::directory_iterator(root_dir + "/3-GroundTruth"))
	{
		const fs::path& filePath = entry.path();
		if (fs::is_regular_file(filePath) && filePath.extension() == ".txt")
		{
			files_transmatrix.push_back(filePath.string());
		}
	}

	// Determine which two stations need to be aligned with each other.
	std::vector<Eigen::Matrix4f> RTs_GroundTruth;
	std::vector<std::pair<int, int>> pair_regis_idx;
	for (const auto& path : files_transmatrix)
	{
		std::ifstream file(path);
		if (!file.is_open())
		{
			spdlog::error("{} Can't Load File.", path);
			return 0;
		}

		// Determine the indexes of the two stations to be registered.
		size_t start = path.find("transformation") + 14; // Find the position after "transformation"
		size_t end1 = path.find("_", start); // find the position of the next "_"
		size_t end2 = path.find(".", end1); // find the position of the next "."
		std::string extractedNumber1 = path.substr(start, end1 - start);
		std::string extractedNumber2 = path.substr(end1 + 1, end2 - end1 - 1);
		int idx1 = std::stoi(extractedNumber1); // Convert the extracted string to an integer.
		int idx2 = std::stoi(extractedNumber2); // Convert the extracted string to an integer.
		pair_regis_idx.push_back(std::pair<int, int>(idx1, idx2));

		Eigen::Matrix4f rt_temp;
		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 4; ++j)
			{
				file >> rt_temp(i, j);
			}
		}
		RTs_GroundTruth.push_back(rt_temp);
	}

	spdlog::info("GroundTruth load completed.");

	// Registration
	std::vector<Eigen::Matrix4f> RTs;
	for (const auto& pair : pair_regis_idx)
	{
		// Initial match based moving alignment
		auto start = clock();
		std::vector<std::vector<std::vector<GPSCO::MatchGroup>>> group_table; // The relationship table between plane groups
		if (GPSCO::Registration::Compute_Group_Table(planes_all[pair.first - 1], planes_all[pair.second - 1],
			PlaneGroups_all[pair.first - 1], PlaneGroups_all[pair.second - 1], params.dist_thresh, group_table))
		{
			auto end = clock();
			auto Time = (double)(end - start) / CLOCKS_PER_SEC;
			time_Match += Time;

			// GPSCO determines the transformation matrix
			start = clock();
			Eigen::Matrix4f rt;
			if (GPSCO::Registration::Get_transformation_matrix(planes_all[pair.first - 1], planes_all[pair.second - 1],
				PlaneGroups_all[pair.first - 1], PlaneGroups_all[pair.second - 1], group_table,
				group_three_vector_all[pair.first - 1], group_three_vector_all[pair.second - 1], rt))
			{
				end = clock();
				Time = (double)(end - start) / CLOCKS_PER_SEC;
				time_Verify += Time;

				RTs.push_back(rt);

				// Export the transformation matrix
				if (argc == 3)
				{
					std::ofstream outfile;
					std::string path = std::string(argv[2]) + "/" + std::to_string(pair.first) + "_" + std::to_string(pair.second) + ".txt";
					outfile.open(path);
					outfile << rt << std::endl;
					outfile.close();
				}
			}
		}
	}

	spdlog::info("Registration completed.");

	// calculation accuracy
	spdlog::info("plane extrction Time_Avg: {}", time_plane_extra / float(files_path.size()));
	spdlog::info("plane cluster Time_Avg: {}", time_plane_cluster / float(files_path.size()));
	spdlog::info("Match Time_Avg: {}", time_Match / float(pair_regis_idx.size()));
	spdlog::info("Verify Time_Avg: {}", time_Verify / float(pair_regis_idx.size()));

	// SRR, R_rmse, T_rmse
	// Rotation threshold is 5°, translation threshold is 1m
	float R_thresh = 5;
	float T_thresh = 1.0;
	int num_rtmatch = 0; // correct match
	std::vector<float> R_error;
	float R_error_avg = 0.0f;
	std::vector<float> T_error;
	float T_error_avg = 0.0f;
	for (int i = 0; i < RTs.size(); i++)
	{
		if (RTs.size() != RTs_GroundTruth.size())
		{
			spdlog::error("The num of RTs unmatching the num of RTs_GroundTruth!");
			return 0;
		}
		else
		{
			Eigen::Matrix3f R1 = RTs_GroundTruth[i].block<3, 3>(0, 0);
			Eigen::Vector3f T1 = RTs_GroundTruth[i].block<3, 1>(0, 3);
			Eigen::Matrix3f R2 = RTs[i].block<3, 3>(0, 0);
			Eigen::Vector3f T2 = RTs[i].block<3, 1>(0, 3);

			Eigen::AngleAxisf aaDiff(R1.transpose() * R2);
			float angleDiff = aaDiff.angle();

			// Convert radians to angles
			float angleDiffDegrees = angleDiff * 180.0 / M_PI;

			if (angleDiffDegrees < R_thresh)
			{
				// Calculate the distance difference
				if ((T1 - T2).norm() < T_thresh)
				{
					num_rtmatch++;
					R_error.push_back(angleDiffDegrees);
					R_error_avg += pow(angleDiffDegrees, 2);
					T_error.push_back((T1 - T2).norm());
					T_error_avg += pow((T1 - T2).norm(), 2);
				}
				else
					spdlog::warn("{} and {} Registration failed!", pair_regis_idx[i].first, pair_regis_idx[i].second);
			}
			else
				spdlog::warn("{} and {} Registration failed!", pair_regis_idx[i].first, pair_regis_idx[i].second);
		}
	}

	float srr = float(num_rtmatch) / float(RTs_GroundTruth.size());
	spdlog::info("Successful Registration Rate(SRR): {}.", srr);

	auto time_end = clock();
	auto Time = (double)(time_end - time_begin) / CLOCKS_PER_SEC;
	Time += (2 * pair_regis_idx.size() - files_path.size()) * ((time_plane_extra + time_plane_cluster) / float(files_path.size()));
	Time /= pair_regis_idx.size();
	spdlog::info("Time: {}.", Time);

	if (srr > 0)
	{
		std::sort(R_error.begin(), R_error.end());
		std::sort(T_error.begin(), T_error.end());

		spdlog::info("Rotation error min: {}.", R_error.front());
		spdlog::info("Rotation error max: {}.", R_error.back());
		spdlog::info("Rotation error RMSE: {}.", sqrt(R_error_avg / R_error.size()));

		spdlog::info("Translation error min: {}.", T_error.front());
		spdlog::info("Translation error max: {}.", T_error.back());
		spdlog::info("Translation error RMSE: {}.", sqrt(T_error_avg / T_error.size()));
	}
	else
	{
		spdlog::info("Rotation error min: NULL.");
		spdlog::info("Rotation error max: NULL.");
		spdlog::info("Rotation error RMSE: NULL.");

		spdlog::info("Translation error min: NULL.");
		spdlog::info("Translation error max: NULL.");
		spdlog::info("Translation error RMSE: NULL.");
	}

	return 0;
}




