/**
 *  \file   exam_regis.cpp
 *  \brief  Registration between two stations
 *  \author fog
 *  \email  luochengwen22@mails.ucas.ac.cn
 *  \date   2023/12/4
 *  \note
 */
//

// spdlog
#include <spdlog/spdlog.h>
// pcl
#include <pcl/io/pcd_io.h>
#include <pcl/io/ply_io.h>
// local
#include "registration.h"

int main(int argc, char** argv)
{
	if (argc != 4)
	{
		spdlog::error("Parameter input error!");
		return -1;
	}

	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_src(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_tgt(new pcl::PointCloud<pcl::PointXYZ>);

	if (pcl::io::loadPLYFile<pcl::PointXYZ>("D:\\Benchmark_HS\\HS_2\\1-RawPointCloud-Sampled\\8_sampled.ply", *cloud_src) == -1 ||
		pcl::io::loadPLYFile<pcl::PointXYZ>("D:\\Benchmark_HS\\HS_2\\1-RawPointCloud-Sampled\\13_sampled.ply", *cloud_tgt) == -1)
	{
		spdlog::error("Tata load failure!");
		return 0;
	}
	else
	{
		spdlog::info("data load success! the num of the source point cloud points is {}; "
					 "the num of the source point cloud points is {}", cloud_src->size(), cloud_tgt->size());
	}

	float density_src, density_tgt;
	GPSCO::Compute_density(cloud_src, density_src);
	GPSCO::Compute_density(cloud_tgt, density_tgt);

	GPSCO::Registration::Options options;
	options.IsRepeat = true;
	options.min_support_points = 800;
	options.max_plane_num = 0;
	options.SmoothnessThreshold = 2.0;
	options.CurvatureThreshold = 1.0;
	options.SegSize = 0.2;
	options.parallel_thresh = 5.0;
	options.coplanar_thresh = 5.0;
	options.e_pl2pldist = 0.05;
	options.min_pl2pldist = 0.3;
	options.max_dist_inspect = 1.0;
	options.overlap_dist = 0.05;
	options.max_dist_evaluate = 0.3;

	GPSCO::Registration regis_(options);
	regis_.SetCloud(cloud_src, cloud_tgt);
	if (regis_.Regis())
	{
		std::ofstream outfile(argv[3]);
		outfile << regis_.GetRT() << std::endl;
		outfile.close();

		spdlog::info("Registration Successful!");
		std::cout << regis_.GetRT() << std::endl;
		spdlog::info("Completed.");
	}
	else
	{
		spdlog::info("Registration failure!");
		spdlog::info("Completed.");
	}

	spdlog::info("plane extrction Time: {}", regis_.time_plane_extra);
	spdlog::info("plane cluster Time: {}", regis_.time_plane_cluster);
	spdlog::info("Match Time: {}", regis_.time_Match);
	spdlog::info("Verify Time: {}", regis_.time_Verify);
	spdlog::info("Total Time: {}", regis_.time_Total);

	return 0;
}




