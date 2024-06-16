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

	if (pcl::io::loadPLYFile<pcl::PointXYZ>(argv[1], *cloud_src) == -1 ||
		pcl::io::loadPLYFile<pcl::PointXYZ>(argv[2], *cloud_tgt) == -1)
	{
		spdlog::error("Tata load failure!");
		return 0;
	}
	else
	{
		spdlog::info("data load success! the num of the source point cloud points is {}; "
					 "the num of the source point cloud points is {}", cloud_src->size(), cloud_tgt->size());
	}

	// If the point cloud quantity is too large (>1e6), downsampling is necessary.
	// Typically, several hundred thousand points per scan are sufficient.

	GPSCO::Registration::Options options;
	options.scale = 0.25; // smaller values(Indoor,such as 0.01); larger values(outdoor, such as 0.25)
	options.wait_time = 50;
	options.min_support_points = 500;
	options.max_plane_num = 30;
	options.SmoothnessThreshold = 5.0;
	options.CurvatureThreshold = 2.0;
	options.parallel_thresh = 5.0;
	options.coplanar_thresh = 2.0;
	options.e_pl2pldist = 0.05;

	GPSCO::Registration regis_(options);
	regis_.SetCloud(cloud_src, cloud_tgt);
	if (regis_.Regis())
	{
		std::ofstream outfile(argv[3]);
		outfile << regis_.GetRT() << std::endl;
		outfile.close();

		spdlog::info("Registration completed!");
		std::cout << regis_.GetRT() << std::endl;
	}
	else
		spdlog::info("Registration failure!");

	spdlog::info("plane extrction Time: {}", regis_.time_plane_extra);
	spdlog::info("plane cluster Time: {}", regis_.time_plane_cluster);
	spdlog::info("Match Time: {}", regis_.time_Match);
	spdlog::info("Verify Time: {}", regis_.time_Verify);
	spdlog::info("Total Time: {}", regis_.time_Total);

	return 0;
}

