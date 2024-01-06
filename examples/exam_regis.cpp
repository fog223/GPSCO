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

	Eigen::Matrix4f RT;
	GPSCO::Params params;
	params.min_support_points = 900;
	params.SmoothnessThreshold = 2.0;
	params.CurvatureThreshold = 1.0;
	params.parallel_thresh = 5.0;
	params.coplanar_thresh = 2.0;
	params.dist_thresh = 0.05;
	if (GPSCO::Registration::Regis(cloud_src, cloud_tgt, params, RT))
	{
		std::ofstream outfile(argv[3]);
		outfile << RT << std::endl;
		outfile.close();

		spdlog::info("Registration complete.");
	}

	return 0;
}




