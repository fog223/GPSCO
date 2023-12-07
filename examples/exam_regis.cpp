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

int main()
{
//	std::string file_src = "F:\\Benchmark\\HS_1\\1-RawPointCloud\\Scan_545_33w.ply";
//	std::string file_tgt = "F:\\Benchmark\\HS_1\\1-RawPointCloud\\Scan_546_50w.ply";
	std::string file_src = "F:\\Data\\Data\\Target_Ball\\Data\\Cloud_852_3DHoPD.ply";
	std::string file_tgt = "F:\\Data\\Data\\Target_Ball\\Data\\Cloud_853_3DHoPD.ply";
	GPSCO::cloudptr cloud_src(new GPSCO::cloud);
	GPSCO::cloudptr cloud_tgt(new GPSCO::cloud);

	if (pcl::io::loadPLYFile<pcl::PointXYZ>(file_src, *cloud_src) == -1 ||
		pcl::io::loadPLYFile<pcl::PointXYZ>(file_tgt, *cloud_tgt) == -1)
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
	params.min_support_points = 1000;
	params.SmoothnessThreshold = 2.0;
	params.CurvatureThreshold = 1.0;
	params.parallel_thresh = 5.0;
	params.coplanar_thresh = 2.0;
	params.dist_thresh = 0.5;
	if (GPSCO::Registration::Regis(cloud_src, cloud_tgt, params, RT))
	{
		spdlog::info("Point cloud registration is successful!!! The transformation matrix is");
		return 0;
	}
	else
	{
		spdlog::error("Point cloud registration is failure!!!");
		return 0;
	}
}




