/**
 *  \file   exam_plane_extraction.cpp
 *  \brief  Testing and validation of planar extraction algorithm
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
#include "plane_extraction.h"

int main()
{
	std::string file = "D:\\Benchmark_HS\\HS_1\\1-RawPointCloud\\Scan_545_33w.ply";
	GPSCO::cloudptr cloud(new GPSCO::cloud);

	if (pcl::io::loadPLYFile<pcl::PointXYZ>(file, *cloud) == -1)
	{
		spdlog::error("{} load failure!");
		return 0;
	}
	else
	{
		spdlog::info("{} load success! the num of the points is {}", file, cloud->size());
	}

	int min_support_points = 1000;
	float SmoothnessThreshold = 2.0;
	float CurvatureThreshold = 1.0;
	std::vector<GPSCO::PLANE> outPlanes;
	if(GPSCO::PLANE_Extraction::PLANE_Tetect_RegionGrow(
		cloud, min_support_points, SmoothnessThreshold, CurvatureThreshold, outPlanes))
	{
		GPSCO::PLANE_Extraction::Export_plane(outPlanes,"D:\\Code\\CLion\\GPSCO\\results\\");
	}

	return 0;
}




