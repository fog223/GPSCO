/**
 *  \file   plane_extraction.h
 *  \brief  Planar extraction based on region growing algorithm
 *  and saving each plane separately using unused colours
 *  \author fog
 *  \email  luochengwen22@mails.ucas.ac.cn
 *  \date   2023/12/4
 *  \note
 */
//

#pragma once

// system
#include <string>
#include <vector>
#include <algorithm>
// local
#include "common.h"

namespace GPSCO
{
	namespace PLANE_Extraction
	{
		/// Planar extraction based on region growing algorithm
		/// \param cloud original point cloud
		/// \param min_support_points
		/// \param max_plane_num
		/// \param SmoothnessThreshold
		/// \param CurvatureThreshold
		/// \param SegSize
		/// \param outPlanes
		/// \return Whether the extraction is successful or not
		bool PLANE_Tetect_RegionGrow(
			cloudptr cloud,
			int min_support_points,
			int max_plane_num,
			float SmoothnessThreshold,
			float CurvatureThreshold,
			float SegSize,
			std::vector<GPSCO::PLANE>& outPlanes);

		/// Save the planes separately, with colour
		/// \param Planes Plane to be exported
		/// \param outpath File Save Path
		/// \return .txt
		bool Export_plane(
			std::vector<GPSCO::PLANE>& Planes,
			std::string outpath);
	}
}