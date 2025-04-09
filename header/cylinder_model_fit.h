#pragma once

#include "../header/utility.h"

namespace Forestry_SLAM
{

int FittingCylinderModel(/*const*/ std::vector<PointType> &point_cloud, std::vector<double>& cylinder_parameters, double & final_error);

}