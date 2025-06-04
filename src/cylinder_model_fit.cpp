//#pragma once
#include "../header/cylinder_model_fit.h"
#include "../header/nanoflann.hpp"
#include "../header/cylinder_estimator.h"

using namespace nanoflann;

namespace Forestry_SLAM
{
	
namespace{

	/**
 * @brief KD-tree adaptor struct for nanoflann to work on std::vector<PointType>.
 */

struct NanoFlannPointCloud
{
	const std::vector<PointType> *cloud_ptr;

	/// Returns the number of data points
	inline size_t kdtree_get_point_count() const { return cloud_ptr->size(); }

	/// Computes L2 squared distance between point and data index
	inline float kdtree_distance(const float *p1, const size_t idx_p2,size_t size) const
	{
		float d0 = p1[0] - cloud_ptr->at(idx_p2).x;
		float d1 = p1[1] - cloud_ptr->at(idx_p2).y;
		float d2 = p1[2] - cloud_ptr->at(idx_p2).z;
		return d0*d0+d1*d1+d2*d2;
	}

	/// Returns the dim'th coordinate of the idx'th point
	inline float kdtree_get_pt(const size_t idx, int dim) const
	{
		return cloud_ptr->at(idx).getVector3fMap()(dim);
	}

	/// Optional bounding-box computation (not needed here)
	template <class BBOX>
	bool kdtree_get_bbox(BBOX &bb) const { return false; }

};

/// Alias for 3D KD-tree using nanoflann
typedef nanoflann::KDTreeSingleIndexAdaptor<
	nanoflann::L2_Simple_Adaptor<float, NanoFlannPointCloud > ,
	NanoFlannPointCloud,
	3 /* dimension */,
	int
> NanoFlannKdTree;

/**
 * @brief Computes the normal vector of a set of 3D points using PCA.
 *
 * @param[in] data Vector of 3D points.
 * @param[out] parameters Output vector containing the computed normal [nx, ny, nz].
 * @return true if the normal is successfully computed, false otherwise (e.g., less than 3 points).
 */

bool ComputePointsNormal(const std::vector<PointType>& data, std::vector<double> &parameters)
{
	if (data.size() < 3)
	{
		return false;
	}
	Eigen::Matrix3d covariance_matrix = Eigen::Matrix3d::Zero();
	Eigen::Matrix<double, 1, 9, Eigen::RowMajor> inter_matrix = Eigen::Matrix<double, 1, 9, Eigen::RowMajor>::Zero();
	size_t point_count = data.size();
	for (int i = 0; i < point_count; ++i)
	{
		inter_matrix[0] += data[i].x * data[i].x;
		inter_matrix[1] += data[i].x * data[i].y;
		inter_matrix[2] += data[i].x * data[i].z;
		inter_matrix[3] += data[i].y * data[i].y;
		inter_matrix[4] += data[i].y * data[i].z;
		inter_matrix[5] += data[i].z * data[i].z;
		inter_matrix[6] += data[i].x;
		inter_matrix[7] += data[i].y;
		inter_matrix[8] += data[i].z;
	}

	inter_matrix /= static_cast<double> (point_count);
	covariance_matrix.coeffRef(0) = inter_matrix[0] - inter_matrix[6] * inter_matrix[6];
	covariance_matrix.coeffRef(1) = inter_matrix[1] - inter_matrix[6] * inter_matrix[7];
	covariance_matrix.coeffRef(2) = inter_matrix[2] - inter_matrix[6] * inter_matrix[8];
	covariance_matrix.coeffRef(4) = inter_matrix[3] - inter_matrix[7] * inter_matrix[7];
	covariance_matrix.coeffRef(5) = inter_matrix[4] - inter_matrix[7] * inter_matrix[8];
	covariance_matrix.coeffRef(8) = inter_matrix[5] - inter_matrix[8] * inter_matrix[8];
	covariance_matrix.coeffRef(3) = covariance_matrix.coeff(1);
	covariance_matrix.coeffRef(6) = covariance_matrix.coeff(2);
	covariance_matrix.coeffRef(7) = covariance_matrix.coeff(5);

	// Eigen decomposition
	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::Matrix3d> solver(covariance_matrix, Eigen::Matrix3d::Identity(), Eigen::ComputeEigenvectors | Eigen::Ax_lBx);

	parameters.clear();
	parameters.resize(3);
	parameters[0] = solver.eigenvectors().col(0)[0];
	parameters[1] = solver.eigenvectors().col(0)[1];
	parameters[2] = solver.eigenvectors().col(0)[2];
	return true;
}

/**
 * @brief Compute surface normal for a specific point given its neighbor indices.
 */
bool ComputePointNormal (const std::vector<PointType>& point_cloud, size_t points_number,int*indices, float &nx, float &ny, float  &nz)
{
	if (points_number < 3) return false;
	std::vector<PointType> points;
	points.reserve(points_number);
	for (size_t i = 0; i < points_number; ++i)
	{
		points.push_back(point_cloud[indices[i]]);
	}
	std::vector<double> parameters;
	ComputePointsNormal(points, parameters);
	nx = parameters[0];
	ny = parameters[1];
	nz = parameters[2];
	return true;
}

/**
 * @brief Estimate 3D plane parameters (normal vector) from point cloud using PCA.
 */

void Estimate3DPlaneParameters(std::vector<PointType> &data, std::vector<double> &parameters){
	 ComputePointsNormal(data, parameters);
}

}

/**
 * @brief Fit a cylinder model to a 3D point cloud using normal estimation and nonlinear optimization.
 *
 * @param[in] point_cloud The raw point cloud.
 * @param[out] cylinder_parameters Estimated parameters [x, y, z, dx, dy, dz, radius].
 * @param[out] final_error RMS error of the cylinder fitting.
 * @return 1 if successful, negative codes for failure modes:
 *         -1: Not enough points, -2: Normal estimation failed, -3: Center shift too large.
 */

int FittingCylinderModel(/*const*/ std::vector<PointType> &point_cloud, std::vector<double>& cylinder_parameters,double & final_error)
{
	if (point_cloud.size() < 3) return -1;
	int neighbor_number = 80;
	int point_number = point_cloud.size();
	if (point_number < neighbor_number) neighbor_number = point_number;

	// Build KD-tree for normal estimation
	NanoFlannPointCloud cloud;
	cloud.cloud_ptr = &point_cloud;
	NanoFlannKdTree kdtree(3, cloud, KDTreeSingleIndexAdaptorParams(10));
	kdtree.buildIndex();

	// Compute centroid
	double center_x = 0, center_y = 0, center_z = 0;
	for (int ii = 0; ii < point_number; ++ii) {
		center_x += point_cloud[ii].x;
		center_y += point_cloud[ii].y;
		center_z += point_cloud[ii].z;
	}
	center_x /= point_number;
	center_y /= point_number;
	center_z /= point_number;

	// Estimate surface normals at each point
	int *pnIndex = NULL;

	int *indices = new int[neighbor_number];
	float *point_distance = new float[neighbor_number];

	double nx, ny, nz;
	std::vector<PointType> point_normal_vector(point_number);

	for (int ii = 0; ii < point_number; ++ii) {
		float query_pt[3] = {
			point_cloud[ii].x,
			point_cloud[ii].y,
			point_cloud[ii].z
		};
		int result_number = kdtree.knnSearch(query_pt, neighbor_number, indices, point_distance);
		if (result_number > 0)
		{
			ComputePointNormal(point_cloud, result_number, indices, point_normal_vector[ii].x, 
				point_normal_vector[ii].y, point_normal_vector[ii].z);
		}

	}
	if (indices)delete[]indices;
	if (point_distance) delete[]point_distance;

	// Estimate cylinder axis direction from normals
	std::vector<double> parameters;
	Estimate3DPlaneParameters(point_normal_vector, parameters);
	if (parameters.size() == 0)
		return -2;
	
	// Build initial model
	std::vector<double> model_parameters;
	CylinderEstimator cylinder_estimator;
	std::vector<PointType> points;
	for (int i = 0; i < point_number; ++i)
		points.push_back(point_cloud[i]);
	model_parameters.resize(7);
	model_parameters[0] = center_x;
	model_parameters[1] = center_y;
	model_parameters[2] = center_z;
	model_parameters[3] = parameters[0];
	model_parameters[4] = parameters[1];
	model_parameters[5] = parameters[2];

	// Estimate initial radius from a sample point
	Eigen::Vector3d point(point_cloud[0].x, point_cloud[0].y, point_cloud[0].z);
	model_parameters[6] = Cylinder::PointToLineDistance(model_parameters, point);

	// Perform nonlinear cylinder fitting
	cylinder_estimator.Estimate(points, model_parameters, final_error);

	Eigen::Vector3d pt_max, pt_min;

	// Compute cylinder's bounding box projection
	Cylinder::BoundingBoxOnProjectionPlane(points, model_parameters, pt_max, pt_min);

	// Final cylinder parameters from fitted result and bounding box center
	cylinder_parameters.resize(7);
	cylinder_parameters[0] = (pt_max(0) + pt_min(0))*0.5;
	cylinder_parameters[1] = (pt_max(1) + pt_min(1))*0.5;
	cylinder_parameters[2] = (pt_max(2) + pt_min(2))*0.5;
	cylinder_parameters[3] = model_parameters[3];
	cylinder_parameters[4] = model_parameters[4];
	cylinder_parameters[5] = model_parameters[5];
	cylinder_parameters[6] = model_parameters[6];

	// Sanity check: if cylinder center shifted too far from point cloud center, it's likely bad fit
	double difference_x = cylinder_parameters[0] - center_x;
	double difference_y = cylinder_parameters[1] - center_y;

	if (std::sqrt(difference_x*difference_x + difference_y*difference_y) > 8 * cylinder_parameters[6])
	{
		return -3;
	}
	return 1;
}

}