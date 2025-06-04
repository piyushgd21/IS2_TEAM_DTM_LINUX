#include "../header/cylinder_estimator.h"
#include <unsupported/Eigen/NonLinearOptimization>

using namespace Eigen;

/**
 * @brief Estimates cylinder parameters from a given 3D point cloud using Levenberg-Marquardt optimization.
 *
 * The model parameters being estimated are:
 * - [0,1,2] : A point on the cylinder axis (x, y, z)
 * - [3,4,5] : Direction vector of the axis (dx, dy, dz) — normalized
 * - [6]     : Radius of the cylinder (r)
 *
 * This method performs non-linear least squares optimization on the input point cloud.
 * It uses Eigen's Levenberg-Marquardt optimizer with numerical differentiation.
 *
 * @param[in] data Input vector of 3D points (PointType).
 * @param[in,out] parameters Initial estimate of the cylinder model (7 values). Will be updated.
 * @param[out] error Final root mean square error of the cylinder fit.
 */

void CylinderEstimator::Estimate(std::vector<PointType> &data, std::vector<double> &parameters,double& error){
	// Cylinder fitting requires at least 4 points
	if (data.size()< 4)
         return;
	
	// Wrap initial guess into Eigen's vector format
	Eigen::VectorXd model_parameters(7);
	model_parameters<<parameters[0],parameters[1],parameters[2], parameters[3],parameters[4],parameters[5],parameters[6];

	// Store input point cloud for internal optimization use
	point_cloud_ = data;
	int data_size = data.size();

	// Construct optimization functor
	OptimizationFunctor functor (data_size, this);
	Eigen::NumericalDiff<OptimizationFunctor > numerical_diff (functor);
	Eigen::LevenbergMarquardt<Eigen::NumericalDiff<OptimizationFunctor>, double> levenberg_marquardt(numerical_diff);

	// Run minimization
	int info = levenberg_marquardt.minimize (model_parameters);

	// Compute root mean square error of the fitting
	error = sqrt(levenberg_marquardt.fvec.squaredNorm()/static_cast<double> (data.size ()-3));

	// Normalize axis direction vector
	Eigen::Vector3d normal (model_parameters[3], model_parameters[4], model_parameters[5]);
	normal.normalize ();

	// Copy optimized values back to parameter vector
	parameters[0] = model_parameters[0];
	parameters[1] = model_parameters[1];
	parameters[2] = model_parameters[2];

	parameters[3] = normal[0];
	parameters[4] = normal[1];
	parameters[5] = normal[2];
	parameters[6] = abs(model_parameters[6]); // Ensure positive radius
}

