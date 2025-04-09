#include "../header/cylinder_estimator.h"
#include <unsupported/Eigen/NonLinearOptimization>

using namespace Eigen;

void CylinderEstimator::Estimate(std::vector<PointType> &data, std::vector<double> &parameters,double& error){
	if (data.size()< 4)
         return;
	Eigen::VectorXd model_parameters(7);
	model_parameters<<parameters[0],parameters[1],parameters[2], parameters[3],parameters[4],parameters[5],parameters[6];
	point_cloud_ = data;
	int data_size = data.size();
	OptimizationFunctor functor (data_size, this);
	Eigen::NumericalDiff<OptimizationFunctor > numerical_diff (functor);
	Eigen::LevenbergMarquardt<Eigen::NumericalDiff<OptimizationFunctor>, double> levenberg_marquardt(numerical_diff);
	int info = levenberg_marquardt.minimize (model_parameters);

	error = sqrt(levenberg_marquardt.fvec.squaredNorm()/static_cast<double> (data.size ()-3));
	Eigen::Vector3d normal (model_parameters[3], model_parameters[4], model_parameters[5]);
	normal.normalize ();
	parameters[0] = model_parameters[0];
	parameters[1] = model_parameters[1];
	parameters[2] = model_parameters[2];

	parameters[3] = normal[0];
	parameters[4] = normal[1];
	parameters[5] = normal[2];
	parameters[6] = abs(model_parameters[6]);
}

