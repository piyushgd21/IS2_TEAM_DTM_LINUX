#pragma once

#include "utility.h"
template<typename _Scalar, int NX=Eigen::Dynamic, int NY=Eigen::Dynamic>
struct Functor
{
  typedef _Scalar Scalar;
  enum 
  {
    InputsAtCompileTime = NX,
    ValuesAtCompileTime = NY
  };

  typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
  typedef Eigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
  typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;

  Functor () : data_size_(ValuesAtCompileTime) {}

  Functor (int data_size) : data_size_(data_size) {}

  virtual ~Functor () {}

  int values () const { return (data_size_); }

  private:
    const int data_size_;
};
namespace  Cylinder
{
	inline double PointToLineDistanceSquare (const Eigen::Vector3d &pt, const Eigen::Vector3d &line_pt, const Eigen::Vector3d &line_direction)
	{
		return (line_direction.cross(line_pt - pt)).squaredNorm () / line_direction.squaredNorm ();
	}

	inline double PointToLineDistance (std::vector<double> &model_parameters,const Eigen::Vector3d &pt)
	{
		Eigen::Vector3d line_pt  (model_parameters[0], model_parameters[1], model_parameters[2]);
		Eigen::Vector3d line_direction (model_parameters[3], model_parameters[4], model_parameters[5]);
		return sqrt(Cylinder::PointToLineDistanceSquare (pt, line_pt, line_direction));
	}

	inline void ProjectPointToLine (const Eigen::Vector3d &pt, const Eigen::Vector3d &line_pt, const Eigen::Vector3d &line_direction,Eigen::Vector3d &projection_point)
	{
		double factor = (pt.dot(line_direction) - line_pt.dot(line_direction)) / line_direction.dot(line_direction);
		projection_point = line_pt + factor * line_direction;
	}


	inline void BoundingBoxOnProjectionPlane(std::vector<PointType> points,std::vector<double> &model_parameters,Eigen::Vector3d &box_max,Eigen::Vector3d &box_min)
	{
		Eigen::Vector3d line_pt(model_parameters[0], model_parameters[1], model_parameters[2]);
		Eigen::Vector3d line_direction(model_parameters[3], model_parameters[4], model_parameters[5]);
		double max_factor = -FLT_MAX;
		double min_factor = FLT_MAX;
		for (int ii = 0; ii < points.size(); ++ii)
		{
			Eigen::Vector3d pt(points[ii].x, points[ii].y, points[ii].z);
			double factor = (pt.dot(line_direction) - line_pt.dot(line_direction)) / line_direction.dot(line_direction);
			if (factor > max_factor)max_factor = factor;
			else if (factor < min_factor)min_factor = factor;
		}
		box_max = line_pt + max_factor * line_direction;
		box_min = line_pt + min_factor * line_direction;
	}
};

class CylinderEstimator
{
public:
	CylinderEstimator(){}
	void Estimate(std::vector<PointType> &data, std::vector<double> &parameters, double& error);

public:
	std::vector<PointType> point_cloud_;
		
};

struct OptimizationFunctor : Functor<double>
{
	OptimizationFunctor (int data_points, CylinderEstimator *model):
    Functor<double> (data_points), model_ (model) {}
	int operator() (const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const
	{
		Eigen::Vector3d line_pt(x[0], x[1], x[2]);
		Eigen::Vector3d line_normal(x[3], x[4], x[5]);

		for (int i = 0; i < values(); ++i)
		{
			Eigen::Vector3d pt(model_->point_cloud_[i].x,
				model_->point_cloud_[i].y,
				model_->point_cloud_[i].z);

			fvec[i] = Cylinder::PointToLineDistanceSquare(pt, line_pt, line_normal) - x[6] * x[6];
		}
		return (0);
	}

	CylinderEstimator *model_;
};
