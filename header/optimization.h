#pragma once
#include "../header/utility.h"
#include <ceres/ceres.h>
#include <ceres/rotation.h>
// #include <eigen3/Eigen/Dense>
// #include <pcl-1.10/pcl/point_cloud.h>
// #include <pcl-1.10/pcl/point_types.h>
// #include <pcl/kdtree/kdtree_flann.h>
// #include <pcl_conversions/pcl_conversions.h>
#include <ceres/jet.h>
//namespace here {
//	template <typename T, int N>
//	inline bool signbit(const ceres::Jet<T, N>& f) {
//		return std::signbit(f.a);
//	}
//}
//
template <typename T>
T ChangeResidual(T prev_residual)
{
	T current_residual;
	if (!ceres::signbit(prev_residual))
	{
		current_residual = prev_residual;
	}
	else
	{
		current_residual = prev_residual*prev_residual*prev_residual;
	}
	return current_residual;
}

template <typename T>
T IncreaseResidualIfDistanceNegative(T prev_residual)
{
	T current_residual;
	T increase_std= T(0.05);
	if (prev_residual>T(0.0))
	{
		current_residual = prev_residual/increase_std;
	}
	else
	{
		current_residual = prev_residual;
	}
	return current_residual;
}


// Point to cylinder residual: T of a lidar point is a function of two epochs
struct CylinderFactorTwoEpoch
{
	Eigen::Vector3d curr_point;
	double s;
	double std;

	CylinderFactorTwoEpoch(Eigen::Vector3d curr_point_, double s_, double std_)
		: curr_point{curr_point_}, s(s_), std(std_) {}

	template <typename T>
	bool operator()(const T *t1, const T *q1, const T *t2, const T *q2, const T *cylinder, T *residual) const
	{

		Eigen::Matrix<T, 3, 1> cp{T(curr_point.x()), T(curr_point.y()), T(curr_point.z())};

		// point coordinates in mapiing
		Eigen::Quaternion<T> q_prev{q1[3], q1[0], q1[1], q1[2]};
		Eigen::Quaternion<T> q_cur_end{q2[3], q2[0], q2[1], q2[2]};
		Eigen::Quaternion<T> q_cur; // intepolated

		T Ts = T(s);
		T T_1_s = T(1 - s);
		q_cur = q_prev.slerp(Ts, q_cur_end);
		Eigen::Matrix<T, 3, 1> t_curr{T_1_s * t1[0] + Ts * t2[0], T_1_s * t1[1] + Ts * t2[1], T_1_s * t1[2] + Ts * t2[2]};

		Eigen::Matrix<T, 3, 1> lp;
		lp = q_cur * cp + t_curr;

		// line representation
		Eigen::Matrix<T, 3, 1> p1{cylinder[0], cylinder[1], cylinder[2]}; // center
		Eigen::Matrix<T, 3, 1> normal{cylinder[3], cylinder[4], cylinder[5]};
		Eigen::Matrix<T, 3, 1> p2 = p1 + normal;
		T rT = cylinder[6];

		//
		T distance = ((lp - p1).cross(lp - p2)).norm() / normal.norm();
		T stdT = T(std);
		T res = (distance - rT) / stdT;
		// res = IncreaseResidualIfDistanceNegative<T>(res);
		residual[0] = res;
		return true;
	}

	static ceres::CostFunction *Create(const Eigen::Vector3d curr_point_, const double s_, const double std_)
	{
		return (new ceres::AutoDiffCostFunction<
				CylinderFactorTwoEpoch, 1, PosBlockSize, OriBlockSize, PosBlockSize, OriBlockSize, CylinderBlockSize>(
			new CylinderFactorTwoEpoch(curr_point_, s_, std_)));
	}
};

// Point to cylinder residual: T of a lidar point is a function of two epochs
struct CylinderFactorOneEpoch
{
	Eigen::Vector3d curr_point;
	double std;

	CylinderFactorOneEpoch(Eigen::Vector3d curr_point_, double std_)
		: curr_point{curr_point_}, std(std_) {}

	template <typename T>
	bool operator()(const T *t, const T *q, const T *cylinder, T *residual) const
	{

		Eigen::Matrix<T, 3, 1> cp{T(curr_point.x()), T(curr_point.y()), T(curr_point.z())};

		// point coordinates in mapiing
		Eigen::Quaternion<T> q_cur{q[3], q[0], q[1], q[2]};
		Eigen::Matrix<T, 3, 1> t_curr{t[0], t[1], t[2]};

		Eigen::Matrix<T, 3, 1> lp;
		lp = q_cur * cp + t_curr;

		// line representation
		Eigen::Matrix<T, 3, 1> p1{cylinder[0], cylinder[1], cylinder[2]}; // center
		Eigen::Matrix<T, 3, 1> normal{cylinder[3], cylinder[4], cylinder[5]};
		Eigen::Matrix<T, 3, 1> p2 = p1 + normal;
		T rT = cylinder[6];

		//
		T distance = ((lp - p1).cross(lp - p2)).norm() / normal.norm();
		T stdT = T(std);
		T test_extra_radius = T(0.0);//0.05
		T res = (distance - rT - test_extra_radius) / stdT;
		 //res = IncreaseResidualIfDistanceNegative<T>(res);
		residual[0] = res;
		//residual[0] = (distance - rT) / stdT;

		return true;
	}

	static ceres::CostFunction *Create(const Eigen::Vector3d curr_point_, const double std_)
	{
		return (new ceres::AutoDiffCostFunction<
				CylinderFactorOneEpoch, 1, PosBlockSize, OriBlockSize, CylinderBlockSize>(
			new CylinderFactorOneEpoch(curr_point_, std_)));
	}
};

// Object point  to cylinder residual: T of a lidar point is a function of two epochs
struct CylinderFactorObjectPoint
{
	Eigen::Vector3d curr_point;
	double std;

	CylinderFactorObjectPoint(Eigen::Vector3d curr_point_, double std_)
		: curr_point{curr_point_}, std(std_) {}

	template <typename T>
	bool operator()(const T *cylinder, T *residual) const
	{
		Eigen::Matrix<T, 3, 1> lp{T(curr_point.x()), T(curr_point.y()), T(curr_point.z())};

		// line representation
		Eigen::Matrix<T, 3, 1> p1{cylinder[0], cylinder[1], cylinder[2]}; // center
		Eigen::Matrix<T, 3, 1> normal{cylinder[3], cylinder[4], cylinder[5]};
		Eigen::Matrix<T, 3, 1> p2 = p1 + normal;
		T rT = cylinder[6];

		T distance = ((lp - p1).cross(lp - p2)).norm() / normal.norm();
		T stdT = T(std);
		T test_extra_radius = T(0.0);//0.05
		T res = (distance - rT - test_extra_radius) / stdT;
		 //res = IncreaseResidualIfDistanceNegative<T>(res);
		residual[0] = res;
		//residual[0] = (distance - rT) / stdT;

		return true;
	}

	static ceres::CostFunction *Create(const Eigen::Vector3d curr_point_, const double std_)
	{
		return (new ceres::AutoDiffCostFunction<
				CylinderFactorObjectPoint, 1, CylinderBlockSize>(
			new CylinderFactorObjectPoint(curr_point_, std_)));
	}
};

struct CylinderPriorConstraints
{
	Eigen::Vector3d tree_location_;
	double std_;

	CylinderPriorConstraints(Eigen::Vector3d tree_location, double std)
		: tree_location_{tree_location}, std_(std) {}

	template <typename T>
	bool operator()(const T *cylinder, T *residuals_ptr) const
	{
		Eigen::Matrix<T, 2, 1> tree_loc_2d{cylinder[0], cylinder[1]}; // center
		Eigen::Map<Eigen::Matrix<T, 2, 1>> residuals(residuals_ptr);
		residuals = (tree_loc_2d - tree_location_.head(2).template cast<T>())/T(std_);

		return true;
	}

	static ceres::CostFunction *Create(const Eigen::Vector3d tree_location, const double std)
	{
		return (new ceres::AutoDiffCostFunction<
				CylinderPriorConstraints, 2, CylinderBlockSize>(
			new CylinderPriorConstraints(tree_location, std)));
	}
};



struct CylinderPriorConstraintsPointToLine
{
	Eigen::Vector3d tree_location_;
	double std_;

	CylinderPriorConstraintsPointToLine(Eigen::Vector3d tree_location, double std)
		: tree_location_{tree_location}, std_(std) {}

	template <typename T>
	bool operator()(const T *cylinder, T *residual) const
	{

		Eigen::Matrix<T, 3, 1> cp{T(tree_location_.x()), T(tree_location_.y()), T(tree_location_.z())};

		// line representation
		Eigen::Matrix<T, 3, 1> p1{cylinder[0], cylinder[1], cylinder[2]}; // center
		Eigen::Matrix<T, 3, 1> normal{cylinder[3], cylinder[4], cylinder[5]};
		Eigen::Matrix<T, 3, 1> p2 = p1 + normal;

		T distance = ((cp - p1).cross(cp - p2)).norm() / normal.norm();
		T stdT = T(std_);

		residual[0] = distance / stdT;
		return true;
	}

	static ceres::CostFunction *Create(const Eigen::Vector3d tree_location, const double std)
	{
		return (new ceres::AutoDiffCostFunction<
				CylinderPriorConstraintsPointToLine, 1, CylinderBlockSize>(
			new CylinderPriorConstraintsPointToLine(tree_location, std)));
	}
};



// Point to plane residual: T of a lidar point is a function of two epochs
struct PlaneFactorTwoEpoch
{
	Eigen::Vector3d curr_point;
	double s;
	double std;

	PlaneFactorTwoEpoch(Eigen::Vector3d curr_point_, double s_, double std_)
		: curr_point(curr_point_), s(s_), std(std_) {}

	template <typename T>
	bool operator()(const T *t1, const T *q1, const T *t2, const T *q2, const T *plane, T *residual) const
	{

		Eigen::Matrix<T, 3, 1> cp{T(curr_point.x()), T(curr_point.y()), T(curr_point.z())};

		// point coordinates in mapiing
		Eigen::Quaternion<T> q_prev{q1[3], q1[0], q1[1], q1[2]};
		Eigen::Quaternion<T> q_cur_end{q2[3], q2[0], q2[1], q2[2]};
		Eigen::Quaternion<T> q_cur; // intepolated

		T Ts = T(s);
		T T_1_s = T(1 - s);
		q_cur = q_prev.slerp(Ts, q_cur_end);
		Eigen::Matrix<T, 3, 1> t_curr{T_1_s * t1[0] + Ts * t2[0], T_1_s * t1[1] + Ts * t2[1], T_1_s * t1[2] + Ts * t2[2]};

		Eigen::Matrix<T, 3, 1> lp;
		lp = q_cur * cp + t_curr;

		// plane representation
		Eigen::Matrix<T, 3, 1> normal{plane[0], plane[1], plane[2]};

		//
		T distance = (plane[0] * lp[0] + plane[1] * lp[1] + plane[2] * lp[2] + plane[3]) / normal.norm();
		T stdT = T(std);
		residual[0] = distance / stdT;

		return true;
	}

	static ceres::CostFunction *Create(const Eigen::Vector3d curr_point_, const double s_, const double std_)
	{
		return (new ceres::AutoDiffCostFunction<
				PlaneFactorTwoEpoch, 1, PosBlockSize, OriBlockSize, PosBlockSize, OriBlockSize, PlaneBlockSize>(
			new PlaneFactorTwoEpoch(curr_point_, s_, std_)));
	}
};

struct PlaneFactorOneEpoch
{
	Eigen::Vector3d curr_point;
	double std;

	PlaneFactorOneEpoch(Eigen::Vector3d curr_point_, double std_)
		: curr_point(curr_point_), std(std_) {}

	template <typename T>
	bool operator()(const T *t, const T *q, const T *plane, T *residual) const
	{

		Eigen::Matrix<T, 3, 1> cp{T(curr_point.x()), T(curr_point.y()), T(curr_point.z())};

		// point coordinates in mapiing
		Eigen::Quaternion<T> q_cur{q[3], q[0], q[1], q[2]};
		Eigen::Matrix<T, 3, 1> t_curr{t[0], t[1], t[2]};

		Eigen::Matrix<T, 3, 1> lp;
		lp = q_cur * cp + t_curr;

		// plane representation
		Eigen::Matrix<T, 3, 1> normal{plane[0], plane[1], plane[2]};

		//
		T distance = (plane[0] * lp[0] + plane[1] * lp[1] + plane[2] * lp[2] + plane[3]) / normal.norm();
		T stdT = T(std);
		residual[0] = distance / stdT;

		return true;
	}

	static ceres::CostFunction *Create(const Eigen::Vector3d curr_point_, const double std_)
	{
		return (new ceres::AutoDiffCostFunction<
				PlaneFactorOneEpoch, 1, PosBlockSize, OriBlockSize, PlaneBlockSize>(
			new PlaneFactorOneEpoch(curr_point_, std_)));
	}
};

struct PlaneFactorObjectPoint
{
	Eigen::Vector3d curr_point;
	double std;

	PlaneFactorObjectPoint(Eigen::Vector3d curr_point_, double std_)
		: curr_point(curr_point_), std(std_) {}

	template <typename T>
	bool operator()(const T *plane, T *residual) const
	{

		Eigen::Matrix<T, 3, 1> lp{T(curr_point.x()), T(curr_point.y()), T(curr_point.z())};

		// plane representation
		Eigen::Matrix<T, 3, 1> normal{plane[0], plane[1], plane[2]};

		//
		T distance = (plane[0] * lp[0] + plane[1] * lp[1] + plane[2] * lp[2] + plane[3]) / normal.norm();
		T stdT = T(std);
		residual[0] = distance / stdT;

		return true;
	}

	static ceres::CostFunction *Create(const Eigen::Vector3d curr_point_, const double std_)
	{
		return (new ceres::AutoDiffCostFunction<
				PlaneFactorObjectPoint, 1, PlaneBlockSize>(
			new PlaneFactorObjectPoint(curr_point_, std_)));
	}
};

// Point to a fixed plane distance, only unknowns are q and t .
struct PlaneFactorOneEpochFixPlane
{
	Eigen::Vector3d curr_point;
	Eigen::Vector4d plane_param;
	double std;
	PlaneFactorOneEpochFixPlane(Eigen::Vector3d curr_point_, Eigen::Vector4d plane_param_, double std_)
		: curr_point(curr_point_), plane_param(plane_param_), std(std_) {}

	template <typename T>
	bool operator()(const T *t, const T *q, T *residual) const
	{
		Eigen::Matrix<T, 3, 1> cp{T(curr_point.x()), T(curr_point.y()), T(curr_point.z())};

		// point coordinates in mapiing
		Eigen::Quaternion<T> q_cur{q[3], q[0], q[1], q[2]};
		Eigen::Matrix<T, 3, 1> t_curr{t[0], t[1], t[2]};

		Eigen::Matrix<T, 3, 1> lp;
		lp = q_cur * cp + t_curr;

		// plane representation
		T nx = T(plane_param(0));
		T ny = T(plane_param(1));
		T nz = T(plane_param(2));
		T nd = T(plane_param(3));

		Eigen::Matrix<T, 3, 1> normal{nx, ny, nz};

		//
		T distance = (nx * lp[0] + ny * lp[1] + nz * lp[2] + nd) / normal.norm();
		T stdT = T(std);
		residual[0] = distance / stdT;

		return true;
	}

	static ceres::CostFunction *Create(const Eigen::Vector3d curr_point_, const Eigen::Vector4d plane_param_, const double std_)
	{
		return (new ceres::AutoDiffCostFunction<
				PlaneFactorOneEpochFixPlane, 1, PosBlockSize, OriBlockSize>(
			new PlaneFactorOneEpochFixPlane(curr_point_, plane_param_, std_)));
	}
};

struct LidarEdgeFactor
{
	LidarEdgeFactor(Eigen::Vector3d curr_point_, Eigen::Vector3d last_point_a_,
					Eigen::Vector3d last_point_b_, double s_)
		: curr_point(curr_point_), last_point_a(last_point_a_), last_point_b(last_point_b_), s(s_) {}

	template <typename T>
	bool operator()(const T *q, const T *t, T *residual) const
	{

		Eigen::Matrix<T, 3, 1> cp{T(curr_point.x()), T(curr_point.y()), T(curr_point.z())};
		Eigen::Matrix<T, 3, 1> lpa{T(last_point_a.x()), T(last_point_a.y()), T(last_point_a.z())};
		Eigen::Matrix<T, 3, 1> lpb{T(last_point_b.x()), T(last_point_b.y()), T(last_point_b.z())};

		// Eigen::Quaternion<T> q_last_curr{q[3], T(s) * q[0], T(s) * q[1], T(s) * q[2]};
		Eigen::Quaternion<T> q_last_curr{q[3], q[0], q[1], q[2]};
		Eigen::Quaternion<T> q_identity{T(1), T(0), T(0), T(0)};
		q_last_curr = q_identity.slerp(T(s), q_last_curr);
		Eigen::Matrix<T, 3, 1> t_last_curr{T(s) * t[0], T(s) * t[1], T(s) * t[2]};

		Eigen::Matrix<T, 3, 1> lp;
		lp = q_last_curr * cp + t_last_curr;

		Eigen::Matrix<T, 3, 1> nu = (lp - lpa).cross(lp - lpb);
		Eigen::Matrix<T, 3, 1> de = lpa - lpb;

		residual[0] = nu.x() / de.norm();
		residual[1] = nu.y() / de.norm();
		residual[2] = nu.z() / de.norm();

		return true;
	}

	static ceres::CostFunction *Create(const Eigen::Vector3d curr_point_, const Eigen::Vector3d last_point_a_,
									   const Eigen::Vector3d last_point_b_, const double s_)
	{
		return (new ceres::AutoDiffCostFunction<
				LidarEdgeFactor, 3, 4, 3>(
			new LidarEdgeFactor(curr_point_, last_point_a_, last_point_b_, s_)));
	}

	Eigen::Vector3d curr_point, last_point_a, last_point_b;
	double s;
};

struct LidarEdgeFactorOneResidual
{
	LidarEdgeFactorOneResidual(Eigen::Vector3d curr_point_, Eigen::Vector3d last_point_a_,
							   Eigen::Vector3d last_point_b_, double s_, double std_)
		: curr_point(curr_point_), last_point_a(last_point_a_), last_point_b(last_point_b_), s(s_), std(std_) {}

	template <typename T>
	bool operator()(const T *q, const T *t, T *residual) const
	{

		Eigen::Matrix<T, 3, 1> cp{T(curr_point.x()), T(curr_point.y()), T(curr_point.z())};
		Eigen::Matrix<T, 3, 1> lpa{T(last_point_a.x()), T(last_point_a.y()), T(last_point_a.z())};
		Eigen::Matrix<T, 3, 1> lpb{T(last_point_b.x()), T(last_point_b.y()), T(last_point_b.z())};

		// Eigen::Quaternion<T> q_last_curr{q[3], T(s) * q[0], T(s) * q[1], T(s) * q[2]};
		Eigen::Quaternion<T> q_last_curr{q[3], q[0], q[1], q[2]};
		Eigen::Quaternion<T> q_identity{T(1), T(0), T(0), T(0)};
		q_last_curr = q_identity.slerp(T(s), q_last_curr);
		Eigen::Matrix<T, 3, 1> t_last_curr{T(s) * t[0], T(s) * t[1], T(s) * t[2]};

		Eigen::Matrix<T, 3, 1> lp;
		lp = q_last_curr * cp + t_last_curr;

		Eigen::Matrix<T, 3, 1> nu = (lp - lpa).cross(lp - lpb);
		Eigen::Matrix<T, 3, 1> de = lpa - lpb;

		residual[0] = nu.norm() / de.norm() / T(std);
		return true;
	}

	static ceres::CostFunction *Create(const Eigen::Vector3d curr_point_, const Eigen::Vector3d last_point_a_,
									   const Eigen::Vector3d last_point_b_, const double s_, const double std_)
	{
		return (new ceres::AutoDiffCostFunction<
				LidarEdgeFactorOneResidual, 1, 4, 3>(
			new LidarEdgeFactorOneResidual(curr_point_, last_point_a_, last_point_b_, s_, std_)));
	}

	Eigen::Vector3d curr_point, last_point_a, last_point_b;
	double s;
	double std;
};

/*Based on two normal vectors*/
struct NormalFactor
{
	NormalFactor(Eigen::Vector3d curr_normal_, Eigen::Vector3d last_normal_)
		: curr_normal(curr_normal_), last_normal(last_normal_) {}

	template <typename T>
	bool operator()(const T *q, T *residual) const
	{

		Eigen::Matrix<T, 3, 1> n1{T(curr_normal.x()), T(curr_normal.y()), T(curr_normal.z())};
		Eigen::Matrix<T, 3, 1> n2{T(last_normal.x()), T(last_normal.y()), T(last_normal.z())};
		// Eigen::Matrix<T, 3, 1> lpb{T(last_point_b.x()), T(last_point_b.y()), T(last_point_b.z())};

		// Eigen::Quaternion<T> q_last_curr{q[3], T(s) * q[0], T(s) * q[1], T(s) * q[2]};
		Eigen::Quaternion<T> q_last_curr{q[3], q[0], q[1], q[2]};

		Eigen::Matrix<T, 3, 1> n1_trans;
		n1_trans = q_last_curr * n1;

		Eigen::Matrix<T, 3, 1> nu = (n1_trans - n2).cross(n1_trans);
		// Eigen::Matrix<T, 3, 1> de = n2;

		residual[0] = nu.x(); //* T(3.0);//*10000.0;// / de.norm();
		residual[1] = nu.y(); //* T(3.0);//*10000.0;// / de.norm();
		residual[2] = nu.z(); //* T(3.0);//*10000.0;// / de.norm();

		return true;
	}

	static ceres::CostFunction *Create(const Eigen::Vector3d curr_normal_, const Eigen::Vector3d last_normal_)
	{
		return (new ceres::AutoDiffCostFunction<NormalFactor, 3, 4>(new NormalFactor(curr_normal_, last_normal_)));
	}

	Eigen::Vector3d curr_normal, last_normal;
	// double s;
};

/*Based on two normal vectors*/
struct NormalFactorOneResidual
{
	NormalFactorOneResidual(Eigen::Vector3d curr_normal_, Eigen::Vector3d last_normal_, double std_)
		: curr_normal(curr_normal_), last_normal(last_normal_), std(std_) {}

	template <typename T>
	bool operator()(const T *q, T *residual) const
	{

		Eigen::Matrix<T, 3, 1> n1{T(curr_normal.x()), T(curr_normal.y()), T(curr_normal.z())};
		Eigen::Matrix<T, 3, 1> n2{T(last_normal.x()), T(last_normal.y()), T(last_normal.z())};
		// Eigen::Matrix<T, 3, 1> lpb{T(last_point_b.x()), T(last_point_b.y()), T(last_point_b.z())};

		// Eigen::Quaternion<T> q_last_curr{q[3], T(s) * q[0], T(s) * q[1], T(s) * q[2]};
		Eigen::Quaternion<T> q_last_curr{q[3], q[0], q[1], q[2]};

		Eigen::Matrix<T, 3, 1> n1_trans;
		n1_trans = q_last_curr * n1;

		Eigen::Matrix<T, 3, 1> nu = (n1_trans - n2);
		// Eigen::Matrix<T, 3, 1> de = n2;

		residual[0] = asin(nu.norm() / T(2.0)) * T(2.0) / T(std); // radius of 5 deg* T(3.0);//*10000.0;// / de.norm();

		return true;
	}

	static ceres::CostFunction *Create(const Eigen::Vector3d curr_normal_, const Eigen::Vector3d last_normal_, const double std_)
	{
		return (new ceres::AutoDiffCostFunction<NormalFactorOneResidual, 1, 4>(new NormalFactorOneResidual(curr_normal_, last_normal_, std_)));
	}

	Eigen::Vector3d curr_normal, last_normal;
	double std;
};

struct LidarPlaneFactor
{
	LidarPlaneFactor(Eigen::Vector3d curr_point_, Eigen::Vector3d last_point_j_,
					 Eigen::Vector3d last_point_l_, Eigen::Vector3d last_point_m_, double s_)
		: curr_point(curr_point_), last_point_j(last_point_j_), last_point_l(last_point_l_),
		  last_point_m(last_point_m_), s(s_)
	{
		ljm_norm = (last_point_j - last_point_l).cross(last_point_j - last_point_m);
		ljm_norm.normalize();
	}

	template <typename T>
	bool operator()(const T *q, const T *t, T *residual) const
	{

		Eigen::Matrix<T, 3, 1> cp{T(curr_point.x()), T(curr_point.y()), T(curr_point.z())};
		Eigen::Matrix<T, 3, 1> lpj{T(last_point_j.x()), T(last_point_j.y()), T(last_point_j.z())};
		// Eigen::Matrix<T, 3, 1> lpl{T(last_point_l.x()), T(last_point_l.y()), T(last_point_l.z())};
		// Eigen::Matrix<T, 3, 1> lpm{T(last_point_m.x()), T(last_point_m.y()), T(last_point_m.z())};
		Eigen::Matrix<T, 3, 1> ljm{T(ljm_norm.x()), T(ljm_norm.y()), T(ljm_norm.z())};

		// Eigen::Quaternion<T> q_last_curr{q[3], T(s) * q[0], T(s) * q[1], T(s) * q[2]};
		Eigen::Quaternion<T> q_last_curr{q[3], q[0], q[1], q[2]};
		Eigen::Quaternion<T> q_identity{T(1), T(0), T(0), T(0)};
		q_last_curr = q_identity.slerp(T(s), q_last_curr);
		Eigen::Matrix<T, 3, 1> t_last_curr{T(s) * t[0], T(s) * t[1], T(s) * t[2]};

		Eigen::Matrix<T, 3, 1> lp;
		lp = q_last_curr * cp + t_last_curr;

		residual[0] = (lp - lpj).dot(ljm);

		return true;
	}

	static ceres::CostFunction *Create(const Eigen::Vector3d curr_point_, const Eigen::Vector3d last_point_j_,
									   const Eigen::Vector3d last_point_l_, const Eigen::Vector3d last_point_m_,
									   const double s_)
	{
		return (new ceres::AutoDiffCostFunction<
				LidarPlaneFactor, 1, 4, 3>(
			new LidarPlaneFactor(curr_point_, last_point_j_, last_point_l_, last_point_m_, s_)));
	}

	Eigen::Vector3d curr_point, last_point_j, last_point_l, last_point_m;
	Eigen::Vector3d ljm_norm;
	double s;
};

struct LidarPlaneFactorStd
{
	LidarPlaneFactorStd(Eigen::Vector3d curr_point_, Eigen::Vector3d last_point_j_,
						Eigen::Vector3d last_point_l_, Eigen::Vector3d last_point_m_, double s_, double std_)
		: curr_point(curr_point_), last_point_j(last_point_j_), last_point_l(last_point_l_),
		  last_point_m(last_point_m_), s(s_), std(std_)
	{
		ljm_norm = (last_point_j - last_point_l).cross(last_point_j - last_point_m);
		ljm_norm.normalize();
	}

	template <typename T>
	bool operator()(const T *q, const T *t, T *residual) const
	{

		Eigen::Matrix<T, 3, 1> cp{T(curr_point.x()), T(curr_point.y()), T(curr_point.z())};
		Eigen::Matrix<T, 3, 1> lpj{T(last_point_j.x()), T(last_point_j.y()), T(last_point_j.z())};
		// Eigen::Matrix<T, 3, 1> lpl{T(last_point_l.x()), T(last_point_l.y()), T(last_point_l.z())};
		// Eigen::Matrix<T, 3, 1> lpm{T(last_point_m.x()), T(last_point_m.y()), T(last_point_m.z())};
		Eigen::Matrix<T, 3, 1> ljm{T(ljm_norm.x()), T(ljm_norm.y()), T(ljm_norm.z())};

		// Eigen::Quaternion<T> q_last_curr{q[3], T(s) * q[0], T(s) * q[1], T(s) * q[2]};
		Eigen::Quaternion<T> q_last_curr{q[3], q[0], q[1], q[2]};
		Eigen::Quaternion<T> q_identity{T(1), T(0), T(0), T(0)};
		q_last_curr = q_identity.slerp(T(s), q_last_curr);
		Eigen::Matrix<T, 3, 1> t_last_curr{T(s) * t[0], T(s) * t[1], T(s) * t[2]};

		Eigen::Matrix<T, 3, 1> lp;
		lp = q_last_curr * cp + t_last_curr;

		residual[0] = (lp - lpj).dot(ljm) / T(std);

		return true;
	}

	static ceres::CostFunction *Create(const Eigen::Vector3d curr_point_, const Eigen::Vector3d last_point_j_,
									   const Eigen::Vector3d last_point_l_, const Eigen::Vector3d last_point_m_,
									   const double s_, const double std_)
	{
		return (new ceres::AutoDiffCostFunction<
				LidarPlaneFactorStd, 1, 4, 3>(
			new LidarPlaneFactorStd(curr_point_, last_point_j_, last_point_l_, last_point_m_, s_, std_)));
	}

	Eigen::Vector3d curr_point, last_point_j, last_point_l, last_point_m;
	Eigen::Vector3d ljm_norm;
	double s;
	double std;
};

struct LidarPlaneNormFactor
{
	LidarPlaneNormFactor(Eigen::Vector3d curr_point_, Eigen::Vector3d plane_unit_norm_,
						 double negative_OA_dot_norm_) : curr_point(curr_point_), plane_unit_norm(plane_unit_norm_),
														 negative_OA_dot_norm(negative_OA_dot_norm_) {}

	template <typename T>
	bool operator()(const T *q, const T *t, T *residual) const
	{
		Eigen::Quaternion<T> q_w_curr{q[3], q[0], q[1], q[2]};
		Eigen::Matrix<T, 3, 1> t_w_curr{t[0], t[1], t[2]};
		Eigen::Matrix<T, 3, 1> cp{T(curr_point.x()), T(curr_point.y()), T(curr_point.z())};
		Eigen::Matrix<T, 3, 1> point_w;
		point_w = q_w_curr * cp + t_w_curr;

		Eigen::Matrix<T, 3, 1> norm(T(plane_unit_norm.x()), T(plane_unit_norm.y()), T(plane_unit_norm.z()));
		residual[0] = norm.dot(point_w) + T(negative_OA_dot_norm);
		return true;
	}

	static ceres::CostFunction *Create(const Eigen::Vector3d curr_point_, const Eigen::Vector3d plane_unit_norm_,
									   const double negative_OA_dot_norm_)
	{
		return (new ceres::AutoDiffCostFunction<
				LidarPlaneNormFactor, 1, 4, 3>(
			new LidarPlaneNormFactor(curr_point_, plane_unit_norm_, negative_OA_dot_norm_)));
	}

	Eigen::Vector3d curr_point;
	Eigen::Vector3d plane_unit_norm;
	double negative_OA_dot_norm;
};

struct LidarDistanceFactor
{

	LidarDistanceFactor(Eigen::Vector3d curr_point_, Eigen::Vector3d closed_point_)
		: curr_point(curr_point_), closed_point(closed_point_) {}

	template <typename T>
	bool operator()(const T *q, const T *t, T *residual) const
	{
		Eigen::Quaternion<T> q_w_curr{q[3], q[0], q[1], q[2]};
		Eigen::Matrix<T, 3, 1> t_w_curr{t[0], t[1], t[2]};
		Eigen::Matrix<T, 3, 1> cp{T(curr_point.x()), T(curr_point.y()), T(curr_point.z())};
		Eigen::Matrix<T, 3, 1> point_w;
		point_w = q_w_curr * cp + t_w_curr;

		residual[0] = point_w.x() - T(closed_point.x());
		residual[1] = point_w.y() - T(closed_point.y());
		residual[2] = point_w.z() - T(closed_point.z());
		return true;
	}

	static ceres::CostFunction *Create(const Eigen::Vector3d curr_point_, const Eigen::Vector3d closed_point_)
	{
		return (new ceres::AutoDiffCostFunction<
				LidarDistanceFactor, 3, 4, 3>(
			new LidarDistanceFactor(curr_point_, closed_point_)));
	}

	Eigen::Vector3d curr_point;
	Eigen::Vector3d closed_point;
};


struct DistanceBetweenTwoEpochsConstraints
{

	DistanceBetweenTwoEpochsConstraints(double distance_threshold,double dis_std)
		: distance_threshold_(distance_threshold),dis_std_(dis_std){}

	template <typename T>
	bool operator()(const T *const p_a_ptr,
					const T *const q_a_ptr,
					const T *const p_b_ptr,
					const T *const q_b_ptr,
					T *residual) const
	{
		Eigen::Map<const Eigen::Matrix<T, 3, 1>> p_a(p_a_ptr);
		Eigen::Map<const Eigen::Quaternion<T>> q_a(q_a_ptr);
		Eigen::Map<const Eigen::Matrix<T, 3, 1>> p_b(p_b_ptr);
		Eigen::Map<const Eigen::Quaternion<T>> q_b(q_b_ptr);
		// Compute the relative transformation between the two frames.
		Eigen::Quaternion<T> q_a_inverse = q_a.conjugate();
		Eigen::Quaternion<T> q_ab_estimated = q_a_inverse * q_b;
		// Represent the displacement between the two frames in the A frame.
		Eigen::Matrix<T, 3, 1> p_ab_estimated = q_a_inverse * (p_b - p_a);
		T distance = p_ab_estimated.norm();
		T diff_res = distance - T(distance_threshold_);
		T diff_res_std = T(dis_std_);
		// if (diff_res < T(0.0))
		// {
		// 	residual[0] = T(0.0);
		// }
		// else
		{
			residual[0] = diff_res / diff_res_std;
		}
		return true;
	}

	static ceres::CostFunction *Create(const double distance_threshold, const double dis_std)
	{
		return (new ceres::AutoDiffCostFunction<
				DistanceBetweenTwoEpochsConstraints, 1, 3,4, 3,4>(
			new DistanceBetweenTwoEpochsConstraints(distance_threshold,dis_std)));
	}

	double distance_threshold_;
	double dis_std_;
};





// http://ceres-solver.org/nnls_tutorial.html?highlight=slam%20pose_graph_3d%20pose_graph_3d%20cc#other-examples  -11
// Computes the error term for two poses that have a relative pose measurement
// between them. Let the hat variables be the measurement. We have two poses x_a
// and x_b. Through sensor measurements we can measure the transformation of
// frame B w.r.t frame A denoted as t_ab_hat. We can compute an error metric
// between the current estimate of the poses and the measurement.
//
// In this formulation, we have chosen to represent the rigid transformation as
// a Hamiltonian quaternion, q, and position, p. The quaternion ordering is
// [x, y, z, w].
// The estimated measurement is:
//      t_ab = [ p_ab ]  = [ R(q_a)^T * (p_b - p_a) ]
//             [ q_ab ]    [ q_a^{-1] * q_b         ]
//
// where ^{-1} denotes the inverse and R(q) is the rotation matrix for the
// quaternion. Now we can compute an error metric between the estimated and
// measurement transformation. For the orientation error, we will use the
// standard multiplicative error resulting in:
//
//   error = [ p_ab - \hat{p}_ab                 ]
//           [ 2.0 * Vec(q_ab * \hat{q}_ab^{-1}) ]
//
// where Vec(*) returns the vector (imaginary) part of the quaternion. Since
// the measurement has an uncertainty associated with how accurate it is, we
// will weight the errors by the square root of the measurement information
// matrix:
//
//   residuals = I^{1/2) * error
// where I is the information matrix which is the inverse of the covariance.
class PoseGraph3dErrorTerm
{
public:
	PoseGraph3dErrorTerm(Eigen::Vector3d t_ab_measured, Eigen::Quaterniond q_ab_measured,
						 Eigen::Matrix<double, 6, 6> sqrt_information)
		: t_ab_measured_(std::move(t_ab_measured)), q_ab_measured_(std::move(q_ab_measured)),
		  sqrt_information_(std::move(sqrt_information)) {}
	template <typename T>
	bool operator()(const T *const p_a_ptr,
					const T *const q_a_ptr,
					const T *const p_b_ptr,
					const T *const q_b_ptr,
					T *residuals_ptr) const
	{
		Eigen::Map<const Eigen::Matrix<T, 3, 1>> p_a(p_a_ptr);
		Eigen::Map<const Eigen::Quaternion<T>> q_a(q_a_ptr);
		Eigen::Map<const Eigen::Matrix<T, 3, 1>> p_b(p_b_ptr);
		Eigen::Map<const Eigen::Quaternion<T>> q_b(q_b_ptr);
		// Compute the relative transformation between the two frames.
		Eigen::Quaternion<T> q_a_inverse = q_a.conjugate();
		Eigen::Quaternion<T> q_ab_estimated = q_a_inverse * q_b;
		// Represent the displacement between the two frames in the A frame.
		Eigen::Matrix<T, 3, 1> p_ab_estimated = q_a_inverse * (p_b - p_a);
		// Compute the error between the two orientation estimates.
		Eigen::Quaternion<T> delta_q =
			q_ab_measured_.template cast<T>() * q_ab_estimated.conjugate();
		// Compute the residuals.
		// [ position         ]   [ delta_p          ]
		// [ orientation (3x1)] = [ 2 * delta_q(0:2) ]
		Eigen::Map<Eigen::Matrix<T, 6, 1>> residuals(residuals_ptr);
		residuals.template block<3, 1>(0, 0) =
			p_ab_estimated - t_ab_measured_.template cast<T>();
		residuals.template block<3, 1>(3, 0) = T(2.0) * delta_q.vec();
		// Scale the residuals by the measurement uncertainty.
		residuals.applyOnTheLeft(sqrt_information_.template cast<T>());
		return true;
	}
	static ceres::CostFunction *Create(
		const Eigen::Vector3d t_ab_measured, const Eigen::Quaterniond q_ab_measured,
		const Eigen::Matrix<double, 6, 6> &sqrt_information)
	{
		return new ceres::AutoDiffCostFunction<PoseGraph3dErrorTerm, 6, 3, 4, 3, 4>(
			new PoseGraph3dErrorTerm(t_ab_measured, q_ab_measured, sqrt_information));
	}
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
	// The measurement for the position of B relative to A in the A frame.
	const Eigen::Vector3d t_ab_measured_;
	const Eigen::Quaterniond q_ab_measured_;
	// The square root of the measurement information matrix.
	const Eigen::Matrix<double, 6, 6> sqrt_information_;
};