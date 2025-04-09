//#pragma once
#ifndef _TRAJECTORY_H_
#define _TRAJECTORY_H_

#include "../header/utility.h"

struct EOP
{
	// unit: m and rad
	double XO;
	double YO;
	double ZO;
	double omega;
	double phi;
	double kappa;
};

struct BOP
{
	double time;
	int id;
	EOP pos;
};

struct RefPoint
{
	int trajID; // first element of the trajectory
	double t_mean;
	Eigen::MatrixXd M_Matrix;
	BOP bop;
};

class CTrajectory
{

public:
	CTrajectory();
	~CTrajectory();
	// function
	// void loadTrajConfig(string trajConfigFile);
	void loadTraj(string trajFile);
	void LoadSlamResultTrajectory(std::string slam_traj_path,const std::vector<double>& lidar_times);
	bool bopInterpolation(double tempTime, int startIndex, int &bopID, Eigen::Vector3d &r, Eigen::Matrix3d &R);
	void loadLidarPara(string paraFile);
//	bool bopInterpolation(double tempTime, int startIndex = -1);
	vector<BOP> bopList; // bop list
	vector<BOP> eopList; // eop list

	Eigen::Vector3d r_lu_b;
	Eigen::Matrix3d R_lu_b;


	//-----------------------
	bool flagRef = false; // 0: based on data rate, 1: based on event marker
	int polyOrder;
	double dataRate;
	int qRef;		  // number of reference points for estimaing polynominal
	double std_pos;	  // m
	double std_om_ph; // degree
	double std_kp;	  // degree
	vector<double> refTimeList;
	vector<BOP> refPoints; // for each refpoint
	vector<RefPoint> refSegments;
	vector<EOP> refCorrections;
	void setRefTime(vector<double> timeList);
	void extractRefPoint();
	void computeEop();
	bool eopInterpolation(double tempTime, int startIndex, int &bopID, Eigen::Vector3d &r, Eigen::Matrix3d &R);


private:
};

// void Compute_Rotation(double om, double phi, double kap, Eigen::Matrix3d &R);
// Eigen::Matrix3d Compute_Rotation(double om, double phi, double kap);
// Eigen::Matrix3f Compute_Rotation(float om, float phi, float kap);

// void Find_Rotation(Eigen::Matrix3d R, double &ome, double &phi, double &kap);
// Eigen::Vector3d Find_Rotation(Eigen::Matrix3d R);
// Eigen::Vector3f Find_Rotation(Eigen::Matrix3f R);
void Interpolate_Pos(BOP b1, BOP b2, double t, EOP &e_new);
void Spherical_Linear_Interpolation(BOP b1, BOP b2, double t_ratio, EOP &e_new);
void Rot2Qua(Eigen::Matrix3d R, Eigen::Vector4d &q);
void Qua2Eul(Eigen::Vector4d q, double &ome, double &phi, double &kap);

#endif