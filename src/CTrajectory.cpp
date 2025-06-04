#include "../header/CTrajectory.h"

/**
 * @brief Constructor for CTrajectory class.
 */
CTrajectory::CTrajectory()
{
}

/**
 * @brief Destructor for CTrajectory class.
 */
CTrajectory::~CTrajectory()
{
}
// insert timetag from image for selecting reference points
/**
 * @brief Set reference times (e.g., from image timestamps) used to select keyframes or matching points.
 * 
 * @param timeList Vector of time values (in seconds).
 */

void CTrajectory::setRefTime(vector<double> timeList)
{
	refTimeList = timeList;
}

/**
 * @brief Load the LiDAR-to-body transformation from a parameter file.
 * 
 * The file should contain: dx dy dz omega phi kappa
 * The transformation is stored in `r_lu_b` and `R_lu_b`.
 * 
 * @param paraFile Path to the LiDAR parameter file.
 */
void CTrajectory::loadLidarPara(string paraFile)
{
	ifstream f_para;
	f_para.open(paraFile);

	if (!f_para)
	{
		cout << paraFile << endl;
		throw std::runtime_error("Wrong Lidar para file");
		exit(0);
	}
	double ome, phi, kap;
	f_para >> r_lu_b(0) >> r_lu_b(1) >> r_lu_b(2) >> ome >> phi >> kap;
	Compute_Rotation(deg2rad(ome), deg2rad(phi), deg2rad(kap), R_lu_b);
}

/**
 * @brief Load trajectory from a GNSS/INS file.
 * 
 * Supports multiple formats based on the first line in the file (gnssinsType).
 * 
 * @param trajFile Path to the trajectory file.
 */

void CTrajectory::loadTraj(string trajFile)
{
	// 0 = Time X Y Z ome phi kap
	// 1 = APX15, POSPac SBET File.
	// 2 = Novatel IGM, Inertial Explorer. TIme, X Y Z Z_alter roll pitch heading ome phi kap
	int gnssinsType = 0; //

	// event file
	ifstream f_tra;
	f_tra.open(trajFile);

	if (!f_tra)
	{
		cout << "failed to open" << trajFile << endl;
		throw std::runtime_error("Wrong trajectory file");

		exit(0);
	}
	f_tra >> gnssinsType;

	vector<BOP> tempTraj;
	BOP bop;
	double ome, phi, kap, temp;
	Eigen::Matrix3d R_b_m, R_b_ned, R_ned_enu, R1, R2, R3;
	string line;
	while (getline(f_tra, line))
	{
		if (/*line.find_first_of("!") == std::string::npos &&*/ (line.find_first_not_of("\t\n\v\f\r") != std::string::npos))
		{
			vector<string> words;
			std::istringstream iss(line);
			do
			{
				// Read a word
				string word;
				iss >> word;
				if (!word.empty())
					words.push_back(word);
			} while (iss);

			if (words.size() < 7)
			{
				cout << line;
				cout << "wrong setting format" << words.size() << endl;
			}

			if (gnssinsType == 2) // inertial explorer
			{
				bop.time = stod(words[0]);
				bop.pos.XO = stod(words[1]);
				bop.pos.YO = stod(words[2]);
				bop.pos.ZO = stod(words[3])/*+1*/;
				bop.pos.omega = deg2rad(stod(words[8]));
				bop.pos.phi = deg2rad(stod(words[9]));
				bop.pos.kappa = deg2rad(stod(words[10]));
			}
			else if(gnssinsType == 3) //SLAM enhanced traj
			{
				bop.time = stod(words[0]);
				bop.pos.XO = stod(words[1]);
				bop.pos.YO = stod(words[2]);
				bop.pos.ZO = stod(words[3]);
				bop.pos.omega = deg2rad(stod(words[4]));
				bop.pos.phi = deg2rad(stod(words[5]));
				bop.pos.kappa = deg2rad(stod(words[6]));
			}
			else if (gnssinsType == 1) // pospac
			{
				bop.time = stod(words[0]);
				bop.pos.XO = stod(words[1]);
				bop.pos.YO = stod(words[2]);
				bop.pos.ZO = stod(words[3]);

				ome = deg2rad(stod(words[4]));
				phi = deg2rad(stod(words[5]));
				kap = deg2rad(stod(words[6]));

				Compute_Rotation(0, 0, kap, R1);
				Compute_Rotation(0, phi, 0, R2);
				Compute_Rotation(ome, 0, 0, R3);
				R_b_ned = R1 * R2 * R3;
				Compute_Rotation(deg2rad(180.0), 0.0, deg2rad(-90.0), R_ned_enu);
				R_b_m = R_ned_enu * R_b_ned;
				Find_Rotation(R_b_m, bop.pos.omega, bop.pos.phi, bop.pos.kappa);
			}
			else if (gnssinsType == 0)
			{
				bop.time = stod(words[0]);
				bop.pos.XO = stod(words[1]);
				bop.pos.YO = stod(words[2]);
				bop.pos.ZO = stod(words[3]);
				bop.pos.omega = deg2rad(stod(words[4]));
				bop.pos.phi = deg2rad(stod(words[5]));
				bop.pos.kappa = deg2rad(stod(words[6]));
			}
			tempTraj.push_back(bop);
		}
	}
	bopList = tempTraj;

	cout << std::fixed << setprecision(6) << bopList[0].time << "\t" << bopList[0].pos.XO << "\t" << bopList[0].pos.YO
		 << "\t" << bopList[0].pos.ZO << "\t" << rad2deg(bopList[0].pos.omega) << "\t"
		 << rad2deg(bopList[0].pos.phi) << "\t" << rad2deg(bopList[0].pos.kappa) << endl;

	cout << bopList.back().time << "\t" << bopList.back().pos.XO << "\t" << bopList.back().pos.YO
		 << "\t" << bopList.back().pos.ZO << "\t" << rad2deg(bopList.back().pos.omega) << "\t"
		 << rad2deg(bopList.back().pos.phi) << "\t" << rad2deg(bopList.back().pos.kappa) << endl;
}

/**
 * @brief Load SLAM-enhanced trajectory (post-processed), aligning times using lidar frame timestamps.
 * 
 * @param slam_traj_path Path to SLAM result file.
 * @param lidar_times List of LiDAR timestamps (in ms).
 */

void CTrajectory::LoadSlamResultTrajectory(std::string slam_traj_path,const std::vector<double>& lidar_times)
{
		// event file
	ifstream f_tra;
	f_tra.open(slam_traj_path);

	if (!f_tra)
	{
		cout << "failed to open" << slam_traj_path << endl;
		throw std::runtime_error("Wrong trajectory file");

		exit(0);
	}

	std::vector<BOP> tempTraj;
	BOP bop;
	string line;
	//first line: header
	
	getline(f_tra, line);
	
	while (getline(f_tra, line))
	{
		if (line.find_first_of("!") == std::string::npos && (line.find_first_not_of("\t\n\v\f\r") != std::string::npos))
		{
			vector<string> words;
			std::istringstream iss(line);
			do
			{
				// Read a word
				string word;
				iss >> word;
				if (!word.empty())
					words.push_back(word);
			} while (iss);

			if (words.size() < 8)
			{
				cout << line;
				cout << "wrong setting format" << words.size() << endl;
			}

			if (stoi(words[0]) == 0)
			{
				int scan_index = (stoi(words[4]));
				bop.time = lidar_times[scan_index]/1000;
				bop.pos.XO = stod(words[1]);
				bop.pos.YO = stod(words[2]);
				bop.pos.ZO = stod(words[3]);
				bop.pos.omega = deg2rad(stod(words[5]));
				bop.pos.phi = deg2rad(stod(words[6]));
				bop.pos.kappa = deg2rad(stod(words[7]));
				tempTraj.push_back(bop);

				// scan_pose[scan_index] = bop;
			}
		}
	}
	bopList = tempTraj;
	eopList = bopList;
	//cout << scan_pose.size() << "pose" <<endl;
}
/**
 * @brief Convert body-frame trajectory to LiDAR-unit frame using transformation parameters.
 */
void CTrajectory::computeEop()
{
	eopList = bopList;
	Eigen::Matrix3d R_b_m, R_lu_m;
	Eigen::Vector3d r_b_m, r_lu_m;
	for (int i = 0; i < bopList.size(); i++)
	{
		Compute_Rotation(bopList[i].pos.omega, bopList[i].pos.phi, bopList[i].pos.kappa, R_b_m);
		r_b_m << bopList[i].pos.XO, bopList[i].pos.YO, bopList[i].pos.ZO;

		r_lu_m = r_b_m + R_b_m * r_lu_b;
		R_lu_m = R_b_m * R_lu_b;

		eopList[i].pos.XO = r_lu_m(0);
		eopList[i].pos.YO = r_lu_m(1);
		eopList[i].pos.ZO = r_lu_m(2);

		Find_Rotation(R_lu_m, eopList[i].pos.omega, eopList[i].pos.phi, eopList[i].pos.kappa);
	}
	// Write out transformed trajectory
    std::string outPass =  output_folder+ "eop_info.txt";
    std::ofstream fEop(outPass, std::ifstream::out);
    fEop << fixed << std::setprecision(6);
	for (int i = 0; i < eopList.size(); i++)
	{
		fEop << eopList[i].time << "\t" << eopList[i].pos.XO << "\t" << eopList[i].pos.YO << "\t" << eopList[i].pos.ZO << "\t"
			 << rad2deg(eopList[i].pos.omega) << "\t" << rad2deg(eopList[i].pos.phi) << "\t" << rad2deg(eopList[i].pos.kappa) << endl;
	}
	fEop.close();
}

/**
 * @brief Interpolates position and orientation at a given time from EOP list.
 */

bool CTrajectory::eopInterpolation(double tempTime, int startIndex, int &bopID, Eigen::Vector3d &r, Eigen::Matrix3d &R)
{
	// Look back slightly to ensure match
	startIndex = max(0, startIndex - 5);

	// find the corresponding trajectory time from trajectory file
	int bopIndex = -1;
	for (int idx_bop = startIndex; idx_bop < bopList.size() - 1; idx_bop++)
	{
		if (bopList[idx_bop].time <= tempTime && bopList[idx_bop + 1].time >= tempTime)
		{
			bopIndex = idx_bop;
			break;
		}
	}

	if (bopIndex == -1)
		return false;

	
	// Linear interpolation of position
	// second option
	Eigen::Vector3d r1(eopList[bopIndex].pos.XO, eopList[bopIndex].pos.YO, eopList[bopIndex].pos.ZO);
	Eigen::Vector3d r2(eopList[bopIndex + 1].pos.XO, eopList[bopIndex + 1].pos.YO, eopList[bopIndex + 1].pos.ZO);
	float t_ratio = (tempTime - eopList[bopIndex].time) / (eopList[bopIndex + 1].time - eopList[bopIndex].time);

	// cout << r1.transpose() <<"\t" << r2.transpose() <<"\t"<< t_ratio <<endl;
	Eigen::Vector3d r_new = r1 + (r2 - r1) * t_ratio;

	// SLERP interpolation of rotation
	Eigen::Matrix3d R1, R2, R_new;

	Eigen::AngleAxisd pitchAngle(eopList[bopIndex].pos.omega, Eigen::Vector3d::UnitX());
	Eigen::AngleAxisd yawAngle(eopList[bopIndex].pos.phi, Eigen::Vector3d::UnitY());
	Eigen::AngleAxisd rollAngle(eopList[bopIndex].pos.kappa, Eigen::Vector3d::UnitZ());
	R1 = pitchAngle * yawAngle * rollAngle;
	Eigen::Quaterniond q1 = Eigen::Quaterniond(R1);

	Eigen::AngleAxisd pitchAngle2(eopList[bopIndex + 1].pos.omega, Eigen::Vector3d::UnitX());
	Eigen::AngleAxisd yawAngle2(eopList[bopIndex + 1].pos.phi, Eigen::Vector3d::UnitY());
	Eigen::AngleAxisd rollAngle2(eopList[bopIndex + 1].pos.kappa, Eigen::Vector3d::UnitZ());
	R2 = pitchAngle2 * yawAngle2 * rollAngle2;
	Eigen::Quaterniond q2 = Eigen::Quaterniond(R2);

	Eigen::Quaterniond q_point_last = q1.slerp(t_ratio, q2);

	R_new = q_point_last.toRotationMatrix();
	Eigen::Vector3d euler_angles = R_new.eulerAngles(0, 1, 2);
	r = r_new;
	R = R_new;
	bopID = bopIndex;

	return true;
}


/* To improve the efficiency, start inpex is optional*/
/**
 * @brief Basic linear + quaternion interpolation from BOP trajectory.
 */
bool CTrajectory::bopInterpolation(double tempTime, int startIndex, int &bopID, Eigen::Vector3d &r, Eigen::Matrix3d &R)
{
	startIndex = max(0, startIndex - 5);
	// find the corresponding trajectory time from trajectory file
	int bopIndex = -1;
	for (int idx_bop = startIndex; idx_bop < bopList.size() - 1; idx_bop++)
	{
		if (bopList[idx_bop].time <= tempTime && bopList[idx_bop + 1].time >= tempTime)
		{
			bopIndex = idx_bop;
			break;
		}
	}

	if (bopIndex == -1)
		return false;

	// second option
	Eigen::Vector3d r1(bopList[bopIndex].pos.XO, bopList[bopIndex].pos.YO, bopList[bopIndex].pos.ZO);
	Eigen::Vector3d r2(bopList[bopIndex + 1].pos.XO, bopList[bopIndex + 1].pos.YO, bopList[bopIndex + 1].pos.ZO);
	float t_ratio = (tempTime - bopList[bopIndex].time) / (bopList[bopIndex + 1].time - bopList[bopIndex].time);

	Eigen::Vector3d r_new = r1 + (r2 - r1) * t_ratio;

	Eigen::Matrix3d R1, R2, R_new;

	Eigen::AngleAxisd pitchAngle(bopList[bopIndex].pos.omega, Eigen::Vector3d::UnitX());
	Eigen::AngleAxisd yawAngle(bopList[bopIndex].pos.phi, Eigen::Vector3d::UnitY());
	Eigen::AngleAxisd rollAngle(bopList[bopIndex].pos.kappa, Eigen::Vector3d::UnitZ());
	R1 = pitchAngle * yawAngle * rollAngle;
	Eigen::Quaterniond q1 = Eigen::Quaterniond(R1);

	Eigen::AngleAxisd pitchAngle2(bopList[bopIndex + 1].pos.omega, Eigen::Vector3d::UnitX());
	Eigen::AngleAxisd yawAngle2(bopList[bopIndex + 1].pos.phi, Eigen::Vector3d::UnitY());
	Eigen::AngleAxisd rollAngle2(bopList[bopIndex + 1].pos.kappa, Eigen::Vector3d::UnitZ());
	R2 = pitchAngle2 * yawAngle2 * rollAngle2;
	Eigen::Quaterniond q2 = Eigen::Quaterniond(R2);

	Eigen::Quaterniond q_point_last = q1.slerp(t_ratio, q2);

	R_new = q_point_last.toRotationMatrix();
	Eigen::Vector3d euler_angles = R_new.eulerAngles(0, 1, 2);
	r = r_new;
	R = R_new;
	bopID = bopIndex;
	return true;
}

/**
 * @brief Interpolate position and orientation between two BOPs.
 * 
 * @param b1 First BOP.
 * @param b2 Second BOP.
 * @param t Target time.
 * @param[out] e_new Interpolated EOP.
 */

void Interpolate_Pos(BOP b1, BOP b2, double t, EOP &e_new)
{
	double t_ratio = (t - b1.time) / (b2.time - b1.time);
	e_new.XO = b1.pos.XO + (b2.pos.XO - b1.pos.XO) * t_ratio;
	e_new.YO = b1.pos.YO + (b2.pos.YO - b1.pos.YO) * t_ratio;
	e_new.ZO = b1.pos.ZO + (b2.pos.ZO - b1.pos.ZO) * t_ratio;

	Spherical_Linear_Interpolation(b1, b2, t_ratio, e_new);
}

/**
 * @brief Performs SLERP between two BOP rotations and updates interpolated EOP.
 */

void Spherical_Linear_Interpolation(BOP b1, BOP b2, double t_ratio, EOP &e_new)
{
	Eigen::Matrix3d R1, R2;
	Compute_Rotation(b1.pos.omega, b1.pos.phi, b1.pos.kappa, R1);
	Compute_Rotation(b2.pos.omega, b2.pos.phi, b2.pos.kappa, R2);

	Eigen::Vector4d q1, q2;
	Rot2Qua(R1, q1);
	Rot2Qua(R2, q2);

	if (q1.dot(q2) < 0)
		q2 = q2 * (-1);

	double C1, C2;
	double theta = acos(q1.dot(q2));
	if (theta > 0.000001)
	{
		C1 = sin((1 - t_ratio) * theta) / sin(theta);
		C2 = sin(t_ratio * theta) / sin(theta);
	}
	else
	{
		C1 = 1 - t_ratio;
		C2 = t_ratio;
	}

	Eigen::Vector4d q_new;
	q_new = C1 * q1 + C2 * q2;

	double ome, phi, kap;
	Qua2Eul(q_new, ome, phi, kap);
	e_new.omega = ome;
	e_new.phi = phi;
	e_new.kappa = kap;
}

/**
 * @brief Convert a rotation matrix to quaternion vector (q0, qx, qy, qz).
 */

void Rot2Qua(Eigen::Matrix3d R, Eigen::Vector4d &q)
{
	double q0, qx, qy, qz;
	if (R(0, 0) + R(1, 1) + R(2, 2) + 1.0 > 0)
	{
		q0 = sqrt(R(0, 0) + R(1, 1) + R(2, 2) + 1.0) / 2;
		qx = (R(2, 1) - R(1, 2)) / (4 * q0);
		qy = (R(0, 2) - R(2, 0)) / (4 * q0);
		qz = (R(1, 0) - R(0, 1)) / (4 * q0);
	}
	else if (R(0, 0) - R(1, 1) - R(2, 2) + 1.0 > 0)
	{
		qx = sqrt(R(0, 0) - R(1, 1) - R(2, 2) + 1.0) / 2;
		qy = (R(0, 1) + R(1, 0)) / (4 * qx);
		qz = (R(0, 2) + R(2, 0)) / (4 * qx);
		q0 = (R(2, 1) - R(1, 2)) / (4 * qx);
	}
	else if (R(1, 1) - R(0, 0) - R(2, 2) + 1.0 > 0)
	{
		qy = sqrt(R(1, 1) - R(0, 0) - R(2, 2) + 1.0) / 2;
		qx = (R(0, 1) + R(1, 0)) / (4 * qy);
		qz = (R(1, 2) + R(2, 1)) / (4 * qy);
		q0 = (R(0, 2) - R(2, 0)) / (4 * qy);
	}
	else if (R(2, 2) - R(0, 0) - R(1, 1) + 1.0 > 0)
	{
		qz = sqrt(R(2, 2) - R(0, 0) - R(1, 1) + 1.0) / 2;
		qx = (R(0, 2) + R(2, 0)) / (4 * qz);
		qy = (R(1, 2) + R(2, 1)) / (4 * qz);
		q0 = (R(1, 0) - R(0, 1)) / (4 * qz);
	}

	//	 Normalize quaternion here
	double magnitude = sqrt(q0 * q0 + qx * qx + qy * qy + qz * qz);
	q0 = q0 / magnitude;
	qx = qx / magnitude;
	qy = qy / magnitude;
	qz = qz / magnitude;

	q << q0, qx, qy, qz;
}

/**
 * @brief Convert a quaternion to Euler angles (omega, phi, kappa).
 */
void Qua2Eul(Eigen::Vector4d q, double &ome, double &phi, double &kap)
{
	double q0, qx, qy, qz;
	q0 = q(0);
	qx = q(1);
	qy = q(2);
	qz = q(3);

	double r11, r12, r13, r23, r33;

	r11 = qx * qx + q0 * q0 - qz * qz - qy * qy;
	r12 = 2 * qx * qy - 2 * q0 * qz;
	r13 = 2 * qx * qz + 2 * q0 * qy;
	r23 = 2 * qy * qz - 2 * q0 * qx;
	r33 = qz * qz - qy * qy - qx * qx + q0 * q0;

	phi = asin(r13);
	kap = atan2(-r12 / cos(phi), r11 / cos(phi));
	ome = atan2(-r23 / cos(phi), r33 / cos(phi));
}
