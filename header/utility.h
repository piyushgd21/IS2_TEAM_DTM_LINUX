//#define TEST_PURPOSE
//#define EXPORT_RESULT
//#define TEST
#pragma once
#define MAP
#define EXPORT_TRAJECTORY
#define EXPORT_LOG

//#define DEBUGMAP
//#define MAPTEST
//#define FIT_CYLINDER_TREE 

//for exporting:
//#define LASLIB_IO

#ifndef _UTILITY_LIDAR_ODOMETRY_H_
#define _UTILITY_LIDAR_ODOMETRY_H_

// #include <ros/ros.h>
// #include <sensor_msgs/Imu.h>
// #include <sensor_msgs/PointCloud2.h>
// #include <nav_msgs/Odometry.h>
// #include "cloud_msgs/cloud_info.h"
// #include <tf/transform_broadcaster.h>
// #include <tf/transform_datatypes.h>
//#include <pcl-1.10/pcl_ros/point_cloud.h>
//#include <pcl-1.10/pcl_conversions/pcl_conversions.h>

#include <Eigen/Dense>

//#include <opencv/cv.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>

#include <pcl/range_image/range_image.h>
#include <pcl/filters/filter.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/common/common.h>
#include <pcl/registration/icp.h>
#include <pcl/octree/octree_search.h>

#include <vector>
#include <cmath>
#include <algorithm>
#include <queue>
#include <deque>
#include <iostream>
#include <fstream>
#include <ctime>
#include <cfloat>
#include <iterator>
#include <sstream>
#include <string>
#include <limits>
#include <iomanip>
#include <array>
#include <thread>
#include <mutex>
#include <ctime>
#include <cstdlib>
#include <chrono>
#include <unordered_map>

#include <boost/filesystem.hpp>

//ceres
#include <ceres/ceres.h>
#include <ceres/rotation.h>

#include <pdal/PointView.hpp>
#include <pdal/PointTable.hpp>
#include <pdal/Dimension.hpp>
#include <pdal/Options.hpp>
#include <pdal/StageFactory.hpp>
#include <pdal/io/BufferReader.hpp>
using namespace pdal;

#define OriBlockSize 4
#define PosBlockSize 3
#define CylinderBlockSize 7
#define PlaneBlockSize 4
#define SigmaFixed 1E-10
#define SigmaFree 1E+10



extern std::ofstream fMapping;
extern std::ofstream fMappingTraj;

extern std::ofstream fTrajectoryRef;
extern std::ofstream fTrajectoryRes;

extern std::ofstream fTrajectoryResFb;
extern std::ofstream fTrajectoryResPb;
extern std::ofstream f_odometry_debug;


extern std::ofstream fMappedPoints;
extern std::ofstream fLog;
extern std::ofstream f_mapping_debug;

extern std::ofstream fDebug;
extern std::ofstream fMappedRef;
extern std::string output_folder;
extern std::string input_folder;

extern std::string output_folder_odometry;
extern std::string output_folder_loopclouse;
extern std::string output_folder_iscan;
extern std::string output_folder_iscan_map;

extern double TimeFeature ;
extern double TimePoint;
extern double TimeLoading ;
extern double TimeSegment ;
extern double Tground, Ttree, Topt, Tcheck, Tp_extraction, Tp_opt;

extern std::mutex mTransMutex;
extern double SCAN_DURATION; //ms


#define PI 3.14159265359

using namespace std;

typedef pcl::PointXYZI PointType;

/*compute angle between two vectors*/
inline float computeAngle(Eigen::Vector3f r1, Eigen::Vector3f r2)
{
    float dot = r1(0) * r2(0) + r1(1) * r2(1) + r1(2) * r2(2);
    float cross_x = (r1(1) * r2(2) - r1(2) * r2(1));
    float cross_y = (r1(2) * r2(0) - r1(0) * r2(2));
    float cross_z = (r1(0) * r2(1) - r1(1) * r2(0));

    float det = sqrt(cross_x * cross_x + cross_y * cross_y + cross_z * cross_z);
    return atan2(det, dot);
}

inline double computeAngle(Eigen::Vector3d r1, Eigen::Vector3d r2)
{
    double dot = r1(0) * r2(0) + r1(1) * r2(1) + r1(2) * r2(2);
    double cross_x = (r1(1) * r2(2) - r1(2) * r2(1));
    double cross_y = (r1(2) * r2(0) - r1(0) * r2(2));
    double cross_z = (r1(0) * r2(1) - r1(1) * r2(0));

    double det = sqrt(cross_x * cross_x + cross_y * cross_y + cross_z * cross_z);
    return atan2(det, dot);
}

/*compute point x0, to line defined by x1, x2 distance*/
inline float computePoint2lineDistance(Eigen::Vector3f p0, Eigen::Vector3f p1, Eigen::Vector3f p2 = Eigen::Vector3f::Zero())
{
    return ((p0 - p1).cross(p0 - p2)).norm() / (p2 - p1).norm();
}

inline double computePoint2lineDistance(Eigen::Vector3d p0, Eigen::Vector3d p1, Eigen::Vector3d p2 = Eigen::Vector3d::Zero())
{
    return ((p0 - p1).cross(p0 - p2)).norm() / (p2 - p1).norm();
}

/*compute point to plane distance: p0 to plane defined by p1, p2, and p3 */
inline double computePoint2PlaneDistance(Eigen::Vector3d p0, Eigen::Vector3d p1, Eigen::Vector3d p2, Eigen::Vector3d p3)
{
    Eigen::Vector3d ljm_norm = (p1 - p2).cross(p1 - p3);
    ljm_norm.normalize();
    return ((p0 - p1).dot(ljm_norm));
}

/*compute point to plane distance: p0 to plane defined by p1, p2, and p3 */
inline Eigen::Vector3d computePoint2LineProjection(Eigen::Vector3d p0, Eigen::Vector3d p1, Eigen::Vector3d p2)
{
    return (p1 + (p0-p1).dot(p2-p1)/((p2-p1).dot(p2-p1)) * (p2-p1));

}

inline float rad2deg(float rad_)
{
    return rad_ * 180.0 / PI;
}

inline Eigen::Vector3f rad2deg(Eigen::Vector3f rad_)
{
    return rad_ * 180.0 / PI;
}

inline Eigen::Vector3d rad2deg(Eigen::Vector3d rad_)
{
    return rad_ * 180.0 / PI;
}

inline double rad2deg(double rad_)
{
    return rad_ * 180.0 / PI;
}

inline float deg2rad(float deg_)
{
    return deg_ / 180.0 * PI;
}

inline double deg2rad(double deg_)
{
    return deg_ / 180.0 * PI;
}

inline void norm2angle(const Eigen::Vector3f n, float &alpha, float &beta)
{
    beta = asin(n(2));
    alpha = atan2(n(1), n(0));
}

inline void angle2norm(float alpha, float beta, Eigen::Vector3f &n)
{
    n << cos(alpha) * cos(beta), sin(alpha) * cos(beta), sin(beta);
}

float median(vector<float> &v);
bool lineFitting(const std::vector<Eigen::Vector3f> points, Eigen::Vector3f &normal, Eigen::Vector3f &centerPoint);
bool planeFitting(const std::vector<Eigen::Vector3f> points, Eigen::Vector4f &params, Eigen::Vector3f &centerPoint);
bool planeFitting_outlier_check(const std::vector<Eigen::Vector3f> points, Eigen::Vector4f &params, Eigen::Vector3f &centerPoint, int &fail_type);

int planeParamTrans(Eigen::Vector4d &params);

Eigen::Vector3f compute2dSimilarity(const std::vector<Eigen::Vector3f> points1, const std::vector<Eigen::Vector3f> points2, int p);

class TicToc
{
public:
    TicToc()
    {
        tic();
    }

    void tic()
    {
        start = std::chrono::system_clock::now();
    }

    double toc()
    {
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        return elapsed_seconds.count() * 1000;
    }

private:
    std::chrono::time_point<std::chrono::system_clock> start, end;
};

class LidarScan;

struct SettingPara
{   
    // input related
    bool intensity_flag = false;

    
    //------------ Odometry Thread ----------------
    //odometry processing setting
    bool bFeatureBased = true;
    bool bLevel = false;
    bool bPointBased = true;
    bool bDistortion = false;
    bool bFPSimultanesouly = false;

    //scans range, must initialize 
    int initScan = -1;
    int endScan = -1;

    int nChannel = -1; // must initialize
    float minRange = 1.0;
    float maxRange = 100.0;

    // segment extraction
    int surfaceLength = 20;
    int minSegLength = 2;
    
    // ground cluster
    int min_points_per_ground_cluster = 500;

    // tree extraction
    float groundBufferTree = 2.5; // meter
    float treeAngleThreshold = 10.0;
    float minTreeDistance = 1.0;
    float maxTreeSegDistance = 1.0;
    int minNumSegTree = 3;

    //tree matching
    int minTreePair = 3;

    //ground point downsampling
    double ground_downsample_distance_per_scan = 0.2;

    //whether use trajectory information as odometry result 
    bool odo_from_trajectory_flag = true;
    //flag for export intermediate result

    bool export_intermediate_result_odometry = false;


    //----------- Mapping Thread -------------
    //global map
    bool global_map_flag = false;
    string global_ground_map_pass;
    string global_tree_map_pass;

    bool global_DTM = false;
    bool global_DTM_only = false;
    bool global_tree_locations = false;


    //scan integration
    bool iscan_opt_flag = true;
    bool bOdometry = false;
    int numInitIntegratedScan = 100;
    int numIntegratedScan = 10;
    double minNumTreePointPerScan = 2.0; //minimum number of tree points per scan
    int minNumPtsPerIscanTree = 100; //for final tree object in iscan

    double iscanTreeStd = 0.2;
    double iscanGroundStd = 0.1;
    double odoPositionStd = 0.1; //in meter
    double odoOrientationStd = 0.1; //deg

    //(TODO)before
    // double odoPositionStdTraj = 0.05; //in meter  if odometry result provided by trajectory
    // double odoOrientationStdTraj = 0.02; //deg if odometry result provided by trajectory

    //(TODO) cuurent(0622)
    //  double odoPositionStdTraj = 0.03; //in meter  if odometry result provided by trajectory
    //  double odoOrientationStdTraj = 0.01; //deg if odometry result provided by trajectory

    //(TODO) for downsampling
     double odoPositionStdTraj = 0.05; //in meter  if odometry result provided by trajectory
     double odoOrientationStdTraj = 0.02; //deg if odometry result provided by trajectory

    //(TODO)after 5.31
    //double odoPositionStdTraj = 1e-6; //in meter  if odometry result provided by trajectory
    //double odoOrientationStdTraj = 1e-6; //deg if odometry result provided by trajectory


    double ground_downsample_distance_scan_integration = 0.1;
    
    int minimum_tree_iscan = 8;


    // scan to map
    double global_tree_loc_std = 0.2; // 0.2 m  for natural
    double mapTreeStd = 0.2;
    double mapGroundStd = 0.1;
    bool bMapTreeRadius = true;
    bool bMapTreeDirec = true;

    //loop closure
    int loop_closure_iscan_num = -1;
    double patch_size_lc = 1.0;
    bool iscan_lc_iterative_flag = false; //conduct iterative lc for iscan level
    bool perform_iscan_lc_at_last_flag = false; //conduct iscan level lc before the last raw scan level LC
    int number_iterations_raw_scan_lc = 1;

    //fetch raw tree points and reconduct optimization
    bool reoptimize_fetch_tree_points = false;
    int reoptimize_iterations = 1;
    double distance_threshold_fetch_points = 0.3;

    //fetch raw points
    bool reoptimize_fetch_raw_points = false;

    //flag for export intermediate result
    bool export_intermediate_result_mapping = false;
    
};

void loadSettingPara(string settingFile, SettingPara &para);

void Checkforconstantparams(ceres::Problem * problem, double * X, double *stdX, const int blocksize, const int nblocks);

struct CylinderPara
{
    Eigen::Vector3d x;
    Eigen::Vector3d n; //not unit, in general nz = 1;
    double r;
};

void Compute_Rotation(double om, double phi, double kap, Eigen::Matrix3d &R);
Eigen::Matrix3d Compute_Rotation(double om, double phi, double kap);
Eigen::Matrix3f Compute_Rotation(float om, float phi, float kap);

void Find_Rotation(Eigen::Matrix3d R, double &ome, double &phi, double &kap);
Eigen::Vector3d Find_Rotation(Eigen::Matrix3d R);
Eigen::Vector3f Find_Rotation(Eigen::Matrix3f R);

/*compute point to cylinder distance*/
inline double computePoint2CylinderDistance(Eigen::Vector3d p0, CylinderPara cylinder)
{
    double d = computePoint2lineDistance(p0, cylinder.x, cylinder.x + cylinder.n);
    return (d - cylinder.r);
}

void downsamplePointCloudDistance(pcl::PointCloud<PointType>::Ptr original_pc, double distance, pcl::PointCloud<PointType> &downsampeld_pc);

pcl::PointCloud<PointType> DownSamplePointCloudBasedOnDistance(pcl::PointCloud<PointType>::Ptr original_pc, double distance);


// inline double computePoint2CylinderDistance(Eigen::Vector3f p0, CylinderPara cylinder)
// {
//     double d = computePoint2lineDistance(p0.cast<double>(), cylinder.x, cylinder.x + cylinder.n);
//     return (d - cylinder.r);
// }

// extern const string pointCloudTopic = "/velodyne_points";
// extern const string imuTopic = "/imu/data";

// // Save pcd
// extern const string fileDirectory = "/tmp/";

// // Using velodyne cloud "ring" channel for image projection (other lidar may have different name for this channel, change "PointXYZIR" below)
// extern const bool useCloudRing = true; // if true, ang_res_y and ang_bottom are not used

// // VLP-16
// extern const int N_SCAN = 16;
// extern const int Horizon_SCAN = 1800;
// extern const float ang_res_x = 0.2;
// extern const float ang_res_y = 2.0;
// extern const float ang_bottom = 15.0 + 0.1;
// extern const int groundScanInd = 7;

// // HDL-32E
// // extern const int N_SCAN = 32;
// // extern const int Horizon_SCAN = 1800;
// // extern const float ang_res_x = 360.0/float(Horizon_SCAN);
// // extern const float ang_res_y = 41.33/float(N_SCAN-1);
// // extern const float ang_bottom = 30.67;
// // extern const int groundScanInd = 20;

// // VLS-128
// // extern const int N_SCAN = 128;
// // extern const int Horizon_SCAN = 1800;
// // extern const float ang_res_x = 0.2;
// // extern const float ang_res_y = 0.3;
// // extern const float ang_bottom = 25.0;
// // extern const int groundScanInd = 10;

// // Ouster users may need to uncomment line 159 in imageProjection.cpp
// // Usage of Ouster imu data is not fully supported yet (LeGO-LOAM needs 9-DOF IMU), please just publish point cloud data
// // Ouster OS1-16
// // extern const int N_SCAN = 16;
// // extern const int Horizon_SCAN = 1024;
// // extern const float ang_res_x = 360.0/float(Horizon_SCAN);
// // extern const float ang_res_y = 33.2/float(N_SCAN-1);
// // extern const float ang_bottom = 16.6+0.1;
// // extern const int groundScanInd = 7;

// // Ouster OS1-64
// // extern const int N_SCAN = 64;
// // extern const int Horizon_SCAN = 1024;
// // extern const float ang_res_x = 360.0/float(Horizon_SCAN);
// // extern const float ang_res_y = 33.2/float(N_SCAN-1);
// // extern const float ang_bottom = 16.6+0.1;
// // extern const int groundScanInd = 15;

// extern const bool loopClosureEnableFlag = false;
// extern const double mappingProcessInterval = 0.3;

// extern const float scanPeriod = 0.1;
// extern const int systemDelay = 0;
// extern const int imuQueLength = 200;

// extern const float sensorMinimumRange = 1.0;
// extern const float sensorMountAngle = 0.0;
// extern const float segmentTheta = 60.0 / 180.0 * M_PI; // decrese this value may improve accuracy
// extern const int segmentValidPointNum = 5;
// extern const int segmentValidLineNum = 3;
// extern const float segmentAlphaX = ang_res_x / 180.0 * M_PI;
// extern const float segmentAlphaY = ang_res_y / 180.0 * M_PI;

// extern const int edgeFeatureNum = 2;
// extern const int surfFeatureNum = 4;
// extern const int sectionsTotal = 6;
// extern const float edgeThreshold = 0.1;
// extern const float surfThreshold = 0.1;
// extern const float nearestFeatureSearchSqDist = 25;

// // Mapping Params
// extern const float surroundingKeyframeSearchRadius = 50.0; // key frame that is within n meters from current pose will be considerd for scan-to-map optimization (when loop closure disabled)
// extern const int surroundingKeyframeSearchNum = 50;        // submap size (when loop closure enabled)
// // history key frames (history submap for loop closure)
// extern const float historyKeyframeSearchRadius = 7.0; // key frame that is within n meters from current pose will be considerd for loop closure
// extern const int historyKeyframeSearchNum = 25;       // 2n+1 number of hostory key frames will be fused into a submap for loop closure
// extern const float historyKeyframeFitnessScore = 0.3; // the smaller the better alignment

// extern const float globalMapVisualizationSearchRadius = 500.0; // key frames with in n meters will be visualized

// struct smoothness_t
// {
//     float value;
//     size_t ind;
// };

// struct by_value
// {
//     bool operator()(smoothness_t const &left, smoothness_t const &right)
//     {
//         return left.value < right.value;
//     }
// };

// /*
//     * A point cloud type that has "ring" channel

// struct PointXYZIR
// {
//     PCL_ADD_POINT4D
//     PCL_ADD_INTENSITY;
//     uint16_t ring;
//     EIGEN_MAKE_ALIGNED_OPERATOR_NEW
// } EIGEN_ALIGN16;

// POINT_CLOUD_REGISTER_POINT_STRUCT (PointXYZIR,
//                                    (float, x, x) (float, y, y)
//                                    (float, z, z) (float, intensity, intensity)
//                                    (uint16_t, ring, ring)
// )
//  */
// /*
//  * A point cloud type that has 6D pose info ([x,y,z,roll,pitch,yaw] intensity is time stamp)
//  */
// struct PointXYZIRPYT
// {
//     PCL_ADD_POINT4D
//     PCL_ADD_INTENSITY;
//     float roll;
//     float pitch;
//     float yaw;
//     double time;
//     EIGEN_MAKE_ALIGNED_OPERATOR_NEW
// } EIGEN_ALIGN16;

// POINT_CLOUD_REGISTER_POINT_STRUCT(PointXYZIRPYT,
//                                   (float, x, x)(float, y, y)(float, z, z)(float, intensity, intensity)(float, roll, roll)(float, pitch, pitch)(float, yaw, yaw)(double, time, time))

// typedef PointXYZIRPYT PointTypePose;

#endif
