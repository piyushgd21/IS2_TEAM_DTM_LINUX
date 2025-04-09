#ifndef _ISCAN_H_
#define _ISCAN_H_

#include "../header/utility.h"
#include "../header/CTrajectory.h"

class ScanInfo
{
public:
    //general information related to scan
    int scanID;
    int partID;
    double mTimeInit; // ms
    double mTimeEnd;
    double mTimeTracking;

    //pose information
    Eigen::Matrix3d R_lu_lup; //leveling rotation from the odometry
    Eigen::Matrix3d R_lut2_lut1; //odometry based estimation
    Eigen::Vector3d r_lut2_lut1; //odometry based estimation

    //if valid feature
    bool valid_feature_flag = false;
    //tree information
    int num_tree = 0;
    std::vector<vector<PointType>> vTreeRawPoints; // points
    std::vector<Eigen::Vector3d> vTreeLoc;  // in lu frame, height range is used to filter some points

    //ground information
    std::vector<PointType> vGroundRawPoints; //all ground points (Can Delete)
    std::vector<PointType> vGroundPlanarRawPoints; //all ground feature points under certain smoothness threshold(downsampled) (will be used in mapping)
    double gH; //ground height after leveling

    // used for the integration step
    Eigen::Vector3d r_lu_local_ini = Eigen::Vector3d::Zero();  //initial, from odometry
    Eigen::Matrix3d R_lu_local_ini = Eigen::Matrix3d::Identity();
    Eigen::Vector3d r_lu_local_updated = Eigen::Vector3d::Zero(); //refined
    Eigen::Matrix3d R_lu_local_updated = Eigen::Matrix3d::Identity();
    Eigen::Vector3d r_lu_local_ref = Eigen::Vector3d::Zero(); //if trajectory is available
    Eigen::Matrix3d R_lu_local_ref = Eigen::Matrix3d::Identity();

    // remove unnecassary info before saving to Iscan
    void cleanRedundantInfo()
    {
        num_tree = vTreeLoc.size();
        vector<PointType>().swap(vGroundRawPoints);           
        vector<PointType>().swap(vGroundPlanarRawPoints);//saved in mvGroundPointsRaw for each iscan
        vector<vector<PointType>> ().swap(vTreeRawPoints); //saved in mvTreePointsRaw for each iscan
        std::vector<Eigen::Vector3d>().swap(vTreeLoc);  // in lu frame, height range is used to filter some points
        std::vector<PointType>().swap(non_ground_raw_points_);
    }

    //non ground points
    vector<PointType> non_ground_raw_points_;

    //backup
    // EOP eop;
    // Eigen::Matrix3d R_lu_m;
    // Eigen::Vector3d r_lu_m;
    //std::vector<pair<int, int>> treePairs;         //<prev, cur> order of current
    //pcl::PointCloud<PointType>::Ptr pGroundPoints;//extracted ground plane points (feature points)
};

//class for the combined scan
class IntegratedScan
{

public:
    //general information related to scan
    double TimeInit; //t(ini) of first scan
    double TimeEnd; //t(end) of last scan
    int numScan;
    int index; //index of Iscan
    // update flag using init or refined
    enum eUpdateFlag{
        INIT=0,
        REFINED=1
    };

    bool flag_iscan_to_map = true;

    //lu frame to the mapping frame at TimeInit
    Eigen::Vector3d r_lu_m_ini;
    Eigen::Matrix3d R_lu_m_ini;
    //transformation from local -> mapping. Refined in ISCAN to map
    Eigen::Matrix3d R_local_m_ini;
    Eigen::Vector3d r_local_m_ini;
    Eigen::Matrix3d R_local_m_updated;
    Eigen::Vector3d r_local_m_updated;
    Eigen::Matrix3d R_local_m_temp; //for computing residual only
    Eigen::Vector3d r_local_m_temp;

    Eigen::Matrix3d R_local_t2_local_t1; //relative transformation between current Iscan and previosu Iscan. Used in LoopClosure() part.  if odo_from_trajectory_flag, this is provided by traj
    Eigen::Vector3d r_local_t2_local_t1;
    Eigen::Matrix3d R_local_t2_local_t1_updated; //relative transformation between current Iscan and previosu Iscan. Not used
    Eigen::Vector3d r_local_t2_local_t1_updated;

    //trajectory information. nScan + 1.
    vector<Eigen::Vector3d> v_r_local; //refined value from scan integration
    vector<Eigen::Matrix3d> v_R_local;
    vector<Eigen::Vector3d> v_r_mapping; // lu to mapping, updated after estimated ISCAN to map
    vector<Eigen::Matrix3d> v_R_mapping;
    vector<Eigen::Vector3d> v_r_mapping_ref;
    vector<Eigen::Matrix3d> v_R_mapping_ref;
    // vector<Eigen::Vector3d> v_r_global_ref;
    // vector<Eigen::Matrix3d> v_R_global_ref;

    //information of individual scans
    vector<ScanInfo> indScans; //individual scans (Can delete)
    
    //tree info
    vector<vector<PointType>> vTreePointRaw; // w.r.t each individual scans, raw
    vector<vector<PointType>> vTreePointLocal; // rilu. w.r.t the t_init of integrated scan
    vector<vector<PointType>> vTreePointMapping; // in the mapping frame
    vector<CylinderPara> vTreeParamLocal; //tree parameters w.r.t to local(t_init)
    vector<CylinderPara> vTreeParamMapping; //tree parameters w.r.t to mapping
    //vector<bool> vValidTree; //whether included in the optimization in the integration
    vector<int> vTreeStatus; //related to the map tree. 0: not matched with map tree. 1: matched with map tree

    //updated, tree points per scan
    vector<vector<vector<PointType>>> raw_tree_points;//scan id,tree id, pt id
    vector<vector<vector<PointType>>> raw_tree_points_mapping;//scan id,tree id, pt id
    vector<vector<int>> matched_map_tree_id;
    vector<int> corresponding_map_tree_ids; //correspond to vTreePointRaw after updating here

    //updated, non ground points per scan
    vector<vector<PointType>> raw_non_ground_points;//scan id, pt id
    vector<vector<PointType>> raw_non_ground_points_mapping;//scan id,tree id, pt id


    void update_vTreePointRaw();



    //ground points
    vector<PointType> vGroundPointRaw; // w.r.t each individual scans, raw
    vector<PointType> vGroundPointLocal; //
    vector<PointType> vGroundPointMapping;
    pcl::PointCloud<PointType>::Ptr pGroundPointLocalDs;//extracted ground plane points after downsampling. intensity == current Iscan index
    pcl::PointCloud<PointType>::Ptr pGroundPointMappingDs;//extracted ground plane points.
    //double groundDsSize = 0.1;

    //functions
    IntegratedScan(){};
    ~IntegratedScan(){
        pGroundPointLocalDs.reset();
        pGroundPointMappingDs.reset();
        vector<vector<PointType>>().swap(vTreePointRaw);     
        vector<vector<PointType>> ().swap(vTreePointLocal);   
        vector<vector<PointType>>().swap( vTreePointMapping); 
        vector<PointType> ().swap(vGroundPointRaw);           
        vector<PointType> ().swap(vGroundPointLocal);
        vector<vector<vector<PointType>>>().swap(raw_tree_points);
        vector<vector<vector<PointType>>>().swap(raw_tree_points_mapping);

         
    };

    void transform2Start(double downsample_distance,int option = 1); //from Raw to Local
    void updateTreeParamLocal();//for trees that are not included in the scan integration
    void computeIndividualPoseMapping(); //trajectory from local to mapping

    //compute point/parameter from local to mapping
    void computeMapGroundPoints(const eUpdateFlag flag);
    void computeMapTreePoints(const eUpdateFlag flag);
    void computeMapTreeParameter(const eUpdateFlag flag);

    //compute point from raw to mapping
    void computeRawGroundPointToMapping();
    void computeRawTreePointToMapping();
    void computeBackupTreePointToMapping();
    void computeNonGroundPointsToMapping();

    
    //test
    void show_number_PointType();
};

#endif