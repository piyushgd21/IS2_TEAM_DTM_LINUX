#ifndef _MAPPING_H_
#define _MAPPING_H_

#include "../header/utility.h"
#include "../header/CTrajectory.h"
#include "../header/optimization.h"
#include "../header/CMapTree.h"

struct RadiusInfo
{
	double radius = 0;
	double error = 0;
	bool success = false;
};

struct CloudInfo {
	//double X, Y, Z;
	//double D;   // distance from perspective center
	//uint32_t r, g, b;
	uint16_t intensity;
	double time;
	uint16_t point_source_id, return_number, no_retrun;
	int8_t scan_angle_rank;
	uint8_t classification;
	uint8_t ground_filtering_type = 0;  // 0: unclassified; 1: terrain; 2: above ground 
	double echo_range;
	uint64_t point_id;
	uint64_t origin_id;
	double gps_time;
	//int isRoadSurface_buffer = 0;
	uint8_t scan_direction_flag;
	uint16_t red; // new
	uint16_t green; // new
	uint16_t blue; // new

	int64_t tree_ID;
	double time_ratio = 0;
};

class Mapping
{   

    //control
    int Min_Num_Tree_ = 8; 
    int Max_Iter_LC = 5;
    int Maximum_Merge_Tree_LC = 5;

    //lc:
    int Min_Tree_Points_Per_Iscan_LC = 200;
    int Max_Point_Per_Tree_Iscan_LC = 1000;
    int Max_Point_Per_Tree_Raw_LC = 3000;
    
    bool Flag_Scan_Inter = true;

    bool Flag_Pose_to_Map = true;

    bool intermediate_result_flag;
    int numTotalScan = 0;
public:
	//vector<vector<pair<int, int>>>& GetmvIntegratedTreeIds() {
	//	return mvIntegratedTreeIds;
	//}
	//vector<Eigen::Vector3d>& GetmvIntegratedTree()
	//{
	//	return mvIntegratedTree;
	//}
	using IntergratedTreeIdsUnit = vector<vector<pair<int, int>>>;
	using TreeLocationUnit = vector<Eigen::Vector3d>;
	std::map<int, int> tree_id_in_iscan_to_location_index_;
	std::map<int, int> tree_id_in_iscan_to_vector_index_;
	IntergratedTreeIdsUnit integrated_tree_Ids_;
	TreeLocationUnit integrated_tree_locations_;

	std::vector<TreeLocationUnit> mvIntegratedTreeVector; //locations of the tree in local frame
	std::vector<IntergratedTreeIdsUnit> mvIntegratedTreeIdsVector;
	int scan_count = 0;
	bool scan_init = true;
	int local_scan_id = 0;
	int running_iscan_index = 0;
	int debug_for_iscan_points_number = 0;
    int last_cut_scan_id_ = 0;
	int number_integrated_scan_ = 0;
    //(TODO:chunxi modification)
    bool DetermineLoopClosureIfExist(Eigen::Vector3d scan_location,int iscan_index,const vector<Eigen::Vector3d>& current_scan_tree_locations,
    Eigen::Matrix3d& transformation,Eigen::Vector3d& translation);

    //(NEW:6.18)
    void OptimizeRawLevelLCUsingSurfaceElements();
	std::vector<RadiusInfo> CalculateMostPortionRadius();
	std::vector<Eigen::Vector3d> tree_portion_points_;
	std::vector<CloudInfo> tree_portion_points_info_;
	bool radius_stategy_ = false;


	int mNumIntegratedScan = 10;//1 sescond.
	int mNumInitialScan = 100; //10 seconds of data @@@TEST
	SettingPara mPara; //parameter from input file

private:
    //***************************************************
    //--------------- from Input ----------------//
    //SettingPara mPara; //parameter from input file
	//bool is2team_traj_input = false;
    CTrajectory *mpTraj; // reference trajectory
    bool mbTraj = false; // bool if trajectory is available
    std::ofstream fMapLog;
    bool mbInit = false; //false: not initialized , tree: general datarate
    int mnCount = 0;    //count number of integration/keyframe
    //int mNumIntegratedScan = 10;//1 sescond.
    //int mNumInitialScan = 100; //10 seconds of data @@@TEST
    int mNumScan;   //num of scans to be integrated
    double mMaxIntegratedTime;
    queue<ScanInfo> mqScan; //saving the scans from odometry

    enum eExportOption{    //setting for output option
        MAPPING_FRAME=0,
        GLOBAL_FRAME=1
    };

    //***************************************************
    //--------------- Scan Integration ----------------//
    bool mbOdoIntegration = false;  //if including odometry result in scan integration
    int mNumIntegrationScan; //number of scans in integration
    int mNumValidTree = 0;
    float MaxHeightTreePoint = 3.5; //height range: 1.5 - 4 m... for branches.. 0.5 - 3.5
    float MinHeightTreePoint = 0.5; //we extract the height range
    double MinNumTreePointPerScan = 2.0; //minimum number of tree points per scan
    double ValidInlierRatio = 0.6; //minimum inlier ratio for a valid tree in optimization
    double MinDistanceIsolatedTree = 2.0; //distance to find isolated trees

    vector<ScanInfo> mvTempScan; //scan to be combined
    //tree information
    vector<vector<PointType>> mvTreePointsRaw; // rIlu. Filtered tree points. Will be updated in the optimization
    vector<Eigen::Vector3d> mvIntegratedTree; //locations of the tree in local frame
    vector<vector<pair<int, int>>> mvIntegratedTreeIds; //scan id + tree id of each tree. Corresponding to scanInfo.treeRawPoints
    vector<Eigen::Vector3f> mvTreeCentersLocal; //center_local: first computed using init T_lu_local
    vector<CylinderPara> mvTreeParamsLocal; //tree parameters: first initialized using center_local, r = 0.0. The parameters for valid trees will be updated in scan integration
    vector<bool> mvValidTrees; //whether included in the optimization, may deactivate in the optimization

    //ground information
    int mNumPlanarPatch = 0;
    double MinNumGroundPoints = 2.0; //minimum number of tree points per scan
    std::vector<PointType> mvGroundPointsRaw; // r_I_lu
    std::vector<PointType> mvGroundPointsLocal; //r_I_local
    std::vector<pair<double, double>> mvGridCenters; //center of each grid in local frame
    vector<vector<int>> mvGridPointIndex; //index of points in mvGroundPointsRaw
    vector<Eigen::Vector4d> mvPlaneParamLocal;

    //used for intepolation of the first scan for exporting the results. the pose for t_ini(0) in integration.
    Eigen::Vector3d mr_lu_local_t0 = Eigen::Vector3d::Zero();
    Eigen::Matrix3d mR_lu_local_t0 = Eigen::Matrix3d::Identity();
    //result from the integration optimization, numScan + 1 poses 
    vector<Eigen::Vector3d> mv_r_lu_local_ini; //from odometry reulst
    vector<Eigen::Matrix3d> mv_R_lu_local_ini;
    vector<Eigen::Vector3d> mv_r_lu_local; //refined in the odometry
    vector<Eigen::Matrix3d> mv_R_lu_local;
    vector<Eigen::Vector3d> mv_r_lu_local_ref; //reference trajectory
    vector<Eigen::Matrix3d> mv_R_lu_local_ref;

    //related to optimization
    bool estimate_tree_normal_ = false; //default
    bool estimate_tree_radius_ = false;
    double tree_outlier_multiplier_ = 1.0;

    //functions
    void integrateScan(); //main function for scan integration
    void derivePoseLocal();
    void integrateTree();
    void integrateGroundPoints();
    void optimizeIntegratedScan(bool fix_traj_flag = false);
    void UpdateMergeTree();
    void broadcastIScan();
    void resetIntegration();

    //test
    void backupTreePoints();
    vector<vector<vector<PointType>>> raw_tree_points_per_iscan;

    Eigen::Quaterniond computeAverageLevelR();

    //***************************************************
    //--------------- Iscan to Map --------------------//
    IntegratedScan *pCurrentIScan;    //all info for the current Iscan
    vector<pair<int, int>> mvPairs; // map tree id & scan tree id
    vector<CylinderPara> mvValidTreeParamsUpdated; //corresponding to previous
    int iter_count = 0;
    //initilization
    void Localize();

    void LocalizeDTMOnly();

    //match tree features to map
    void matchTreeScantoMap(const vector<Eigen::Vector3d> *mappingTree, const vector<Eigen::Vector3d> *preScanTree, const vector<double> *threshold, std::vector<pair<int, int>> &pairs);
    void matchTreeScantoMap(const vector<Eigen::Vector3d> *preScanTree, const vector<double> *threshold, std::vector<pair<int, int>> &pairs); //defaultly match to mapTree
    void compute2dTransScantoMap(const vector<Eigen::Vector3d> *mappingTree, const vector<Eigen::Vector3d> *scanTree, const vector<pair<int, int>> *pairs, std::vector<bool> &validTreePair,
                                 Eigen::Matrix3d &R, Eigen::Vector3d &t, Eigen::Vector4d initResult);
    //compute the pose of the integrated scan to mapping
    void computePoseToMap();
    bool optimizeIScantoMap();
    void addIscanTreetoMapTree();
    void addFeaturetoMap();
    void computeRefTraj();

    void ExportFeatureIscanToMap(const std::string outPass);



    //***************************************************
    //--------------- Map information --------------------//
    vector<MapTree *> mpMapTree;   //tree features
    pcl::PointCloud<PointType>::Ptr mpMapGroundPoint; //ground points
    pcl::octree::OctreePointCloudSearch<PointType> octreeGroundPointsFromMap;
    //function
    void computeMapTreeInfo();
    double computeTreeToTreeResidual(int refId, int candId);
    double computeTreeToTreeResidual(int refId, int candId, vector<double> &residual_distribution); 

    int mergeMapTree();
    void CountMapTreeType();

    //--------------- Global map --------------------//
    bool global_map_flag_ = false;
    bool DTM_only_flag_ =false;
    double ref_ground_std_ = 0.02;
    Eigen::Vector3d global_const_shift_ = Eigen::Vector3d::Zero();
    bool LoadGlobalGroundMap();
    bool LoadGlobalTreeMap();


    //***************************************************
    //--------------- Loop closure --------------------//
    int iscan_index_prev_loop_closure_ = -1;
    int iscan_index_start_ = 0; //refine the pose
    int iscan_index_end_ = 0;   //refine the pose
    int iscan_index_feature_start_ = 0; //include the feature for planar
    int count_lc = 0;
    bool raw_point_flag_ = false;//false;
    //planar patch
    vector<vector<Eigen::Vector3d>> point_ref_per_grid_;
    vector<vector<pair<int,int>>> point_index_per_grid_;
    vector<Eigen::Vector4d> plane_params_;
    //tree info
    vector<int> used_trees_;
    string prefix;
    void ConductIntermediateLC();
    void ConductFinalLC(int num_iterations);
    void LoopClosure();
    void DeriveIndexLC();
    void DerivePlanarPatchLC();
    void OptimizeLC();
    void OptimizeRawLevelLC();
    void UpdatePosePointsLC();
    void UpdateMapTreeLC();
    void ExportFeatureLC(const std::string outPass);
    void ExportTrajectoryLC(const std::string outPass);
    void ExportTrajectoryRefLC(const std::string outPass);


    //***************************************************
    //--------------- Establish Tree features again --------------------//
    void FindRawTreePointsPerMapTree(double distance_threshold = 0.3); 
    void Update_Map_Tree_Obs(); //update the observation of each map tree

    void FindRawPointsPerMapTree(double distance_threshold = 0.3);


    //***************************************************
    //--------------- Evaluation & export --------------------//
    Eigen::Vector3d r_m_global = Eigen::Vector3d::Zero();
    Eigen::Matrix3d R_m_global = Eigen::Matrix3d::Identity();
    Eigen::Vector3d r_m_global_ref_traj = Eigen::Vector3d::Zero();
    Eigen::Matrix3d R_m_global_ref_traj = Eigen::Matrix3d::Identity();


    //scan integration
    void exportIntegratedTrajectory(const std::string outPass);
    void exportIntegratedOptFeatures(const std::string outPass, int flag, bool bGroundPoint);
    void exportIntegratedFeatures(const std::string outPass, int flag, bool bGroundPoint);
    void exportIscanFeatureBroadcast(const std::string outPass, int ds_flag);

    //for mapping
    void exportMapTreeParam(const std::string outPass, eExportOption exportOption = MAPPING_FRAME);
    void exportMapTreeModel(const std::string outPass, eExportOption exportOption = MAPPING_FRAME); // Line representing the tree model
    void exportFeatureMapping(const std::string outPass, eExportOption exportOption = MAPPING_FRAME);
    void exportFeatureMappingRef(const std::string outPass, eExportOption exportOption = MAPPING_FRAME);
    void exportGroundMap(const std::string outPass, eExportOption exportOption = MAPPING_FRAME);
    void exportTrajGlobalCheck(const std::string outPass);
    void exportTrajGlobal(const std::string outPass,pcl::PointCloud<pcl::PointXYZ>* slam_trajectory =nullptr, pcl::PointCloud<pcl::PointXYZ>* gnss_trajectory =nullptr);
	void UpdateMappingToGlobalTransformation(const pcl::PointCloud<pcl::PointXYZ>& slam_trajectory, const pcl::PointCloud<pcl::PointXYZ>& gnss_trajectory);
    void ExportFeatureMappingLas(const std::string outPass, eExportOption exportOption = MAPPING_FRAME);
    void ExportBackupTreeFeatureMappingLas(const std::string outPass, eExportOption exportOption = MAPPING_FRAME);
	Eigen::Vector3d offset_trans_traj_;
    //***************************************************
public:
    //initialization
    Mapping(SettingPara setting, CTrajectory *pTraj = NULL) : mPara(setting), mpTraj(pTraj), mpMapGroundPoint(new pcl::PointCloud<PointType>()), octreeGroundPointsFromMap(0.2)
    {
        //create Map Log file
#ifdef EXPORT_LOG
        std::string outPass = output_folder + "_Map_LOG.txt";
        fMapLog.open(outPass, std::ifstream::out);
        fMapLog << fixed << std::setprecision(4);

#endif
        if (mpTraj)
            mbTraj = true;

        if (mPara.bOdometry)
            mbOdoIntegration = true;

        if(mPara.global_map_flag)
        {
            global_map_flag_ = true;
            ref_ground_std_ = 0.02;
        }

        if(mPara.global_DTM_only)
        {
            DTM_only_flag_ = true;
            ref_ground_std_ = 0.0005;//0.005
        }

        intermediate_result_flag = mPara.export_intermediate_result_mapping;
        mNumIntegratedScan = mPara.numIntegratedScan;
        mNumInitialScan = mPara.numInitIntegratedScan;
        MinNumTreePointPerScan = mPara.minNumTreePointPerScan;
        Min_Num_Tree_ = mPara.minimum_tree_iscan;

        fMapLog << "Minimum number of extracted tree/matches for iscan:" << Min_Num_Tree_ <<endl;

        output_folder_iscan = output_folder + "scan_integration/";
        boost::filesystem::create_directories(output_folder_iscan);

        if(intermediate_result_flag)
        {
            output_folder_iscan_map = output_folder + "scan_to_map/";
            boost::filesystem::create_directories(output_folder_iscan_map);
        }
    };
    ~Mapping(){};

    //--------------- Map information --------------------//
    vector<IntegratedScan *> mvpIScans; //all Iscan in the mapping
    void buildMap();

    //--------------- Receive scan from odo -------------//
    void insertScan(ScanInfo scan)
    {
        numTotalScan++;
        mqScan.push(scan);
    }


    //************************* check *****************************
    void exportFeatures(const std::string outPass, int flag, bool bpoint);
    void exportFeaturesDis(const std::string outPass, int flag, bool bpoint);
    void exportFeaturePoints(const std::string outPass, int flag, bool bpoint);
    void exportFeaturePointsDis(const std::string outPass, int flag, bool bpoint);
    void exportTreeModel(const std::string outPass);
    void exportIntegratedTreeIni(const std::string outPass);
    void exportIntegratedTree(const std::string outPass);
    void exportIntegratedTreeRef(const std::string outPass);
    void exportIntegratedTreeCenter(const std::string outPass);
    void exportFinalTree(const std::string outPass);
    void exportGroundPoints(const std::string outPass);

    void freeSpace()
    {
        for (int i = 0; i < mvpIScans.size(); i++)
        {
        delete mvpIScans[i];

        }
    }

};

#endif