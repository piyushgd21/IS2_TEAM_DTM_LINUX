/**
 * @file featureExtraction.cpp
 * @brief Core logic for extracting LiDAR features (corners, surfaces) from point clouds.
 * 
 * This file handles sharp corner detection, flat surface detection, downsampling,
 * filtering, and exporting features for SLAM or mapping purposes.
 */
#include "../header/featureExtraction.h"
#include "../header/optimization.h"
#include <filesystem>

// ========== Global Variables ==========

// Pointer to curvature buffer used in sorting
float *mpCurv;

// Sorting comparator for curvature values
bool comp(int i, int j) { return (mpCurv[i] < mpCurv[j]); }

// Global scan and file counters
int ScanCount = 0;
int SubFile = 1;

// Buffers for rotation/translation states from previous scans
vector<Eigen::Matrix3d> lastR;
vector<Eigen::Vector3d> lastr;
vector<Eigen::Matrix3d> lastRefR;
vector<Eigen::Vector3d> lastRefr;
vector<double> refT;

// Flag and value for continuity in ground height tracking
bool previous_ground_info_flag = false;
float ground_height_previous;
Eigen::Matrix3f R_lu_lup_previous;

// Counter for failed odometry tracking attempts
int successive_odometry_failure_count  = 0;

/* Load data from current scan, fill mpLaserCloud
Inputs:
    1. file name for point cloud, 2. scan ID, 3. pointer to tracked scan (Null or tracked), 4.ID for tracking sub maps, 5. initTime, 6, duration of a scan,
    7. setting struct, 8. pointer to reference trajectory (Null), 9. pointer to map thread (Null)
 */
// ========== Destructor ==========

/**
 * @brief Destructor for LidarScan - releases memory and resets global buffers.
 */
LidarScan::~LidarScan() {
	mpSurfPointsFlatGround.reset();
	pSurfPointsLessFlatGround_ds_.reset();
	pSurfPointsLessFlatGround_.reset();
	mpCornerPointsSharp.reset();
	mpCornerPointsLessSharp.reset();
	mpSurfPointsFlat.reset();
	mpSurfPointsLessFlat.reset();
	refT.clear();
	lastR.clear();
	lastr.clear();
	lastRefR.clear();
	lastRefr.clear();
	previous_ground_info_flag = false;
	ScanCount = 0;
	SubFile = 1;
	mpCurv = nullptr;
	successive_odometry_failure_count = 0;
};

/**
 * @brief Constructor for LidarScan that handles scan data loading, initialization, and odometry estimation.
 *
 * @param lidar_data_path       Path to raw binary scan file.
 * @param index                 Scan ID.
 * @param pPrev                 Pointer to previous scan (nullptr for first).
 * @param partID                Submap ID.
 * @param timetag               Scan start time.
 * @param duration              Scan duration.
 * @param setting               User-specified SLAM and extraction settings.
 * @param single_scan_tree_file Optional: Tree feature file.
 * @param single_scan_ground_file Optional: Ground feature file.
 * @param points_number         Total number of LiDAR points.
 * @param pTraj                 Pointer to ground truth trajectory.
 * @param pMap                  Pointer to mapping thread logic.
 */
LidarScan::LidarScan(const std::string lidar_data_path, int index, LidarScan *pPrev, int partID, double timetag, double duration, SettingPara setting,
	std::string single_scan_tree_file, std::string single_scan_ground_file, int points_number, CTrajectory *pTraj, Mapping *pMap)
    : id(index), mPara(setting), mpPrev(pPrev), mPartID(partID), mTimeInit(timetag), mScanDuration(duration),
	single_scan_points_number_(points_number), mpTraj(pTraj), mpMap(pMap)
{
    if (!mpPrev)
        mbInit = true;

    if (mpTraj)
        mbTraj = true;

    if (mpMap)
        mbMap = true;
    
    flag_export_intermediate_result = setting.export_intermediate_result_odometry;

    mTimeEnd = mTimeInit + mScanDuration;

    // Trajectory interpolation for initial pose
    int gnss_idx = -1;
    if (!(mpTraj->eopInterpolation(timetag / 1000.0, -1, gnss_idx, r_lu_m_t_ini_, R_lu_m_t_ini_)))
    {
        throw std::runtime_error("Wrong Time");
    }

#ifdef EXPORT_LOG
    fLog << "-----------------------------" << endl;
    fLog << "Scan " << id << "(" << mTimeInit << "-" << mTimeEnd << ")"
         << " w.r.t ";
    if (mpPrev)
        fLog << mpPrev->id << "(" << mpPrev->mTimeInit << "-" << mpPrev->mTimeEnd << ")" << endl;
    else
        fLog << "Null" << endl;
#endif

    /*Compute expected translation/rotation using the results from previous scan*/
    if (mpPrev)
    {
        // t_end - t_end_prev
        mTimeTracking = mTimeEnd - mpPrev->mTimeEnd;
        float s = mTimeTracking / mpPrev->mTimeTracking;

        //(TODO:slerp is only for interpolation)
        Eigen::Quaterniond qPrev(mpPrev->R_lut2_lut1);

        // Estimate rotation and translation via SLERP
        q_ini = Eigen::Quaterniond::Identity().slerp(double(s), qPrev);
        t_ini = double(s) * mpPrev->r_lut2_lut1;

#ifdef EXPORT_LOG
        fLog << "Tracking time: " << mTimeTracking << "v.s." << mScanDuration << endl;
        fLog << "Initial from last " << t_ini.transpose() << "\t" << rad2deg(Find_Rotation(q_ini.toRotationMatrix())).transpose() << endl;
#endif
    }
    else
    {
        mTimeTracking = mScanDuration;
    }

    cout << endl;
    cout << "------------ Scan " << id << " --------------" << endl;

    //*********** 1. Load data ****************
	//(MODI)
    TicToc tLoading;
    //load_raw_data(lidar_data_path);
    TimeLoading += tLoading.toc();

#ifdef TEST_PURPOSE
    std::string outPass(output_folder + to_string(index) + "_info.txt");
    fInfo.open(outPass, std::ifstream::out);
    fInfo << fixed << std::setprecision(4);
#endif

    //*********** 2. Processing ****************
    //------ 2.1 Initialize & segmentation
    TicToc tSeg;
    init();
	//(MODI)
    //computeAttribute();
    TimeSegment += tSeg.toc();

    //------ 2.2 Feature based odometry -------------------
	//(MODI)
	valid_feature_based_odo_flag_ = true;
    //if (mbFeatureBased)
    //{
    //    TicToc tFb;
    //    if(conduct_feature_based_odo())
    //    {
    //        valid_feature_based_odo_flag_ = true;
    //    }
    //    else
    //    {
    //        //firsr scan and valid features are extracted
    //        if(mbInit && valid_feature_flag_)
    //        {
    //            valid_feature_based_odo_flag_ = true;
    //        }
    //        else
    //            valid_feature_based_odo_flag_ = false;
    //    }
    //    TimeFeature += tFb.toc();
    //}

    //------ 2.3 Point based odometry -------------------
	//(MODI)
	mbPointBased = false;
    if (mbPointBased)
    {
        // if the point.feature based will conduct simultaneously or the feature-based fail
        if (mbSimultaneously || !valid_feature_based_odo_flag_)
        {
            TicToc tPb;

            if (conduct_point_based_odo())
            {
                valid_point_based_odo_flag_ = true;
            }
            else
            {
                // firsr scan and valid points are extracted
                if (mbInit && valid_point_flag_)
                {
                    valid_point_based_odo_flag_ = true;
                }
                else
                    valid_feature_based_odo_flag_ = false;
            }
            TimePoint += tPb.toc();
        }
    }


    //export extracted features in the feature-based odometry
    if(false/*flag_export_intermediate_result*/)
    {
        std::string out_point_attribute(output_folder_odometry + to_string(id) + "_point_attribute.txt");
        export_scan_point(out_point_attribute);
        std::string out_ground_cluster(output_folder_odometry + to_string(id) + "_segments.txt");
        export_scan_segment(out_ground_cluster);

        std::string out_segment_leveled(output_folder_odometry + to_string(id) + "_segments_final_leveled.txt");
        export_scan_segment(out_segment_leveled, true);
        //system("read -p 'Press Enter to continue...' var");
    }


    //------ 3.1 Final result for the odometry-------------------
    // flag for success tracking or not
    valid_odo_flag_ = valid_feature_based_odo_flag_ || valid_point_based_odo_flag_;


    // compute trajectory in mapping
    if (valid_odo_flag_)
    {
#ifdef EXPORT_LOG
        fLog << "-- Valid! Final odometry result from ";
        if (valid_feature_based_odo_flag_)
            fLog << "feature-based" << endl;
        else
            fLog << "point-based" << endl;
        fLog << "\tEstimated" << r_lut2_lut1.transpose() << "\t" << rad2deg(Find_Rotation(R_lut2_lut1)).transpose() << endl;
        if (mpTraj)
        {
            fLog << "\tReference" << r_lut2_lut1_ref.transpose() << "\t" << rad2deg(Find_Rotation(R_lut2_lut1_ref)).transpose() << endl;
        }
#endif
		//(MODI)
        //computeTransformation();
        successive_odometry_failure_count = 0;
    }
    else
    {
        successive_odometry_failure_count ++;
        mPartID = -1;

        if(mpPrev)
        {
#ifdef EXPORT_LOG
            fLog << "Error: the odometry is failed, use odometry result from previous scan. Current successive scans: " << successive_odometry_failure_count << endl;
            f_odometry_debug << id << "\tThe odometry is failed, use odometry result from previous scan. Current successive scans: " << successive_odometry_failure_count << endl;
#endif      
            //use the odometry result from previous scan
            R_lut2_lut1 = mpPrev->R_lut2_lut1;
            r_lut2_lut1 = mpPrev->r_lut2_lut1;
            valid_odo_flag_ = true;
        }
        else
            throw std::runtime_error("Current implementation does not support fail scan in mapping thread");

    }

    //------ 3.2 Send results to map thread -------------------
    if (mbMap)
    {
#ifdef EXPORT_LOG
        fLog << "Send scan to map: " << endl;
#endif
        //AddScantoMap();
		//std::string single_file, ground_file;
		AddScanToMap(single_scan_tree_file, single_scan_ground_file, index);
    }

    if(flag_export_intermediate_result&&valid_feature_flag_)
    {
        std::string out_final_ground_planar(output_folder_odometry + to_string(id) + "_ground_planar.txt");
        export_final_planar_points(out_final_ground_planar);    
        //system("read -p 'Press Enter to continue...' var");
    }


    //------ 3.3 Export residuals for the odometry-------------------
    if (mpTraj && !mbInit)
    {
        // residuals in the odometry
#ifdef EXPORT_TRAJECTORY
        if (valid_feature_based_odo_flag_)
        {
            Eigen::Matrix3d rot_dif = R_lut2_lut1_gt * R_lut2_lut1_ref.inverse();
            Eigen::Vector3d pos_dif = r_lut2_lut1_gt - r_lut2_lut1_ref;
            Eigen::Vector3d angles_dif = Find_Rotation(rot_dif);
            fTrajectoryResFb << id << "\t" << std::setw(8) << pos_dif(0) << "\t" << std::setw(8) << pos_dif(1) << "\t" << std::setw(8) << pos_dif(2) << "\t"
                             << std::setw(8) << rad2deg(angles_dif(0)) << "\t" << std::setw(8) << rad2deg(angles_dif(1)) << "\t" << std::setw(8) << rad2deg(angles_dif(2)) << "\t" << mPartID << endl;
        }

        if (valid_point_based_odo_flag_)
        {
            Eigen::Matrix3d rot_dif = R_lut2_lut1_pb * R_lut2_lut1_ref.inverse();
            Eigen::Vector3d pos_dif = r_lut2_lut1_pb - r_lut2_lut1_ref;
            Eigen::Vector3d angles_dif = Find_Rotation(rot_dif);
            fTrajectoryResPb << id << "\t" << std::setw(8) << pos_dif(0) << "\t" << std::setw(8) << pos_dif(1) << "\t" << std::setw(8) << pos_dif(2) << "\t"
                             << std::setw(8) << rad2deg(angles_dif(0)) << "\t" << std::setw(8) << rad2deg(angles_dif(1)) << "\t" << std::setw(8) << rad2deg(angles_dif(2)) << "\t" << mPartID << endl;
        }
#endif
    }

}

/*Main function for the feature based odometry */
/**
 * @brief Performs feature-based odometry estimation.
 *
 * This function handles both leveled and non-leveled cases. It extracts ground and tree features,
 * estimates relative transformation using geometric constraints, and matches with previous scan.
 *
 * @return true if successful, false otherwise.
 */

bool LidarScan::conduct_feature_based_odo()
{
    fLog << "-- Conduct feature-based odometry" <<endl;
    // Case 1. Estimate leveling rotation given ground is flat and conduct matching in 2d
    if (!mbLevel)
    {
        // ------------------ 1. extract ground information --------------------
        //Basic assumption: Ground is flat and relatively leveled, so the leveling rotation and ground height can be use
        TicToc tGround;
        // 1.1 find ground clusters to get initial ground representation
        if (!extractGroundCluster(false))
        {
            //if fail in extraction, use the info from previous
            if(previous_ground_info_flag)
            {
                fetch_previous_ground_info();   
            }
            else if(mbTraj)
            {
                // ground info from mounting parameter
                fLog << "Get ground info from trajectory & system info" <<endl;
                R_lu_lup = R_lu_b_ref.cast<float>();
                mGroundH = ground_height_ref;
            }
            else
            {
                return false;
            }
        }
        
        // export extracted features in the feature-based odometry
        if (flag_export_intermediate_result)
        {
            std::string out_segment_leveled(output_folder_odometry + to_string(id) + "_segments_ground_cluster_leveled.txt");
            export_scan_segment(out_segment_leveled, true);
            // system("read -p 'Press Enter to continue...' var");
        }

        // 1.2 compute the centroid of segment in the leveled lu
        computeLevelCenters();

        // 1.3 extract ground segments and update plane parameters
        if (!extract_ground_segments(true))
        {
            return false;
        }

        // 1.4 for the mapping part, extract ground planar points
        if(!extract_planar_points_from_ground_segment())
        {
            return false;
        }

        //save the previous ground info
        update_previous_ground_info();
        valid_ground_info_flag_ = true;
        Tground += tGround.toc();

        // ------------------ 2. extract tree clusters ------------------------
        //only if the ground is extracted correctly -> leveling rotation is correct
        TicToc tTree;
        if(!extractTreeCluster())
        {
            return false;
        }
        valid_feature_flag_ = true;
        Ttree += tTree.toc();

        //------------------- 3. optimization ---------------------------------
        //conduct optimization in case previous scan has valid features
        if (mpPrev && mpPrev->valid_feature_flag_)
        {
            // 3.1 get pose prediction from the init value and R_level
            TicToc tOpt;
            PredictLevelInit();
            
            // 3.2 match trees using 2D
            if(!matchTrees2d(true))
            {
                return false;
            }

            // 3.3 conduct pose estimation
            if(!compute2dSimilarityLeveled())
            {
                return false;
            }
            Topt += tOpt.toc();
        }
        else
        {
            if (mpPrev) //previous scan is available but not enough feature
            {
#ifdef EXPORT_LOG
                fLog << "Error: No enough feature from previous steps for matching  " << endl;
                f_odometry_debug << id << "\tNo enough feature from previous steps for matching  " << endl;
#endif
            }
            return false;
        }

        //feature based odometry is successful
        return true;
    }

    // Case 2. Cannot estimate leveling rotation, using nominal leveling rotation
    // TODO: modify 
    else
    {
        if(!mbTraj)
            throw std::runtime_error("Mounting parameters are not available");

        // ------------------ 1. extract ground information --------------------
        //Basic assumption: leveling rotation is valid, the ground cannot be modeled as a flat horizontal terrain (ground height no valid)
        // leveling info from mounting parameter
        R_lu_lup = R_lu_b_ref.cast<float>();
        TicToc tGround;

        // 1.1 compute the centroid of segment in the leveled lu
        computeLevelCenters();

        // 1.2 extract the bottom layer and extract corresponding ground lidar points
        if(!extract_ground_segments(false))
        {
            return false;
        }
        
        // 1.3 extract ground planar points for both odometry and mapping thread
        if(!extract_planar_points_from_ground_segment())
        {
            return false;
        }
        
        valid_ground_info_flag_ = true;
        Tground += tGround.toc();
        
        // ------------------ 2. extract tree clusters ------------------------
        TicToc tTree;
        if (!extractTreeCluster())
        {
            return false;
        }
        valid_feature_flag_ = true;
        Ttree += tTree.toc();

        //------------------- 3. optimization ---------------------------------
        //conduct optimization in case previous scan has valid features
        if (mpPrev && mpPrev->valid_feature_flag_)
        {
            // 3.1 get pose prediction from the init value and R_level
            TicToc tOpt;
            PredictLevelInit();

            // 3.2 match trees using 3D
            if(!matchTrees3d(true))
            {
                return false;
            }

            // 3.3 conduct pose estimation
            if(!computeOdometryLeveled())
            {
                return false;
            }
            Topt += tOpt.toc();
        }
        else
        {
            if (mpPrev) // previous scan is available but not enough feature
            {
#ifdef EXPORT_LOG
                fLog << "Error: No enough feature from previous steps for matching  " << endl;
                f_odometry_debug << id << "\tNo enough feature from previous steps for matching  " << endl;
#endif
            }
            return false;
        }

        // feature based odometry is successful
        return true;
    }
}


/*Main function for the point based odometry */
/**
 * @brief Performs point-based odometry estimation using planar/edge features.
 * 
 * Used when feature-based odometry fails or is not reliable.
 * 
 * @return true if odometry is computed successfully, false otherwise.
 */
bool LidarScan::conduct_point_based_odo()
{
    fLog << "-- Conduct Point-based odometry:" << endl;
    
    // 1. Extract points from scan
    TicToc tPointExtraction;
    if(!extract_edge_planar_points())
    {
        return false;
    }
    
    // 2. Check if enough points are extracted from previous scan
    if (mpPrev)
    {   
        //if point has been extracted from previous scan
        if(!mpPrev->point_extraction_flag_)
        {
            //extract 
            fLog << "Extract points from previous scan:" <<endl;
            mpPrev->extract_edge_planar_points();
        }

        if(!mpPrev->valid_point_flag_)
        {
#ifdef EXPORT_LOG
            fLog << "Error: No enough point from previous scans for odometry  " << endl;
            f_odometry_debug << id << "\tError: No enough point from previous scans for odometry   " << endl;
#endif
            return false;
        }
    }
    else
    {
        return false;
    }
    Tp_extraction += tPointExtraction.toc();

    // 3. compute odometry 
    TicToc tOpt;
    if(!compute_point_based_odometry())
    {
        return false;
    }
    
    Tp_opt += tOpt.toc();
    return true;

}


/* Load binary data for the raw scan*/
/**
 * @brief Loads raw LiDAR data from binary .bin file into mpLaserCloud.
 *
 * File must contain points in [x, y, z, time, (optional intensity)] format.
 *
 * @param pass Path to the binary file.
 */
void LidarScan::load_raw_data(const string pass)
{
    std::ifstream lidar_data_file(pass, std::ifstream::in | std::ifstream::binary);
    if (!lidar_data_file)
        throw std::runtime_error("Wrong file name");

    lidar_data_file.seekg(0, std::ios::end);
    const size_t num_elements = lidar_data_file.tellg() / sizeof(float);
    lidar_data_file.seekg(0, std::ios::beg);
    cout << num_elements << endl;
    std::vector<float> lidar_data_buffer(num_elements);
    lidar_data_file.read(reinterpret_cast<char *>(&lidar_data_buffer[0]), num_elements * sizeof(float));
    mpLaserCloud.reset(new pcl::PointCloud<PointType>());

    // save file to mpLaserCloud
    int numInvalidPoint = 0;
    int num_float_per_point = mPara.intensity_flag ? 5 : 4;

    for (std::size_t i = 0; i < lidar_data_buffer.size(); i += num_float_per_point)
    {
        PointType point;
        point.x = lidar_data_buffer[i];
        point.y = lidar_data_buffer[i + 1];
        point.z = lidar_data_buffer[i + 2];
        float timeRatio = lidar_data_buffer[i + 3];

        // current time: t = timeInit + timeRatio*mScanDuration.
        // ratio w.r.t odometry: (t - mprev->TimeEnd)/ mTimeTracking . Theoretically,  timeInit == mprev->TimeEnd,  mTimeTracking = mScanDuration
        point.intensity = 1.0 - (1.0 - timeRatio) * mScanDuration / mTimeTracking;
        if (point.intensity > 1.0)
        {
            throw std::runtime_error("Wrong time ratio for lidar point");
        }
        mpLaserCloud->push_back(point);
    }

    size_t numPoints = mpLaserCloud->size();
    std::vector<int> indices;
    pcl::removeNaNFromPointCloud(*mpLaserCloud, *mpLaserCloud, indices);
    if (numPoints != mpLaserCloud->size())
        throw std::runtime_error("Nan in input");
}

/* Add tree to map*/
// void LidarScan::AddTreetoMap()
// {
//     if (!mbTrackFeature)
//         return;

//     vector<bool> vecValid(mvClusterTree.size(), false);
//     for (int nPair = 0; nPair < mVecTreePairs.size(); nPair++)
//     {
//         int PrevId = mVecTreePairs[nPair].first;
//         int CurId = mVecTreePairs[nPair].second;
//         vecValid[CurId] = true;
//         mpMap->addTree(mpPrev->mvClusterTree[PrevId].featureId, id, CurId, mvClusterTree[CurId].mvRawPoints);
//         mvClusterTree[CurId].featureId = mpPrev->mvClusterTree[PrevId].featureId;
//     }

//     for (int nCur = 0; nCur < mvClusterTree.size(); nCur++)
//     {
//         if (vecValid[nCur])
//             continue;

//         mvClusterTree[nCur].featureId = mpMap->initTree(id, nCur, mvClusterTree[nCur].mvRawPoints);
//     }
// }


/* Add tree to map*/

/**
 * @brief Adds current scan’s ground and tree features to the map module.
 * 
 * Prepares all necessary information like tree clusters, ground points, planar surfaces
 * and trajectory estimates, and inserts them into `mpMap`.
 */

void LidarScan::AddScantoMap()
{
    vector<Eigen::Vector3d> vTreeLoc;
    vector<vector<PointType>> vTreeRawPoints;
    vector<PointType> vGroundRawPoints; //all points in ground segments
    vector<PointType> vGroundPlanarPointsDs; //all ground planar points
    //vector<PointType> vGroundPlanarPoints;

    // save tree info
    if (valid_feature_flag_)
    {

        for (int nT = 0; nT < mvClusterTree.size(); nT++)
        {

            // remove clusters with large height information in the leveled coordinates
            if (mvClusterTree[nT].vCenterTrans(2) < 10.0)
            {
                vTreeLoc.push_back((R_lu_lup.inverse() * mvClusterTree[nT].vCenterTrans).cast<double>()); // in lidar unit frame
                vTreeRawPoints.push_back(mvClusterTree[nT].mvRawPoints);
            }
        }
        fLog << "\tNumber of tree clusters " << vTreeLoc.size() << endl;
        // save ground info
        for (int nG = 0; nG < ground_cluster_.size(); nG++)
        {
            vGroundRawPoints.insert(vGroundRawPoints.end(), ground_cluster_[nG].mvRawPoints.begin(), ground_cluster_[nG].mvRawPoints.end());
        }

        fLog << "\tNumber of points in ground segments (vGroundRawPoints) " << vGroundRawPoints.size() << endl;
        fLog << "\tNumber of ground planar points (vGroundPlanarPointsDs) " << pSurfPointsLessFlatGround_ds_->size() << endl;

        for (int nGroundP = 0; nGroundP < pSurfPointsLessFlatGround_ds_->size(); nGroundP++)
        {
            vGroundPlanarPointsDs.push_back(pSurfPointsLessFlatGround_ds_->points[nGroundP]);
        }
    }
    else
    {
        fLog << "\tFeature extraction is failed" <<endl;
    }
    // for(int nGroundP = 0; nGroundP < mpSurfPointsLessFlatGround->size(); nGroundP++)
    // {
    //     vGroundPlanarPoints.push_back(mpSurfPointsLessFlatGroundNonDs->points[nGroundP]);
    // }

    // std::string output3(output_folder + to_string(id) + "ground_points.txt");
    // std::ofstream fSegFile(output3, std::ifstream::out);
    // fSegFile << fixed << std::setprecision(8);

    // for (int i = 0; i < vGroundRawPoints.size(); i++)
    // {
    //     fSegFile << 0 << "\t" << vGroundRawPoints[i].x << "\t" << vGroundRawPoints[i].y << "\t" << vGroundRawPoints[i].z << "\t"
    //              << vGroundRawPoints[i].intensity << endl;
    // }
    // for (int i = 0; i < mpSurfPointsLessFlatGround->size(); i++)
    // {
    //     fSegFile << 1 << "\t" << mpSurfPointsLessFlatGround->points[i].x << "\t" << mpSurfPointsLessFlatGround->points[i].y << 
    //     "\t" << mpSurfPointsLessFlatGround->points[i].z << "\t"
    //              << mpSurfPointsLessFlatGround->points[i].intensity << endl;
    // }
    //     for (int i = 2; i < mpSurfPointsLessFlatGroundNonDs->size(); i++)
    // {
    //     fSegFile << 2 << "\t" << mpSurfPointsLessFlatGroundNonDs->points[i].x << "\t" << mpSurfPointsLessFlatGroundNonDs->points[i].y << 
    //     "\t" << mpSurfPointsLessFlatGroundNonDs->points[i].z << "\t"
    //              << mpSurfPointsLessFlatGroundNonDs->points[i].intensity << endl;
    // }

    // fSegFile.close();
    // system("read -p 'Press Enter to continue...' var");

	ScanInfo tem;
	tem.scanID = id; tem.partID = mPartID; tem.mTimeInit = mTimeInit; tem.mTimeEnd = mTimeEnd; tem.mTimeTracking = mTimeTracking;
	tem.R_lu_lup = R_lu_lup.cast<double>();
    tem.valid_feature_flag = valid_feature_flag_;

    //if use trajectory as odometry result
    if (mPara.odo_from_trajectory_flag)
    {
        if (!mpTraj)
        {
            throw std::runtime_error("Cannot use trajectory information when it is unavailable ");
        }
        else
        {
            tem.R_lut2_lut1 = R_lut2_lut1_ref;
            tem.r_lut2_lut1 = r_lut2_lut1_ref;
        }
    }
    else
    {
        tem.R_lut2_lut1 = R_lut2_lut1;
        tem.r_lut2_lut1 = r_lut2_lut1;
    }

	tem.vTreeRawPoints = vTreeRawPoints; tem.vTreeLoc = vTreeLoc; tem.vGroundRawPoints = vGroundRawPoints;
	tem.vGroundPlanarRawPoints = vGroundPlanarPointsDs; tem.gH = mGroundH;

    //mvPoints
    if (mPara.reoptimize_fetch_raw_points)
    {
        std::vector<PointType> non_ground_points;
        for (int np = 0; np < mvPoints.size(); np++)
        {
            if (mvPoints[np].bValid && mvPoints[np].clusterType != 1)
            {
                Eigen::Vector3f r_I_lu(mvPoints[np].point.x, mvPoints[np].point.y, mvPoints[np].point.z);
                Eigen::Vector3f r_I_lup;
                r_I_lup = R_lu_lup * r_I_lu;
                if (r_I_lup(2) < 8.0)
                    non_ground_points.push_back(mvPoints[np].point);
            }
        }
        tem.non_ground_raw_points_ = non_ground_points;
    }

    mpMap->insertScan(tem);
}

/**
 * @brief Checks if a given file exists using Boost filesystem.
 * 
 * @param filename The path of the file to check.
 * @return true if the file exists, false otherwise.
 */

bool fileExists(const std::string& filename) {
	//std::ifstream file(filename);
	return boost::filesystem::exists(filename);
}

/**
 * @brief Loads ground points from a saved feature file (e.g., text file) and transforms them into global map frame.
 *
 * @param single_scan_ground_features_file Path to the ground features file (e.g., txt).
 * @param R_lu_lup Rotation matrix from LiDAR to leveled unit plane.
 * @param ground_height Output: average height of ground points after leveling.
 * @return std::vector<PointType> containing all transformed ground points.
 */

std::vector<PointType> LidarScan::GetGroundFeaturesLoadOnce(std::string single_scan_ground_features_file, 
	Eigen::Matrix3f R_lu_lup, double& ground_height)
{
	std::vector<PointType> ground;

	// 打开文件并读入全部内容
	std::ifstream file(single_scan_ground_features_file, std::ios::in | std::ios::binary | std::ios::ate);
	if (!file.is_open())
	{
		std::cout << "Check file path please." << std::endl;
		throw std::runtime_error("Wrong path file");
		return ground;
	}

	std::streamsize size = file.tellg();
	file.seekg(0, std::ios::beg);

	std::vector<char> buffer(size);
	if (!file.read(buffer.data(), size))
	{
		throw std::runtime_error("Error reading file");
	}

	// 将整个文件的内容转为字符串
	std::string fileContent(buffer.begin(), buffer.end());

	int attribute_number = 4;
	double sum_height = 0;

	// 使用stringstream分割每一行
	std::istringstream fstr(fileContent);
	std::string line;
	getline(fstr, line); // 跳过文件的第一行
	int start_gnss_idx = -1;
	while (getline(fstr, line))
	{
		std::istringstream string_stream(line);
		std::string attribute;
		std::vector<std::string> attributes;
		while (getline(string_stream, attribute, ' '))
		{
			//if(point_coordinate.compare(" ")) continue;
			attributes.push_back(attribute);
		}
		std::vector<std::string> filtered_attributes;
		for (auto attri : attributes)
		{
			if (!attri.empty())
			{
				filtered_attributes.push_back(attri);
			}
		}
		if (filtered_attributes.size() != attribute_number)
		{
			std::cout << "Check file please." << std::endl;
			throw std::runtime_error("Wrong path file");
			return ground;
		}
		PointType point;
		point.x = std::stod(filtered_attributes[0]);
		point.y = std::stod(filtered_attributes[1]);
		point.z = std::stod(filtered_attributes[2]);
		point.intensity = std::stod(filtered_attributes[3]);
		double echo_range = std::sqrt(point.x*point.x + point.y*point.y + point.z*point.z);
		if (echo_range<mRangeTreshold.first || echo_range>mRangeTreshold.second)
		{
			continue;
		}
		Eigen::Vector3d result_p = R_lu_lup.cast<double>()*Eigen::Vector3d(point.x, point.y, point.z);
		sum_height += result_p(2);

		///////////////////////////////////////////////////////
		Eigen::Vector3d r_lu_m_t_current;
		Eigen::Matrix3d R_lu_m_t_current;
		int gnss_idx = 0;
		double current_t = mTimeInit + point.intensity * mScanDuration;
		//if (!(mpTraj->eopInterpolation(current_t / 1000.0, start_gnss_idx, gnss_idx, r_lu_m_t_current, R_lu_m_t_current)))
		//{
		//	throw std::runtime_error("Wrong Time");
		//}

        if (!(mpTraj->eopInterpolation(current_t / 1000.0, start_gnss_idx, gnss_idx, r_lu_m_t_current, R_lu_m_t_current)))
        {
            throw std::runtime_error("Wrong Time");
        }
		start_gnss_idx = gnss_idx;

		Eigen::Matrix3d R_lut2_lut1_gnss = R_lu_m_t_ini_.inverse() * R_lu_m_t_current;
		Eigen::Vector3d r_lut2_lut1_gnss = R_lu_m_t_ini_.inverse() * (r_lu_m_t_current - r_lu_m_t_ini_);

		// Eigen::Matrix3d R_lu_m_t_current_slam = R_lu_m_t1 * R_lut2_lut1_gnss;
		// Eigen::Vector3d r_lu_m_t_current_slam = r_lu_m_t1 + R_lu_m_t1 * r_lut2_lut1_gnss;

		Eigen::Vector3d r_I_m = R_lut2_lut1_gnss * point.getVector3fMap().cast<double>() + r_lut2_lut1_gnss;
		point.getVector3fMap() = r_I_m.cast<float>();
		///////////////////////////////////////////////////////
		ground.push_back(point);
	}
	ground_height = sum_height / ground.size();
	file.close();
	return ground;
}

/**
 * @brief Alternative ground feature loader (line-by-line stream version).
 *
 * This version uses standard line-wise reading instead of full-buffer reading.
 * Each point is transformed into the global frame using trajectory-based pose interpolation.
 *
 * @param single_scan_ground_features_file Path to the feature file.
 * @param R_lu_lup Rotation matrix for leveling.
 * @param ground_height Output: average height of all valid ground points.
 * @return std::vector<PointType> Ground points in global frame.
 */
std::vector<PointType> LidarScan::GetGroundFeatures(std::string single_scan_ground_features_file, Eigen::Matrix3f R_lu_lup, double& ground_height)
{
	std::vector<PointType> ground;

	//(TODO:ground points should be there)
	//if (!fileExists(single_scan_ground_features_file))
	//{
	//	std::cout << single_scan_ground_features_file << " does not exist." << std::endl;
	//	return ground;
	//}
	int attribute_number = 4;
	std::ifstream file_stream(single_scan_ground_features_file);
	if (!file_stream)
	{
		std::cout << "Check file path please." << std::endl;
		throw std::runtime_error("Wrong path file");
		return ground;
	}
	std::string line;
	int scan_id;
	double sum_height = 0;
	getline(file_stream, line);
	int start_gnss_idx = -1;
	while (getline(file_stream, line))
	{
		std::istringstream string_stream(line);
		std::string attribute;
		std::vector<std::string> attributes;
		while (getline(string_stream, attribute, ' '))
		{
			//if(point_coordinate.compare(" ")) continue;
			attributes.push_back(attribute);
		}
		std::vector<std::string> filtered_attributes;
		for (auto attri : attributes)
		{
			if (!attri.empty())
			{
				filtered_attributes.push_back(attri);
			}
		}
		if (filtered_attributes.size() != attribute_number)
		{
			std::cout << "Check file please." << std::endl;
			throw std::runtime_error("Wrong path file");
			return ground;
		}
		PointType point;
		point.x = std::stod(filtered_attributes[0]);
		point.y = std::stod(filtered_attributes[1]);
		point.z = std::stod(filtered_attributes[2]);
		point.intensity = std::stod(filtered_attributes[3]);
		double echo_range = std::sqrt(point.x*point.x + point.y*point.y + point.z*point.z);
		if (echo_range<mRangeTreshold.first || echo_range>mRangeTreshold.second)
		{
			continue;
		}
		Eigen::Vector3d result_p = R_lu_lup.cast<double>()*Eigen::Vector3d(point.x, point.y, point.z);
		sum_height += result_p(2);

        ///////////////////////////////////////////////////////
        Eigen::Vector3d r_lu_m_t_current;
        Eigen::Matrix3d R_lu_m_t_current;
        int gnss_idx = 0;
        double current_t =mTimeInit + point.intensity * mScanDuration;
        if (!(mpTraj->eopInterpolation(current_t / 1000.0, start_gnss_idx, gnss_idx, r_lu_m_t_current, R_lu_m_t_current)))
        {
            throw std::runtime_error("Wrong Time");
        }
         start_gnss_idx=gnss_idx;

        Eigen::Matrix3d R_lut2_lut1_gnss = R_lu_m_t_ini_.inverse() * R_lu_m_t_current;
        Eigen::Vector3d r_lut2_lut1_gnss = R_lu_m_t_ini_.inverse() * (r_lu_m_t_current - r_lu_m_t_ini_);

        // Eigen::Matrix3d R_lu_m_t_current_slam = R_lu_m_t1 * R_lut2_lut1_gnss;
        // Eigen::Vector3d r_lu_m_t_current_slam = r_lu_m_t1 + R_lu_m_t1 * r_lut2_lut1_gnss;

		Eigen::Vector3d r_I_m = R_lut2_lut1_gnss * point.getVector3fMap().cast<double>() + r_lut2_lut1_gnss;
        point.getVector3fMap() = r_I_m.cast<float>();
        ///////////////////////////////////////////////////////
		ground.push_back(point);
	}
	ground_height = sum_height / ground.size();
	file_stream.close();
	return ground;
}

/**
 * @brief Loads tree and ground feature files from a single scan and inserts them into the map.
 *
 * This function:
 * 1. Reads tree cluster data and associates them with unique IDs.
 * 2. Transforms all points using GNSS + interpolation.
 * 3. Computes averaged locations and updates integrated tree memory.
 * 4. Loads ground planar points and updates ScanInfo metadata.
 * 
 * @param single_scan_features_file Path to tree features (with tree IDs).
 * @param single_scan_ground_file Path to ground feature points.
 * @param single_scan_id Index of the scan.
 */
void LidarScan::AddScanToMap(std::string single_scan_features_file, std::string single_scan_ground_file,int single_scan_id)
{
#if 1
	int iscan_size = mpMap->mNumIntegratedScan;
    mpMap->scan_count++;
	//if (mpMap->scan_init)
	//{
	//	iscan_size = mpMap->mNumInitialScan;
	//}
	int read_scan_number = (single_scan_id - mpMap->mPara.initScan);
	if (/*mpMap->scan_count*/read_scan_number%iscan_size == 0 && read_scan_number != 0)
	{
		std::cout << "number of iscan points: " << mpMap->debug_for_iscan_points_number << std::endl;
		std::cout << "mpMap->integrated_tree_locations_ size: " << mpMap->integrated_tree_locations_.size() << std::endl;
		std::cout << "mpMap->integrated_tree_Ids_ size: " << mpMap->integrated_tree_Ids_.size() << std::endl;
#ifdef EXPORT_LOG
		fLog << "mpMap->integrated_tree_locations_ size: " << mpMap->integrated_tree_locations_.size() << std::endl;
		fLog << "mpMap->integrated_tree_Ids_ size: " << mpMap->integrated_tree_Ids_.size() << std::endl;
#endif
		mpMap->mvIntegratedTreeVector.push_back(mpMap->integrated_tree_locations_);
		mpMap->mvIntegratedTreeIdsVector.push_back(mpMap->integrated_tree_Ids_);
		mpMap->integrated_tree_locations_.clear();
		mpMap->integrated_tree_Ids_.clear();
		mpMap->tree_id_in_iscan_to_location_index_.clear();
		mpMap->tree_id_in_iscan_to_vector_index_.clear();
		mpMap->local_scan_id = 0;
		mpMap->scan_init = false;
		std::cout<<"iscan finished"<<std::endl;
		mpMap->debug_for_iscan_points_number = 0;
	}
	bool tree_feature_file_exist = fileExists(single_scan_features_file);
	bool ground_feature_file_exist = fileExists(single_scan_ground_file);

	if (!(tree_feature_file_exist&&ground_feature_file_exist))
	{
		ScanInfo tem_empty;
		tem_empty.scanID = single_scan_id; tem_empty.partID = mPartID; tem_empty.mTimeInit = mTimeInit; tem_empty.mTimeEnd = mTimeEnd; tem_empty.mTimeTracking = mTimeTracking;
		tem_empty.R_lu_lup = R_lu_b_ref.cast<double>();
		tem_empty.valid_feature_flag = valid_feature_flag_;

		tem_empty.R_lut2_lut1 = R_lut2_lut1_ref;
		tem_empty.r_lut2_lut1 = r_lut2_lut1_ref;

		tem_empty.gH = ground_height_ref;
		std::vector<std::vector<PointType>> vTreeRawPoints_empty;
		vector<Eigen::Vector3d> vTreeLoc_empty;
		tem_empty.vTreeRawPoints = vTreeRawPoints_empty; 
		tem_empty.vTreeLoc = vTreeLoc_empty;
		std::vector<PointType> ground_empty;
		tem_empty.vGroundPlanarRawPoints = ground_empty;
		mpMap->insertScan(tem_empty);
		mpMap->local_scan_id++;
		return;
	}

	if (!fileExists(single_scan_features_file))
	{
		std::cout << single_scan_features_file << " does not exist." << std::endl;
		return;
	}
	if (!fileExists(single_scan_ground_file))
	{
		std::cout << single_scan_ground_file << " does not exist." << std::endl;
		return;
	}
	
	vector<Eigen::Vector3d> vIntegratedTree_= mpMap->integrated_tree_locations_;           // locations of the tree (vector3d)
	vector<vector<pair<int, int>>> vIntegratedTreeIds_ = mpMap->integrated_tree_Ids_; // scan id + tree id of each tree
	std::map<int, int> tree_id_in_iscan_to_location_index = mpMap->tree_id_in_iscan_to_location_index_;

	std::map<int, int> tree_id_in_iscan_to_vector_index = mpMap->tree_id_in_iscan_to_vector_index_;

	int attribute_number = 5;
	std::ifstream file_stream(single_scan_features_file);
	if (!file_stream)
	{
		std::cout << "Check file path please." << std::endl;
		return;
	}
	std::string line;
	std::vector<std::vector<PointType>> vTreeRawPoints;
	vector<Eigen::Vector3d> vTreeLoc;
	int scan_id= mpMap->local_scan_id;
	std::map<int, std::pair<int, int>> iscan_tree_id_to_raw_scan_id_and_raw_scan_tree_id;
	std::map<int, int> raw_scan_tree_id_to_vector_index;
	int raw_scan_tree_id = 0;
	std::map<int,int> tree_id_in_iscan_to_id_in_raw_scan;
	getline(file_stream, line);
	int start_gnss_idx = -1;
	while (getline(file_stream, line))
	{
		std::istringstream string_stream(line);
		std::string attribute;
		std::vector<std::string> attributes;
		while (getline(string_stream, attribute, ' '))
		{
			//if(point_coordinate.compare(" ")) continue;
			attributes.push_back(attribute);
		}
		std::vector<std::string> filtered_attributes;
		for (auto attri : attributes)
		{
			if (!attri.empty())
			{
				filtered_attributes.push_back(attri);
			}
		}
		if (filtered_attributes.size() != attribute_number)
		{
			std::cout << "Check file please." << std::endl;
			return;
		}
		PointType point;
		point.x = std::stod(filtered_attributes[0]);
		point.y = std::stod(filtered_attributes[1]);
		point.z = std::stod(filtered_attributes[2]);
		point.intensity = std::stod(filtered_attributes[3]);
		// (ADD 2.5.2025)
		double echo_range = std::sqrt(point.x*point.x + point.y*point.y + point.z*point.z);
		if (echo_range<mRangeTreshold.first || echo_range>mRangeTreshold.second)
		{
			continue;
		}
        //(NEW 6.18)
        ///////////////////////////////////////////////////////
        Eigen::Vector3d r_lu_m_t_current;
        Eigen::Matrix3d R_lu_m_t_current;
        int gnss_idx = 0;
        double current_t =mTimeInit + point.intensity * mScanDuration;
        if (!(mpTraj->eopInterpolation(current_t / 1000.0, start_gnss_idx, gnss_idx, r_lu_m_t_current, R_lu_m_t_current)))
        {
            throw std::runtime_error("Wrong Time");
        }
         start_gnss_idx=gnss_idx;

        Eigen::Matrix3d R_lut2_lut1_gnss = R_lu_m_t_ini_.inverse() * R_lu_m_t_current;
        Eigen::Vector3d r_lut2_lut1_gnss = R_lu_m_t_ini_.inverse() * (r_lu_m_t_current - r_lu_m_t_ini_);

        // Eigen::Matrix3d R_lu_m_t_current_slam = R_lu_m_t1 * R_lut2_lut1_gnss;
        // Eigen::Vector3d r_lu_m_t_current_slam = r_lu_m_t1 + R_lu_m_t1 * r_lut2_lut1_gnss;

		Eigen::Vector3d r_I_m = R_lut2_lut1_gnss * point.getVector3fMap().cast<double>() + r_lut2_lut1_gnss;
        point.getVector3fMap() = r_I_m.cast<float>();
        ///////////////////////////////////////////////////////

		mpMap->debug_for_iscan_points_number++;

		bool new_tree_in_raw_scan = false;
		int iscan_tree_id = std::stoi(filtered_attributes[4]); //(CHANGE)
		if (tree_id_in_iscan_to_id_in_raw_scan.find(iscan_tree_id) == tree_id_in_iscan_to_id_in_raw_scan.end())
		{
			tree_id_in_iscan_to_id_in_raw_scan[iscan_tree_id] = raw_scan_tree_id;
			std::vector<PointType> raw_tree_points;
			raw_tree_points.push_back(point);
			vTreeRawPoints.push_back(raw_tree_points);
			raw_scan_tree_id++;
			new_tree_in_raw_scan = true;
		}
		else
		{
			vTreeRawPoints[tree_id_in_iscan_to_id_in_raw_scan[iscan_tree_id]].push_back(point);
		}
		std::pair<int, int> scan_id_and_tree_id_in_raw_scan = std::make_pair(scan_id, tree_id_in_iscan_to_id_in_raw_scan[iscan_tree_id]);
		if (tree_id_in_iscan_to_vector_index.find(iscan_tree_id) == tree_id_in_iscan_to_vector_index.end())
		{
			/*vTreeRawPoints[raw_scan_tree_id_to_vector_index[raw_scan_tree_id]].push_back(point);*/
			std::vector<std::pair<int, int>> scan_id_and_tree_id_in_raw_scan_vector;
			scan_id_and_tree_id_in_raw_scan_vector.push_back(scan_id_and_tree_id_in_raw_scan);
			tree_id_in_iscan_to_vector_index[iscan_tree_id] = vIntegratedTreeIds_.size();
			vIntegratedTreeIds_.push_back(scan_id_and_tree_id_in_raw_scan_vector);
		}
		else
		{
			if (new_tree_in_raw_scan)
			{
				vIntegratedTreeIds_[tree_id_in_iscan_to_vector_index[iscan_tree_id]].push_back(scan_id_and_tree_id_in_raw_scan);
			}
		}

		if (tree_id_in_iscan_to_location_index.find(iscan_tree_id) == tree_id_in_iscan_to_location_index.end())
		{

			tree_id_in_iscan_to_location_index[iscan_tree_id] = vIntegratedTree_.size();
			Eigen::Vector3d tree_location;
			tree_location.setZero();
			vIntegratedTree_.push_back(tree_location);
		}
	}
	file_stream.close();

	std::map<int, int> tree_id_raw_scan_to_id_in_iscan;
	for (const auto& ele : tree_id_in_iscan_to_id_in_raw_scan)
	{
		tree_id_raw_scan_to_id_in_iscan[ele.second] = ele.first;
	}

	R_lu_lup = R_lu_b_ref.cast<float>();
	mGroundH = ground_height_ref;

	ScanInfo tem;
	tem.scanID = single_scan_id; tem.partID = mPartID; tem.mTimeInit = mTimeInit; tem.mTimeEnd = mTimeEnd; tem.mTimeTracking = mTimeTracking;
	tem.R_lu_lup = R_lu_lup.cast<double>(); 
	tem.valid_feature_flag = valid_feature_flag_;
	int treeid_in_raw_scan=0;
	for (const auto& single_tree_points : vTreeRawPoints)
	{
		Eigen::Vector3f sum_result;
		sum_result.setZero();
		Eigen::Vector3f offset;
		int count = 0;
		for(const auto& tree_point: single_tree_points)
		{
			if (count == 0)
			{
				offset = tree_point.getVector3fMap();
			}
			else
			{
				sum_result += (tree_point.getVector3fMap() - offset);
			}
			count++;
		}
		Eigen::Vector3f tree_location = (sum_result / count) + offset;
		vTreeLoc.push_back(tree_location.cast<double>());
		auto & loc = vIntegratedTree_[tree_id_in_iscan_to_location_index[tree_id_raw_scan_to_id_in_iscan[treeid_in_raw_scan]]];
		Eigen::Vector3d leveled_tree_location = R_lu_lup.cast<double>()*tree_location.cast<double>();
		if (fabs(loc(0)) < 1e-4&& fabs(loc(1)) < 1e-4&& fabs(loc(2)) < 1e-4)
		{
			loc = leveled_tree_location.cast<double>();
		}
		else
		{
			loc = (loc + leveled_tree_location.cast<double>()) / 2;
		}
		treeid_in_raw_scan++;
	}
	tem.R_lut2_lut1 = R_lut2_lut1_ref;
	tem.r_lut2_lut1 = r_lut2_lut1_ref;
	tem.vTreeRawPoints = vTreeRawPoints; tem.vTreeLoc = vTreeLoc; /*tem.vGroundRawPoints = vGroundRawPoints;*/
	std::string ground_file = single_scan_ground_file;
	double ground_height;
	const auto ground = GetGroundFeaturesLoadOnce(ground_file, R_lu_lup, ground_height);
	tem.vGroundPlanarRawPoints = ground; tem.gH = ground_height;
	//(TODO:) another tem.non_ground_raw_points_ = non_ground_points;

	mpMap->insertScan(tem);
	mpMap->integrated_tree_locations_.swap(vIntegratedTree_);
	mpMap->integrated_tree_Ids_.swap(vIntegratedTreeIds_);
	mpMap->tree_id_in_iscan_to_location_index_ = tree_id_in_iscan_to_location_index;
	mpMap->tree_id_in_iscan_to_vector_index_ = tree_id_in_iscan_to_vector_index;
	mpMap->local_scan_id++;

	if (mpMap->mPara.endScan == single_scan_id/*(mpMap->mPara.initScan + mpMap->scan_count-1)*/)
	{
		mpMap->mvIntegratedTreeVector.push_back(mpMap->integrated_tree_locations_);
		mpMap->mvIntegratedTreeIdsVector.push_back(mpMap->integrated_tree_Ids_);
		mpMap->integrated_tree_locations_.clear();
		mpMap->integrated_tree_Ids_.clear();
		mpMap->tree_id_in_iscan_to_location_index_.clear();
		mpMap->tree_id_in_iscan_to_vector_index_.clear();
		mpMap->local_scan_id = 0;
		mpMap->scan_init = false;
		std::cout << "iscan finished" << std::endl;
	}

	// mpMap->scan_count++;
#endif
}

/************************************
Compute Odometry based on the extracted planar points and tree clusters from leveled point cloud
Requirement: mbTrackFeature & !mbInit
Notes: distortion not used here
Output:
    Odometry result, r_lut2_lut1, R_lut2_lut1
mbTrackFeature: fail when valid number of tree pairs < MinTreePair
************************************/
/**
 * @brief Computes the 6-DOF odometry between two scans using leveled planar ground points and tree cluster features.
 * 
 * @return true if a valid transformation is estimated; false otherwise.
 * 
 * Requirements: Ground must be reasonably flat and enough tree matches must exist.
 * Uses a joint optimization over ground-plane residuals and cylindrical tree features.
 */

bool LidarScan::computeOdometryLeveled()
{
    fLog <<"Compute transformation through ground planar points and tree clusters"<<endl;

    int numPlanar = mpSurfPointsFlatGround->points.size();
    int numPair = mVecTreePairs.size();
    cout << "Planar points: " << mpPrev->mpSurfPointsFlatGround->size() << " with " << mpSurfPointsFlatGround->size() << endl;

    // build tree
    pcl::KdTreeFLANN<pcl::PointXYZI>::Ptr kdtreeSurfLast(new pcl::KdTreeFLANN<pcl::PointXYZI>());
    kdtreeSurfLast->setInputCloud(mpPrev->mpSurfPointsFlatGround);

    // initialize parameters
    double para_q[4] = {q_ini_level.x(), q_ini_level.y(), q_ini_level.z(), q_ini_level.w()};
    double para_t[3] = {t_ini_level[0], t_ini_level[1], t_ini_level[2]};

    Eigen::Map<Eigen::Quaterniond> q_last_curr(para_q); // requires array in [x, y, z, w]
    Eigen::Map<Eigen::Vector3d> t_last_curr(para_t);

    cout << "Initial " << t_last_curr.transpose() << "\t" << rad2deg(Find_Rotation(q_last_curr.toRotationMatrix())).transpose() << endl;

#ifdef EXPORT_LOG
    fLog << "Initial levled after tree matching " << t_last_curr.transpose() << "\t" << rad2deg(Find_Rotation(q_last_curr.toRotationMatrix())).transpose() << endl;
#endif
    // current estimation of q and t
    qTem = q_last_curr;
    tTem = t_last_curr;

    int iter = 0;
    int maxIter = 2;
    vector<bool> validTreePair(numPair, true);
    while (iter < maxIter)
    {
        iter++;
        // build ceres
        ceres::LossFunction *loss_function = new ceres::HuberLoss(0.1);
        ceres::LocalParameterization *q_parameterization = new ceres::EigenQuaternionParameterization(); // x, y, z, w
        ceres::Problem::Options problem_options;
        ceres::Problem problem(problem_options);

        problem.AddParameterBlock(para_q, 4, q_parameterization);
        problem.AddParameterBlock(para_t, 3);

        // find correspondence
        int numPlanarCorres = 0;
        PointType pointSel;
        std::vector<int> pointSearchInd;
        std::vector<float> pointSearchSqDis;

        vector<Eigen::Vector3d> planarPoints;
        double res_planar = 0.0;
        for (int i = 0; i < numPlanar; ++i)
        {
            // TransformToStart(&(surfPointsFlat->points[i]), &pointSel);
            TransformToStart(&(mpSurfPointsFlatGround->points[i]), &pointSel);
            kdtreeSurfLast->nearestKSearch(pointSel, 1, pointSearchInd, pointSearchSqDis);

            int closestPointInd = -1, minPointInd2 = -1, minPointInd3 = -1;
            if (pointSearchSqDis[0] < DISTANCE_SQ_THRESHOLD)
            {
                closestPointInd = pointSearchInd[0];
                // get closest point's scan ID
                int closestPointScanID = int(mpPrev->mpSurfPointsFlatGround->points[closestPointInd].intensity);
                double minPointSqDis2 = DISTANCE_SQ_THRESHOLD, minPointSqDis3 = DISTANCE_SQ_THRESHOLD;

                // search in the direction of increasing scan line
                for (int j = closestPointInd + 1; j < (int)mpPrev->mpSurfPointsFlatGround->points.size(); ++j)
                {
                    // if not in nearby scans, end the loop
                    if (int(mpPrev->mpSurfPointsFlatGround->points[j].intensity) > (closestPointScanID + NEARBY_SCAN))
                        break;

                    double pointSqDis = (mpPrev->mpSurfPointsFlatGround->points[j].x - pointSel.x) * (mpPrev->mpSurfPointsFlatGround->points[j].x - pointSel.x) +
                                        (mpPrev->mpSurfPointsFlatGround->points[j].y - pointSel.y) * (mpPrev->mpSurfPointsFlatGround->points[j].y - pointSel.y) +
                                        (mpPrev->mpSurfPointsFlatGround->points[j].z - pointSel.z) * (mpPrev->mpSurfPointsFlatGround->points[j].z - pointSel.z);

                    // if in the same or lower scan line
                    if (int(mpPrev->mpSurfPointsFlatGround->points[j].intensity) <= closestPointScanID && pointSqDis < minPointSqDis2)
                    {
                        minPointSqDis2 = pointSqDis;
                        minPointInd2 = j;
                    }
                    // if in the higher scan line
                    else if (int(mpPrev->mpSurfPointsFlatGround->points[j].intensity) > closestPointScanID && pointSqDis < minPointSqDis3)
                    {
                        minPointSqDis3 = pointSqDis;
                        minPointInd3 = j;
                    }
                }
                // search in the direction of decreasing scan line
                for (int j = closestPointInd - 1; j >= 0; --j)
                {
                    // if not in nearby scans, end the loop
                    if (int(mpPrev->mpSurfPointsFlatGround->points[j].intensity) < (closestPointScanID - NEARBY_SCAN))
                        break;

                    double pointSqDis = (mpPrev->mpSurfPointsFlatGround->points[j].x - pointSel.x) * (mpPrev->mpSurfPointsFlatGround->points[j].x - pointSel.x) +
                                        (mpPrev->mpSurfPointsFlatGround->points[j].y - pointSel.y) * (mpPrev->mpSurfPointsFlatGround->points[j].y - pointSel.y) +
                                        (mpPrev->mpSurfPointsFlatGround->points[j].z - pointSel.z) * (mpPrev->mpSurfPointsFlatGround->points[j].z - pointSel.z);

                    // if in the same or higher scan line
                    if (int(mpPrev->mpSurfPointsFlatGround->points[j].intensity) >= closestPointScanID && pointSqDis < minPointSqDis2)
                    {
                        minPointSqDis2 = pointSqDis;
                        minPointInd2 = j;
                    }
                    else if (int(mpPrev->mpSurfPointsFlatGround->points[j].intensity) < closestPointScanID && pointSqDis < minPointSqDis3)
                    {
                        // find nearer point
                        minPointSqDis3 = pointSqDis;
                        minPointInd3 = j;
                    }
                }

                if (minPointInd2 >= 0 && minPointInd3 >= 0)
                {

                    Eigen::Vector3d curr_point(mpSurfPointsFlatGround->points[i].x,
                                               mpSurfPointsFlatGround->points[i].y,
                                               mpSurfPointsFlatGround->points[i].z);
                    Eigen::Vector3d last_point_a(mpPrev->mpSurfPointsFlatGround->points[closestPointInd].x,
                                                 mpPrev->mpSurfPointsFlatGround->points[closestPointInd].y,
                                                 mpPrev->mpSurfPointsFlatGround->points[closestPointInd].z);
                    Eigen::Vector3d last_point_b(mpPrev->mpSurfPointsFlatGround->points[minPointInd2].x,
                                                 mpPrev->mpSurfPointsFlatGround->points[minPointInd2].y,
                                                 mpPrev->mpSurfPointsFlatGround->points[minPointInd2].z);
                    Eigen::Vector3d last_point_c(mpPrev->mpSurfPointsFlatGround->points[minPointInd3].x,
                                                 mpPrev->mpSurfPointsFlatGround->points[minPointInd3].y,
                                                 mpPrev->mpSurfPointsFlatGround->points[minPointInd3].z);

                    double s;
                    // if (mLidarDistortion)
                    //     s = (mpSurfPointsFlatGround->points[i].intensity - int(mpSurfPointsFlatGround->points[i].intensity));
                    // else
                    s = 1.0;

                    if (s > 1.0)
                    {
                        cout << s << endl;
                        throw std::runtime_error("Wrong scale of time");
                    }

                    // ceres::CostFunction *cost_function = LidarPlaneFactor::Create(curr_point, last_point_a, last_point_b, last_point_c, s);
                    ceres::CostFunction *cost_function = LidarPlaneFactorStd::Create(curr_point, last_point_a, last_point_b, last_point_c, s, 0.8);
                    problem.AddResidualBlock(cost_function, loss_function, para_q, para_t);
                    numPlanarCorres++;
                    planarPoints.push_back(curr_point);
                    planarPoints.push_back(last_point_a);
                    planarPoints.push_back(last_point_b);
                    planarPoints.push_back(last_point_c);
                    double d = computePoint2PlaneDistance(q_last_curr * curr_point + t_last_curr, last_point_a, last_point_b, last_point_c);
                    res_planar = res_planar + d * d;
                }
            }
        }

        float res_dis = 0.0, res_angle = 0.0;
        int numValidTree = 0;
        for (int i = 0; i < numPair; ++i)
        {
            if (!validTreePair[i])
                continue;
            Eigen::Vector3d norm_prev = mpPrev->mvClusterTree[mVecTreePairs[i].first].vNorm.cast<double>();
            Eigen::Vector3d center_prev = mpPrev->mvClusterTree[mVecTreePairs[i].first].vCenterTrans.cast<double>();

            Eigen::Vector3d norm_cur = mvClusterTree[mVecTreePairs[i].second].vNorm.cast<double>();
            Eigen::Vector3d center_cur = mvClusterTree[mVecTreePairs[i].second].vCenterTrans.cast<double>();

            float d = computePoint2lineDistance(q_last_curr * center_cur + t_last_curr, center_prev, center_prev + norm_prev);
            float angle = computeAngle(norm_prev, q_last_curr * norm_cur);
            // cout << d << " " << rad2deg(angle) << endl;
            res_dis = res_dis + d * d;
            res_angle = res_angle + angle * angle;
            double s;
            // if (mLidarDistortion)
            //     s = (mpSurfPointsFlat->points[i].intensity - int(mpSurfPointsFlat->points[i].intensity));
            // else
            s = 1.0;

            ceres::CostFunction *cost_function = LidarEdgeFactorOneResidual::Create(center_cur, center_prev, center_prev + norm_prev, s, 0.8);
            problem.AddResidualBlock(cost_function, loss_function, para_q, para_t);

            ceres::CostFunction *cost_function2 = NormalFactorOneResidual::Create(norm_cur, norm_prev, deg2rad(3.0));
            problem.AddResidualBlock(cost_function2, loss_function, para_q);
            numValidTree++;
        }

        TicToc t_solver;
        ceres::Solver::Options options;
        options.linear_solver_type = ceres::DENSE_QR;
        options.max_num_iterations = 4;
        options.minimizer_progress_to_stdout = false;
        ceres::Solver::Summary summary;
        ceres::Solve(options, &problem, &summary);

        qTem = q_last_curr;
        tTem = t_last_curr;

        cout << "Refined " << tTem.transpose() << "\t" << rad2deg(Find_Rotation(qTem.toRotationMatrix())).transpose() << endl;

#ifdef EXPORT_LOG
        fLog << "P:"
             << " iter-" << to_string(iter)
             << "\t" << numPlanarCorres << "/" << numValidTree << "\t" << endl;
        fLog << "\t" << tTem.transpose() << "\t" << rad2deg(Find_Rotation(qTem.toRotationMatrix())).transpose() << endl;
#endif

        double res_planar_after = 0.0;
        for (int i = 0; i < numPlanarCorres; ++i)
        {
            double d = computePoint2PlaneDistance(q_last_curr * planarPoints[i * 4] + t_last_curr, planarPoints[i * 4 + 1], planarPoints[i * 4 + 2], planarPoints[i * 4 + 3]);
            res_planar_after = res_planar_after + d * d;
        }

        float res_dis_after = 0.0, res_angle_after = 0.0;
        int numValidTreePairs = 0;
        for (int i = 0; i < numPair; ++i)
        {
            if (!validTreePair[i])
                continue;
            Eigen::Vector3d norm_prev = mpPrev->mvClusterTree[mVecTreePairs[i].first].vNorm.cast<double>();
            Eigen::Vector3d center_prev = mpPrev->mvClusterTree[mVecTreePairs[i].first].vCenterTrans.cast<double>();

            Eigen::Vector3d norm_cur = mvClusterTree[mVecTreePairs[i].second].vNorm.cast<double>();
            Eigen::Vector3d center_cur = mvClusterTree[mVecTreePairs[i].second].vCenterTrans.cast<double>();

            float d = computePoint2lineDistance(q_last_curr * center_cur + t_last_curr, center_prev, center_prev + norm_prev);
            float angle = computeAngle(norm_prev, q_last_curr * norm_cur);

            if (d > 0.5 || angle > deg2rad(10.0))
                // if (d > 0.5)
                validTreePair[i] = false;
            else
                numValidTreePairs++;

            res_dis_after = res_dis_after + d * d;
            res_angle_after = res_angle_after + angle * angle;

            // cout << mVecTreePairs[i].first << " to " << mVecTreePairs[i].second << " : " << d << " " << rad2deg(angle) << endl;
            // double s;
        }

        cout << "Plane res: " << sqrt(res_planar / float(numPlanarCorres)) << " -> " << sqrt(res_planar_after / float(numPlanarCorres))
             << ", distance res: " << sqrt(res_dis / float(numPair)) << " -> " << sqrt(res_dis_after / float(numPair)) << ", angle res: "
             << rad2deg(sqrt(res_angle / float(numPair))) << " -> " << rad2deg(sqrt(res_angle_after / float(numPair))) << endl;
#ifdef EXPORT_LOG
        fLog << "\tPlane res: " << sqrt(res_planar / float(numPlanarCorres)) << " -> " << sqrt(res_planar_after / float(numPlanarCorres))
             << "\tDistance res: " << sqrt(res_dis / float(numPair)) << " -> " << sqrt(res_dis_after / float(numPair))
             << "\tAngle res: " << rad2deg(sqrt(res_angle / float(numPair))) << " -> " << rad2deg(sqrt(res_angle_after / float(numPair))) << endl;
#endif
        if (numValidTreePairs < MinTreePair) // quit the function
        {

#ifdef EXPORT_LOG
            f_odometry_debug << id << "\tFailed, number of tree matches in optimization at iter " << iter << ": " << numValidTreePairs << endl;
            fLog << "Error: remaining tree matches in too few: " << numValidTreePairs << endl;
#endif
            return false; 
        }
    }

    // save valid tree pairs
    std::vector<pair<int, int>>
        vecTreePairs;
    for (int i = 0; i < numPair; ++i)
    {
        if (!validTreePair[i])
            continue;

        vecTreePairs.push_back(mVecTreePairs[i]);
    }
    mVecTreePairs = vecTreePairs;

    // the rotation/translation after applying the leveling r_lu_lup
    Eigen::Vector3d r_lut2_lut1_level = t_last_curr;
    Eigen::Matrix3d R_lut2_lut1_level = q_last_curr.toRotationMatrix();

    // compute transformation in original frame
    Eigen::Matrix3d R_lu_lup_t1 = mpPrev->R_lu_lup.cast<double>();
    Eigen::Matrix3d R_lu_lup_t2 = R_lu_lup.cast<double>();

    r_lut2_lut1_gt = R_lu_lup_t1.inverse() * r_lut2_lut1_level;
    R_lut2_lut1_gt = R_lu_lup_t1.inverse() * R_lut2_lut1_level * R_lu_lup_t2;

    if (mpTraj) // reference trajectory
    {
        Eigen::Vector3d r_lut2_lut1_level_ref = R_lu_lup_t1 * r_lut2_lut1_ref;
        Eigen::Matrix3d R_lut2_lut1_level_ref = R_lu_lup_t1 * R_lut2_lut1_ref * R_lu_lup_t2.inverse();

        cout << "Reference level: " << r_lut2_lut1_level_ref.transpose() << "\t" << rad2deg(Find_Rotation(R_lut2_lut1_level_ref)).transpose() << endl;

#ifdef EXPORT_LOG
        fLog << "Reference leveled: " << endl;
        fLog << "\t" << r_lut2_lut1_level_ref.transpose() << "\t" << rad2deg(Find_Rotation(R_lut2_lut1_level_ref)).transpose() << endl;
#endif

#ifdef EXPORT_TRAJECTORY

        Eigen::Matrix3d rot_dif = R_lut2_lut1_level * R_lut2_lut1_level_ref.inverse();
        Eigen::Vector3d pos_dif = r_lut2_lut1_level - r_lut2_lut1_level_ref;
        Eigen::Vector3d angles_dif = Find_Rotation(rot_dif);
        fTrajectoryRes << id << "\t" << std::setw(8) << pos_dif(0) << "\t" << std::setw(8) << pos_dif(1) << "\t" << std::setw(8) << pos_dif(2) << "\t"
                       << std::setw(8) << rad2deg(angles_dif(0)) << "\t" << std::setw(8) << rad2deg(angles_dif(1)) << "\t" << std::setw(8) << rad2deg(angles_dif(2)) << "\t" << mPartID << endl;

#endif
    }

    // based on the magnitude of result
    double distance = r_lut2_lut1_gt.norm();
    double angular_mag = Find_Rotation(R_lut2_lut1_gt).norm();

    if (distance > 0.5 || angular_mag > deg2rad(30.0))
    {
#ifdef EXPORT_LOG
        f_odometry_debug << id << "\tErroneous value for odometry" << distance << " " << angular_mag << endl;
        fLog << "Error: erroneous value for odometry" << distance << " " << angular_mag << endl;
#endif
        return false;
    }

    r_lut2_lut1 = r_lut2_lut1_gt;
    R_lut2_lut1 = R_lut2_lut1_gt;
    return true;
}

/************************************
Compute Odometry based on the ground information and 2d
Requirement: mbTrackFeature & !mbInit
Notes: distortion not used here
Output:
    Odometry result, r_lut2_lut1, R_lut2_lut1
mbTrackFeature: fail when valid number of tree pairs < MinTreePair
************************************/
/**
 * @brief Estimates 2D transformation (rotation + translation) between scans based on tree cluster centroids.
 * 
 * @return true if similarity transform is successfully estimated; false if too few valid tree matches remain.
 * 
 * Assumes leveled terrain and relies on least-squares optimization over 2D projection of tree centers.
 */

bool LidarScan::compute2dSimilarityLeveled()
{

    fLog <<"Compute transformation through 2d similarity "<<endl;
    int numPair = mVecTreePairs.size();
    vector<bool> validTreePair(numPair, true);

    double scale_, theta_, tx_, ty_;
    int iter = 0;
    int maxIter = 4;
    int numValidPair = numPair;
    Eigen::Vector4d tempResult(1, 0, 0, 0);
    double DistanceThreshold = 0.5;
    while (iter < maxIter)
    {
        const int n = 2 * numValidPair;

        Eigen::MatrixXd AA(n, 4); // cos, sin, tx, ty. x' = x*cos + (-y)*sin + tx   y' = x*sin + y*cos+ty;
        Eigen::MatrixXd yy(n, 1);

        int count = 0;
        for (int i = 0; i < numPair; i++)
        {
            if (!validTreePair[i])
                continue;

            Eigen::Vector3d center_prev = mpPrev->mvClusterTree[mVecTreePairs[i].first].vCenterTrans.cast<double>();
            Eigen::Vector3d center_cur = mvClusterTree[mVecTreePairs[i].second].vCenterTrans.cast<double>();

            AA(count * 2, 0) = center_cur(0);
            AA(count * 2, 1) = -center_cur(1);
            AA(count * 2, 2) = 1;
            AA(count * 2, 3) = 0;
            AA(count * 2 + 1, 0) = center_cur(1);
            AA(count * 2 + 1, 1) = center_cur(0);
            AA(count * 2 + 1, 2) = 0;
            AA(count * 2 + 1, 3) = 1;

            yy(count * 2) = center_prev(0);
            yy(count * 2 + 1) = center_prev(1);

            count++;
        }
        Eigen::VectorXd res1, res2;
        res1 = yy - AA * tempResult;
        Eigen::Vector4d result_ = AA.colPivHouseholderQr().solve(yy);
        res2 = yy - AA * result_;

        tempResult = result_;
        double s1, s2;
        s1 = res1.transpose() * res1;
        s1 = sqrt(s1 / double(numValidPair));
        s2 = res2.transpose() * res2;
        s2 = sqrt(s2 / double(numValidPair));
        cout << s1 << " -> " << s2 << endl;
        iter++;

#ifdef EXPORT_LOG
        fLog << "\t2d iter " << iter << ",\tPair: " << numValidPair << "\tResidual: " << s1 << " -> " << s2 << endl;
#endif
        count = 0;
        numValidPair = 0;
        bool bStop = true; // flag for whether tree pair removed

        for (int i = 0; i < numPair; ++i)
        {
            if (!validTreePair[i])
                continue;

            double d_sqr = res2(count * 2) * res2(count * 2) + res2(count * 2 + 1) * res2(count * 2 + 1);
            count++;

            if (d_sqr > DistanceThreshold * DistanceThreshold)
            {
                validTreePair[i] = false;
                bStop = false;
            }
            else
                numValidPair++;
        }

        //CONTROL
        if (numValidPair < MinTreePair)
        {

#ifdef EXPORT_LOG
            f_odometry_debug << id << "\tFailed tree matches in 2d optimization at iter " << iter << ": " << numValidPair << endl;
            fLog << "Error: num of remaining tree matches is too few: " << numValidPair << endl;
#endif
            return false;
        }

        scale_ = sqrt(result_(0) * result_(0) + result_(1) * result_(1));
        theta_ = atan2(result_(1), result_(0));
        tx_ = result_(2);
        ty_ = result_(3);
        cout << "Result: " << scale_ << "\t" << rad2deg(theta_) << "\t" << tx_ << "\t" << ty_ << endl;

        DistanceThreshold = 0.9 * DistanceThreshold;

        if (bStop)
            break;
    }

    // save valid tree pairs
    std::vector<pair<int, int>> vecTreePairs;
    for (int i = 0; i < numPair; ++i)
    {
        if (!validTreePair[i])
            continue;

        vecTreePairs.push_back(mVecTreePairs[i]);
    }
    mVecTreePairs = vecTreePairs;

    // compute the odo based on previous
    double tz_ = double(mpPrev->mGroundH - mGroundH);
    Eigen::AngleAxisd rollAngle_(theta_, Eigen::Vector3d::UnitZ());
    Eigen::Matrix3d R_kap_ = rollAngle_.matrix();

    // compute transformation in original frame
    Eigen::Matrix3d R_lu_lup_t1 = mpPrev->R_lu_lup.cast<double>();
    Eigen::Matrix3d R_lu_lup_t2 = R_lu_lup.cast<double>();

    Eigen::Matrix3d R_lut2_lut1_level = R_kap_;
    Eigen::Vector3d r_lut2_lut1_level = R_lut2_lut1_level * Eigen::Vector3d(tx_, ty_, 0) + Eigen::Vector3d(0, 0, tz_);

#ifdef EXPORT_LOG
    fLog << "\tEstimated 2d similarity: " << scale_ << "\t" << rad2deg(theta_) << "\t" << tx_ << "\t" << ty_ << endl;
    fLog << "\tEstimated transformation in leveled: " << r_lut2_lut1_level.transpose() << "\t" << rad2deg(Find_Rotation(R_lut2_lut1_level)).transpose() << endl;
#endif

    R_lut2_lut1_gt = R_lu_lup_t1.inverse() * R_lut2_lut1_level * R_lu_lup_t2;
    r_lut2_lut1_gt = R_lu_lup_t1.inverse() * r_lut2_lut1_level;

    if (mpTraj) // reference trajectory
    {
        Eigen::Vector3d r_lut2_lut1_level_ref = R_lu_lup_t1 * r_lut2_lut1_ref;
        Eigen::Matrix3d R_lut2_lut1_level_ref = R_lu_lup_t1 * R_lut2_lut1_ref * R_lu_lup_t2.inverse();

        cout << "Reference level: " << r_lut2_lut1_level_ref.transpose() << "\t" << rad2deg(Find_Rotation(R_lut2_lut1_level_ref)).transpose() << endl;

#ifdef EXPORT_LOG
        fLog << "\tReference transformation in leveled: "<< r_lut2_lut1_level_ref.transpose() << "\t" << rad2deg(Find_Rotation(R_lut2_lut1_level_ref)).transpose() << endl;
#endif

#ifdef EXPORT_TRAJECTORY

        Eigen::Matrix3d rot_dif = R_lut2_lut1_level * R_lut2_lut1_level_ref.inverse();
        Eigen::Vector3d pos_dif = r_lut2_lut1_level - r_lut2_lut1_level_ref;
        Eigen::Vector3d angles_dif = Find_Rotation(rot_dif);
        fTrajectoryRes << id << "\t" << std::setw(8) << pos_dif(0) << "\t" << std::setw(8) << pos_dif(1) << "\t" << std::setw(8) << pos_dif(2) << "\t"
                       << std::setw(8) << rad2deg(angles_dif(0)) << "\t" << std::setw(8) << rad2deg(angles_dif(1)) << "\t" << std::setw(8) << rad2deg(angles_dif(2)) << "\t" << mPartID << endl;
#endif
    }

    double distance = r_lut2_lut1_gt.norm();
    double angular_mag = Find_Rotation(R_lut2_lut1_gt).norm();
    if (distance > 0.5 || angular_mag > deg2rad(30.0))
    {
#ifdef EXPORT_LOG
        f_odometry_debug << id << "\tUnreasonable value for feature-based odometry: distance - " << distance << " angles " << rad2deg(angular_mag) << endl;
        fLog << "Error: Unreasonable value for feature-based odometry : distance - " << distance << " angles " << rad2deg(angular_mag) << endl;
#endif
        return false;
    }

    r_lut2_lut1 = r_lut2_lut1_gt;
    R_lut2_lut1 = R_lut2_lut1_gt;

    return true;
}

/*Get ground info from previous valid scans*/
/**
 * @brief Restores ground height and leveling rotation from the previous scan.
 * 
 * Used when current scan cannot extract new ground features reliably.
 */

void LidarScan::fetch_previous_ground_info()
{
    fLog << "Get ground info from previous scan" <<endl;
    mGroundH = ground_height_previous;
    R_lu_lup =  R_lu_lup_previous;
}

/*UPdate previous ground info*/
/**
 * @brief Stores the ground height and leveling rotation for use in subsequent scans if needed.
 */

void LidarScan::update_previous_ground_info()
{
    previous_ground_info_flag = true;
    ground_height_previous = mGroundH;
    R_lu_lup_previous = R_lu_lup;
}

/*Compute odometry based on planar and edge points*/
/**
 * @brief Estimates pose transformation using edge (corner) and planar features via Ceres optimization.
 * 
 * @return true if pose estimation converges and passes threshold checks; false otherwise.
 * 
 * Applies both geometric constraints (lines, planes) and optionally distortion correction.
 */

bool LidarScan::compute_point_based_odometry()
{
    fLog << "Compute transformation " <<endl;
    // number of features
    int numCorner = mpCornerPointsSharp->points.size();
    int numPlanar = mpSurfPointsFlat->points.size();

    cout << mpPrev->mpCornerPointsSharp->size() << "/" << mpPrev->mpSurfPointsFlat->size() << " with " << mpCornerPointsSharp->size() << "/" << mpSurfPointsFlat->size() << endl;

    // build tree
    pcl::KdTreeFLANN<pcl::PointXYZI>::Ptr kdtreeCornerLast(new pcl::KdTreeFLANN<pcl::PointXYZI>());
    pcl::KdTreeFLANN<pcl::PointXYZI>::Ptr kdtreeSurfLast(new pcl::KdTreeFLANN<pcl::PointXYZI>());
    kdtreeCornerLast->setInputCloud(mpPrev->mpCornerPointsSharp);
    kdtreeSurfLast->setInputCloud(mpPrev->mpSurfPointsFlat);

    // initialize parameters
    double para_q[4] = {q_ini.x(), q_ini.y(), q_ini.z(), q_ini.w()};
    double para_t[3] = {t_ini[0], t_ini[1], t_ini[2]};

    Eigen::Map<Eigen::Quaterniond> q_last_curr(para_q); // requires array in [x, y, z, w]
    Eigen::Map<Eigen::Vector3d> t_last_curr(para_t);

    cout << "Initial " << t_last_curr.transpose() << "\t" << rad2deg(Find_Rotation(q_last_curr.toRotationMatrix())).transpose() << endl;

    // current estimation of q and t
    qTem = q_last_curr;
    tTem = t_last_curr;

    int iter = 0;
    int maxIter = 2;

    while (iter < maxIter)
    {
        iter++;

        // build ceres
        ceres::LossFunction *loss_function = new ceres::HuberLoss(0.1);
        ceres::LocalParameterization *q_parameterization = new ceres::EigenQuaternionParameterization(); // x, y, z, w
        ceres::Problem::Options problem_options;
        ceres::Problem problem(problem_options);

        problem.AddParameterBlock(para_q, 4, q_parameterization);
        problem.AddParameterBlock(para_t, 3);

        // find correspondence
        int numCornerCorres = 0;
        int numPlanarCorres = 0;
        PointType pointSel;
        std::vector<int> pointSearchInd;
        std::vector<float> pointSearchSqDis;

        // string outPass2 = output_folder + "matches_" + to_string(id) + ".txt";
        // ofstream fMatch;
        // // fMatch.open(outPass2, std::ifstream::out);
        // fMatch << fixed << std::setprecision(6);

        // find correspondence for corner features
        for (int i = 0; i < numCorner; ++i)
        {
            TransformToStart(&(mpCornerPointsSharp->points[i]), &pointSel);
            // pointSel = mpCornerPointsSharp->points[i];
            kdtreeCornerLast->nearestKSearch(pointSel, 1, pointSearchInd, pointSearchSqDis);

            int closestPointInd = -1, minPointInd2 = -1;
            if (pointSearchSqDis[0] < DISTANCE_SQ_THRESHOLD)
            {
                closestPointInd = pointSearchInd[0];
                int closestPointScanID = int(mpPrev->mpCornerPointsSharp->points[closestPointInd].intensity);

                double minPointSqDis2 = DISTANCE_SQ_THRESHOLD;
                // search in the direction of increasing scan line
                for (int j = closestPointInd + 1; j < (int)mpPrev->mpCornerPointsSharp->points.size(); ++j)
                {
                    // if in the same scan line, continue
                    if (int(mpPrev->mpCornerPointsSharp->points[j].intensity) <= closestPointScanID)
                        continue;

                    // if not in nearby scans, end the loop
                    if (int(mpPrev->mpCornerPointsSharp->points[j].intensity) > (closestPointScanID + NEARBY_SCAN))
                        break;

                    double pointSqDis = (mpPrev->mpCornerPointsSharp->points[j].x - pointSel.x) * (mpPrev->mpCornerPointsSharp->points[j].x - pointSel.x) +
                                        (mpPrev->mpCornerPointsSharp->points[j].y - pointSel.y) * (mpPrev->mpCornerPointsSharp->points[j].y - pointSel.y) +
                                        (mpPrev->mpCornerPointsSharp->points[j].z - pointSel.z) * (mpPrev->mpCornerPointsSharp->points[j].z - pointSel.z);

                    if (pointSqDis < minPointSqDis2)
                    {
                        // find nearer point
                        minPointSqDis2 = pointSqDis;
                        minPointInd2 = j;
                    }
                }

                // search in the direction of decreasing scan line
                for (int j = closestPointInd - 1; j >= 0; --j)
                {
                    // if in the same scan line, continue
                    if (int(mpPrev->mpCornerPointsSharp->points[j].intensity) >= closestPointScanID)
                        continue;

                    // if not in nearby scans, end the loop
                    if (int(mpPrev->mpCornerPointsSharp->points[j].intensity) < (closestPointScanID - NEARBY_SCAN))
                        break;

                    double pointSqDis = (mpPrev->mpCornerPointsSharp->points[j].x - pointSel.x) * (mpPrev->mpCornerPointsSharp->points[j].x - pointSel.x) +
                                        (mpPrev->mpCornerPointsSharp->points[j].y - pointSel.y) * (mpPrev->mpCornerPointsSharp->points[j].y - pointSel.y) +
                                        (mpPrev->mpCornerPointsSharp->points[j].z - pointSel.z) * (mpPrev->mpCornerPointsSharp->points[j].z - pointSel.z);

                    if (pointSqDis < minPointSqDis2)
                    {
                        // find nearer point
                        minPointSqDis2 = pointSqDis;
                        minPointInd2 = j;
                    }
                }
            }
            if (minPointInd2 >= 0) // both closestPointInd and minPointInd2 is valid
            {
                Eigen::Vector3d curr_point(mpCornerPointsSharp->points[i].x,
                                           mpCornerPointsSharp->points[i].y,
                                           mpCornerPointsSharp->points[i].z);
                Eigen::Vector3d last_point_a(mpPrev->mpCornerPointsSharp->points[closestPointInd].x,
                                             mpPrev->mpCornerPointsSharp->points[closestPointInd].y,
                                             mpPrev->mpCornerPointsSharp->points[closestPointInd].z);
                Eigen::Vector3d last_point_b(mpPrev->mpCornerPointsSharp->points[minPointInd2].x,
                                             mpPrev->mpCornerPointsSharp->points[minPointInd2].y,
                                             mpPrev->mpCornerPointsSharp->points[minPointInd2].z);

                double s;
                if (mLidarDistortion)
                    s = (mpCornerPointsSharp->points[i].intensity - int(mpCornerPointsSharp->points[i].intensity));
                else
                    s = 1.0;

                if (s > 1.0)
                {
                    cout << s << endl;
                    throw std::runtime_error("Wrong scale of time");
                }
                ceres::CostFunction *cost_function = LidarEdgeFactor::Create(curr_point, last_point_a, last_point_b, s);
                problem.AddResidualBlock(cost_function, loss_function, para_q, para_t);
                numCornerCorres++;

                // fMatch << numCornerCorres << "\t" << curr_point.transpose() << "\t " << mpCornerPointsSharp->points[i].intensity << endl;
                // fMatch << numCornerCorres << "\t" << last_point_a.transpose() << "\t " << mpPrev->mpCornerPointsSharp->points[closestPointInd].intensity << endl;
                // fMatch << numCornerCorres << "\t" << last_point_b.transpose() << "\t " << mpPrev->mpCornerPointsSharp->points[minPointInd2].intensity << endl;
            }
        }
        cout << "Number of corner correspondence " << numCornerCorres << endl;

        // find correspondence for plane features
        for (int i = 0; i < numPlanar; ++i)
        {
            // TransformToStart(&(surfPointsFlat->points[i]), &pointSel);
            TransformToStart(&(mpSurfPointsFlat->points[i]), &pointSel);
            // pointSel = mpSurfPointsFlat->points[i];
            kdtreeSurfLast->nearestKSearch(pointSel, 1, pointSearchInd, pointSearchSqDis);

            int closestPointInd = -1, minPointInd2 = -1, minPointInd3 = -1;
            if (pointSearchSqDis[0] < DISTANCE_SQ_THRESHOLD)
            {
                closestPointInd = pointSearchInd[0];

                // get closest point's scan ID
                int closestPointScanID = int(mpPrev->mpSurfPointsFlat->points[closestPointInd].intensity);
                double minPointSqDis2 = DISTANCE_SQ_THRESHOLD, minPointSqDis3 = DISTANCE_SQ_THRESHOLD;

                // search in the direction of increasing scan line
                for (int j = closestPointInd + 1; j < (int)mpPrev->mpSurfPointsFlat->points.size(); ++j)
                {
                    // if not in nearby scans, end the loop
                    if (int(mpPrev->mpSurfPointsFlat->points[j].intensity) > (closestPointScanID + NEARBY_SCAN))
                        break;

                    double pointSqDis = (mpPrev->mpSurfPointsFlat->points[j].x - pointSel.x) * (mpPrev->mpSurfPointsFlat->points[j].x - pointSel.x) +
                                        (mpPrev->mpSurfPointsFlat->points[j].y - pointSel.y) * (mpPrev->mpSurfPointsFlat->points[j].y - pointSel.y) +
                                        (mpPrev->mpSurfPointsFlat->points[j].z - pointSel.z) * (mpPrev->mpSurfPointsFlat->points[j].z - pointSel.z);

                    // if in the same or lower scan line
                    if (int(mpPrev->mpSurfPointsFlat->points[j].intensity) <= closestPointScanID && pointSqDis < minPointSqDis2)
                    {
                        minPointSqDis2 = pointSqDis;
                        minPointInd2 = j;
                    }
                    // if in the higher scan line
                    else if (int(mpPrev->mpSurfPointsFlat->points[j].intensity) > closestPointScanID && pointSqDis < minPointSqDis3)
                    {
                        minPointSqDis3 = pointSqDis;
                        minPointInd3 = j;
                    }
                }

                // search in the direction of decreasing scan line
                for (int j = closestPointInd - 1; j >= 0; --j)
                {
                    // if not in nearby scans, end the loop
                    if (int(mpPrev->mpSurfPointsFlat->points[j].intensity) < (closestPointScanID - NEARBY_SCAN))
                        break;

                    double pointSqDis = (mpPrev->mpSurfPointsFlat->points[j].x - pointSel.x) * (mpPrev->mpSurfPointsFlat->points[j].x - pointSel.x) +
                                        (mpPrev->mpSurfPointsFlat->points[j].y - pointSel.y) * (mpPrev->mpSurfPointsFlat->points[j].y - pointSel.y) +
                                        (mpPrev->mpSurfPointsFlat->points[j].z - pointSel.z) * (mpPrev->mpSurfPointsFlat->points[j].z - pointSel.z);

                    // if in the same or higher scan line
                    if (int(mpPrev->mpSurfPointsFlat->points[j].intensity) >= closestPointScanID && pointSqDis < minPointSqDis2)
                    {
                        minPointSqDis2 = pointSqDis;
                        minPointInd2 = j;
                    }
                    else if (int(mpPrev->mpSurfPointsFlat->points[j].intensity) < closestPointScanID && pointSqDis < minPointSqDis3)
                    {
                        // find nearer point
                        minPointSqDis3 = pointSqDis;
                        minPointInd3 = j;
                    }
                }

                if (minPointInd2 >= 0 && minPointInd3 >= 0)
                {

                    Eigen::Vector3d curr_point(mpSurfPointsFlat->points[i].x,
                                               mpSurfPointsFlat->points[i].y,
                                               mpSurfPointsFlat->points[i].z);
                    Eigen::Vector3d last_point_a(mpPrev->mpSurfPointsFlat->points[closestPointInd].x,
                                                 mpPrev->mpSurfPointsFlat->points[closestPointInd].y,
                                                 mpPrev->mpSurfPointsFlat->points[closestPointInd].z);
                    Eigen::Vector3d last_point_b(mpPrev->mpSurfPointsFlat->points[minPointInd2].x,
                                                 mpPrev->mpSurfPointsFlat->points[minPointInd2].y,
                                                 mpPrev->mpSurfPointsFlat->points[minPointInd2].z);
                    Eigen::Vector3d last_point_c(mpPrev->mpSurfPointsFlat->points[minPointInd3].x,
                                                 mpPrev->mpSurfPointsFlat->points[minPointInd3].y,
                                                 mpPrev->mpSurfPointsFlat->points[minPointInd3].z);

                    double s;
                    if (mLidarDistortion)
                        s = (mpSurfPointsFlat->points[i].intensity - int(mpSurfPointsFlat->points[i].intensity));
                    else
                        s = 1.0;

                    if (s > 1.0)
                    {
                        cout << s << endl;
                        throw std::runtime_error("Wrong scale of time");
                    }
                    ceres::CostFunction *cost_function = LidarPlaneFactor::Create(curr_point, last_point_a, last_point_b, last_point_c, s);
                    problem.AddResidualBlock(cost_function, loss_function, para_q, para_t);
                    numPlanarCorres++;

                    // fMatch << numPlanarCorres << "\t" << curr_point.transpose() << "\t " << mpSurfPointsFlat->points[i].intensity << endl;
                    // fMatch << numPlanarCorres << "\t" << last_point_a.transpose() << "\t " << mpPrev->mpSurfPointsFlat->points[closestPointInd].intensity << endl;
                    // fMatch << numPlanarCorres << "\t" << last_point_b.transpose() << "\t " << mpPrev->mpSurfPointsFlat->points[minPointInd2].intensity << endl;
                    // fMatch << numPlanarCorres << "\t" << last_point_c.transpose() << "\t " << mpPrev->mpSurfPointsFlat->points[minPointInd3].intensity << endl;
                }
            }
        }
        // fMatch.close();
        cout << "Number of corner/planar correspondence " << numCornerCorres << "/" << numPlanarCorres << endl;
        fLog << "\tIter " << iter <<" Number of corner/planar correspondence: " << numCornerCorres << "/" << numPlanarCorres << endl;

        //CONTROL:
        if (numCornerCorres + numPlanarCorres < 20)
        {
#ifdef EXPORT_LOG
            fLog << "\tError: Point matches is too few - " << numCornerCorres + numPlanarCorres << " at iter " << iter << endl;
            f_odometry_debug << id << "\tPoint matches is too few - " << numCornerCorres + numPlanarCorres << " at iter " << iter << endl;
#endif
            return false;
        }

        TicToc t_solver;
        ceres::Solver::Options options;
        options.linear_solver_type = ceres::DENSE_QR;
        options.max_num_iterations = 4;
        options.minimizer_progress_to_stdout = false;
        ceres::Solver::Summary summary;
        ceres::Solve(options, &problem, &summary);

        qTem = q_last_curr;
        tTem = t_last_curr;

        cout << "Refined " << tTem.transpose() << "\t" << rad2deg(Find_Rotation(qTem.toRotationMatrix())).transpose() << endl;

#ifdef EXPORT_LOG
        //fLog << "P:\t" << numCornerCorres << "/" << numPlanarCorres << "\t" << endl;
        fLog << "\tEstimated transformation: " << tTem.transpose() << "\t" << rad2deg(Find_Rotation(qTem.toRotationMatrix())).transpose() << endl;
#endif
    }

    r_lut2_lut1_pb = t_last_curr; // float(t_last_curr(0), t_last_curr(1), t_last_curr(2);
    R_lut2_lut1_pb = q_last_curr.toRotationMatrix();

    // CONTROL: based on the magnitude of result
    double distance = r_lut2_lut1_pb.norm();
    double angular_mag = Find_Rotation(R_lut2_lut1_pb).norm();
    if (distance > 0.5 || angular_mag > deg2rad(30.0))
    {
#ifdef EXPORT_LOG
        f_odometry_debug << id << "\tUnreasonable value for point-based odometry: distance - " << distance << " angles " << rad2deg(angular_mag) << endl;
        fLog << "Error: Unreasonable value for point-based odometry : distance - " << distance << " angles " << rad2deg(angular_mag) << endl;
#endif
        return false;
    }

    // //@@@ Here, it will be affected by the wrong value computed from this step
    // if (mLidarDistortion)
    // {
    //     for (int i = 0; i < numCorner; i++)
    //         TransformToEnd(&mpCornerPointsSharp->points[i], &mpCornerPointsSharp->points[i]);

    //     for (int i = 0; i < numPlanar; i++)
    //         TransformToEnd(&mpSurfPointsFlat->points[i], &mpSurfPointsFlat->points[i]);
    // }

    // if the feature-based approach is not valid
    if (!valid_feature_based_odo_flag_)
    {
        r_lut2_lut1 = r_lut2_lut1_pb;
        R_lut2_lut1 = R_lut2_lut1_pb;
    }
    return true;
}

/* Based on the setting, initialize required parameters:
mRangeTreshold, mNScan, mNFiring, mTolerateGap */
/**
 * @brief Initialize LiDAR scan parameters and precompute continuity thresholds.
 * 
 * Sets up sensor-specific configuration such as range thresholds, scan structure,
 * and computes acceptable range continuity using angular resolution.
 */

void LidarScan::init()
{
    // initial parameters
    // mTimeScan = 0.1;                                            //@ 0.1 second
    mRangeTreshold = make_pair(mPara.minRange, mPara.maxRange); // 1-100
    mNScan = mPara.nChannel;                                    //@ 32
    TreeAngThreshold = deg2rad(mPara.treeAngleThreshold);       // 10.0 - 20.0
    mGroundBufferTreeInit = mPara.groundBufferTree;
    mMinSegLength = mPara.minSegLength; // 2-3
    mSegLength = mPara.surfaceLength;   // 20 for classifying
    MaxTreeSegDistance = mPara.maxTreeSegDistance;
    MinTreeDistance = mPara.minTreeDistance;
    MinNumSegTree = mPara.minNumSegTree;
    MinTreePair = mPara.minTreePair;

    mbSimultaneously = mPara.bFPSimultanesouly;
    mbFeatureBased = mPara.bFeatureBased;
    mbLevel = mPara.bLevel;
    mbPointBased = mPara.bPointBased;
    mLidarDistortion = mPara.bDistortion;

    mTolerateGap = 3;

    mNPoints = single_scan_points_number_;
    mNFiring = mNPoints / mNScan;

    // horizontal angular res
    double horiAngularRes = 2.0 * PI / double(mNFiring);
    // cout << rad2deg(horiAngularRes) << endl;
    double minInterAnge = 10.0 * PI / 180.0;

    /*compute the triangle: angle between light ray and surface. Let the range of first as 1.
    Case 1: range 2 -> 180 - minInterAnge, range 1 -> 180 - horiAngularRes - (180 - minInterAnge)
    Case 2: range 2 -> minInterAnge, range 1 -> 180 - horiAngularRes -  minInterAnge
    */
    for (size_t nGap = 1; nGap <= mTolerateGap; nGap++)
    {
        double range2_upper = 1.0 / (sin(PI - horiAngularRes * nGap - (PI - minInterAnge))) * sin(PI - minInterAnge);
        double range2_lower = 1.0 / (sin(PI - horiAngularRes * nGap - minInterAnge)) * sin(minInterAnge);
        // cout << range2_lower << "\t" << range2_upper << endl;
        mvContRange.push_back(make_pair(range2_lower, range2_upper));
    }

    /*Compute reference trajectory information */
    if (mpTraj)
    {
        r_lu_b_ref = mpTraj->r_lu_b;
        R_lu_b_ref = mpTraj->R_lu_b;
        ground_height_ref = -1.4;

        //relative transformation from previous t_end to current t_end
        if (mpPrev)
        {
            Eigen::Vector3d r_b_m_t1, r_b_m_t2, r_lu_m_t1, r_lu_m_t2;
            Eigen::Matrix3d R_b_m_t1, R_b_m_t2, R_lu_m_t1, R_lu_m_t2;
            int idx;
            if (!(mpTraj->bopInterpolation((mpPrev->mTimeEnd) / 1000.0, -1, idx, r_b_m_t1, R_b_m_t1) && mpTraj->bopInterpolation(mTimeEnd / 1000.0, -1, idx, r_b_m_t2, R_b_m_t2)))
                throw std::runtime_error("Wrong Time");

            r_lu_m_t1 = r_b_m_t1 + R_b_m_t1 * r_lu_b_ref;
            R_lu_m_t1 = R_b_m_t1 * R_lu_b_ref;

            r_lu_m_t2 = r_b_m_t2 + R_b_m_t2 * r_lu_b_ref;
            R_lu_m_t2 = R_b_m_t2 * R_lu_b_ref;

            R_lut2_lut1_ref = R_lu_m_t1.inverse() * R_lu_m_t2;
            r_lut2_lut1_ref = R_lu_m_t1.inverse() * (r_lu_m_t2 - r_lu_m_t1);

            R_lu_m_ref = mpPrev->R_lu_m_ref * R_lut2_lut1_ref;
            r_lu_m_ref = mpPrev->r_lu_m_ref + mpPrev->R_lu_m_ref * r_lut2_lut1_ref;
        }
        //from t_ini to t_end
        else
        {
            Eigen::Vector3d r_b_m_t1, r_b_m_t2, r_lu_m_t1, r_lu_m_t2;
            Eigen::Matrix3d R_b_m_t1, R_b_m_t2, R_lu_m_t1, R_lu_m_t2;
            int idx;
            if (!(mpTraj->bopInterpolation(mTimeInit / 1000.0, -1, idx, r_b_m_t1, R_b_m_t1) && mpTraj->bopInterpolation(mTimeEnd / 1000.0, -1, idx, r_b_m_t2, R_b_m_t2)))
                throw std::runtime_error("Wrong Time");

            r_lu_m_t1 = r_b_m_t1 + R_b_m_t1 * r_lu_b_ref;
            R_lu_m_t1 = R_b_m_t1 * R_lu_b_ref;

            r_lu_m_t2 = r_b_m_t2 + R_b_m_t2 * r_lu_b_ref;
            R_lu_m_t2 = R_b_m_t2 * R_lu_b_ref;

            R_lut2_lut1_ref = R_lu_m_t1.inverse() * R_lu_m_t2;
            r_lut2_lut1_ref = R_lu_m_t1.inverse() * (r_lu_m_t2 - r_lu_m_t1);

            Eigen::Matrix3d R_temp = R_lu_m_ref * R_lut2_lut1_ref;
            R_lu_m_ref = R_temp;
            Eigen::Vector3d r_temp = r_lu_m_ref + R_lu_m_ref * r_lut2_lut1_ref;
            r_lu_m_ref = r_temp;
        }
    }
}

/* Compute attribute of the mvPoints.
for each point, prev and next continuous point
segments of each channle.
b*/
/**
 * @brief Compute attributes such as continuity, segmentation, and smoothness for each point.
 * 
 * Connects neighboring valid points, assigns segment IDs, and calculates smoothness values
 * to assist in edge/surface classification.
 */

void LidarScan::computeAttribute()
{
    int numInvalidPoint = 0;
    for (std::size_t i = 0; i < mNPoints; i++)
    {
        size_t rowIdn, columnIdn;
        rowIdn = i % mNScan;
        columnIdn = int(i / mNScan);
        mpLaserCloud->points[i].intensity = mpLaserCloud->points[i].intensity + float(rowIdn);
        LidarPoint thisPoint(mpLaserCloud->points[i], rowIdn, columnIdn);

        // cout << thisPoint.range << "\t" << thisPoint.point.x << endl;
        // thisPoint.point.x = 10.0;
        // cout << mpLaserCloud->points[i].x << "\t" << thisPoint.point.x << endl;

        if (thisPoint.range < mRangeTreshold.first || thisPoint.range > mRangeTreshold.second)
            thisPoint.bValid = false;
        else
        {
            thisPoint.bValid = true;
            numInvalidPoint++;
        }
        mvPoints.push_back(thisPoint);
    }

    cout << "Number of valid point " << numInvalidPoint << " out of " << mNPoints << endl;
    fLog << "Preprocess:" <<endl;
    fLog << "\tNumber of point with valid range information " << numInvalidPoint << " out of " << mNPoints << endl;

    // compute continuity along each channel
    for (std::size_t nFiring = 0; nFiring < mNFiring; nFiring++)
    {
        for (std::size_t nChannel = 0; nChannel < mNScan; nChannel++)
        {
            // pointer to the points
            LidarPoint *pThisPoint;
            LidarPoint *pNextPoint;

            int index = nFiring * mNScan + nChannel;
            pThisPoint = &(mvPoints[index]);
            // if the range is not valid
            if (!pThisPoint->bValid)
                continue;

            int nGap = 1;
            while (nGap <= mTolerateGap)
            {
                int nextFiring = (nFiring + nGap) % mNFiring; // firing ID
                int indexNext = nextFiring * mNScan + nChannel;
                pNextPoint = &(mvPoints[indexNext]);

                if (!pNextPoint->bValid) // range not valid
                {
                    nGap++;
                    continue;
                }

                // compute range ratio
                float ratio = pThisPoint->range / pNextPoint->range;

                // if range are continous
                if (ratio > mvContRange[nGap - 1].first && ratio < mvContRange[nGap - 1].second)
                {
                    if (pNextPoint->neighborPrevious != -1) // have already been linked. Delete the link for the previsous. As a result, there is no ambiguity
                    {
                        mvPoints[pNextPoint->neighborPrevious].neighborNext = -1;
                    }

                    pThisPoint->neighborNext = indexNext;
                    pNextPoint->neighborPrevious = index;
                    Eigen::Vector3f vec(pNextPoint->point.x - pThisPoint->point.x, pNextPoint->point.y - pThisPoint->point.y,
                                        pNextPoint->point.z - pThisPoint->point.z);
                    pThisPoint->vecNext = vec;
                    pNextPoint->vecPrev = vec;
                    break;
                }
                else
                {
                    nGap++;
                    continue;
                }
            }
        }
    }

    // compute segment ID
    for (std::size_t nChannel = 0; nChannel < mNScan; nChannel++)
    {
        int segId = 0;
        std::vector<LidarSegment> vecSegs; // vector of a segment in a scan line
        for (std::size_t nFiring = 0; nFiring < mNFiring; nFiring++)
        {
            LidarPoint *pThisPoint;
            int index = nFiring * mNScan + nChannel;
            pThisPoint = &(mvPoints[index]);

            int indexPrev = pThisPoint->neighborPrevious;
            int indexNext = pThisPoint->neighborNext;

            if (indexPrev == -1 && indexNext != -1) // starting point
            {
                // initialize a lidar segment
                LidarSegment lidarSeg(nChannel, mMinSegLength, mSegLength, this);
                lidarSeg.vecIndex.push_back(index);

                // fill next point
                do
                {
                    lidarSeg.vecIndex.push_back(indexNext);
                    indexNext = mvPoints[indexNext].neighborNext;

                } while (indexNext != -1);

                // compute attribute of segment. classify
                lidarSeg.computeSegAttribute();
                if (lidarSeg.valid())
                {
                    // lidarSeg.computeNormDirection();
                    vecSegs.push_back(lidarSeg);

                    // assign segment ID and type to the mvPoints
                    for (std::size_t np = 0; np < lidarSeg.vecIndex.size(); np++)
                    {
                        mvPoints[lidarSeg.vecIndex[np]].segID = segId;
                        mvPoints[lidarSeg.vecIndex[np]].seg_init_type = lidarSeg.initType;
                    }

                    segId++;
                }
            }
        }

        mvSegments.push_back(vecSegs);
    }

    // compute smoothness
    for (int nFiring = 0; nFiring < mNFiring; nFiring++)
    {
        for (int nChannel = 0; nChannel < mNScan; nChannel++)
        {
            //  cout << nFiring << "-" << nChannel << endl;
            // pointer to the points
            LidarPoint *pThisPoint;
            LidarPoint *pNextPoint;

            int index = nFiring * mNScan + nChannel;
            pThisPoint = &(mvPoints[index]);
            // if the range is not valid
            if (!pThisPoint->bValid)
                continue;

            float diffX = -float(N_Smooth * 2) * pThisPoint->point.x,
                  diffY = -float(N_Smooth * 2) * pThisPoint->point.y,
                  diffZ = -float(N_Smooth * 2) * pThisPoint->point.z;

            bool valid = true;
            int n = 0;
            int df = 0;
            // find the neighboring valid measurments for computing smoothness
            // if the points is too away from the current point (df too large), fail
            while (n < N_Smooth)
            {
                df++;
                if (df > 20)
                {
                    valid = false;
                    break;
                }
                int tempFiring = (nFiring + df) < mNFiring ? nFiring + df : nFiring + df - mNFiring;
                int tempIndex = tempFiring * mNScan + nChannel;
                //     cout << "A Cur firing:" << nFiring << "\t"
                //                     << "df " << df << "\t" << tempIndex << "/" << mvPoints.size() << "\t" << tempFiring << "\t" << nChannel << endl;
                pNextPoint = &(mvPoints[tempIndex]);
                if (pNextPoint->bValid)
                {
                    n++;
                    diffX = diffX + pNextPoint->point.x;
                    diffY = diffY + pNextPoint->point.y;
                    diffZ = diffZ + pNextPoint->point.z;
                }
            }
            n = 0;
            df = 0;
            while (n < N_Smooth)
            {
                df++;
                if (df > 20)
                {
                    valid = false;
                    break;
                }
                int tempFiring = (nFiring - df) < 0 ? nFiring - df + mNFiring : nFiring - df;
                int tempIndex = tempFiring * mNScan + nChannel;
                //    cout << "B Cur firing:" << nFiring << "\t" << tempIndex << "/" << mvPoints.size() << "\t" << tempFiring << "\t" << nChannel << endl;

                pNextPoint = &(mvPoints[tempIndex]);

                if (pNextPoint->bValid)
                {
                    n++;
                    diffX = diffX + pNextPoint->point.x;
                    diffY = diffY + pNextPoint->point.y;
                    diffZ = diffZ + pNextPoint->point.z;
                }
            }

            if (valid)
                pThisPoint->smooth = diffX * diffX + diffY * diffY + diffZ * diffZ;
            else
                pThisPoint->smooth = -2.0;
        }
    }

    for (int nChannel = 0; nChannel < mNScan; nChannel++)
    {
        for (std::size_t nSeg = 0; nSeg < mvSegments[nChannel].size(); nSeg++)
        {
            mvSegments[nChannel][nSeg].computePlaneAttribute(Smooth_Threshold); //@@@@
        }
    }

    // system("read -p 'Press Enter to continue...' var");
}

/*Point-based approach: extract planar and edge points
True: Enough number of points are extracted
False: else
*/
/**
 * @brief Extract sharp corner and flat surface features from the point cloud.
 * 
 * Segments are further classified into sharp, less sharp, flat, and less flat features 
 * based on curvature. Features are stored in respective point clouds for downstream odometry.
 * 
 * @return true if a sufficient number of features are found.
 */

bool LidarScan::extract_edge_planar_points()
{
    // cout << "Compute Smooth " << endl;
    fLog << "Extract planar / edge points: " <<endl;
    mpCurv = new float[mNPoints];

    vector<int> mvValidIndex; // index in the mvPoints

    int *cloudSortInd = new int[mNPoints];        // index for used points
    int *cloudNeighborPicked = new int[mNPoints]; // curvature array of all point.
    int *cloudLabel = new int[mNPoints];          // curvature array of all point.

    vector<int> startIndex;
    vector<int> endIndex;

    int count = 0;
    // for each channel, find valid point and save the start and end index
    for (int nChannel = 0; nChannel < mNScan; nChannel++)
    {
        startIndex.push_back(count + N_Smooth);
        for (int nFiring = 0; nFiring < mNFiring; nFiring++)
        {
            int index = nFiring * mNScan + nChannel;
            if (mvPoints[index].bValid && mvPoints[index].segID >= 0)
            {
                mvValidIndex.push_back(index);
                count++;
            }
        }
        endIndex.push_back(max(0, count - N_Smooth - 1));
    }

    // compute the curvature, only valid for points within the start and end index
    for (int nP = N_Smooth; nP < count - N_Smooth; nP++)
    {
        LidarPoint *pThisPoint;
        LidarPoint *pNextPoint1;
        LidarPoint *pNextPoint2;

        pThisPoint = &(mvPoints[mvValidIndex[nP]]);

        float diffX = -float(N_Smooth * 2) * pThisPoint->point.x,
              diffY = -float(N_Smooth * 2) * pThisPoint->point.y,
              diffZ = -float(N_Smooth * 2) * pThisPoint->point.z;
        for (int nN = 1; nN <= N_Smooth; nN++)
        {
            pNextPoint1 = &(mvPoints[mvValidIndex[nP - nN]]);
            pNextPoint2 = &(mvPoints[mvValidIndex[nP + nN]]);

            diffX = diffX + pNextPoint1->point.x + pNextPoint2->point.x;
            diffY = diffY + pNextPoint1->point.y + pNextPoint2->point.y;
            diffZ = diffZ + pNextPoint1->point.z + pNextPoint2->point.z;
        }

        mpCurv[nP] = diffX * diffX + diffY * diffY + diffZ * diffZ;
        cloudNeighborPicked[nP] = 0;
        cloudSortInd[nP] = nP;
        pThisPoint->smooth2 = diffX * diffX + diffY * diffY + diffZ * diffZ;
    }

    // find planar point and edge point
    // cout << "Find POints " << endl;

    mpCornerPointsSharp.reset(new pcl::PointCloud<PointType>());
    mpCornerPointsLessSharp.reset(new pcl::PointCloud<PointType>());
    mpSurfPointsFlat.reset(new pcl::PointCloud<PointType>());
    mpSurfPointsLessFlat.reset(new pcl::PointCloud<PointType>());

    // pcl::PointCloud<PointType> cornerPointsSharp;
    // pcl::PointCloud<PointType> cornerPointsLessSharp;
    // pcl::PointCloud<PointType> surfPointsFlat;
    // pcl::PointCloud<PointType> surfPointsLessFlat;

    for (int nChannel = 0; nChannel < mNScan; nChannel++)
    {
        // if number of points are too few
        int numChannelPoint = endIndex[nChannel] - startIndex[nChannel] + 1;
        if (numChannelPoint < 40)
            continue;
        pcl::PointCloud<PointType>::Ptr surfPointsLessFlatScan(new pcl::PointCloud<PointType>);

        int numSubArea = max(std::min(numChannelPoint / MinSubPoints, MaxSubArea), 1);
        // cout << numSubArea << endl;
        // cout << startIndex[nChannel] << " to " << endIndex[nChannel] << endl;
        for (int nArea = 0; nArea < numSubArea; nArea++)
        {
            int sp = startIndex[nChannel] + numChannelPoint * nArea / numSubArea;
            int ep = startIndex[nChannel] + numChannelPoint * (nArea + 1) / numSubArea - 1;

            // cout << sp << "-" << ep << endl;
            //  for (int k = sp; k <= ep; k++)
            //      cout << cloudSortInd[k] << " ";
            //                  cout << endl;

            // for (int k = sp; k <= ep; k++)
            //     cout << mpCurv[k] << " ";
            // cout << endl;

            // for (int k = sp; k <= ep; k++)
            //  cout << cloudSortInd[k] << " ";
            //  cout << endl;
            //  system("read -p 'Press Enter to continue...' var");

            // sort index from minimum to maximum
            std::sort(cloudSortInd + sp, cloudSortInd + ep + 1, comp);
            // cout << "sort" << endl;

            // edge points
            int largestPickedNum = 0;
            for (int k = ep; k >= sp; k--)
            {
                int ind = cloudSortInd[k];
                // cout << ind << "-" << count << "\t" << mvValidIndex.size() << "-" << mvValidIndex[ind]<< endl;
                if (cloudNeighborPicked[ind] == 0 && mpCurv[ind] > Smooth_Threshold)
                {
                    // cout << 1 << endl;

                    largestPickedNum++;
                    if (largestPickedNum <= NumSubEdge)
                    {
                        cloudLabel[ind] = 2;
                        mvPoints[mvValidIndex[ind]].point_type_id = 2;
                        // cornerPointsSharp.push_back(mvPoints[mvValidIndex[ind]].point);
                        // cornerPointsLessSharp.push_back(mvPoints[mvValidIndex[ind]].point);

                        mpCornerPointsSharp->push_back(mvPoints[mvValidIndex[ind]].point);
                        mpCornerPointsLessSharp->push_back(mvPoints[mvValidIndex[ind]].point);
                    }
                    else if (largestPickedNum <= NumSubEdgeLess)
                    {
                        cloudLabel[ind] = 1;
                        mvPoints[mvValidIndex[ind]].point_type_id = 1;
                        // cornerPointsLessSharp.push_back(mvPoints[mvValidIndex[ind]].point);
                        mpCornerPointsLessSharp->push_back(mvPoints[mvValidIndex[ind]].point);
                    }
                    else
                        break;

                    cloudNeighborPicked[ind] = 1;

                    // cout << 2 << endl;

                    for (int l = 1; l <= N_Smooth; l++)
                    {
                        float diffX = mvPoints[mvValidIndex[ind + l]].point.x - mvPoints[mvValidIndex[ind + l - 1]].point.x;
                        float diffY = mvPoints[mvValidIndex[ind + l]].point.y - mvPoints[mvValidIndex[ind + l - 1]].point.y;
                        float diffZ = mvPoints[mvValidIndex[ind + l]].point.z - mvPoints[mvValidIndex[ind + l - 1]].point.z;

                        if (diffX * diffX + diffY * diffY + diffZ * diffZ > 0.05)
                            break;

                        cloudNeighborPicked[ind + l] = 1;
                    }
                    //     cout << 3 << endl;

                    for (int l = -1; l >= -N_Smooth; l--)
                    {
                        float diffX = mvPoints[mvValidIndex[ind + l]].point.x - mvPoints[mvValidIndex[ind + l + 1]].point.x;
                        float diffY = mvPoints[mvValidIndex[ind + l]].point.y - mvPoints[mvValidIndex[ind + l + 1]].point.y;
                        float diffZ = mvPoints[mvValidIndex[ind + l]].point.z - mvPoints[mvValidIndex[ind + l + 1]].point.z;
                        if (diffX * diffX + diffY * diffY + diffZ * diffZ > 0.05)
                            break;

                        cloudNeighborPicked[ind + l] = 1;
                    }
                }
            }
            // cout << "edge" << endl;

            // planar points
            int smallestPickedNum = 0;
            for (int k = sp; k <= ep; k++)
            {
                int ind = cloudSortInd[k];
                if (cloudNeighborPicked[ind] == 0 && mpCurv[ind] < Smooth_Threshold)
                {
                    cloudLabel[ind] = -1;
                    mvPoints[mvValidIndex[ind]].point_type_id = 3;
                    //                    surfPointsFlat.push_back(mvPoints[mvValidIndex[ind]].point);
                    mpSurfPointsFlat->push_back(mvPoints[mvValidIndex[ind]].point);

                    smallestPickedNum++;
                    if (smallestPickedNum >= 4)
                        break;

                    cloudNeighborPicked[ind] = 1;
                    for (int l = 1; l <= N_Smooth; l++)
                    {
                        float diffX = mvPoints[mvValidIndex[ind + l]].point.x - mvPoints[mvValidIndex[ind + l - 1]].point.x;
                        float diffY = mvPoints[mvValidIndex[ind + l]].point.y - mvPoints[mvValidIndex[ind + l - 1]].point.y;
                        float diffZ = mvPoints[mvValidIndex[ind + l]].point.z - mvPoints[mvValidIndex[ind + l - 1]].point.z;
                        if (diffX * diffX + diffY * diffY + diffZ * diffZ > 0.05)
                            break;

                        cloudNeighborPicked[ind + l] = 1;
                    }
                    for (int l = -1; l >= -N_Smooth; l--)
                    {
                        float diffX = mvPoints[mvValidIndex[ind + l]].point.x - mvPoints[mvValidIndex[ind + l + 1]].point.x;
                        float diffY = mvPoints[mvValidIndex[ind + l]].point.y - mvPoints[mvValidIndex[ind + l + 1]].point.y;
                        float diffZ = mvPoints[mvValidIndex[ind + l]].point.z - mvPoints[mvValidIndex[ind + l + 1]].point.z;
                        if (diffX * diffX + diffY * diffY + diffZ * diffZ > 0.05)
                            break;

                        cloudNeighborPicked[ind + l] = 1;
                    }
                }
            }
            // cout << "planar" << endl;

            // flat points for map
            for (int k = sp; k <= ep; k++)
            {
                int ind = cloudSortInd[k];
                if (mpCurv[ind] < Smooth_Threshold)
                {
                    surfPointsLessFlatScan->push_back(mvPoints[mvValidIndex[ind]].point);
                    if (mvPoints[mvValidIndex[ind]].point_type_id != 3)
                        mvPoints[mvValidIndex[ind]].point_type_id = 4;
                }
                else
                    break;
            }
        }

        pcl::PointCloud<PointType> surfPointsLessFlatScanDS;
        pcl::VoxelGrid<PointType> downSizeFilter;
        downSizeFilter.setInputCloud(surfPointsLessFlatScan);
        downSizeFilter.setLeafSize(0.2, 0.2, 0.2);
        downSizeFilter.filter(surfPointsLessFlatScanDS);
        *mpSurfPointsLessFlat += surfPointsLessFlatScanDS;
    }

    cout << "Valid point for smooth: " << count << "/" << mNPoints << endl;
    cout << "Number of edge points: " << mpCornerPointsSharp->size() << "/" << mpCornerPointsLessSharp->size() << endl;
    cout << "Number of planar points: " << mpSurfPointsFlat->size() << "/" << mpSurfPointsLessFlat->size() << endl;

    fLog << "\tValid point for smooth: " << count << "/" << mNPoints << endl;
    fLog << "\tNumber of edge points: " << mpCornerPointsSharp->size() << "/" << mpCornerPointsLessSharp->size() << endl;
    fLog << "\tNumber of planar points: " << mpSurfPointsFlat->size() << "/" << mpSurfPointsLessFlat->size() << endl;
    delete []mpCurv;

    //set the flag to true
    point_extraction_flag_ = true;

    if (mpCornerPointsSharp->size() > 20 && mpSurfPointsFlat->size() > 100)
    {
        valid_point_flag_ = true;
        return true;
    }
    else
    {
#ifdef EXPORT_LOG
        f_odometry_debug << id << "\tInsufficient number of points are extracted" << endl;
        fLog << "Error: Insufficient number of points are extracte" << endl;
#endif
        valid_point_flag_ = false;
        return false;
    }
    // for (int i = 0; i < startIndex.size(); i++)
    //     cout << startIndex[i] << "-" << endIndex[i] << "\t";
    // cout << endl;
    // system("read -p 'Press Enter to continue...' var");
}

/************************************
Extract ground clusters based on ground candidate segments. Three cases:
1. The derived cluster is reliable, ground information is estimated
2. No cluster or incompatible ground info, use the ground info from previous scan
3. If previous scan with valid ground is not valid, return false;
Input:
    bFinal: save the extracted ground cluster as final extraction or not.
Output:
    mGroundH: ground height
    R_lu_lup: leveling rotation
Return:
    true: ground information valid (either by itself or from previous scan)
    false: no ground information valid
mbTrackFeature: fail when false
mbGroundPlane: fail when false
************************************/
/**
 * @brief Extracts ground clusters by segmenting lower scan lines and fitting planes.
 * 
 * Uses geometric relationships and PCA fitting to isolate ground points.
 * Computes leveling rotation and filters out outlier segments.
 * 
 * @param bFinal If true, assigns ground labels to the final feature map.
 * @return true if valid ground information is extracted or retrieved from the previous scan.
 */

bool LidarScan::extractGroundCluster(bool bFinal)
{

    fLog << "Extract ground clusters:"<<endl;
    std::vector<ClusterSegment> candCluster;
    int countClusterID = 0;

    //in this code, the scan line with largest id is closest to the ground
    //1. find all clusters for the long segments: candCluster, segment.clusterID is changed
    for (std::size_t nChannel = mNScan - 1; nChannel > 0; nChannel--)
    {
        int numSeg1 = mvSegments[nChannel].size();
        int numSeg2 = mvSegments[nChannel - 1].size();

        int nSeg1 = 0, nSeg2 = 0;
        while (nSeg1 < numSeg1 && nSeg2 < numSeg2)
        {
            int s1 = mvSegments[nChannel][nSeg1].minFiring;
            int e1 = mvSegments[nChannel][nSeg1].maxFiring;
            int l1 = mvSegments[nChannel][nSeg1].length;
            int type1 = mvSegments[nChannel][nSeg1].initType;

            //if not the long segment, skip
            if (type1 != 1)
            {
                nSeg1++;
                continue;
            }

            // if e1 < s1, line seg pass the 0, so e1 = e1 + mnFiring
            e1 = e1 < s1 ? e1 + mNFiring : e1;

            int s2 = mvSegments[nChannel - 1][nSeg2].minFiring;
            int e2 = mvSegments[nChannel - 1][nSeg2].maxFiring;
            int l2 = mvSegments[nChannel - 1][nSeg2].length;
            int type2 = mvSegments[nChannel - 1][nSeg2].initType;
            if (type2 != 1) //!= 1
            {
                nSeg2++;
                continue;
            }
            e2 = e2 < s2 ? e2 + mNFiring : e2;

            if (e1 < s2) // seg1 lower than s2
            {
                nSeg1++;
                continue;
            }
            else if (e2 < s1)
            {
                nSeg2++;
                continue;
            }
            else // common area
            {
                if (type1 != type2) // not same type
                {
                    throw std::runtime_error("Wrong type");
                    e1 < e2 ? nSeg1++ : nSeg2++;
                    continue;
                }

                int overlap = min(e1, e2) - max(s1, s2);
                if (float(overlap) / float(min(l1, l2)) > 0.3) // enough overlap
                {
                    if (mvSegments[nChannel][nSeg1].ground_cluster_id == -1) // new cluster
                    {
                        // initialize a new cluster
                        ClusterSegment cluster(1);
                        cluster.vSegments.push_back(make_pair(nChannel, nSeg1));
                        cluster.vSegments.push_back(make_pair(nChannel - 1, nSeg2));
                        cluster.featureId = countClusterID;
                        candCluster.push_back(cluster);

                        mvSegments[nChannel][nSeg1].ground_cluster_id = countClusterID;
                        mvSegments[nChannel - 1][nSeg2].ground_cluster_id = countClusterID;
                        countClusterID++;
                    }
                    else
                    {
                        mvSegments[nChannel - 1][nSeg2].ground_cluster_id = mvSegments[nChannel][nSeg1].ground_cluster_id;
                        candCluster[mvSegments[nChannel][nSeg1].ground_cluster_id].vSegments.push_back(make_pair(nChannel - 1, nSeg2));
                    }
                    // if found the matches. the candidate seg + 1
                    // e1 < e2 ? nSeg1++ : nSeg2++;
                    nSeg2++;
                    continue;
                }
                else // not enough overlap
                {
                    //    e1 < e2 ? nSeg1++ : nSeg2++;
                    nSeg2++;
                    continue;
                }
            }
        }
    }
    
    //2. remove clusters with few numbers of segments: ground_cluster_
    int i = 0;
    for (size_t nC = 0; nC < candCluster.size(); nC++)
    {
        if (candCluster[nC].vSegments.size() >= MinNumSegGround)
        {
            ground_cluster_.push_back(candCluster[nC]);
            ground_cluster_[i].len = ground_cluster_[i].vSegments.size();

            for(int n_seg = 0; n_seg < ground_cluster_[i].len; ++n_seg)
            {
                mvSegments[ground_cluster_[i].vSegments[n_seg].first][ground_cluster_[i].vSegments[n_seg].second].outlier_removal_flag = 10;
            }
            i++;
        }
    }
    cout << "Number of initial ground-based clusters and the ones with enough segments " << candCluster.size() << "/" << ground_cluster_.size() <<endl;
    fLog << "\tNumber of ground clusters - initial/with enough segments " << candCluster.size() << "/" << ground_cluster_.size() <<endl;

    //3. remove clusters with number of lidar points < 1000, and the clusters who doesnt constitude a plane.
    // if plane exists, save the normal vector
    vector<float> vAlpha;
    vector<float> vBeta;
    vector<Eigen::Vector3f> cand_normal_vec;
    std::vector<ClusterSegment> vTempCluster;
    std::vector<int> vNumberPoints;
    int num_failed_clusters_point_num = 0, num_failed_clusters_fitting = 0, num_failed_clusters_fitting_pca = 0;
    for (size_t nC = 0; nC < ground_cluster_.size(); nC++)
    {
        std::vector<Eigen::Vector3f> allPoints;
        int numPoints = 0;
        for (size_t nSeg = 0; nSeg < ground_cluster_[nC].len; nSeg++)
        {
            int channel = ground_cluster_[nC].vSegments[nSeg].first;
            int segID = ground_cluster_[nC].vSegments[nSeg].second;
            allPoints.insert(allPoints.end(), mvSegments[channel][segID].vecPlanePoints.begin(), mvSegments[channel][segID].vecPlanePoints.end());
            numPoints += mvSegments[channel][segID].vecPlanePoints.size();
        }
        if (allPoints.size() < mPara.min_points_per_ground_cluster) // size limitation
        {
            num_failed_clusters_point_num++;
            for(int n_seg = 0; n_seg < ground_cluster_[nC].len; ++n_seg)
            {
                mvSegments[ground_cluster_[nC].vSegments[n_seg].first][ground_cluster_[nC].vSegments[n_seg].second].outlier_removal_flag = 1;
            }
            continue;
        }

        Eigen::Vector4f params;
        int fail_type = 0;
        if (planeFitting_outlier_check(allPoints, params, ground_cluster_[nC].vCenterPlane,fail_type ))
        {
            Eigen::Vector3f normVec(params(0), params(1), params(2));
            float alpha, beta;
            norm2angle(normVec, alpha, beta);
            vAlpha.push_back(alpha);
            vBeta.push_back(beta);
            vTempCluster.push_back(ground_cluster_[nC]);
            vNumberPoints.push_back(numPoints);
            cand_normal_vec.push_back(normVec);
        }
        else
        {
            for(int n_seg = 0; n_seg < ground_cluster_[nC].len; ++n_seg)
            {
                //fail in pca: 1, fail in residuals: 2
                mvSegments[ground_cluster_[nC].vSegments[n_seg].first][ground_cluster_[nC].vSegments[n_seg].second].outlier_removal_flag = fail_type + 1;
            }
            num_failed_clusters_fitting++;
            if(fail_type == 1) num_failed_clusters_fitting++;
        }
    }
    ground_cluster_ = vTempCluster;
    cout << "Number of valid ground clusters after plane fitting " << ground_cluster_.size() << endl;
    fLog << "\tFailed clusters in number of points/plane fitting (pca): " << num_failed_clusters_point_num << "/" << num_failed_clusters_fitting
         << "(" << num_failed_clusters_fitting << "), number of remaining: " << ground_cluster_.size() << endl;

    // CONTROL: if no valid ground cluster is extracted
    if (ground_cluster_.size() < 1)
    {
        // if no ground information from previous scan
        if (!mpPrev || !mpPrev->valid_ground_info_flag_)
        {

#ifdef EXPORT_LOG
            f_odometry_debug << id << "\tExtracting ground plane information failed: " << ground_cluster_.size() << endl;
            fLog << "Error: no valid plane info" <<endl;
#endif
            return false;
        }
        else
        {
            R_lu_lup = mpPrev->R_lu_lup;
            mGroundH = mpPrev->mGroundH;

#ifdef EXPORT_LOG
            fLog << "\tNo valid ground cluster, used ground info from previous scan " << endl;
            f_odometry_debug << id << "\tNo valid ground cluster, used ground info from previous scan " << endl;
#endif
            return true;
        }
    }

    // 4. Derive the rotation matrix using normal vector consuses
    // sufficient number of clusters, take the medium. Else, use the one with most number of points
    float phi, ome;
    float ang_threshold = deg2rad(20.0);
    int best_index = -1;
    int best_conseseus = -1;
    // for each cluster, compute the angle between current and remaining
    for (int i = 0; i < cand_normal_vec.size(); ++i)
    {
        int num_consesus = 0;
        for (int j = 0; j < cand_normal_vec.size(); ++j)
        {
            if (i == j)
                continue;

            float angle = computeAngle(cand_normal_vec[i], cand_normal_vec[j]);
            // f_odometry_debug <<id <<": "<<endl;
            // f_odometry_debug << ground_cluster_[i].featureId << " - " << ground_cluster_[j].featureId <<": " <<rad2deg(angle)<< "\t";

            if (angle < ang_threshold)
            {
                num_consesus++;
            }
        }
        // f_odometry_debug <<"consesus " << num_consesus <<endl;
        if (num_consesus > best_conseseus ||
            (num_consesus == best_conseseus && vNumberPoints[i] > vNumberPoints[best_index])) // if same consesus, the one with the more number of points
        {
            best_conseseus = num_consesus;
            best_index = i;
        }
        // f_odometry_debug<< "Best cluster: " << ground_cluster_[best_index].featureId << " " << best_conseseus <<endl;
    }
    Eigen::Vector3f best_normal_vector = cand_normal_vec[best_index];
    phi = asin(best_normal_vector(0));
    ome = atan2(-best_normal_vector(1), best_normal_vector(2));
    R_lu_lup = Compute_Rotation(ome, phi, float(0.0));
    R_lu_lup.transposeInPlace();

    // 5. remove clusters with incompatible normal vector
    vTempCluster.clear();
    for (size_t nC = 0; nC < ground_cluster_.size(); nC++)
    {
        if (computeAngle(cand_normal_vec[nC], best_normal_vector) < ang_threshold)
        {
            vTempCluster.push_back(ground_cluster_[nC]);
        }
        else
        {
            for (int n_seg = 0; n_seg < ground_cluster_[nC].len; ++n_seg)
            {
                mvSegments[ground_cluster_[nC].vSegments[n_seg].first][ground_cluster_[nC].vSegments[n_seg].second].outlier_removal_flag = 4;
            }
        }
    }
    ground_cluster_ = vTempCluster;
    cout << "Number of valid ground-based clusters after normal check " << vTempCluster.size() << endl;
    fLog << "\tNumber of valid ground-based clusters after normal check " << vTempCluster.size() << endl;

    // 6. derive ground height and remove clusters based on ground height
    vector<float> vH;
    for (size_t nC = 0; nC < ground_cluster_.size(); nC++)
    {
        ground_cluster_[nC].vCenterPlaneTrans = R_lu_lup * ground_cluster_[nC].vCenterPlane;
        vH.push_back(ground_cluster_[nC].vCenterPlaneTrans(2));
    }
    mGroundH = *min_element(vH.begin(), vH.end());
    vTempCluster.clear();
    for (size_t nC = 0; nC < ground_cluster_.size(); nC++)
    {
        if ((ground_cluster_[nC].vCenterPlaneTrans(2) - mGroundH) < mGroundBuffer)
        {
            vTempCluster.push_back(ground_cluster_[nC]);
        }
        else
        {
            for(int n_seg = 0; n_seg < ground_cluster_[nC].len; ++n_seg)
            {
                mvSegments[ground_cluster_[nC].vSegments[n_seg].first][ground_cluster_[nC].vSegments[n_seg].second].outlier_removal_flag = 5;
            }
        }
    }

    ground_cluster_ = vTempCluster;
    cout << "Number of final valid ground-based clusters " << vTempCluster.size() << endl;
    fLog << "\tNumber of valid ground-based clusters with reasonale height " << vTempCluster.size() << endl;

#ifdef TEST_PURPOSE
    fInfo << "Groud height: " << mGroundH << " , ome/phi: " << rad2deg(ome) << " " << rad2deg(phi) << endl;
    //  cout << "Ground buffer: " << mGroundH << " " << mGroundH + mGroundBuffer << endl;
    // system("read -p 'Press Enter to continue...' var");
#endif

#ifdef EXPORT_LOG
    fLog << "\tGround info for best cluster:\t" << ground_cluster_.size() << "\t" << mGroundH << "\t" << rad2deg(ome) << "\t" << rad2deg(phi) << endl;
#endif

    // assign attribute to mvPoints
    if (bFinal)
    {
        for (size_t nC = 0; nC < ground_cluster_.size(); nC++)
        {
            for (size_t nSeg = 0; nSeg < ground_cluster_[nC].len; nSeg++)
            {
                int channelID = ground_cluster_[nC].vSegments[nSeg].first;
                int segID = ground_cluster_[nC].vSegments[nSeg].second;

                mvSegments[channelID][segID].final_feature_id = nC;
                mvSegments[channelID][segID].final_feature_type = 1;

                for (std::size_t nP = 0; nP < mvSegments[channelID][segID].length; nP++)
                {
                    int index = mvSegments[channelID][segID].vecIndex[nP];
                    mvPoints[index].clusterID = nC;
                    mvPoints[index].clusterType = 1;
                }
            }
        }
    }

    //------------------------------------------------
    // another round of plane fitting
    std::vector<Eigen::Vector3f> allPointsClusters;
    for (size_t nC = 0; nC < ground_cluster_.size(); nC++)
    {
        for (size_t nSeg = 0; nSeg < ground_cluster_[nC].len; nSeg++)
        {
            int channel = ground_cluster_[nC].vSegments[nSeg].first;
            int segID = ground_cluster_[nC].vSegments[nSeg].second;
            allPointsClusters.insert(allPointsClusters.end(), mvSegments[channel][segID].vecPlanePoints.begin(), mvSegments[channel][segID].vecPlanePoints.end());
        }
    }
    Eigen::Vector4f paramsFinal;
    Eigen::Vector3f centerFinal;
    planeFitting(allPointsClusters, paramsFinal, centerFinal);

    Eigen::Vector3f normVecFinal(paramsFinal(0), paramsFinal(1), paramsFinal(2));
    float phi2 = asin(normVecFinal(0));
    float ome2 = atan2(-normVecFinal(1), normVecFinal(2));
    Eigen::Matrix3f mRlevel;
    mRlevel = Compute_Rotation(ome2, phi2, float(0.0));
    mRlevel.transposeInPlace();
    Eigen::Vector3f p(0.0, 0.0, -paramsFinal(3) / paramsFinal(2));
    Eigen::Vector3f p_trans = mRlevel * p;

#ifdef EXPORT_LOG
    fLog << "\tGround info using all clusters:\t" << ground_cluster_.size() << "\t" << p_trans(2) << "\t" << rad2deg(ome2) << "\t" << rad2deg(phi2) << endl;
#endif
    mGroundH = p_trans(2);
    R_lu_lup = mRlevel;

    mOme = ome2;
    mPhi = phi2;

    // CONTROL: if ground information from previous scan is valid, check consistency
    if (mpPrev && mpPrev->valid_ground_info_flag_)
    {
        // ground information very different from previous scan
        if (abs(mOme - mpPrev->mOme) + abs(mPhi - mpPrev->mPhi) > deg2rad(20.0) || abs(mGroundH - mpPrev->mGroundH) > 0.5)
        {
            R_lu_lup = mpPrev->R_lu_lup;
            mGroundH = mpPrev->mGroundH;
#ifdef EXPORT_LOG
            fLog << "\tFind inconsistency in ground info, use ground info from previous scan " << endl;
            f_odometry_debug << id << "\tFind inconsistency in ground info, used ground info from previous scan " << endl;

#endif
        }
    }
    return true;
}

/************************************
Extract ground using segments only, with the assumption of valid ground plane
Input:
    ground_plane_model_flag: if the ground is modeled as a plane
    R_lu_lup: leveling rotation
    mGroundH: ground height (only used when ground_plane_model_flag)
Output:
    ground segments.
    if ground_plane_model_flag, follwing will be refined
        mGroundH: ground height
        R_lu_lup: leveling rotation
return false:
    number of points from ground segments is too few
    if ground_plane_model_flag: ground info is incompatible with the previous scan
************************************/
/**
 * @brief Extracts ground segments from LiDAR scan using relative positioning and leveling information.
 * 
 * @param ground_plane_model_flag If true, uses an existing ground plane model for consistency checks and plane refinement.
 * @return true If valid ground segments are found and optionally refined into a consistent ground plane.
 * @return false If too few ground points are extracted or ground info is inconsistent with the previous frame.
 */

bool LidarScan::extract_ground_segments(bool ground_plane_model_flag)
{

    fLog<< "Extract all ground segments:" << ground_plane_model_flag  <<endl;
    // leveled point cloud
    Eigen::Vector3f surf_normal(0.0, 0.0, 1.0);
    vector<Eigen::Vector3f> centers;
    vector<bool> validCenter;
    vector<pair<int, int>> validSegments;
    int totalSegments = 0, segment_close_to_ground = 0;
    int numGroundSeg = 0;

    //Step 1, add ground segments based on relative position
    ClusterSegment cluster(1);
    for (int nChannel = mNScan - 1; nChannel >= 0; nChannel--)
    {
        for (std::size_t nSeg = 0; nSeg < mvSegments[nChannel].size(); nSeg++)
        {
            LidarSegment *thisSegment = &(mvSegments[nChannel][nSeg]);
            Eigen::Vector3f center_cur = thisSegment->centerPointTrans;

            totalSegments++;

            //1. if the plane model is valid, remove segments that are not close to the plane
            if (ground_plane_model_flag)
            {
                // with the assumption of valid surface, height extracted from previsou step can be used
                float dz = center_cur(2) - mGroundH;
                float angle = atan2(center_cur(2) - mGroundH, sqrt(center_cur(1) * center_cur(1) + center_cur(0) + center_cur(0)));
                if (abs(dz) > 0.2 && abs(angle) > deg2rad(10.0))
                    continue;
            }
            segment_close_to_ground++;


            //2. go through the segments that have been classified as ground, check the relative postion
            bool valid = true;
            for (std::size_t i_valid = 0; i_valid < centers.size(); i_valid++)
            {
                // if this reference ground segment has been removed, continue
                if (!validCenter[i_valid])
                    continue;

                Eigen::Vector3f dif = center_cur - centers[i_valid];
                float dxy_sqr = dif(0) * dif(0) + dif(1) * dif(1);

                //angular threshold in degree, if valid plane, the threshold is smaller
                float angular_threshold = ground_plane_model_flag ? 30.0 : 45.0;

                // Case 1. current segment is lower than the reference ground segment 
                if (dif(2) < 0.0)
                {
                    // if the current segment is much lower, remove the reference ground segment, continue;
                    dif = -dif;
                    if (atan(dif(2) / sqrt(dxy_sqr)) > deg2rad(angular_threshold)) 
                    {
                        validCenter[i_valid] = false;
                    }
                    continue;
                }
                
                // Case 2. current segment is higher than the reference ground segment 
                //else
                {
                    // if the current segment not much higher than the reference ground segment, continue;
                    if (atan(dif(2) / sqrt(dxy_sqr)) < deg2rad(angular_threshold))
                    {
                        continue;
                    }
                    //else, current segment is not ground segment, break
                    else
                    {
                        valid = false;
                        break;
                    }
                }
            }
            
            //3. if current segment is valid, include it in the ground segments
            if (valid)
            {
                centers.push_back(center_cur);
                validCenter.push_back(true);
                validSegments.push_back(make_pair(nChannel, nSeg));
            }
        }
    }

    //Step 2, remove ground segments based on the direction of segment - the ones that perpendicular segments
    int num_valid_seg_step1 = 0, num_valid_seg_step2 = 0;
    for (int nSeg = 0; nSeg < validSegments.size(); nSeg++)
    {
        if (!validCenter[nSeg])
            continue;
        
        num_valid_seg_step1++;

        int channelId = validSegments[nSeg].first;
        int segId = validSegments[nSeg].second;
        LidarSegment *thisSegment = &(mvSegments[channelId][segId]);
        Eigen::Vector3f center_cur = thisSegment->centerPointTrans;

        Eigen::Vector3f dire_trans = R_lu_lup * thisSegment->direction;
        dire_trans = dire_trans(2) > 0 ? dire_trans : dire_trans * -1.0;
        float angle = computeAngle(dire_trans, surf_normal);

        //angular threshold in degree, between the direction of a segment and horizontal axis
        //if valid plane, the threshold is smaller
        float angular_threshold = ground_plane_model_flag ? 30.0 : 45.0;

        if ((angle < deg2rad(90.0 - angular_threshold ))) // computed angle is the direction to the vertical direction 
            validCenter[nSeg] = false;
        else
            num_valid_seg_step2++;
    }
    fLog << "\tValid segments: all, close to ground, ground segments based on relative position, after removal based on segment direction: " << totalSegments << "/" << segment_close_to_ground << "/"
          << num_valid_seg_step1 << "/" << num_valid_seg_step2 << endl;
    //Step 3, if ground_plane_model_flag, we will use the ground segments to derive plane information again
    if(ground_plane_model_flag)
    {
        // another round of plane fitting. iterative. Only use the points with small smootheness vecPlanePoints
        int iterPlane = 0, iterPlaneMax = 2;
        Eigen::Vector4f paramsFinal;
        Eigen::Vector3f centerFinal;
        while (iterPlane < iterPlaneMax)
        {
            numGroundSeg = 0;
            std::vector<Eigen::Vector3f> allPointsClusters;
            allPointsClusters.clear();
            int num_all_points = 0;

            for (int nSeg = 0; nSeg < validSegments.size(); nSeg++)
            {
                if (!validCenter[nSeg])
                    continue;
                numGroundSeg++;
                int channelId = validSegments[nSeg].first;
                int segId = validSegments[nSeg].second;
                allPointsClusters.insert(allPointsClusters.end(), mvSegments[channelId][segId].vecPlanePoints.begin(), mvSegments[channelId][segId].vecPlanePoints.end());

                num_all_points += mvSegments[channelId][segId].length;
            }

            // CONTROL: if number of points too few for ground
            if (allPointsClusters.size() < 500)
            {

#ifdef EXPORT_LOG
                f_odometry_debug << id << "\tInsufficient ground points- segments (ground/all): " << numGroundSeg << "/" << totalSegments << "  points in ground segment (smooth/all): " << allPointsClusters.size() << " " << num_all_points << endl;
                fLog << "\tInsufficient ground points- segments (ground/all): " << numGroundSeg << "/" << totalSegments << "  points in ground segment (smooth/all): " << allPointsClusters.size() << " " << num_all_points << endl;
#endif
                return false;
            }

            planeFitting(allPointsClusters, paramsFinal, centerFinal);
            iterPlane++;

#ifdef EXPORT_LOG
            fLog << "\tGround segments - iter " << iterPlane << ": " << numGroundSeg << "/" << totalSegments << "  points: " << allPointsClusters.size() << " " << num_all_points << endl;
#endif

            cout << "Ground segments - " << iterPlane << ": " << numGroundSeg << "/" << totalSegments << "  points: " << allPointsClusters.size() << " " << num_all_points << endl;
            for (int nSeg = 0; nSeg < validSegments.size(); nSeg++)
            {
                if (!validCenter[nSeg])
                    continue;

                int channelId = validSegments[nSeg].first;
                int segId = validSegments[nSeg].second;
                LidarSegment *thisSegment = &(mvSegments[channelId][segId]);
                Eigen::Vector3f center_cur = thisSegment->centerPoint;
                Eigen::Vector4f center_cur_;
                center_cur_ << center_cur, 1.0;
                float d = center_cur_.transpose() * paramsFinal;
                if (abs(d) > 0.5)
                    validCenter[nSeg] = false;
            }
        }

        //Update the ground information
        Eigen::Vector3f normVecFinal(paramsFinal(0), paramsFinal(1), paramsFinal(2));
        float phi2 = asin(normVecFinal(0));
        float ome2 = atan2(-normVecFinal(1), normVecFinal(2));
        Eigen::Matrix3f mRlevel = Compute_Rotation(ome2, phi2, float(0.0));
        Eigen::Vector3f p(0.0, 0.0, -paramsFinal(3) / paramsFinal(2));
        mRlevel.transposeInPlace();
        Eigen::Vector3f p_trans = mRlevel * p;
        mGroundH = p_trans(2);
        R_lu_lup = mRlevel;
        mOme = ome2;
        mPhi = phi2;

        #ifdef EXPORT_LOG
            fLog << "\tGround info from all segments:\t" << p_trans(2) << "\t" << rad2deg(ome2) << "\t" << rad2deg(phi2) << endl;
        #endif


        // update level centers
        computeLevelCenters();
    }

    //Step 4, add valid segment to the cluster
    numGroundSeg = 0;
    int numPoints = 0;
    for (int nSeg = 0; nSeg < validSegments.size(); nSeg++)
    {
        if (!validCenter[nSeg])
            continue;
        numGroundSeg++;

        int channelId = validSegments[nSeg].first;
        int segId = validSegments[nSeg].second;

        LidarSegment *thisSegment = &(mvSegments[channelId][segId]);
        Eigen::Vector3f center_cur = thisSegment->centerPoint;

        cluster.addSegment(make_pair(channelId, segId), center_cur);
        numPoints += thisSegment->length;
    }
    cout << "Ground segments: " << numGroundSeg << "/" << totalSegments << "  points: " << numPoints << endl;
    fLog << "\tFinal ground segments: " << numGroundSeg << "/" << totalSegments << "  points: " << numPoints << endl;

    ground_cluster_.clear();
    ground_cluster_.push_back(cluster);

    // save the raw point, assign value to mvPoints
    for (size_t nC = 0; nC < ground_cluster_.size(); nC++)
    {
        for (size_t nSeg = 0; nSeg < ground_cluster_[nC].len; nSeg++)
        {
            int channelID = ground_cluster_[nC].vSegments[nSeg].first;
            int segID = ground_cluster_[nC].vSegments[nSeg].second;
            mvSegments[channelID][segID].final_feature_id = nC;
            mvSegments[channelID][segID].final_feature_type = 1;

            for (std::size_t nP = 0; nP < mvSegments[channelID][segID].length; nP++)
            {
                int index = mvSegments[channelID][segID].vecIndex[nP];
                ground_cluster_[nC].mvRawPoints.push_back(mvPoints[index].point);
                mvPoints[index].clusterID = nC;
                mvPoints[index].clusterType = 1;
            }
        }
    }


    // CONTROL: if number of points too few for ground
    if (numPoints < 500)
    {
#ifdef EXPORT_LOG
        f_odometry_debug << id << "\tInsufficient number of points in ground segments " << numPoints << endl;
        fLog << "Error: Insufficient number of points in ground segments " << numPoints << endl;
#endif
        return false;
    }

    // CONTROL: if the ground info compared to previous is very different
    if (ground_plane_model_flag)
    {
        if (mpPrev && mpPrev->valid_ground_info_flag_)
        {
            if (abs(mOme - mpPrev->mOme) + abs(mPhi - mpPrev->mPhi) > deg2rad(20.0) || abs(mGroundH - mpPrev->mGroundH) > 0.5) // consider as fail
            {
#ifdef EXPORT_LOG
                f_odometry_debug << id << "\tIncompatible final ground infor - differences in ome/phi " << rad2deg(abs(mOme - mpPrev->mOme) + abs(mPhi - mpPrev->mPhi))
                                 << " or height: " << abs(mGroundH - mpPrev->mGroundH) << endl;
                fLog << "\tIncompatible final ground infor - differences in ome/phi " << rad2deg(abs(mOme - mpPrev->mOme) + abs(mPhi - mpPrev->mPhi))
                     << " or height: " << abs(mGroundH - mpPrev->mGroundH) << endl;
#endif
                return false;
            }
        }
    }
    return true;
}



/************************************
Extract ground using segments only, while valid ground plane is not available
Output:
    ground segments.
************************************/
/*void LidarScan::extractGround()
{
    fLog<< "Extract all ground segments:" <<endl;
    // leveled point cloud
    Eigen::Vector3f surf_normal(0.0, 0.0, 1.0);
    vector<Eigen::Vector3f> centers;
    vector<bool> validCenter;
    vector<pair<int, int>> validSegments;
    int totalSegments = 0;
    int numGroundSeg = 0;
    ClusterSegment cluster(1);

    for (int nChannel = mNScan - 1; nChannel >= 0; nChannel--)
    {
        for (std::size_t nSeg = 0; nSeg < mvSegments[nChannel].size(); nSeg++)
        {
            LidarSegment *thisSegment = &(mvSegments[nChannel][nSeg]);
            Eigen::Vector3f center_cur = thisSegment->centerPointTrans;

            totalSegments++;

            bool valid = true;
            for (std::size_t i_valid = 0; i_valid < centers.size(); i_valid++)
            {
                if (!validCenter[i_valid])
                    continue;

                Eigen::Vector3f dif = center_cur - centers[i_valid];
                float dxy_sqr = dif(0) * dif(0) + dif(1) * dif(1);

                float threshold;
                if (dxy_sqr > 2.0 * 2.0)
                    threshold = 45.0;
                else
                    threshold = 45.0;

                // current seg is lower
                if (dif(2) < 0.0)
                {
                    dif = -dif;
                    if (atan(dif(2) / sqrt(dxy_sqr)) > deg2rad(threshold)) // current center not valid
                        validCenter[i_valid] = false;

                    continue;
                }

                if (atan(dif(2) / sqrt(dxy_sqr)) < deg2rad(threshold))
                    continue;
                else
                {
                    valid = false;
                    break;
                }
            }

            if (valid)
            {
                centers.push_back(center_cur);
                validCenter.push_back(true);
                validSegments.push_back(make_pair(nChannel, nSeg));
            }
        }
    }

    for (int nSeg = 0; nSeg < validSegments.size(); nSeg++)
    {
        if (!validCenter[nSeg])
            continue;

        int channelId = validSegments[nSeg].first;
        int segId = validSegments[nSeg].second;
        LidarSegment *thisSegment = &(mvSegments[channelId][segId]);
        Eigen::Vector3f center_cur = thisSegment->centerPointTrans;

        Eigen::Vector3f dire_trans = R_lu_lup * thisSegment->direction;
        dire_trans = dire_trans(2) > 0 ? dire_trans : dire_trans * -1.0;
        float angle = computeAngle(dire_trans, surf_normal);

        if (!(angle < deg2rad(45.0))) // 60 if the segment is vertical
        {
            numGroundSeg++;

            // thisSegment->type = 1;

            cluster.addSegment(make_pair(channelId, segId), center_cur);
            // for (std::size_t nP = 0; nP < thisSegment->length; nP++)
            // {
            //    int index = thisSegment->vecIndex[nP];
            //                    mvClusterTree[nC].mvRawPoints.push_back(mvPoints[index].point);
            //    mvPoints[index].clusterID = 0;
            // mvPoints[index].clusterType = 1;
            //   }
        }
    }

    ground_cluster_.push_back(cluster);

    // save the raw point, assign value to mvPoints
    for (size_t nC = 0; nC < ground_cluster_.size(); nC++)
    {
        for (size_t nSeg = 0; nSeg < ground_cluster_[nC].len; nSeg++)
        {
            int channelID = ground_cluster_[nC].vSegments[nSeg].first;
            int segID = ground_cluster_[nC].vSegments[nSeg].second;
            mvSegments[channelID][segID].final_feature_id = nC;
            mvSegments[channelID][segID].final_feature_type = 1;

            for (std::size_t nP = 0; nP < mvSegments[channelID][segID].length; nP++)
            {
                int index = mvSegments[channelID][segID].vecIndex[nP];
                ground_cluster_[nC].mvRawPoints.push_back(mvPoints[index].point);
                mvPoints[index].clusterID = nC;
                mvPoints[index].clusterType = 1;
            }
        }
    }
    cout << "Ground segments: " << numGroundSeg << "/" << totalSegments << endl;
}
*/



/************************************
Extract planar based on ground segments
Output:
    mpSurfPointsFlatGround: ground planar points
Return:
    true: number of mpSurfPointsFlatGround > 50
    else, false
mbTrackFeature & mbExtractedFeature: fail when false
************************************/
/**
 * @brief Extracts planar points from previously identified ground segments.
 * 
 * This function computes smoothness for each point and selects low-curvature points
 * as planar features. It outputs the leveled point cloud for ground surfaces.
 * 
 * @return true If at least 50 planar points are extracted successfully.
 * @return false Otherwise.
 */

bool LidarScan::extract_planar_points_from_ground_segment()
{
    mpSurfPointsFlatGround.reset(new pcl::PointCloud<PointType>());
    pSurfPointsLessFlatGround_ds_.reset(new pcl::PointCloud<PointType>());
    pSurfPointsLessFlatGround_.reset(new pcl::PointCloud<PointType>());


    fLog << "Extract ground planar points from ground segments:"<<endl;
    // cout << "Compute Smooth " << endl;
    mpCurv = new float[mNPoints];

    vector<int> mvValidIndex; // index in the mvPoints

    int *cloudSortInd = new int[mNPoints];        // index for used points
    int *cloudNeighborPicked = new int[mNPoints]; // curvature array of all point.
    int *cloudLabel = new int[mNPoints];          // curvature array of all point.

    vector<int> startIndex;
    vector<int> endIndex;

    int count = 0;
    // for each channel, find valid point and save the start and end index
    for (int nChannel = 0; nChannel < mNScan; nChannel++)
    {
        startIndex.push_back(count + N_Smooth);
        for (int nFiring = 0; nFiring < mNFiring; nFiring++)
        {
            int index = nFiring * mNScan + nChannel;

            if (mvPoints[index].clusterID >= 0 && mvPoints[index].clusterType == 1) //@@@Tian
            {
                mvValidIndex.push_back(index);
                count++;
            }
        }
        endIndex.push_back(max(0, count - N_Smooth - 1));
    }
    cout << count << endl;

    // compute the curvature, only valid for points within the start and end index
    for (int nP = N_Smooth; nP < count - N_Smooth; nP++)
    {
        LidarPoint *pThisPoint;
        LidarPoint *pNextPoint1;
        LidarPoint *pNextPoint2;

        pThisPoint = &(mvPoints[mvValidIndex[nP]]);

        float diffX = -float(N_Smooth * 2) * pThisPoint->point.x,
              diffY = -float(N_Smooth * 2) * pThisPoint->point.y,
              diffZ = -float(N_Smooth * 2) * pThisPoint->point.z;
        for (int nN = 1; nN <= N_Smooth; nN++)
        {
            pNextPoint1 = &(mvPoints[mvValidIndex[nP - nN]]);
            pNextPoint2 = &(mvPoints[mvValidIndex[nP + nN]]);

            diffX = diffX + pNextPoint1->point.x + pNextPoint2->point.x;
            diffY = diffY + pNextPoint1->point.y + pNextPoint2->point.y;
            diffZ = diffZ + pNextPoint1->point.z + pNextPoint2->point.z;
        }

        mpCurv[nP] = diffX * diffX + diffY * diffY + diffZ * diffZ;
        cloudNeighborPicked[nP] = 0;
        cloudSortInd[nP] = nP;
        pThisPoint->smooth2 = diffX * diffX + diffY * diffY + diffZ * diffZ;
    }

    //double ds_size_max = 1.0, ds_size_min = 0.2;
    for (int nChannel = 0; nChannel < mNScan; nChannel++)
    {
        // if number of points are too few
        int numChannelPoint = endIndex[nChannel] - startIndex[nChannel] + 1;
        if (numChannelPoint < 40)
            continue;
        pcl::PointCloud<PointType>::Ptr surfPointsLessFlatScan(new pcl::PointCloud<PointType>);

        int numSubArea = max(std::min(numChannelPoint / MinSubPoints, MaxSubArea), 1);

        for (int nArea = 0; nArea < numSubArea; nArea++)
        {
            int sp = startIndex[nChannel] + numChannelPoint * nArea / numSubArea;
            int ep = startIndex[nChannel] + numChannelPoint * (nArea + 1) / numSubArea - 1;

            // sort index from minimum to maximum
            std::sort(cloudSortInd + sp, cloudSortInd + ep + 1, comp);

            // planar points
            int smallestPickedNum = 0;
            for (int k = sp; k <= ep; k++)
            {
                int ind = cloudSortInd[k];
                if (cloudNeighborPicked[ind] == 0 && mpCurv[ind] < Smooth_Threshold)
                {
                    cloudLabel[ind] = -1;
                    mvPoints[mvValidIndex[ind]].point_type_id = 5;

                    PointType pointLevel;
                    TransformToLevel(&(mvPoints[mvValidIndex[ind]].point), &pointLevel);
                    mpSurfPointsFlatGround->push_back(pointLevel);

                    smallestPickedNum++;
                    if (smallestPickedNum >= 4)
                        break;

                    cloudNeighborPicked[ind] = 1;
                    for (int l = 1; l <= N_Smooth; l++)
                    {
                        float diffX = mvPoints[mvValidIndex[ind + l]].point.x - mvPoints[mvValidIndex[ind + l - 1]].point.x;
                        float diffY = mvPoints[mvValidIndex[ind + l]].point.y - mvPoints[mvValidIndex[ind + l - 1]].point.y;
                        float diffZ = mvPoints[mvValidIndex[ind + l]].point.z - mvPoints[mvValidIndex[ind + l - 1]].point.z;
                        if (diffX * diffX + diffY * diffY + diffZ * diffZ > 0.05)
                            break;

                        cloudNeighborPicked[ind + l] = 1;
                    }
                    for (int l = -1; l >= -N_Smooth; l--)
                    {
                        float diffX = mvPoints[mvValidIndex[ind + l]].point.x - mvPoints[mvValidIndex[ind + l + 1]].point.x;
                        float diffY = mvPoints[mvValidIndex[ind + l]].point.y - mvPoints[mvValidIndex[ind + l + 1]].point.y;
                        float diffZ = mvPoints[mvValidIndex[ind + l]].point.z - mvPoints[mvValidIndex[ind + l + 1]].point.z;
                        if (diffX * diffX + diffY * diffY + diffZ * diffZ > 0.05)
                            break;

                        cloudNeighborPicked[ind + l] = 1;
                    }
                }
            }
            // cout << "planar" << endl;

            // flat points for map
            for (int k = sp; k <= ep; k++)
            {
                int ind = cloudSortInd[k];
                if (mpCurv[ind] < Smooth_Threshold)
                {
                    // PointType pointLevel;
                    // TransformToLevel(&(mvPoints[mvValidIndex[ind]].point), &pointLevel);
                    // surfPointsLessFlatScan->push_back(pointLevel);
                    surfPointsLessFlatScan->push_back(mvPoints[mvValidIndex[ind]].point);
                    if (mvPoints[mvValidIndex[ind]].point_type_id != 5)
                        mvPoints[mvValidIndex[ind]].point_type_id = 6;
                }
                else
                    break;
            }
        }

        //double ds_size = ds_size_max - (ds_size_max - ds_size_min)/double(mNScan-1) * double(nChannel);

        pcl::PointCloud<PointType> surfPointsLessFlatScanDS;

        //test. distance based downsample
        double ds_size = mPara.ground_downsample_distance_per_scan;
        if (0) // voxel will change the intensity value (ratio)
        {
            pcl::VoxelGrid<PointType> downSizeFilter;
            downSizeFilter.setInputCloud(surfPointsLessFlatScan);
            downSizeFilter.setLeafSize(ds_size, ds_size, ds_size);
            downSizeFilter.filter(surfPointsLessFlatScanDS);
        }
        else //intensity wont change
        {
            downsamplePointCloudDistance(surfPointsLessFlatScan, ds_size, surfPointsLessFlatScanDS);
        }

        *pSurfPointsLessFlatGround_ds_ += surfPointsLessFlatScanDS;
        *pSurfPointsLessFlatGround_ += *(surfPointsLessFlatScan);
    }

    delete []cloudSortInd;        
    delete []cloudNeighborPicked; 
    delete []cloudLabel ;         
    delete []mpCurv;

    cout << "Valid point for smooth: " << count << "/" << mNPoints << endl;
    cout << "Number of plane points from ground cluster only: " << mpSurfPointsFlatGround->size() << "/" << pSurfPointsLessFlatGround_ds_->size() << endl;

#ifdef EXPORT_LOG
    fLog << "\tValid point for smooth: " << count << "/" << mNPoints << endl;
    fLog << "\tNumber of ground planar points (fixed number in each subregion): " << mpSurfPointsFlatGround->size() << endl;
    fLog << "\tNumber of ground planar points below threshold, before/after ds: " << pSurfPointsLessFlatGround_->size() << "\t"
         << pSurfPointsLessFlatGround_ds_->size() << endl;
#endif


    //CONTROL:if number of reliable ground planar points is too few
    if (mpSurfPointsFlatGround->size() < 50)
    {

#ifdef EXPORT_LOG
        f_odometry_debug << id << "\tNumber of planar points: " << mpSurfPointsFlatGround->size() << endl;

#endif
        return false;
    }
    else
        return true;
}

/*for each segment, compute centerPointTrans based on computed R_lu_lup*/
/**
 * @brief Computes the leveled center point for each segment using the R_lu_lup rotation matrix.
 * 
 * This is typically used after ground leveling to transform segment centers into a normalized frame.
 */

void LidarScan::computeLevelCenters()
{
    for (int nChannel = 0; nChannel < mNScan; nChannel++)
    {
        for (std::size_t nSeg = 0; nSeg < mvSegments[nChannel].size(); nSeg++)
        {
            LidarSegment *thisSegment = &(mvSegments[nChannel][nSeg]);
            thisSegment->centerPointTrans = R_lu_lup * thisSegment->centerPoint;
        }
    }
}

/************************************
Extract tree cluster using non-ground segments
Input:
    mGroundH: ground height (optional!)
    R_lu_lup: leveling rotation
Output:
    mvTreeClusters: tree clusters
Return:
    true: >= MinTreePair
    else, false
mbTrackFeature & mbExtractedFeature: fail when number of trees < MinTreePair
************************************/
/**
 * @brief Extracts tree clusters from non-ground LiDAR segments.
 * 
 * Clusters are formed by aggregating vertically aligned segments, applying geometric and directional filters,
 * and finally merging clusters that are spatially and directionally similar.
 * 
 * @return true If the number of valid tree clusters is greater than or equal to MinTreePair.
 * @return false Otherwise (e.g., too few trees).
 */

bool LidarScan::extractTreeCluster()
{

    fLog << "Extract tree clusters:" <<endl;
    for (int nChannel = mNScan - 1; nChannel >= 0; nChannel--)
    {
        for (std::size_t nSeg = 0; nSeg < mvSegments[nChannel].size(); nSeg++)
        {
            // cout << nChannel << "\t" << nSeg << "\t" << mvSegments[nChannel][nSeg].type << endl;
            LidarSegment *thisSegment = &(mvSegments[nChannel][nSeg]);

            float hCur = thisSegment->centerPointTrans(2);

            // segment used in defining the ground
            if (thisSegment->final_feature_type == 1)
                continue;

            // // height limitation needed criteria @@
            // if (!mbLevel)
            // {
            //     if (hCur <= (mGroundH + mGroundBuffer))
            //         continue;
            // }

            int bestTC = -1;
            float minScore = 100.0;
            float bestDistance = 0.0;
            // go through all available clusters, find the best compare to the normal vector defined by ground
            for (std::size_t nTC = 0; nTC < mvClusterTree.size(); nTC++)
            {
                // info of the last segment of a cluster
                float hPrev = mvClusterTree[nTC].segCenters.back()(2);
                int channelPrev = mvClusterTree[nTC].vSegments.back().first;

                float dist = (thisSegment->centerPointTrans - mvClusterTree[nTC].segCenters.back()).norm();
                int dChannel = channelPrev - nChannel;
                if (dChannel == 0)
                    continue;
                if (dist > MaxTreeSegDistance)
                    continue;
                dist = dist / float(dChannel);
                if (dChannel < 1)
                {
                    cout << dChannel << endl;
                    throw std::runtime_error("Wrong dChannel");
                }

                if (hCur > hPrev && dChannel <= ChannelContinous)
                {
                    // check if the distance is very different from previous
                    if (mvClusterTree[nTC].len > 1 && (dist > mvClusterTree[nTC].uDistance * 1.5 || dist < mvClusterTree[nTC].uDistance / 1.5))
                        continue;

                    // compute the angle for add to available tree trunck
                    // //1. based on a seg
                    // int indexRef = max(0, mvClusterTree[nTC].len - 2);                                            // not choosing the closest segment, the 2nd closest
                    // Eigen::Vector3f v1 = thisSegment->centerPointTrans - mvClusterTree[nTC].segCenters[indexRef]; // test.back();

                    // 2. based on the mean seg
                    Eigen::Vector3f v1 = thisSegment->centerPointTrans - mvClusterTree[nTC].vCenterTrans;

                    // compute angles. Points are transformed!
                    float angle = computeAngle(v1, mvClusterTree[nTC].vNorm); //@surf_normal);
                    if (angle < TreeAngThreshold)
                    {
                        float score = angle * dist;
                        if (score < minScore)
                        {
                            minScore = score;
                            bestTC = nTC;
                            bestDistance = dist;
                        }
                    }
                }
            }

            // add to the best cluster
            if (bestTC != -1)
                mvClusterTree[bestTC].addSegment(make_pair(nChannel, nSeg), thisSegment->centerPointTrans, bestDistance);
            else // create a new cluster
            {
                // initialize a new cluster
                //Height cretiria valid only if the ground is flat and horizontal (with valid mGroundH)
                if (mbLevel || (!mbLevel && hCur < mGroundH + mGroundBufferTreeInit))
                {
                    ClusterSegment cluster(2);
                    cluster.addSegment(make_pair(nChannel, nSeg), thisSegment->centerPointTrans);
                    mvClusterTree.push_back(cluster);
                }
            }
        }
    }
    //save cluster id to segments
    for (size_t nC = 0; nC < mvClusterTree.size(); nC++)
    {
        for (size_t nSeg = 0; nSeg < mvClusterTree[nC].len; nSeg++)
        {
            int channelID = mvClusterTree[nC].vSegments[nSeg].first;
            int segID = mvClusterTree[nC].vSegments[nSeg].second;
            mvSegments[channelID][segID].tree_cluster_id = nC;
           // mvSegments[channelID][segID].tree_outlier_removal_flag = 10;

        }
    }
    cout << "Number of tree clusters " << mvClusterTree.size() << endl;
    fLog << "\tNumber of initial tree clusters " << mvClusterTree.size() << endl;

    // check whether a cluster is valid or not
    std::vector<ClusterSegment> vTempCluster;
    int num_fail_num_segment = 0, num_fail_fit_line = 0;
    for (size_t nC = 0; nC < mvClusterTree.size(); nC++)
    {
        // creteria 1: number of segment
        if (mvClusterTree[nC].len >= MinNumSegTree)
        {
            // creteria 2: fit the line
            if (mvClusterTree[nC].fitLine())
            {
                vTempCluster.push_back(mvClusterTree[nC]);
            }
            else
            {
                num_fail_fit_line++;
                for (size_t nSeg = 0; nSeg < mvClusterTree[nC].len; nSeg++)
                {
                    int channelID = mvClusterTree[nC].vSegments[nSeg].first;
                    int segID = mvClusterTree[nC].vSegments[nSeg].second;
                    mvSegments[channelID][segID].tree_outlier_removal_flag = 1;
                }
                
            }
        }
        else
        {
            num_fail_num_segment++;
            for (size_t nSeg = 0; nSeg < mvClusterTree[nC].len; nSeg++)
            {
                int channelID = mvClusterTree[nC].vSegments[nSeg].first;
                int segID = mvClusterTree[nC].vSegments[nSeg].second;
                mvSegments[channelID][segID].tree_outlier_removal_flag = 0;
            }
        }
    }
    mvClusterTree = vTempCluster;
    vTempCluster.clear();
    cout << "Number of tree based clusters with enough segments and survive from fit line " << mvClusterTree.size() << endl;
    fLog << "\tFailed clusters in number of segments/line fitting: " << num_fail_num_segment << "/" << num_fail_fit_line << ", number of remaining "
         << mvClusterTree.size() << endl;

    // merge tree clusters that belongs to same clusters iteratively
    while (1)
    {
        for (size_t nC = 0; nC < mvClusterTree.size(); nC++)
        {
            if (!mvClusterTree[nC].valid) // invalid in the previous
                continue;

            Eigen::Vector3f center1 = mvClusterTree[nC].vCenterTrans;
            int bestTC = -1;
            float minAngle = PI;

            for (size_t nC2 = nC + 1; nC2 < mvClusterTree.size(); nC2++)
            {
                Eigen::Vector3f center2 = mvClusterTree[nC2].vCenterTrans;
                float d_sqr = (center2(0) - center1(0)) * (center2(0) - center1(0)) +
                              (center2(1) - center1(1)) * (center2(1) - center1(1));

                if (d_sqr > MinTreeDistance * MinTreeDistance) // large enough
                {
                    continue;
                }
                else // two cluster close to each other
                {
                    float angle1 = computeAngle(mvClusterTree[nC].vNorm, mvClusterTree[nC2].vNorm);
                    float angle2 = computeAngle((center2 - center1), mvClusterTree[nC2].vNorm);
                    // angle2 = angle2 < PI ? angle2 : 2.0 * PI - angle2;
                    // cout << angle1 << " " << angle2 << endl;
                    if (angle1 < deg2rad(20.0) && angle2 < deg2rad(20.0)) // merge
                    {
                        if (angle2 < minAngle)
                        {
                            minAngle = angle2;
                            bestTC = nC2;
                        }
                    }
                }
            }
            // find a valid tree to merge
            if (bestTC >= 0)
            {
                ClusterSegment temp = mvClusterTree[nC];
                temp.mergeCluster(mvClusterTree[bestTC]);
                if (temp.fitLine())
                {
                    // cout << "Merge: " << bestTC << " to " << nC << endl;
                    mvClusterTree[nC] = temp;

                    // deactivate the tree being merged
                    mvClusterTree[bestTC].valid = false;
                }
                else
                {
                    // cout << "Merge but fail: " << MinTreeDistance << endl;
                }
            }

            if (mvClusterTree[nC].valid)
            {
                vTempCluster.push_back(mvClusterTree[nC]);
            }
        }
        int prev_tree_num = mvClusterTree.size();
        mvClusterTree = vTempCluster;
        vTempCluster.clear();
        if (prev_tree_num == mvClusterTree.size())
            break;
    }
    cout << "Number of valid tree clusters after merging " << mvClusterTree.size() << endl;
    fLog << "\tNumber of tree clusters after merging " << mvClusterTree.size() << endl;


    // creteria 3: height range of the segment should be large enough
    for (size_t nC = 0; nC < mvClusterTree.size(); nC++)
    {
        Eigen::Vector3f v1 = mvClusterTree[nC].segCenters[0];
        Eigen::Vector3f v2 = mvClusterTree[nC].segCenters.back();
        // creteria height threshold
        if (mvClusterTree[nC].compute_height_range() >= 1.0)
        {
            vTempCluster.push_back(mvClusterTree[nC]);
        }
        else
        {
            for (size_t nSeg = 0; nSeg < mvClusterTree[nC].len; nSeg++)
            {
                int channelID = mvClusterTree[nC].vSegments[nSeg].first;
                int segID = mvClusterTree[nC].vSegments[nSeg].second;
                mvSegments[channelID][segID].tree_outlier_removal_flag = 3;
            }
        }
    }
    mvClusterTree = vTempCluster;
    cout << "Number of tree based clusters with enough height difference " << mvClusterTree.size() << endl;
    fLog << "\tNumber of tree based clusters with enough height difference " << mvClusterTree.size() << endl;

    // save the raw point, assign cluster and cluster id value to mvPoints
    for (size_t nC = 0; nC < mvClusterTree.size(); nC++)
    {
        for (size_t nSeg = 0; nSeg < mvClusterTree[nC].len; nSeg++)
        {
            int channelID = mvClusterTree[nC].vSegments[nSeg].first;
            int segID = mvClusterTree[nC].vSegments[nSeg].second;
            mvSegments[channelID][segID].final_feature_id = nC;
            mvSegments[channelID][segID].final_feature_type = 2;
            mvSegments[channelID][segID].tree_outlier_removal_flag = 10;

            for (std::size_t nP = 0; nP < mvSegments[channelID][segID].length; nP++)
            {
                int index = mvSegments[channelID][segID].vecIndex[nP];
                mvClusterTree[nC].mvRawPoints.push_back(mvPoints[index].point);
                mvPoints[index].clusterID = nC;
                mvPoints[index].clusterType = 2;
            }
        }
    }

#ifdef EXPORT_LOG
    fLog << "\tExtracted trees: " << mvClusterTree.size() << endl;
#endif

    if (mvClusterTree.size() < MinTreePair)
    {
#ifdef EXPORT_LOG
        fLog << "Error: not enough trees" << endl;
        f_odometry_debug << id << "\tTree clusters: " << mvClusterTree.size() << endl;
#endif
        return false;
    }

    return true;
}

/*Compute incremental rotation and translation, and EOP in mapping frame */
/**
 * @brief Computes the transformation from local (LiDAR) frame to global map frame.
 * 
 * This function updates the pose `T_lu_m`, using either feature-based estimation or trajectory reference,
 * and handles initialization and propagation across multiple scans.
 */

void LidarScan::computeTransformation()
{
    if (mbInit) // first scan
    {
        if (lastR.empty())
        {
            if (valid_feature_based_odo_flag_)
            {
                r_lu_m = Eigen::Vector3d::Zero();
                R_lu_m = R_lu_lup.cast<double>(); // Eigen::Matrix3f::Identity();
            }
        }
        else
        {
            r_lu_m = lastr.back();
            R_lu_m = lastR.back();
        }

        T_lu_m.rotate(R_lu_m);
        T_lu_m.pretranslate(r_lu_m);

        if (mpTraj)
        {
            if (lastR.empty())
            {
                R_lu_m_ref = R_lu_m;
                r_lu_m_ref = r_lu_m;
            }
            else
            {

                Eigen::Vector3d r_b_m_t1, r_b_m_t2, r_lu_m_t1, r_lu_m_t2, r_lut2_lut1_rr;
                Eigen::Matrix3d R_b_m_t1, R_b_m_t2, R_lu_m_t1, R_lu_m_t2, R_lut2_lut1_rr;
                int idx;
                if (!(mpTraj->bopInterpolation(refT.back() / 1000.0, -1, idx, r_b_m_t1, R_b_m_t1) && mpTraj->bopInterpolation(mTimeEnd / 1000.0, -1, idx, r_b_m_t2, R_b_m_t2)))
                    throw std::runtime_error("Wrong Time");

                r_lu_m_t1 = r_b_m_t1 + R_b_m_t1 * r_lu_b_ref;
                R_lu_m_t1 = R_b_m_t1 * R_lu_b_ref;

                r_lu_m_t2 = r_b_m_t2 + R_b_m_t2 * r_lu_b_ref;
                R_lu_m_t2 = R_b_m_t2 * R_lu_b_ref;

                R_lut2_lut1_rr = R_lu_m_t1.inverse() * R_lu_m_t2;
                r_lut2_lut1_rr = R_lu_m_t1.inverse() * (r_lu_m_t2 - r_lu_m_t1);

                R_lu_m_ref = lastRefR.back() * R_lut2_lut1_rr;
                r_lu_m_ref = lastRefr.back() + lastRefR.back() * r_lut2_lut1_rr;
            }
        }
        cout << "init" << endl;
        //  system("read -p 'Press Enter to continue...' var");
    }
    else
    {

        R_lu_m = mpPrev->R_lu_m * R_lut2_lut1;
        r_lu_m = mpPrev->r_lu_m + mpPrev->R_lu_m * r_lut2_lut1;

        T_lu_m.rotate(R_lu_m);
        T_lu_m.pretranslate(r_lu_m);
    }

    lastR.clear();
    lastR.push_back(R_lu_m);
    lastr.clear();
    lastr.push_back(r_lu_m);
    if (mbTraj)
    {
        lastRefR.clear();
        lastRefR.push_back(R_lu_m_ref);
        lastRefr.clear();
        lastRefr.push_back(r_lu_m_ref);
        refT.clear();
        refT.push_back(mTimeEnd);
    }
    Eigen::Vector3d euler_angles = Find_Rotation(R_lu_m);
    mEop.XO = r_lu_m(0);
    mEop.YO = r_lu_m(1);
    mEop.ZO = r_lu_m(2);
    mEop.omega = euler_angles(0);
    mEop.phi = euler_angles(1);
    mEop.kappa = euler_angles(2);

#ifdef EXPORT_TRAJECTORY
    if (mbTraj)
    {
        Eigen::Vector3d angles = Find_Rotation(R_lu_m_ref);
        fTrajectoryRef << id << "\t" << r_lu_m_ref(0) << "\t" << r_lu_m_ref(1) << "\t" << r_lu_m_ref(2) << "\t"
                       << rad2deg(angles(0)) << "\t" << rad2deg(angles(1)) << "\t" << rad2deg(angles(2)) << "\t" << mPartID << endl;
    }
    int type = 2; // point based
    if (mbInit)
    {
        type = 0;
        // cout << "type: " << type << endl;
        //    system("read -p 'Press Enter to continue...' var");
    }
    else
    {
        if (valid_feature_based_odo_flag_)
            type = 1; // feature based
    }
    fMapping << id << "\t" << r_lu_m(0) << "\t" << r_lu_m(1) << "\t" << r_lu_m(2) << "\t"
             << rad2deg(euler_angles(0)) << "\t" << rad2deg(euler_angles(1)) << "\t" << rad2deg(euler_angles(2)) << "\t" << mvClusterTree.size()
             << "\t" << type << endl;
#endif
}

/************************************
Match tree clusters between current and previous trees based on 3D info
Requirement: mbTrackFeature & !mbInit  & mpPrev->mbExtractedFeature
Input:
    bSeach: true, search different kappa value, false, don't use
Output:
    mvTreeClusters: tree clusters
Return:
    true: number of mpSurfPointsFlatGround > 50
    else, false
mbTrackFeature: fail when number of pairing trees < MinTreePair
************************************/
/**
 * @brief Matches tree clusters between current and previous LiDAR scans using 3D alignment (center and orientation).
 * 
 * @param bSearch Whether to iterate over a range of kappa (yaw) angles for best alignment.
 * @return true If enough matching tree pairs (≥ MinTreePair) are found.
 * @return false Otherwise.
 */

bool LidarScan::matchTrees3d(bool bSearch)
{
    fLog << "Match trees between scans in 3d" <<endl;
    int numPrev = mpPrev->mvClusterTree.size();
    int numCur = mvClusterTree.size();

    Eigen::Vector3d t = t_ini_level;
    Eigen::Vector3d angles = Find_Rotation(q_ini_level.toRotationMatrix());

#ifdef EXPORT_LOG
    fLog << "\tInitial transformation values " << t.transpose() << "\t" << rad2deg(angles).transpose() << endl;
#endif

    // best kap value for tree matching
    int maxPairs = -1;
    int bestKap = -100;
    Eigen::Quaterniond q_best;

    vector<double> vecDistanceThreshold;
    for (int nCur = 0; nCur < numCur; nCur++) // current scan
        vecDistanceThreshold.push_back(max(1.0, 0.05 * mvClusterTree[nCur].vCenterTrans.norm()));

    vector<Eigen::Vector3d> vecPrevCenter;
    vector<Eigen::Vector3d> vecPrevNorm;
    for (int nPrev = 0; nPrev < numPrev; nPrev++) // prev scan
    {
        Eigen::Vector3d norm_prev = mpPrev->mvClusterTree[nPrev].vNorm.cast<double>();
        Eigen::Vector3d center_prev = mpPrev->mvClusterTree[nPrev].vCenterTrans.cast<double>();
        vecPrevCenter.push_back(center_prev);
        vecPrevNorm.push_back(norm_prev);
    }

    int scan_dif = bSearch ? id - mpPrev->id : 0;
    for (int i = -20 * scan_dif; i <= 20 * scan_dif; i++)
    {
        Eigen::Quaterniond q;
        if (bSearch)
        {
            double kap = deg2rad(double(i));
            Eigen::Matrix3d candRot = Compute_Rotation(angles(0), angles(1), kap);
            q = candRot;
        }
        else
            q = q_ini_level;
        int numPairs;
        std::vector<pair<int, int>> pointPairs;

        // compute estimated center and norm
        vector<Eigen::Vector3d> vecEstimateCenter;
        vector<Eigen::Vector3d> vecEstimateNorm;

        for (int nCur = 0; nCur < numCur; nCur++) // current scan
        {
            Eigen::Vector3d norm_cur = mvClusterTree[nCur].vNorm.cast<double>();
            Eigen::Vector3d center_cur = mvClusterTree[nCur].vCenterTrans.cast<double>();
            vecEstimateCenter.push_back(q * center_cur + t);
            vecEstimateNorm.push_back(q * norm_cur);
        }

        // matching trees
        std::vector<int> pair1;                   // cur -> prev
        for (int nCur = 0; nCur < numCur; nCur++) // current scan
        {
            double score_best = 100.0;
            int id_best = -1;
            for (int nPrev = 0; nPrev < numPrev; nPrev++) // prev scan
            {
                TicToc tTemp;

                // Eigen::Vector3d norm_prev = mpPrev->mvClusterTree[nPrev].vNorm.cast<double>();
                //  Eigen::Vector3d center_prev = mpPrev->mvClusterTree[nPrev].vCenterTrans.cast<double>();

                // Eigen::Vector3d norm_cur = mvClusterTree[nCur].vNorm.cast<double>();
                // Eigen::Vector3d center_cur = mvClusterTree[nCur].vCenterTrans.cast<double>();

                double d = computePoint2lineDistance(vecEstimateCenter[nCur], vecPrevCenter[nPrev], vecPrevCenter[nPrev] + vecPrevNorm[nPrev]);
                double angle = computeAngle(vecPrevNorm[nPrev], vecEstimateNorm[nCur]);
                double threshold = vecDistanceThreshold[nCur];
                Tcheck += tTemp.toc();

                if (d < threshold && angle < deg2rad(20.0))
                {
                    double score = d * angle;
                    if (score < score_best)
                    {
                        score_best = score;
                        id_best = nPrev;
                    }
                }
            }
            pair1.push_back(id_best);
        }

        std::vector<int> pair2;                       // prev -> cur
        for (int nPrev = 0; nPrev < numPrev; nPrev++) // prev scan
        {
            double score_best = 100.0;
            int id_best = -1;
            // Eigen::Vector3d norm_prev = mpPrev->mvClusterTree[nPrev].vNorm.cast<double>();
            // Eigen::Vector3d center_prev = mpPrev->mvClusterTree[nPrev].vCenterTrans.cast<double>();

            for (int nCur = 0; nCur < numCur; nCur++) // current scan
            {
                // Eigen::Vector3d norm_cur = mvClusterTree[nCur].vNorm.cast<double>();
                // Eigen::Vector3d center_cur = mvClusterTree[nCur].vCenterTrans.cast<double>();
                TicToc tTemp;

                double d = computePoint2lineDistance(vecEstimateCenter[nCur], vecPrevCenter[nPrev], vecPrevCenter[nPrev] + vecPrevNorm[nPrev]);
                double angle = computeAngle(vecPrevNorm[nPrev], vecEstimateNorm[nCur]);
                double threshold = vecDistanceThreshold[nCur];
                Tcheck += tTemp.toc();

                if (d < threshold && angle < deg2rad(20.0))
                {
                    double score = d * angle;
                    if (score < score_best)
                    {
                        score_best = score;
                        id_best = nCur;
                    }
                }
            }
            pair2.push_back(id_best);
        }

        // checkpair
        for (int nCur = 0; nCur < numCur; nCur++) // go through current scan
        {
            bool bMatch = false;
            int PrevId = pair1[nCur];
            if (PrevId != -1)
            {
                if (pair2[PrevId] == nCur)
                {
                    pointPairs.push_back(make_pair(PrevId, nCur));
                    bMatch = true;
                }
            }
        }

        numPairs = pointPairs.size();

        if (numPairs > maxPairs)
        {
            maxPairs = numPairs;
            bestKap = i;
            mVecTreePairs = pointPairs;
            q_best = q;
        }

    } // end of go through different kap

    // assign new rotation for the refinement
    q_ini_level = q_best;

    cout << "Pairs: " << maxPairs << " at kap = " << bestKap << endl;
#ifdef EXPORT_LOG
    fLog << "\tTree pairs: " << maxPairs << " at kap = " << bestKap << endl;
#endif
    if (maxPairs < MinTreePair) // fail
    {

#ifdef EXPORT_LOG
        f_odometry_debug << id << "\tNo enough valid tree matches: " << maxPairs << endl;
        fLog << "Error: No enough valid tree matches: " << maxPairs << endl;

#endif
        return false;
    }
    else
        return true;
}

/************************************
Match tree clusters between current and previous trees based on 2D locations only
Requirement: mbTrackFeature & !mbInit  & mpPrev->mbExtractedFeature
Input:
    bSeach: true, search different kappa value, false, don't use
Output:
    mvTreeClusters: tree clusters
Return:
    true: number of mpSurfPointsFlatGround > 50
    else, false
mbTrackFeature: fail when number of pairing trees < MinTreePair
************************************/
/**
 * @brief Matches tree clusters between current and previous LiDAR scans using 2D positional alignment (XY-plane).
 * 
 * @param bSearch Whether to iterate over kappa (yaw) values to find best alignment.
 * @return true If enough matching tree pairs are found.
 * @return false Otherwise.
 */

bool LidarScan::matchTrees2d(bool bSearch)
{    
    fLog << "Match trees between scans in 2d" <<endl;
    int numPrev = mpPrev->mvClusterTree.size();
    int numCur = mvClusterTree.size();

    Eigen::Vector3d t = t_ini_level;
    Eigen::Vector3d angles = Find_Rotation(q_ini_level.toRotationMatrix());

#ifdef EXPORT_LOG
    fLog << "\tInitial transformation values " << t.transpose() << "\t" << rad2deg(angles).transpose() << endl;
#endif

    // best kap value for tree matching
    int maxPairs = -1;
    int bestKap = -100;
    Eigen::Quaterniond q_best;

    int scan_dif = bSearch ? id - mpPrev->id : 0;

    for (int i = -20 * scan_dif; i <= 20 * scan_dif; i++)
    {
        Eigen::Quaterniond q;
        if (bSearch)
        {
            double kap = deg2rad(double(i));
            Eigen::Matrix3d candRot = Compute_Rotation(angles(0), angles(1), kap);
            q = candRot;
        }
        else
        {
            // q = q_ini_level;
            t = Eigen::Vector3d::Zero();
            q = Eigen::Quaterniond::Identity();
        }
        int numPairs;
        std::vector<pair<int, int>> pointPairs;

        // matching trees
        std::vector<int> pair1;                   // cur -> prev
        for (int nCur = 0; nCur < numCur; nCur++) // current scan
        {
            double score_best = 100.0;
            int id_best = -1;
            for (int nPrev = 0; nPrev < numPrev; nPrev++) // prev scan
            {
                Eigen::Vector3d center_prev = mpPrev->mvClusterTree[nPrev].vCenterTrans.cast<double>();
                Eigen::Vector3d center_cur = mvClusterTree[nCur].vCenterTrans.cast<double>();
                Eigen::Vector3d center_cur_predicted = q * center_cur + t;
                double d = (center_cur_predicted(0) - center_prev(0)) * (center_cur_predicted(0) - center_prev(0)) + (center_cur_predicted(1) - center_prev(1)) * (center_cur_predicted(1) - center_prev(1));
                double threshold = max(1.0, 0.05 * center_cur.norm());

                if (d < threshold)
                {
                    double score = d;
                    if (score < score_best)
                    {
                        score_best = score;
                        id_best = nPrev;
                    }
                }
            }
            pair1.push_back(id_best);
        }

        std::vector<int> pair2;                       // prev -> cur
        for (int nPrev = 0; nPrev < numPrev; nPrev++) // prev scan
        {
            double score_best = 100.0;
            int id_best = -1;
            for (int nCur = 0; nCur < numCur; nCur++) // current scan
            {
                Eigen::Vector3d center_prev = mpPrev->mvClusterTree[nPrev].vCenterTrans.cast<double>();
                Eigen::Vector3d center_cur = mvClusterTree[nCur].vCenterTrans.cast<double>();

                Eigen::Vector3d center_cur_predicted = q * center_cur + t;
                double d = (center_cur_predicted(0) - center_prev(0)) * (center_cur_predicted(0) - center_prev(0)) + (center_cur_predicted(1) - center_prev(1)) * (center_cur_predicted(1) - center_prev(1));
                double threshold = max(1.0, 0.05 * center_cur.norm());

                if (d < threshold)
                {
                    double score = d;
                    if (score < score_best)
                    {
                        score_best = score;
                        id_best = nCur;
                    }
                }
            }
            pair2.push_back(id_best);
        }

        // checkpair
        for (int nCur = 0; nCur < numCur; nCur++) // go through current scan
        {
            bool bMatch = false;
            int PrevId = pair1[nCur];
            if (PrevId != -1)
            {
                if (pair2[PrevId] == nCur)
                {
                    pointPairs.push_back(make_pair(PrevId, nCur));
                    bMatch = true;
                }
            }
        }

        numPairs = pointPairs.size();

        if (numPairs > maxPairs)
        {
            maxPairs = numPairs;
            bestKap = i;
            mVecTreePairs = pointPairs;
            q_best = q;
        }

    } // end of go through different kap

    // assign new rotation for the refinement
    q_ini_level = q_best;

    cout << "Pairs: " << maxPairs << " at kap = " << bestKap << endl;
#ifdef EXPORT_LOG
    fLog << "\tTree pairs: " << maxPairs << " at kap = " << bestKap << endl;
#endif
    if (maxPairs < MinTreePair) // fail
    {

#ifdef EXPORT_LOG
        f_odometry_debug << id << "\tNo enough valid tree matches: " << maxPairs << endl;
        fLog << "Error: No enough valid tree matches: " << maxPairs << endl;

#endif
        return false;
    }
    else
        return true;
}

/* Compute attribute of the lidar seg:
length, vecPoints, centerPoint, vecAngle, uAngle, maxfiring, minfiring*/
/**
 * @brief Computes geometric and angular attributes for a LiDAR segment.
 * 
 * This includes segment length, average direction vector, average angular smoothness, and center point.
 * The segment is classified based on these features (e.g., ground, vertical).
 */

void LidarSegment::computeSegAttribute()
{
    length = vecIndex.size();
    float ux = 0, uy = 0, uz = 0;
    for (std::size_t i = 0; i < length; i++)
    {
        int index = vecIndex[i];
        Eigen::Vector3f thisPoint(pLidarscan->mvPoints[index].point.x, pLidarscan->mvPoints[index].point.y, pLidarscan->mvPoints[index].point.z);

        ux += thisPoint(0);
        uy += thisPoint(1);
        uz += thisPoint(2);
        vecPoints.push_back(thisPoint);
    }

    for (std::size_t i = 1; i < length - 1; i++)
    {
        int index = vecIndex[i];
        vecAngle.push_back(computeAngle(pLidarscan->mvPoints[index].vecPrev, pLidarscan->mvPoints[index].vecNext));
    }
    uAngle = std::accumulate(vecAngle.begin(), vecAngle.end(), 0.0) / float(length - 2);

    direction = vecPoints.back() - vecPoints[0];

    ux /= float(length);
    uy /= float(length);
    uz /= float(length);
    centerPoint << ux, uy, uz;

    // firing range
    minFiring = pLidarscan->mvPoints[vecIndex[0]].firing;
    maxFiring = pLidarscan->mvPoints[vecIndex.back()].firing;

    // classify based on the length and angles
    classify();

    // //
    // pCurv = new float[length];
    // for (std::size_t i = 0; i < length; i++)
    //     pCurv[i] = -1.0;
    // if (type == 1)
    // {
    //     for (int i = N_Smooth; i < length - N_Smooth; i++)
    //     {
    //         float diffX = -float(N_Smooth * 2) * vecPoints[i](0),
    //               diffY = -float(N_Smooth * 2) * vecPoints[i](1),
    //               diffZ = -float(N_Smooth * 2) * vecPoints[i](2);

    //         for (int j = 1; j <= N_Smooth; j++)
    //         {
    //             diffX = diffX + vecPoints[i - j](0) + vecPoints[i + j](0);
    //             diffY = diffY + vecPoints[i - j](1) + vecPoints[i + j](1);
    //             diffZ = diffZ + vecPoints[i - j](2) + vecPoints[i + j](2);
    //         }
    //         pCurv[i] = diffX * diffX + diffY * diffY + diffZ * diffZ;
    //     }

    //     // find subsegment based on the curvature
    //     int s = 0;
    //     int e = -1;
    //     int i = 0;
    //     while (i < length - N_Smooth)
    //     {
    //         if (pCurv[i] > Smooth_Threshold) // edge point
    //         {
    //             e = i;
    //             if (e - s >= 5)
    //             {
    //                 SubSegment sub{.s = s, .e = e, .len = e - s};
    //                 vecSubsegment.push_back(sub);
    //             }
    //             s = i;
    //         }
    //         i++;
    //     }
    //     e = length - 1;
    //     if (e - s >= 5)
    //     {
    //         SubSegment sub{.s = s, .e = e, .len = e - s};
    //         vecSubsegment.push_back(sub);
    //     }
    // }
}

/* Attribute of Lidar segment has been calculated, based on the smooth derive attribute using points belong to plane*/
/**
 * @brief Computes plane-related attributes of a segment using smooth points.
 * 
 * Extracts low-smoothness points and computes their mean position to determine the segment’s planar properties.
 * 
 * @param thSmooth Smoothness threshold to filter points used in plane extraction.
 */

void LidarSegment::computePlaneAttribute(float thSmooth)
{
    float ux = 0, uy = 0, uz = 0;
    int nP = 0;
    for (std::size_t i = 0; i < length; i++)
    {
        int index = vecIndex[i];

        if (pLidarscan->mvPoints[index].smooth >= 0 && pLidarscan->mvPoints[index].smooth < thSmooth)
        {
            ux += vecPoints[i](0);
            uy += vecPoints[i](1);
            uz += vecPoints[i](2);

            nP++;
            vecPlanePoints.push_back(vecPoints[i]);
        }
    }
    planeLength = nP;
    ux /= float(planeLength);
    uy /= float(planeLength);
    uz /= float(planeLength);
    planeCenterPoint << ux, uy, uz;
}

// Exporting-------------------------------------------------------

// Export functions

//Export points in the scan with its attributes
/**
 * @brief Exports LiDAR scan points with attributes to a text file.
 * 
 * @param outPass File path for the exported data.
 */

void LidarScan::export_scan_point(const std::string outPass)
{
    std::ofstream fPointFile(outPass, std::ifstream::out);
    fPointFile << fixed << std::setprecision(4);
    fPointFile << "X_lu" << "\t" << "Y_lu" << "\t" << "Z_lu" <<"\t" << "range_valid" <<"\t" << "laser_beam_id"  <<"\t" 
            << "Firing_id"  <<"\t" << "smoothness" <<"\t" << "seg_id" <<"\t" << "seg_type" <<"\t" << "point_type" <<endl;
    for (size_t i = 0; i < mNPoints; ++i)
    {
        fPointFile << mvPoints[i].point.x << "\t" << mvPoints[i].point.y << "\t" << mvPoints[i].point.z << "\t"
                   << mvPoints[i].bValid << "\t" << mvPoints[i].channel << "\t" << mvPoints[i].firing << "\t"
                   << mvPoints[i].smooth << "\t" << mvPoints[i].segID << "\t" << mvPoints[i].seg_init_type << "\t"
                   << mvPoints[i].point_type_id << "\t" << endl;
        // mvPoints[i].neighborPrevious << "\t" << mvPoints[i].neighborNext << endl;
    }
    fPointFile.close();
}

//Export segments with its attributes
/**
 * @brief Exports segment-level information for all points in a LiDAR scan.
 * 
 * @param outPass File path for the exported segment data.
 * @param flag_level Whether to apply leveling transformation before exporting.
 */

void LidarScan::export_scan_segment(const std::string outPass, bool flag_level)
{
    std::ofstream fSegFile(outPass, std::ifstream::out);
    fSegFile << fixed << std::setprecision(8);
    fSegFile << "X_lu" << "\t" << "Y_lu" << "\t" << "Z_lu" <<"\t" << "laser_beam_id" <<"\t" << "seg_type"  <<"\t" 
            <<"ground_cluster_id" << "\t" << "ground_outlier_flag" <<"\t"
             <<"tree_cluster_id" << "\t" << "tree_outlier_flag" <<"\t"
              << "feature_type" << "\t" << "feature_id" <<endl;
    
    for (std::size_t nChannel = 0; nChannel < mNScan; nChannel++)
    {
        for (std::size_t segID = 0; segID < mvSegments[nChannel].size(); segID++)
        {
            for (std::size_t nP = 0; nP < mvSegments[nChannel][segID].length; nP++)
            {
                int index = mvSegments[nChannel][segID].vecIndex[nP];
                Eigen::Vector3f r_I_lu(mvPoints[index].point.x, mvPoints[index].point.y, mvPoints[index].point.z);
                Eigen::Vector3f r_I_lup;
                if (flag_level)
                    r_I_lup = R_lu_lup * r_I_lu;
                else
                    r_I_lup = r_I_lu;

                //fSegFile << mvPoints[index].point.x << "\t" << mvPoints[index].point.y << "\t" << mvPoints[index].point.z << "\t"
                fSegFile << r_I_lup(0) << "\t" << r_I_lup(1) << "\t" << r_I_lup(2) << "\t"
                         << mvPoints[index].channel << "\t" << mvSegments[nChannel][segID].initType << "\t"
                         << mvSegments[nChannel][segID].ground_cluster_id << "\t" << mvSegments[nChannel][segID].outlier_removal_flag << "\t"
                         << mvSegments[nChannel][segID].tree_cluster_id << "\t" << mvSegments[nChannel][segID].tree_outlier_removal_flag << "\t"
                         << mvSegments[nChannel][segID].final_feature_type << "\t" << mvSegments[nChannel][segID].final_feature_id << endl;
                //<< mvSegments[nChannel][segID].fClusterID << "\t" << mvPoints[index].smooth << "\t"
                // << mvSegments[nChannel][segID].fType << "\t" << mvSegments[nChannel][segID].initType << "\t"
                // << mvSegments[nChannel][segID].initType << "\t" << mvSegments[nChannel][segID].type << endl;
            }
        }
    }
    fSegFile.close();
}

//Export segments with its attributes
/**
 * @brief Exports final planar points (typically ground features) used in map construction.
 * 
 * @param outPass File path for output.
 */

void LidarScan::export_final_planar_points(const std::string outPass)
{
    std::ofstream f_planar_point(outPass, std::ifstream::out);
    f_planar_point << fixed << std::setprecision(8);
    f_planar_point << "X_lu" << "\t" << "Y_lu" << "\t" << "Z_lu" << "\t" << "scan_time" <<endl;
    
       for (int i = 0; i < pSurfPointsLessFlatGround_ds_->size(); i++)
    {
        f_planar_point  << pSurfPointsLessFlatGround_ds_->points[i].x << "\t" << pSurfPointsLessFlatGround_ds_->points[i].y << 
        "\t" << pSurfPointsLessFlatGround_ds_->points[i].z << "\t"
                 << pSurfPointsLessFlatGround_ds_->points[i].intensity << endl;
    }
    f_planar_point.close();
}

/**
 * @brief Exports all points with basic attributes such as coordinates, firing/channel info, and connectivity.
 * 
 * @param outPass File path for export.
 */

void LidarScan::exportPoint(const std::string outPass)
{
    std::ofstream fPointFile(outPass, std::ifstream::out);
    fPointFile << fixed << std::setprecision(4);

    for (size_t i = 0; i < mNPoints; ++i)
    {
        fPointFile << mvPoints[i].point.x << "\t" << mvPoints[i].point.y << "\t" << mvPoints[i].point.z << "\t"
                   << mvPoints[i].channel << "\t" << mvPoints[i].firing << "\t"
                   << mvPoints[i].segID << "\t" << mvPoints[i].neighborPrevious << "\t" << mvPoints[i].neighborNext << endl;
    }

    fPointFile.close();
}

/**
 * @brief Exports all points with basic attributes such as coordinates, firing/channel info, and connectivity.
 * 
 * @param outPass File path for export.
 */

void LidarScan::exportPointChannel(const std::string outPass)
{
    std::ofstream fPointFile(outPass, std::ifstream::out);
    fPointFile << fixed << std::setprecision(4);

    for (std::size_t nChannel = 0; nChannel < mNScan; nChannel++)
    {
        for (std::size_t nFiring = 0; nFiring < mNFiring; nFiring++)
        {
            std::size_t index = nFiring * mNScan + nChannel;
            fPointFile << mvPoints[index].point.x << "\t" << mvPoints[index].point.y << "\t" << mvPoints[index].point.z << "\t"
                       << index << "\t" << mvPoints[index].channel << "\t" << mvPoints[index].firing << "\t"
                       << mvPoints[index].segID << "\t" << mvPoints[index].neighborPrevious << "\t" << mvPoints[index].neighborNext << endl;
        }
    }

    fPointFile.close();
}

/**
 * @brief Exports extracted feature points including corner and planar types.
 * 
 * @param outPass Output path for the features file.
 */

void LidarScan::exportFeature(const std::string outPass)
{
    std::ofstream fPointFile(outPass, std::ifstream::out);
    fPointFile << fixed << std::setprecision(8);
    Eigen::Matrix3f Rot = Eigen::Matrix3f::Identity();
    // Rot = R_lu_lup;
    for (int nP = 0; nP < mpCornerPointsSharp->size(); nP++)
    {
        Eigen::Vector3f rilu(mpCornerPointsSharp->points[nP].x, mpCornerPointsSharp->points[nP].y, mpCornerPointsSharp->points[nP].z);
        Eigen::Vector3f rib2;
        rib2 = Rot * rilu;
        fPointFile << 1 << "\t" << rib2(0) << "\t" << rib2(1) << "\t" << rib2(2) << "\t" << mpCornerPointsSharp->points[nP].intensity << endl;
    }

    for (int nP = 0; nP < mpCornerPointsLessSharp->size(); nP++)
    {
        Eigen::Vector3f rilu(mpCornerPointsLessSharp->points[nP].x, mpCornerPointsLessSharp->points[nP].y, mpCornerPointsLessSharp->points[nP].z);
        Eigen::Vector3f rib2;
        rib2 = Rot * rilu;
        fPointFile << 2 << "\t" << rib2(0) << "\t" << rib2(1) << "\t" << rib2(2) << "\t" << mpCornerPointsLessSharp->points[nP].intensity << endl;
    }

    for (int nP = 0; nP < mpSurfPointsFlat->size(); nP++)
    {
        Eigen::Vector3f rilu(mpSurfPointsFlat->points[nP].x, mpSurfPointsFlat->points[nP].y, mpSurfPointsFlat->points[nP].z);
        Eigen::Vector3f rib2;
        rib2 = Rot * rilu;
        fPointFile << 3 << "\t" << rib2(0) << "\t" << rib2(1) << "\t" << rib2(2) << "\t" << mpSurfPointsFlat->points[nP].intensity << endl;
    }

    for (int nP = 0; nP < mpSurfPointsLessFlat->size(); nP++)
    {
        Eigen::Vector3f rilu(mpSurfPointsLessFlat->points[nP].x, mpSurfPointsLessFlat->points[nP].y, mpSurfPointsLessFlat->points[nP].z);
        Eigen::Vector3f rib2;
        rib2 = Rot * rilu;
        fPointFile << 4 << "\t" << rib2(0) << "\t" << rib2(1) << "\t" << rib2(2) << "\t" << mpSurfPointsLessFlat->points[nP].intensity << endl;
    }
    fPointFile.close();
}

/**
 * @brief Exports transformed LiDAR points (leveled frame) with associated features and segment data.
 * 
 * @param outPass Path for the output file.
 */

void LidarScan::exportPointChannelTransformed(const std::string outPass)
{
    std::ofstream fPointFile(outPass, std::ifstream::out);
    fPointFile << fixed << std::setprecision(8);

    for (std::size_t nChannel = 0; nChannel < mNScan; nChannel++)
    {
        for (std::size_t nFiring = 0; nFiring < mNFiring; nFiring++)
        {
            std::size_t index = nFiring * mNScan + nChannel;
            if (!mvPoints[index].bValid)
                continue;
            Eigen::Vector3f rilu(mvPoints[index].point.x, mvPoints[index].point.y, mvPoints[index].point.z);
            Eigen::Vector3f rib2;
            rib2 = R_lu_lup * rilu;
            // point-based: smooth2, clusterType
            if (mvPoints[index].segID < 0 && mvPoints[index].seg_init_type < 0)
                continue;
            //     if (mvPoints[index].clusterID < 0)
            //        continue;
            int treeid = -1;
            if (mvPoints[index].clusterType == 2)
                treeid = mvClusterTree[mvPoints[index].clusterID].featureId;
            fPointFile << rib2(0) << "\t" << rib2(1) << "\t" << rib2(2) << "\t" << mvPoints[index].clusterType << "\t" //<< mvPoints[index].type << "\t" << mvPoints[index].segID << "\t"
                       << treeid << "\t" << mvPoints[index].point_type_id << "\t"                                          //<< mvPoints[index].segID << "\t"<< mvPoints[index].clusterID
                       << mvPoints[index].clusterID << "\t" << mvPoints[index].channel << "\t" << endl;
        }
    }

    fPointFile.close();
}

/**
 * @brief Exports transformed LiDAR scan points into the global frame for map building.
 */

void LidarScan::exportPointMap()
{

    for (std::size_t nChannel = 0; nChannel < mNScan; nChannel++)
    {
        for (std::size_t nFiring = 0; nFiring < mNFiring; nFiring++)
        {
            std::size_t index = nFiring * mNScan + nChannel;
            if (!mvPoints[index].bValid)
                continue;

            Eigen::Vector3f rilu(mvPoints[index].point.x, mvPoints[index].point.y, mvPoints[index].point.z);
            Eigen::Vector3d rib;
            rib = T_lu_m * rilu.cast<double>();

            int treeid = -1;
            if (mvPoints[index].clusterType == 2)
                treeid = mvClusterTree[mvPoints[index].clusterID].featureId;

            fMappedPoints << id << "\t" << rib(0) << "\t" << rib(1) << "\t" << rib(2) << "\t" << mvPoints[index].clusterType << "\t"
                          // mvPoints[index].type << "\t" << mvPoints[index].segID << "\t"
                          << treeid << "\t" << endl; //<< mvPoints[index].featureId << "\t" //<< mvPoints[index].segID << "\t"<< mvPoints[index].clusterID
                                                     // << mvPoints[index].clusterID << "\t" << mvPoints[index].channel << "\t"
        }
    }
}

/**
 * @brief Exports all points with basic attributes such as coordinates, firing/channel info, and connectivity.
 * 
 * @param outPass File path for export.
 */

void LidarScan::exportPointPerChannel(const std::string folderName)
{

    for (std::size_t nChannel = 0; nChannel < mNScan; nChannel++)
    {
        std::string outPass = folderName + to_string(nChannel) + ".txt";
        std::ofstream fPointFile(outPass, std::ifstream::out);
        fPointFile << fixed << std::setprecision(4);

        for (std::size_t nFiring = 0; nFiring < mNFiring; nFiring++)
        {
            std::size_t index = nFiring * mNScan + nChannel;
            fPointFile << mvPoints[index].point.x << "\t" << mvPoints[index].point.y << "\t" << mvPoints[index].point.z << "\t"
                       << index << "\t" << mvPoints[index].firing << "\t"
                       << mvPoints[index].segID << "\t" << mvPoints[index].neighborPrevious << "\t" << mvPoints[index].neighborNext << endl;
        }
        fPointFile.close();
    }
}

/**
 * @brief Exports segment information such as center, length, angle, and transformed coordinates.
 * 
 * @param folderName Output directory or prefix for file generation.
 */

void LidarScan::exportSegmentPerChannel(const std::string folderName)
{

    for (std::size_t nChannel = 0; nChannel < mNScan; nChannel++)
    {
        std::string outPass = folderName + "Seg" + to_string(nChannel) + ".txt";
        std::ofstream fSegFile(outPass, std::ifstream::out);
        fSegFile << fixed << std::setprecision(4);

        for (std::size_t segID = 0; segID < mvSegments[nChannel].size(); segID++)
        {
            for (std::size_t nP = 0; nP < mvSegments[nChannel][segID].length; nP++)
            {
                // mvPoints[index].neighborPrevious                         << "\t" << mvPoints[index].neighborNext << index << "\t" << mvPoints[index].firing << "\t"
                float angle = (nP == 0 || nP == mvSegments[nChannel][segID].length - 1) ? 0 : mvSegments[nChannel][segID].vecAngle[nP - 1];
                int index = mvSegments[nChannel][segID].vecIndex[nP];
                fSegFile << mvPoints[index].point.x << "\t" << mvPoints[index].point.y << "\t" << mvPoints[index].point.z << "\t"
                         << rad2deg(mvSegments[nChannel][segID].uAngle) << "\t"
                         << mvPoints[index].segID << "\t" << rad2deg(angle) << endl;
            }
        }

        fSegFile.close();
    }
}

/**
 * @brief Exports segment information such as center, length, angle, and transformed coordinates.
 * 
 * @param folderName Output directory or prefix for file generation.
 */

void LidarScan::exportSegmentInfo(const std::string outPass)
{
    std::ofstream fSegFile(outPass, std::ifstream::out);
    fSegFile << fixed << std::setprecision(4);

    for (std::size_t nChannel = 0; nChannel < mNScan; nChannel++)
    {
        for (std::size_t segID = 0; segID < mvSegments[nChannel].size(); segID++)
        {
            // int index = mvSegments[nChannel][segID].vecIndex[nP];
            Eigen::Vector3f r;
            r = R_lu_lup * mvSegments[nChannel][segID].centerPoint;
            fSegFile << r(0) << "\t" << r(1) << "\t"
                     << r(2) << "\t" << segID << "\t" << mvSegments[nChannel][segID].length
                     << endl;

            // fSegFile << mvSegments[nChannel][segID].centerPoint(0) << "\t" << mvSegments[nChannel][segID].centerPoint(1) << "\t"
            //          << mvSegments[nChannel][segID].centerPoint(2) << "\t" << segID << "\t" << mvSegments[nChannel][segID].length
            //          <<  endl;
        }
    }
    fSegFile.close();
}

/**
 * @brief Exports segment information such as center, length, angle, and transformed coordinates.
 * 
 * @param folderName Output directory or prefix for file generation.
 */

void LidarScan::exportSegmentTransInfo(const std::string outPass)
{
    std::ofstream fSegFile(outPass, std::ifstream::out);
    fSegFile << fixed << std::setprecision(4);

    for (std::size_t nChannel = 0; nChannel < mNScan; nChannel++)
    {
        for (std::size_t segID = 0; segID < mvSegments[nChannel].size(); segID++)
        {
            // int index = mvSegments[nChannel][segID].vecIndex[nP];
            fSegFile << mvSegments[nChannel][segID].centerPointTrans(0) << "\t" << mvSegments[nChannel][segID].centerPointTrans(1) << "\t"
                     << mvSegments[nChannel][segID].centerPointTrans(2) << "\t" << segID << "\t" << mvSegments[nChannel][segID].length << endl;
        }
    }
    fSegFile.close();
}

/**
 * @brief Exports segment information such as center, length, angle, and transformed coordinates.
 * 
 * @param folderName Output directory or prefix for file generation.
 */

void LidarScan::exportSegmentInfoPerChannel(const std::string folderName)
{
    for (std::size_t nChannel = 0; nChannel < mNScan; nChannel++)
    {
        std::string outPass = folderName + "InfoSeg" + to_string(nChannel) + ".txt";
        std::ofstream fSegFile(outPass, std::ifstream::out);
        fSegFile << fixed << std::setprecision(4);

        for (std::size_t segID = 0; segID < mvSegments[nChannel].size(); segID++)
        {
            // int index = mvSegments[nChannel][segID].vecIndex[nP];
            fSegFile << mvSegments[nChannel][segID].centerPoint(0) << "\t" << mvSegments[nChannel][segID].centerPoint(1) << "\t"
                     << mvSegments[nChannel][segID].centerPoint(2) << "\t" << segID << "\t" << mvSegments[nChannel][segID].length
                     << "\t" << mvSegments[nChannel][segID].uAngle * 180.0 / PI << endl;
        }

        fSegFile.close();
    }
}

/**
 * @brief Exports final segmented feature points with cluster and feature IDs, optionally in transformed global frame.
 * 
 * @param outPass File path for export.
 */

void LidarScan::exportSegments(const std::string outPass)
{
    std::ofstream fSegFile(outPass, std::ifstream::out);
    fSegFile << fixed << std::setprecision(8);

    for (std::size_t nChannel = 0; nChannel < mNScan; nChannel++)
    {
        for (std::size_t segID = 0; segID < mvSegments[nChannel].size(); segID++)
        {
            for (std::size_t nP = 0; nP < mvSegments[nChannel][segID].length; nP++)
            {
                // float angle = (nP == 0 || nP == mvSegments[nChannel][segID].length - 1) ? 0 : mvSegments[nChannel][segID].vecAngle[nP - 1];

                int index = mvSegments[nChannel][segID].vecIndex[nP];
                // fSegFile << mvPoints[index].point.x << "\t" << mvPoints[index].point.y << "\t" << mvPoints[index].point.z << "\t"
                //          << segID << "\t" << nChannel << "\t" << mvSegments[nChannel][segID].type << endl;
                fSegFile << mvPoints[index].point.x << "\t" << mvPoints[index].point.y << "\t" << mvPoints[index].point.z << "\t"
                         //<< mvSegments[nChannel][segID].fClusterID << "\t" << mvPoints[index].smooth << "\t"
                         // << mvSegments[nChannel][segID].fType << "\t" << mvSegments[nChannel][segID].initType << "\t"
                         << mvSegments[nChannel][segID].initType << endl;
            }
        }
    }
    fSegFile.close();
}

/**
 * @brief Exports final segmented feature points with cluster and feature IDs, optionally in transformed global frame.
 * 
 * @param outPass File path for export.
 */

void LidarScan::exportFinal(const std::string outPass)
{
    std::ofstream fSegFile(outPass, std::ifstream::out);
    fSegFile << fixed << std::setprecision(4);

    for (std::size_t nChannel = 0; nChannel < mNScan; nChannel++)
    {
        for (std::size_t segID = 0; segID < mvSegments[nChannel].size(); segID++)
        {
            for (std::size_t nP = 0; nP < mvSegments[nChannel][segID].length; nP++)
            {

                int index = mvSegments[nChannel][segID].vecIndex[nP];

                Eigen::Vector3f rilu(mvPoints[index].point.x, mvPoints[index].point.y, mvPoints[index].point.z);
                Eigen::Vector3f rib;
                rib = R_lu_lup * rilu;

                fSegFile << rib(0) << "\t" << rib(1) << "\t" << rib(2) << "\t" << mvSegments[nChannel][segID].final_feature_type << "\t"
                         << mvSegments[nChannel][segID].final_feature_id << endl;
            }
        }
    }
    fSegFile.close();
}

/**
 * @brief Exports final segmented feature points with cluster and feature IDs, optionally in transformed global frame.
 * 
 * @param outPass File path for export.
 */

void LidarScan::exportFinalMapping(const std::string outPass)
{
    std::ofstream fSegFile(outPass, std::ifstream::out);
    fSegFile << fixed << std::setprecision(4);

    for (std::size_t nChannel = 0; nChannel < mNScan; nChannel++)
    {
        for (std::size_t segID = 0; segID < mvSegments[nChannel].size(); segID++)
        {
            for (std::size_t nP = 0; nP < mvSegments[nChannel][segID].length; nP++)
            {

                int index = mvSegments[nChannel][segID].vecIndex[nP];

                Eigen::Vector3f rilu(mvPoints[index].point.x, mvPoints[index].point.y, mvPoints[index].point.z);
                Eigen::Vector3d rib;
                //                rib = mT * (mR_lu_lulevel * rilu);
                //  rib = T_lu_m * rilu;
                rib = R_lu_m * rilu.cast<double>() + r_lu_m;
                fSegFile << rib(0) << "\t" << rib(1) << "\t" << rib(2) << "\t" << mvSegments[nChannel][segID].final_feature_type << "\t"
                         << mvSegments[nChannel][segID].final_feature_id << endl;
                // if (mvSegments[nChannel][segID].fType != 0)
                //     fMappedPoints << id << "\t" << rib(0) << "\t" << rib(1) << "\t" << rib(2) << "\t" << mvSegments[nChannel][segID].fType << "\t"
                //                   << mvSegments[nChannel][segID].fClusterID << endl;
            }
        }
    }
    fSegFile.close();
}

/**
 * @brief Exports planar ground points with optional leveling applied.
 * 
 * @param outPass Output file path.
 * @param flagLevel If true, applies leveling rotation to point coordinates.
 */

void LidarScan::exportGroundPoints(const std::string outPass, bool flagLevel)
{
    std::ofstream fGroundPoints(outPass, std::ifstream::out);
    fGroundPoints << fixed << std::setprecision(4);

    for (std::size_t nP = 0; nP < pSurfPointsLessFlatGround_ds_->size(); nP++)
    {
        Eigen::Vector3f rilu(pSurfPointsLessFlatGround_ds_->points[nP].x, pSurfPointsLessFlatGround_ds_->points[nP].y, pSurfPointsLessFlatGround_ds_->points[nP].z);
        Eigen::Vector3f rilup;
        if (flagLevel)
            rilup = R_lu_lup * rilu;
        else
            rilup = rilu;

        fGroundPoints << id << "\t" << rilup(0) << "\t" << rilup(1) << "\t" << rilup(2) << "\t" << pSurfPointsLessFlatGround_ds_->points[nP].intensity << endl;
    }
    fGroundPoints.close();
}

/**
 * @brief Exports planar ground points with optional leveling applied.
 * 
 * @param outPass Output file path.
 * @param flagLevel If true, applies leveling rotation to point coordinates.
 */

void LidarScan::exportGroundPoints2(const std::string outPass, bool flagLevel)
{
    std::ofstream fGroundPoints(outPass, std::ifstream::out);
    fGroundPoints << fixed << std::setprecision(4);

    for (std::size_t nP = 0; nP < pSurfPointsLessFlatGround_->size(); nP++)
    {
        Eigen::Vector3f rilu(pSurfPointsLessFlatGround_->points[nP].x, pSurfPointsLessFlatGround_->points[nP].y, pSurfPointsLessFlatGround_->points[nP].z);
        Eigen::Vector3f rilup;
        if (flagLevel)
            rilup = R_lu_lup * rilu;
        else
            rilup = rilu;

        fGroundPoints << id << "\t" << rilup(0) << "\t" << rilup(1) << "\t" << rilup(2) << "\t" << pSurfPointsLessFlatGround_->points[nP].intensity << endl;
    }
    fGroundPoints.close();
}

/**
 * @brief Exports all final scan data (features only) in global coordinates for mapping.
 */

void LidarScan::exportScan()
{
    for (std::size_t nChannel = 0; nChannel < mNScan; nChannel++)
    {
        for (std::size_t segID = 0; segID < mvSegments[nChannel].size(); segID++)
        {
            for (std::size_t nP = 0; nP < mvSegments[nChannel][segID].length; nP++)
            {

                int index = mvSegments[nChannel][segID].vecIndex[nP];
                Eigen::Vector3f rilu(mvPoints[index].point.x, mvPoints[index].point.y, mvPoints[index].point.z);
                Eigen::Vector3d rib;
                rib = T_lu_m * rilu.cast<double>();
                // rib= rilu;
                if (mvSegments[nChannel][segID].final_feature_type != 0)
                    fMappedPoints << id << "\t" << rib(0) << "\t" << rib(1) << "\t" << rib(2) << "\t" << mvSegments[nChannel][segID].final_feature_type << "\t"
                                  << mvSegments[nChannel][segID].final_feature_id << endl;
            }
        }
    }
}

/**
 * @brief Exports transformed LiDAR scan points into the global frame for map building.
 */

void LidarScan::exportRefScan()
{
    Eigen::Vector3d r_b_m_, r_lu_m_;
    Eigen::Matrix3d R_b_m_, R_lu_m_;
    int startIndex = -1;
    if (!mpTraj->bopInterpolation(mTimeInit / 1000.0, -1, startIndex, r_b_m_, R_b_m_))
        return;
    // cout << "Scan start index" << startIndex << "\t" <<mTime / 1000.0<<endl;

    for (std::size_t nChannel = 0; nChannel < mNScan; nChannel++)
    {
        for (std::size_t segID = 0; segID < mvSegments[nChannel].size(); segID++)
        {
            for (std::size_t nP = 0; nP < mvSegments[nChannel][segID].length; nP++)
            {

                int index = mvSegments[nChannel][segID].vecIndex[nP];

                Eigen::Vector3d rilu(mvPoints[index].point.x, mvPoints[index].point.y, mvPoints[index].point.z);

                double time = mTimeInit / 1000.0 + mScanDuration / 1000.0 * double(mvPoints[index].point.intensity - int(mvPoints[index].point.intensity));
                Eigen::Vector3d r_b_m_, r_lu_m_;
                Eigen::Matrix3d R_b_m_, R_lu_m_;
                int idx;
                if (mpTraj->bopInterpolation(time, startIndex, idx, r_b_m_, R_b_m_))
                {
                    r_lu_m_ = r_b_m_ + R_b_m_ * mpTraj->r_lu_b;
                    R_lu_m_ = R_b_m_ * mpTraj->R_lu_b;
                    Eigen::Vector3d rib;
                    rib = R_lu_m_ * rilu + r_lu_m_;

                    // rib= rilu;
                    //  if (mvSegments[nChannel][segID].fType != 0)
                    fMappedRef << id << "\t" << rib(0) << "\t" << rib(1) << "\t" << rib(2) << "\t" << time << "\t" << mvSegments[nChannel][segID].final_feature_type << "\t" << endl;
                }
                else
                {
                    //     cout << mTime << endl;
                    //     cout << mTimeScan * double(mvPoints[index].point.intensity - int(mvPoints[index].point.intensity)) << endl;
                    //     cout << time << endl;
                    throw std::runtime_error("Wrong Time");
                }
            }
        }
    }
}

/**
 * @brief Exports structural information about identified tree clusters.
 * 
 * @param outPass Output path for writing cluster segment centers and associations.
 */

void LidarScan::exportTreeClusterInfo(const std::string outPass)
{
    std::ofstream fClusterFile(outPass, std::ifstream::out);
    fClusterFile << fixed << std::setprecision(4);
    for (std::size_t nC = 0; nC < mvClusterTree.size(); nC++)
    {
        for (std::size_t nSeg = 0; nSeg < mvClusterTree[nC].len; nSeg++)

            fClusterFile << mvClusterTree[nC].segCenters[nSeg](0) << "\t" << mvClusterTree[nC].segCenters[nSeg](1) << "\t"
                         << mvClusterTree[nC].segCenters[nSeg](2) << "\t" << nC << "\t" << mvClusterTree[nC].vSegments[nSeg].first
                         << "\t" << mvClusterTree[nC].vSegments[nSeg].second << endl;
    }
    fClusterFile.close();
}
