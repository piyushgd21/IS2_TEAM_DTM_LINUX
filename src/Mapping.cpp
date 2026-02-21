#include "../header/Mapping.h"
#include "../header/cylinder_model_fit.h"

constexpr double kDownsamplingSize = 0.01; //0.10 outperforms than 0.15
constexpr double KDownsampleLCGround = 0.01;
constexpr double KDownsampleSizeRawLevelGround = 0.05;
constexpr double KDownsampleSizeRawLevelTree = 0.05;
// constexpr double KRef_ground_std = 0.005; // previous:0.005 
constexpr double KGround_resolution = 1.5; //Previous:1.0 2.0

namespace {

// uint64_t GetkeyIndex(const Eigen::Vector3d& point)
// {
// 	double cell_size = 0.01;
// 	uint64_t length = 65535;
// 	uint64_t x_index = static_cast<uint64_t>(point(0) / cell_size);
// 	uint64_t y_index = static_cast<uint64_t>(point(1) / cell_size);
// 	uint64_t z_index = static_cast<uint64_t>(point(2) / cell_size);
// 	uint64_t x_index64 = static_cast<uint64_t>(x_index%length);
// 	uint64_t y_index64 = static_cast<uint64_t>(y_index%length);
// 	uint64_t z_index64 = static_cast<uint64_t>(z_index%length);
// 	uint64_t key_index = (x_index64 << 32) + (y_index64 << 16) + z_index64;
// 	return key_index;
// }

uint64_t GetVoxelkeyIndex(const Eigen::Vector3d& point,double voxel_size)
{
	double cell_size = voxel_size;
	uint64_t length = 65535;
	uint64_t x_index = static_cast<uint64_t>(point(0) / cell_size);
	uint64_t y_index = static_cast<uint64_t>(point(1) / cell_size);
	uint64_t z_index = static_cast<uint64_t>(point(2) / cell_size);
	uint64_t x_index64 = static_cast<uint64_t>(x_index%length);
	uint64_t y_index64 = static_cast<uint64_t>(y_index%length);
	uint64_t z_index64 = static_cast<uint64_t>(z_index%length);
	uint64_t key_index = (x_index64 << 32) + (y_index64 << 16) + z_index64;
	return key_index;
}

Eigen::Matrix4f Trajectory_ICP(const pcl::PointCloud<pcl::PointXYZ>& comapred_points, const pcl::PointCloud<pcl::PointXYZ>& reference_points)
{
	pcl::IterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ> icp;
	icp.setInputSource(comapred_points.makeShared());
	icp.setInputTarget(reference_points.makeShared());
	icp.setMaxCorrespondenceDistance(10.0);
	icp.setMaximumIterations(100);

	pcl::PointCloud<pcl::PointXYZ> Final;
	icp.align(Final);

	if (!icp.hasConverged())
		throw std::runtime_error("Wrong initialization");

	//fMapLog << icp.getFinalTransformation() << std::endl;
	return icp.getFinalTransformation();
}
}
/*main function of build map*/
/**
 * @brief Main function for building the map from LiDAR scans.
 *
 * This function performs the complete mapping pipeline:
 *  - Integration of raw LiDAR scans
 *  - Feature extraction
 *  - Pose estimation and refinement
 *  - Tree model fitting and long-term optimization
 *  - Optional global map alignment and trajectory export
 *  - Final export of mapping results (ground, trees, trajectories)
 *
 * @note Includes logic for intermediate loop closure, raw-level optimization,
 *       fetching and updating tree observations, and merging map trees.
 *
 * @throws std::runtime_error on trajectory time mismatch.
 */
void Mapping::buildMap()
{
    std::string outPass = output_folder + "TrajectoryMapping.txt";
	last_cut_scan_id_ = mPara.initScan;

    fMappingTraj.open(outPass, std::ifstream::out);
    fMappingTraj << fixed << std::setprecision(4);

    double t_integration_ = 0.0, t_scan_to_map_ = 0.0, t_add_feature_to_map_ = 0.0, t_loop_closure_ = 0.0;

    //--------------- Setting related to the case when global map is available 
    //check if global map is available
    if(global_map_flag_)
    {
        //derive constant shift from trajectory or map (TODO)
        global_const_shift_ << floor(mpTraj->bopList[0].pos.XO),floor(mpTraj->bopList[0].pos.YO), floor(mpTraj->bopList[0].pos.ZO);
        fMapLog << "Applied shift: " << global_const_shift_.transpose() <<endl;

        //load global map
        LoadGlobalGroundMap();

        //load tree map
        if(!DTM_only_flag_)
        {
            LoadGlobalTreeMap();
        }
    }

    int terminateCount = 0;
    mnCount = 0;

    // for initialization, more scan are used.
    if(global_map_flag_)
    {
        mNumScan = mNumIntegratedScan;
    }
    else
    {
        mNumScan = mNumInitialScan;
    }
    mMaxIntegratedTime = double(mNumScan) * 1.5 * SCAN_DURATION; // 1.5 seconds in ms




    TicToc t_Mapping;

    // 20 seconds without receiving
	int max_count = intermediate_result_flag?1000:200;
    while (terminateCount < max_count) // 0.1 s * countMax
    {

        // this does not hold for non-real time cases as time of sending scan depends on processing time only
        // if (mqScan.size() >= mNumScan || (terminateCount >= int(mMaxIntegratedTime / SCAN_DURATION) && !mqScan.empty())) 
        if (mqScan.size() >= mNumScan * 2 || (terminateCount >= int(mMaxIntegratedTime / SCAN_DURATION) && !mqScan.empty())) 
        {
            TicToc t_mapping;

            unique_lock<mutex> lock(mTransMutex);

            // int numSavedScan = min(mNumScan, int(mqScan.size()));
            // double t_ini = mqScan.front().mTimeInit;
            // for (int i = 0; i < numSavedScan; i++)
            // {
            //     if (mqScan.front().mTimeInit - t_ini > mMaxIntegratedTime) // within 1.5 second
            //         break;

            //     mvTempScan.push_back(mqScan.front());
            //     mqScan.pop();
            // }

            double t_ini = mqScan.front().mTimeInit;
            int num_valid_scan = 0;
            while(1)
            {
                //check!
                // if (mqScan.front().mTimeInit - t_ini > mNumScan) * 1.5 * SCAN_DURATION) // within 1.5 * num of valid scan
                //     break;
                
                //if no trajectory info is used, count number of valid scans
                if(mqScan.front().valid_feature_flag || mPara.odo_from_trajectory_flag)
                {
                    num_valid_scan++;
                }
                //mvTempScan.push_back(mqScan.front());
                //mqScan.pop();

				int current_scan_id = mqScan.front().scanID;
				int number_in_scans = current_scan_id - last_cut_scan_id_;
				//bool cut = false;
				if (number_in_scans >= mNumInitialScan )
				{
					//cut = true;
					number_integrated_scan_++;
					last_cut_scan_id_ = number_integrated_scan_*mNumInitialScan+ mPara.initScan;
					std::cout << "last_cut_scan_id_: " << last_cut_scan_id_ << std::endl;
					break;
				} 

				mvTempScan.push_back(mqScan.front());
				mqScan.pop();
                
                //if(num_valid_scan >= mNumScan)
                //{
                //    break;
                //}

                if(mqScan.empty())
                {
                    break;
                }
            }
            lock.unlock();
            int num_scan = mvTempScan.size();

            // 1. integrate different scans into trees % groud patches
            TicToc tIntegration;
            integrateScan();
            fMapLog << "Total time for scan integration " << tIntegration.toc() / 1000.0 << " s " << endl;
            mvpIScans.push_back(pCurrentIScan);
            t_integration_ += tIntegration.toc() / 1000.0;

            // 2. compute transformation from integrated scan results to map
            TicToc tScantoMap;
            computePoseToMap();
            fMapLog << "Total time for compute Iscan to map " << tScantoMap.toc() / 1000.0 << " s " << endl;
            t_scan_to_map_ += tScantoMap.toc() / 1000.0;

            // 3. add feature to map
            TicToc t_add_feature;
            addFeaturetoMap();
            fMapLog << "Total time for add feature to map " << t_add_feature.toc() / 1000.0 << " s " << endl;
            t_add_feature_to_map_ += t_add_feature.toc()/1000.0;

            //4. save current Iscan to mvpIScans
            //pCurrentIScan->show_number_PointType();
            pCurrentIScan = nullptr;

            //5. update mapTree status
            // for (int nPair = 0; nPair < mvPairs.size(); nPair++)
            // {
            //     int mapTreeId = mvPairs[nPair].first;
            //     mpMapTree[mapTreeId]->is_established();
            // }
  
            //6. setting for next Iscan 
            if (!mbInit)
            {
                mbInit = true;
                mNumScan = mNumIntegratedScan;
                mMaxIntegratedTime = double(mNumScan) * 1.5 * SCAN_DURATION; // 1.5 seconds in ms
            }

            terminateCount = 0;
            mnCount++;
            fMapLog << "--------------------------------------------------" << endl;
            fMapLog << "Total time: " <<  t_mapping.toc() / 1000.0 << " s for " << num_scan << " scans in this iteration" << endl;
            fMapLog <<endl;

            //merge map tree
			//(MODI/TODO)
            if(mnCount % 1 == 0) //mnCount % 10 %5
            {
                fMapLog << "*******************************************************" << endl;
                TicToc t_merging;
                fMapLog << "Merge Tree" << endl;
                CountMapTreeType();
                mergeMapTree();
                CountMapTreeType();
                fMapLog << "Time for merging " <<  t_merging.toc() <<endl;
                fMapLog <<endl;

            }

            //long term optimization
            if (mPara.loop_closure_iscan_num > 0)
            {
                if (mnCount % mPara.loop_closure_iscan_num == 0)
                {
                    TicToc t_loopclosure;
                    ConductIntermediateLC();
                    t_loop_closure_ += t_loopclosure.toc() / 1000.0;

                }
            }

            // reset tree parameters for INIT, FITTED and TBD, prevent from fitting in wrong local miminum
            for (int nMapT = 0; nMapT < mpMapTree.size(); nMapT++)
            {
                if (mpMapTree[nMapT]->status == MapTree::FITTED || mpMapTree[nMapT]->status == MapTree::TBD || mpMapTree[nMapT]->status == MapTree::INIT)
                {
                    mpMapTree[nMapT]->ResetPara();
                }
            }
        }
        else
        {
            terminateCount++;
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
    }

    //*************** Conduct final raw level long term optimization ************
    if (mPara.loop_closure_iscan_num > 0)//&& iscan_index_end_ != mnCount - 1)
    {
        if (mPara.perform_iscan_lc_at_last_flag)
        {
            TicToc t_loopclosure;
            ConductIntermediateLC();
            t_loop_closure_ += t_loopclosure.toc() / 1000.0;
        }

        // conduct raw level optimization
        TicToc t_loopclosure;
        prefix = "Final_LC_";
        ConductFinalLC(mPara.number_iterations_raw_scan_lc);
        t_loop_closure_ += t_loopclosure.toc() / 1000.0;
    }

    //*******************************Fetch raw tree points and conduct long term optimization ---
    if(mPara.reoptimize_fetch_tree_points)
    {
        int iter = 0, max_iter = mPara.reoptimize_iterations;
        double distance_threshold = mPara.distance_threshold_fetch_points;
        while(iter < max_iter)
        {   
            fMapLog<<endl;
            fMapLog << "*******************************************************" << endl;
            fMapLog << "* Conduct optimization round after update tree observations using tree points from odometry " << ++iter << endl;
            cout << "Conduct optimization round after update tree observations using tree points from odometry " << iter << endl;
            fMapLog << "*******************************************************" << endl;
            prefix = "Final_LC_Round_" +to_string(iter) +"_";

            // derive the tree points in mapping
            for (int nIscan = 0; nIscan < mvpIScans.size(); nIscan++)
            {
                mvpIScans[nIscan]->computeBackupTreePointToMapping();
            }

            // match points to tree
            FindRawTreePointsPerMapTree(distance_threshold);
            distance_threshold = max(0.3, distance_threshold-0.2);

            //Then, update vTreePointRaw per iscan
            for (int nIscan = 0; nIscan < mvpIScans.size(); nIscan++)
            {
                mvpIScans[nIscan]->update_vTreePointRaw();
            }

            //checking only
            std::string outBackupTreeFeature(output_folder + "MapBackupTreeMapping_Round_" + to_string(iter) + ".las");
            ExportBackupTreeFeatureMappingLas(outBackupTreeFeature);

            //reform tree observations and update tree status
            Update_Map_Tree_Obs();
            CountMapTreeType();

            //conduct optimization
            TicToc t_loopclosure;
            ConductFinalLC(mPara.number_iterations_raw_scan_lc);
            t_loop_closure_ += t_loopclosure.toc() / 1000.0;
        }
    }

    //*******************************Fetch raw tree points and conduct long term optimization ---
    if(mPara.reoptimize_fetch_raw_points)
    {
        int iter = 0, max_iter = mPara.reoptimize_iterations;
        double distance_threshold = mPara.distance_threshold_fetch_points;
        while(iter < max_iter)
        {   
            fMapLog<<endl;
            fMapLog << "*******************************************************" << endl;
            fMapLog << "* Conduct optimization round after update tree observations using all non ground points " << ++iter << endl;
            cout << "* Conduct optimization round after update tree observations using all non ground points " << iter << endl;

            fMapLog << "*******************************************************" << endl;
            prefix = "Final_LC_Round_" +to_string(iter) +"_";

            // derive the tree points in mapping
            for (int nIscan = 0; nIscan < mvpIScans.size(); nIscan++)
            {
                mvpIScans[nIscan]->computeNonGroundPointsToMapping();
            }

            TicToc t_find_raw_points;
            FindRawPointsPerMapTree(distance_threshold);
            fMapLog << "\t" << t_find_raw_points.toc() / 1000.0 << " seconds" <<endl;
            distance_threshold = max(0.3, distance_threshold-0.2);

            // //checking only
            // std::string outBackupTreeFeature(output_folder + "MapBackupTreeMapping_Round_" + to_string(iter) + ".las");
            // ExportBackupTreeFeatureMappingLas(outBackupTreeFeature);

            //reform tree observations and update tree status
            Update_Map_Tree_Obs();
            CountMapTreeType();

            //conduct optimization
            TicToc t_loopclosure;
            ConductFinalLC(mPara.number_iterations_raw_scan_lc);
            t_loop_closure_ += t_loopclosure.toc() / 1000.0;
        }
    }

    //---------------------------------- Exporting the results -----------------------------------
    //compute residual for each tree
    for (int nMapT = 0; nMapT < mpMapTree.size(); nMapT++)
    {
        if (mpMapTree[nMapT]->status >= 0)
        {
            mpMapTree[nMapT]->rmse =computeTreeToTreeResidual(nMapT,nMapT, mpMapTree[nMapT]->residual_distribution);
        }
    }

    // derive backup tree points again
    for (int nIscan = 0; nIscan < mvpIScans.size(); nIscan++)
    {
        mvpIScans[nIscan]->computeBackupTreePointToMapping();
    }
    // match points to tree
    FindRawTreePointsPerMapTree();

    //end of mapping, compute the mapping to global if available
    if (mbTraj)
    {
        int index_for_map_to_global = 5;
      
        Eigen::Vector3d r_b_global_t1, r_lu_global_t1, r_lu_m_t1;
        Eigen::Matrix3d R_b_global_t1, R_lu_global_t1, R_lu_m_t1;

        int idx = 0;

        if (!mpTraj->bopInterpolation((mvpIScans[0]->indScans[index_for_map_to_global].mTimeEnd - mvpIScans[0]->indScans[index_for_map_to_global].mTimeTracking) / 1000.0, -1, idx, r_b_global_t1, R_b_global_t1))
            throw std::runtime_error("Wrong Time");

        r_lu_global_t1 = r_b_global_t1 + R_b_global_t1 * mpTraj->r_lu_b;
        R_lu_global_t1 = R_b_global_t1 * mpTraj->R_lu_b;

        r_lu_m_t1 = mvpIScans[0]->v_r_mapping_ref[index_for_map_to_global];
        R_lu_m_t1 = mvpIScans[0]->v_R_mapping_ref[index_for_map_to_global];        

        R_m_global_ref_traj = R_lu_global_t1 * R_lu_m_t1.transpose() ;
        r_m_global_ref_traj = r_lu_global_t1 - R_m_global_ref_traj * r_lu_m_t1;

        r_lu_m_t1 = mvpIScans[0]->v_r_mapping[index_for_map_to_global];
        R_lu_m_t1 = mvpIScans[0]->v_R_mapping[index_for_map_to_global];      

        R_m_global = R_lu_global_t1 * R_lu_m_t1.transpose() ;
        r_m_global = r_lu_global_t1 - R_m_global * r_lu_m_t1;

        // R_m_global = R_m_global_ref_traj;
        // r_m_global = r_m_global_ref_traj;

        if(global_map_flag_)
        {
            R_m_global = Eigen::Matrix3d::Identity();
            r_m_global = global_const_shift_;
        }
        
        fMapLog<< " Trans from mapping to global: " << r_m_global.transpose() << "\t" << rad2deg(Find_Rotation(R_m_global)).transpose() <<endl;

    }
    fMapLog << t_Mapping.toc()/1000.0 << " s"<< endl;
    fMapLog << "Scan Integration: " << t_integration_ << " s"<< endl;
    fMapLog << "Scan to Map Transformation: " << t_scan_to_map_ << " s"<< endl;
    fMapLog << "Add feature to Map: " << t_add_feature_to_map_ << " s"<< endl;
    fMapLog << "Long term optimization: " << t_loop_closure_ << " s"<< endl;

    std::string outTreeParamCheck(output_folder + "MapTreeParamCheck" + ".txt");
    exportMapTreeParam(outTreeParamCheck);

    computeMapTreeInfo();
    


    // export results
    std::string outTreeModel(output_folder + "MapTreeModel" + ".txt");
    exportMapTreeModel(outTreeModel);
    std::string outTreeParam(output_folder + "MapTreeParam" + ".txt");
    exportMapTreeParam(outTreeParam);

    std::string outGroundMap(output_folder + "MapGround" + ".txt");
    exportGroundMap(outGroundMap);

    // for (int nTree = 0; nTree < mpMapTree.size(); nTree++)
    // {
    //     mpMapTree[nTree]->ComputeCentroid();
    //     mpMapTree[nTree]->compute_surface_ratio();
    // }
    // std::string outTreeParam2(output_folder + "MapTreeParam2222222" + ".txt");
    // exportMapTreeParam(outTreeParam2);

    //test
    // for (int nIscan = 0; nIscan < mvpIScans.size(); nIscan++)
    // {
    //     mvpIScans[nIscan]->computeBackupTreePointToMapping();
    // }

	pcl::PointCloud<pcl::PointXYZ> slam_trajectory;
	pcl::PointCloud<pcl::PointXYZ> gnss_trajectory;
    if (mbTraj)
    {
        //(TODO)(MODI):traj export traj firstly
        std::string outTrajGlobalCheck(output_folder + "TrajectoryGlobalCheck" + ".txt");
        exportTrajGlobalCheck(outTrajGlobalCheck);
        std::string outTrajGlobal(output_folder + "TrajectoryGlobal" + ".txt");
		if (!global_map_flag_)
		{
			exportTrajGlobal(outTrajGlobal,&slam_trajectory,&gnss_trajectory);
		}
		else
		{
			exportTrajGlobal(outTrajGlobal);
		}


        std::string outBackupTreeFeature(output_folder + "MapBackupTreeGlobal" + ".las");
        ExportBackupTreeFeatureMappingLas(outBackupTreeFeature,GLOBAL_FRAME);

        // std::string outFinalFeatureRef(output_folder + "MapFeatureRef" + ".txt");
        // exportFeatureMappingRef(outFinalFeatureRef, GLOBAL_FRAME);

        std::string outTreeModelGlobal(output_folder + "MapTreeModelGlobal" + ".txt");
        exportMapTreeModel(outTreeModelGlobal, GLOBAL_FRAME);
        std::string outTreeParamGlobal(output_folder + "MapTreeParamGlobal" + ".txt");
        exportMapTreeParam(outTreeParamGlobal,GLOBAL_FRAME);

        std::string outFinalFeatureGlobal(output_folder + "MapFeatureGlobal" + ".las");
        ExportFeatureMappingLas(outFinalFeatureGlobal, GLOBAL_FRAME);
        //std::string outFinalFeatureGlobal(output_folder + "MapFeatureGlobal" + ".txt");
        //exportFeatureMapping(outFinalFeatureGlobal, GLOBAL_FRAME);
        // std::string outFinalFeature_las(output_folder + "MapFeaturelas" + ".las");
        // ExportLasTest(outFinalFeature_las);

    }
    else
    {
        std::string outFinalFeature(output_folder + "MapFeature" + ".las");
        ExportFeatureMappingLas(outFinalFeature);
//         std::string outFinalFeature(output_folder + "MapFeature" + ".txt");
//         exportFeatureMapping(outFinalFeature);
    }
	if (!global_map_flag_)
	{
		UpdateMappingToGlobalTransformation(slam_trajectory, gnss_trajectory);
		std::string outTreeModelGlobal(output_folder + "MapTreeModelGlobal_trans" + ".txt");
		exportMapTreeModel(outTreeModelGlobal, GLOBAL_FRAME);
        std::string outTreeParamGlobal(output_folder + "MapTreeParamGlobal_trans" + ".txt");
        exportMapTreeParam(outTreeParamGlobal,GLOBAL_FRAME);

		//std::string outFinalFeatureGlobal(output_folder + "MapFeatureGlobal_trans" + ".las");
		//ExportFeatureMappingLas(outFinalFeatureGlobal, GLOBAL_FRAME);
		std::string outTrajGlobalCheck(output_folder + "TrajectoryGlobalCheck_trans" + ".txt");
		exportTrajGlobalCheck(outTrajGlobalCheck);
		std::string outTrajGlobal(output_folder + "TrajectoryGlobal_trans" + ".txt");
		exportTrajGlobal(outTrajGlobal);
	}

    // check
    fMapLog << "end" << endl;
}

/**
 * @brief Performs intermediate-level loop closure optimization.
 *
 * This function runs loop closure (LC) on ISCAN level.
 * If configured to run iteratively, it continues until tree merging stabilizes.
 */

void Mapping::ConductIntermediateLC()
{
    fMapLog << endl;
    fMapLog << "*******************************************************" << endl;
    raw_point_flag_ = false;
    int num_iter = 0;
    fMapLog << "Intermediate optimization round, ISCAN level "<<endl;
    if(mPara.iscan_lc_iterative_flag)
        fMapLog<<"Conduct LC iteratively until number of merge trees is smaller than: " << Maximum_Merge_Tree_LC<< ", max iterations: " << Max_Iter_LC << endl;

    //iterative long term optimization
    while (1)
    {
        num_iter++;
        fMapLog << "--------------------------------------------------" << endl;
        fMapLog << "Iteration " << num_iter << endl;
        
        //conduct long term optimization
        LoopClosure();

        //merge tree
        fMapLog << "--------------------------------------------------" << endl;
        TicToc t_merging;
        fMapLog << "Merge Tree" << endl;
        int num_merge_tree = mergeMapTree();
        CountMapTreeType();
        fMapLog << "Time for merging " << t_merging.toc() << endl;
        fMapLog << endl;

        // reset tree parameters for INIT, FITTED and TBD, prevent from fitting in wrong local miminum
        for (int nMapT = 0; nMapT < mpMapTree.size(); nMapT++)
        {
            if (mpMapTree[nMapT]->status == MapTree::FITTED || mpMapTree[nMapT]->status == MapTree::TBD || mpMapTree[nMapT]->status == MapTree::INIT)
            {
                mpMapTree[nMapT]->ResetPara();
            }
        }

        //if not conduct lc iteratively, break
        if(!mPara.iscan_lc_iterative_flag)
        {
            break;
        }
        else
        {  
            //if number of merge tree is few or number of iterations is enough
            if (num_merge_tree <= Maximum_Merge_Tree_LC ||  num_iter >= Max_Iter_LC)
            {
                break;
            }
        }

    }
}

/**
 * @brief Conducts final raw-level loop closure optimization.
 *
 * @param num_iterations Number of optimization iterations to run.
 */

void Mapping::ConductFinalLC(int num_iterations)
{
    fMapLog << endl;
    fMapLog << "*******************************************************" << endl;
    raw_point_flag_ = true;
    int num_iter = 0;
    fMapLog << "Final optimization round, raw level: " << num_iterations << " iterations "<<endl;
    
    //iterative long term optimization
    while (1)
    {
        num_iter++;
        fMapLog << "--------------------------------------------------" << endl;
        fMapLog << "Iteration " << num_iter << endl;

        //conduct long term optimization
        LoopClosure();

        //check residual
        fMapLog << "--------------------------------------------------" << endl;
        fMapLog << "Check map tree with big residual, remove the reference tree flag: ";
        int count_deactivate_reference_tree = 0;
        for (int nMapT = 0; nMapT < mpMapTree.size(); nMapT++)
        {
            if (mpMapTree[nMapT]->status >= 0)
            {
                mpMapTree[nMapT]->rmse = computeTreeToTreeResidual(nMapT, nMapT);
                // if residual is too large, check if caused by wrong match to uav map tree
                if (mpMapTree[nMapT]->rmse > 0.5 && mpMapTree[nMapT]->map_tree_flag)
                {
                    mpMapTree[nMapT]->map_tree_flag = false;
                    count_deactivate_reference_tree++;
                    if (mpMapTree[nMapT]->status >= 3) // established or solid
                    {
                        mpMapTree[nMapT]->ResetPara();
                    }
                }
            }
        }
        fMapLog << count_deactivate_reference_tree << endl;
        CountMapTreeType();

        if (num_iter >= num_iterations)
        {
            break;
        }

        fMapLog << "--------------------------------------------------" << endl;
        TicToc t_merging;
        fMapLog << "Merge Tree" << endl;
        mergeMapTree();
        CountMapTreeType();
        fMapLog << "Time for merging " << t_merging.toc() << endl;
        fMapLog << endl;

        // reset tree parameters for INIT, FITTED and TBD, prevent from fitting in wrong local miminum
        for (int nMapT = 0; nMapT < mpMapTree.size(); nMapT++)
        {
            if (mpMapTree[nMapT]->status == MapTree::FITTED || mpMapTree[nMapT]->status == MapTree::TBD || mpMapTree[nMapT]->status == MapTree::INIT)
            {
                mpMapTree[nMapT]->ResetPara();
            }
        }
    } // end of iterative final LC

    // reset tree parameters for INIT and TBD for final output
    for (int nMapT = 0; nMapT < mpMapTree.size(); nMapT++)
    {
        if (mpMapTree[nMapT]->status == MapTree::INIT || mpMapTree[nMapT]->status == MapTree::TBD)
        {
            mpMapTree[nMapT]->ResetPara();
        }
    }
}

/*********************************************************
 * For original tree points per scan, find the corresponding tree points
 **********************************************************/
/**
 * @brief Assigns tree points (per-scan backups) to closest corresponding map trees.
 *
 * @param distance_threshold Maximum RMSE allowed between tree point and map tree.
 */
void Mapping::FindRawTreePointsPerMapTree(double distance_threshold)
{
    fMapLog << "--------------------------------------------------" << endl;
    fMapLog << "Finding corresponding maptree for backup tree points, distance threshold: " <<  distance_threshold<<endl;
    int num_total_points = 0, num_found_points = 0;
    for (int nIscan = 0; nIscan < mvpIScans.size(); nIscan++)
    {
        int iscanId = nIscan;
        //reset match tree information
        vector<vector<int>>().swap( mvpIScans[iscanId]->matched_map_tree_id);
        for (int nScan = 0; nScan < mvpIScans[iscanId]->raw_tree_points_mapping.size(); nScan++) //individual scan
        {
            vector<int> match_tree_id;
            for (int nTree = 0; nTree < mvpIScans[iscanId]->raw_tree_points_mapping[nScan].size(); nTree++) //individual tree
            {
                //for each individual tree from a scan, find the corresponding map tree
                int treeId = nTree;
                int best_tree_status = -1, best_map_tree_id = -1;
                double best_score = 100.0;
                for(int nMapT = 0; nMapT< mpMapTree.size(); nMapT++)
                {
                    //only consider tree >= init (0-4)
                    if(mpMapTree[nMapT]->status < 0) 
                    {
                        continue;
                    }
                    
                    // for each tree point
                    int num_point = 0;
                    double dis_rmse = 0.0;
                    bool stop_flag = false;
                    for (int nP = 0; nP < mvpIScans[iscanId]->raw_tree_points_mapping[nScan][treeId].size(); nP++)
                    {
                        // points in mapping
                        Eigen::Vector3d r_I_m(mvpIScans[iscanId]->raw_tree_points_mapping[nScan][treeId][nP].x,
                                            mvpIScans[iscanId]->raw_tree_points_mapping[nScan][treeId][nP].y,
                                            mvpIScans[iscanId]->raw_tree_points_mapping[nScan][treeId][nP].z);
                        double dis = mpMapTree[nMapT]->compute_dis(r_I_m);
                        
                        //if already too far from the tree model, stop
                        if(dis > 5.0)
                        {
                            stop_flag = true;
                            break;
                        }
                        dis_rmse = dis_rmse + dis * dis;
                        num_point++;
                    }

                    if(stop_flag)
                    {
                        continue;
                    }
                    dis_rmse = sqrt(dis_rmse / double(num_point));

                    //check rmse and tree status
                    if(dis_rmse < distance_threshold)
                    {
                        if(mpMapTree[nMapT]->status > best_tree_status || 
                        (mpMapTree[nMapT]->status == best_tree_status && dis_rmse < best_score))
                        {
                            best_tree_status = mpMapTree[nMapT]->status;
                            best_map_tree_id = nMapT;
                            best_score = dis_rmse;
                        }
                    }
                }
                match_tree_id.push_back(best_map_tree_id);
                num_total_points += mvpIScans[iscanId]->raw_tree_points_mapping[nScan][treeId].size();
                if(best_map_tree_id!= -1)
                {
                    num_found_points+= mvpIScans[iscanId]->raw_tree_points_mapping[nScan][treeId].size();
                }
            }
            mvpIScans[iscanId]->matched_map_tree_id.push_back(match_tree_id);
        }
    }
    fMapLog << "\tTotal number of points and the ones that find map trees: " << num_total_points << " - " << num_found_points<<endl;

}

/*********************************************************
 * For non ground raw points per scan, find the corresponding tree points
 **********************************************************/
/**
 * @brief Assigns non-ground points to the nearest valid map tree.
 *
 * @param distance_threshold Distance threshold to associate raw points with trees.
 */
void Mapping::FindRawPointsPerMapTree(double distance_threshold)
{
    fMapLog << "--------------------------------------------------" << endl;
    fMapLog << "Finding corresponding maptree for backup raw points, distance threshold: " <<  distance_threshold<<endl;
    int num_total_points = 0, num_found_points = 0;

    //clearn the raw point vector

    for (int nIscan = 0; nIscan < mvpIScans.size(); nIscan++)
    {
        int iscanId = nIscan;

        vector<vector<PointType>> tree_points_per_map_tree;
        vector<int> map_tree_ids;
        int num_tree_points = 0;
        int num_all_points = 0;
        for (int nScan = 0; nScan < mvpIScans[iscanId]->raw_non_ground_points.size(); nScan++) // individual scan
        {
            for (int nP = 0; nP < mvpIScans[iscanId]->raw_non_ground_points_mapping[nScan].size(); nP++)
            {
                num_all_points++;
                Eigen::Vector3d r_I_m(mvpIScans[iscanId]->raw_non_ground_points_mapping[nScan][nP].x,
                                      mvpIScans[iscanId]->raw_non_ground_points_mapping[nScan][nP].y,
                                      mvpIScans[iscanId]->raw_non_ground_points_mapping[nScan][nP].z);
                
                //for each individual point from a scan, find the corresponding map tree
                int best_tree_status = -1, best_map_tree_id = -1;
                double best_score = 100.0;
                for (int nMapT = 0; nMapT < mpMapTree.size(); nMapT++)
                {
                    // only consider tree >= init (0-4)
                    if (mpMapTree[nMapT]->status < 0)
                    {
                        continue;
                    }

                    if (r_I_m(2) < mpMapTree[nMapT]->ground_point(2) + MaxHeightTreePoint && r_I_m(2) > mpMapTree[nMapT]->ground_point(2) + MinHeightTreePoint)
                    {
                        double dis = mpMapTree[nMapT]->compute_dis(r_I_m); //if negative, its inside the cylinder, accept
                        //check rmse and tree status
                        if(dis < distance_threshold)
                        {
                            if(mpMapTree[nMapT]->status > best_tree_status || 
                            (mpMapTree[nMapT]->status == best_tree_status && abs(dis) < best_score))
                            {
                                best_tree_status = mpMapTree[nMapT]->status;
                                best_map_tree_id = nMapT;
                                best_score = abs(dis);
                            }
                        }
                    }
                }

                //find the corresponding tree
                num_total_points ++;
                if(best_map_tree_id == -1)
                {
                    continue;
                }

                num_found_points++;
                num_tree_points++;
                int index = -1;
                for (int i = 0; i < map_tree_ids.size(); i++)
                {
                    if (best_map_tree_id == map_tree_ids[i])
                    {
                        index = i;
                        break;
                    }
                }
                // if not, create a new one
                if (index == -1)
                {
                    index = map_tree_ids.size();
                    map_tree_ids.push_back(best_map_tree_id);
                    tree_points_per_map_tree.push_back(vector<PointType>());
                }

                // add the points in raw_tree_points
                PointType p = mvpIScans[iscanId]->raw_non_ground_points[nScan][nP];
                p.intensity += float(100 * nScan);
                tree_points_per_map_tree[index].push_back(p);
                
            }
        }

        mvpIScans[iscanId]->vTreePointRaw = tree_points_per_map_tree;
        mvpIScans[iscanId]->corresponding_map_tree_ids = map_tree_ids;
        mvpIScans[iscanId]->computeRawTreePointToMapping();
        fDebug<<"Iscan " << mvpIScans[iscanId]->index << ": " << num_all_points << " -> " <<  num_tree_points <<" points, " << map_tree_ids.size() << "trees" <<endl;
    }
    fMapLog << "\tTotal number of points and the ones that find map trees: " << num_total_points << " - " << num_found_points<<endl;
}

/**
 * @brief Updates visibility and point associations for all map trees.
 *
 * This clears previous associations and repopulates them based on current observations.
 */

void Mapping::Update_Map_Tree_Obs()
{
    //remove observation for each tree with points (status >= 0)
    for (int nMapT = 0; nMapT < mpMapTree.size(); nMapT++)
    {
        if (mpMapTree[nMapT]->status >= 0)
        {
            mpMapTree[nMapT]->clean_obs();
        }
    }

    for (int nIscan = 0; nIscan < mvpIScans.size(); nIscan++)
    {    
        int iscanId = nIscan;
        for (int nTree = 0; nTree < mvpIScans[iscanId]->corresponding_map_tree_ids.size(); nTree++)
        {
            mpMapTree[mvpIScans[iscanId]->corresponding_map_tree_ids[nTree]]->visibleInfo.push_back(make_pair(nIscan, nTree));
        }
    }
    for (int nMapT = 0; nMapT < mpMapTree.size(); nMapT++)
    {
        int mapTreeId = nMapT;
        if (mpMapTree[mapTreeId]->status >= 0)
        {
            int numPoint = 0;
            for (int nT = 0; nT < mpMapTree[mapTreeId]->visibleInfo.size(); nT++)
            {
                int iscanId = mpMapTree[mapTreeId]->visibleInfo[nT].first;
                int treeId = mpMapTree[mapTreeId]->visibleInfo[nT].second;
                numPoint += mvpIScans[iscanId]->vTreePointRaw[treeId].size();
            }
            fDebug<< "Map Tree " << mapTreeId << ". Points from " <<mpMapTree[mapTreeId]->numPoint << " -> " << numPoint << ". Status from " << mpMapTree[mapTreeId]->status;
            mpMapTree[mapTreeId]->numPoint = numPoint;
            mpMapTree[mapTreeId]->update_status_new_obs();
            fDebug << " -> " << mpMapTree[mapTreeId]->status <<endl;
        }
    }
    // //fMapLog << "Finding corresponding maptree for backup tree points" <<endl;
    // for (int nIscan = 0; nIscan < mvpIScans.size(); nIscan++)
    // {
    //     int iscanId = nIscan;
    //     //reset match tree information
    //     vector<vector<int>>().swap( mvpIScans[iscanId]->matched_map_tree_id);
    //     for (int nScan = 0; nScan < mvpIScans[iscanId]->raw_tree_points_mapping.size(); nScan++) //individual scan
    //     {
    //         vector<int> match_tree_id;
    //         for (int nTree = 0; nTree < mvpIScans[iscanId]->raw_tree_points_mapping[nScan].size(); nTree++) //individual tree
    //         {
    //             //for each individual tree from a scan, find the corresponding map tree
    //             int treeId = nTree;
    //             int best_tree_status = -1, best_map_tree_id = -1;
    //             double best_score = 100.0;
    //             for(int nMapT = 0; nMapT< mpMapTree.size(); nMapT++)
    //             {
    //                 //only consider tree >= init (0-4)
    //                 if(mpMapTree[nMapT]->status < 0) 
    //                 {
    //                     continue;
    //                 }
                    
    //                 // for each tree point
    //                 int num_point = 0;
    //                 double dis_rmse = 0.0;
    //                 bool stop_flag = false;
    //                 for (int nP = 0; nP < mvpIScans[iscanId]->raw_tree_points_mapping[nScan][treeId].size(); nP++)
    //                 {
    //                     // points in mapping
    //                     Eigen::Vector3d r_I_m(mvpIScans[iscanId]->raw_tree_points_mapping[nScan][treeId][nP].x,
    //                                         mvpIScans[iscanId]->raw_tree_points_mapping[nScan][treeId][nP].y,
    //                                         mvpIScans[iscanId]->raw_tree_points_mapping[nScan][treeId][nP].z);
    //                     double dis = mpMapTree[nMapT]->compute_dis(r_I_m);
                        
    //                     //if already too far from the tree model, stop
    //                     if(dis > 5.0)
    //                     {
    //                         stop_flag = true;
    //                         break;
    //                     }
    //                     dis_rmse = dis_rmse + dis * dis;
    //                     num_point++;
    //                 }

    //                 if(stop_flag)
    //                 {
    //                     continue;
    //                 }
    //                 dis_rmse = sqrt(dis_rmse / double(num_point));

    //                 //check rmse and tree status
    //                 if(dis_rmse < 0.3)
    //                 {
    //                     if(mpMapTree[nMapT]->status > best_tree_status || 
    //                     (mpMapTree[nMapT]->status == best_tree_status && dis_rmse < best_score))
    //                     {
    //                         best_tree_status = mpMapTree[nMapT]->status;
    //                         best_map_tree_id = nMapT;
    //                         best_score = dis_rmse;
    //                     }
    //                 }
    //             }
    //             match_tree_id.push_back(best_map_tree_id);
    //         }
    //         mvpIScans[iscanId]->matched_map_tree_id.push_back(match_tree_id);
    //     }
    // }

}

/********************************************************
Main function to integrate different scans to one integral (in a local frame)
Input:
    mvTempScan
Ouput:
    pCurrentIScan
********************************************************/
/**
 * @brief Main function to integrate multiple scans into a single integrated scan in local frame.
 * 
 * @details This function integrates scans stored in mvTempScan. It performs the following steps:
 * - Derives local frame pose from odometry or trajectory
 * - Backs up tree points per scan
 * - Integrates tree and ground points into local map
 * - Optionally performs optimization on integrated data
 * - Updates trajectory and features
 * - Broadcasts final integrated data
 */
void Mapping::integrateScan()
{
    //Flag_Scan_Inter = true;
    Flag_Scan_Inter = mPara.iscan_opt_flag;

    mNumIntegrationScan = mvTempScan.size();
    fMapLog << "*******************************************************" << endl;
    fMapLog <<  mnCount << "\t" << "Scan index: " << mvTempScan[0].scanID << " - " << mvTempScan.back().scanID << ", Number of scan: " << mNumIntegrationScan << endl;
    cout <<endl;
    cout << "*******************************************************" << endl;
    cout << "Scan index: " << mvTempScan[0].scanID << " - " << mvTempScan.back().scanID << ", Number of scan: " << mNumIntegrationScan << endl;

    fMapLog << "--------------------------------------------------" << endl;
    fMapLog << "1. Scan integration" << endl;

    //from the odometry, derive pose in local frame
    derivePoseLocal();

    //backup tree Points
    backupTreePoints();

    // integrate tree from these scans
    integrateTree();

    // extract planar patches from these scans
    integrateGroundPoints();

    if (intermediate_result_flag)
    {
        // Features in Optimization
        std::string outFeatureOptIni(output_folder_iscan + to_string(mnCount) + "_IscanOptFeatureInit_" + to_string(mvTempScan[0].scanID) + ".txt");
        exportIntegratedOptFeatures(outFeatureOptIni, 0, 1);
    }

    // next step. optimize
    if(Flag_Scan_Inter)
    {
        optimizeIntegratedScan();
    }

    //if fail, fix trajectory and redo 
    if(!Flag_Scan_Inter)
    {
        fMapLog << "Fix traj and update feature parameters" <<endl;
        optimizeIntegratedScan(true);
    }

    // intermediate results for paper

    //1. Trajectory 
    std::string outTrajectory(output_folder_iscan + to_string(mnCount) + "_IscanTraj_" + to_string(mvTempScan[0].scanID) + ".txt");
    exportIntegratedTrajectory(outTrajectory);

    if (intermediate_result_flag)
    {
        // //2. Features in Optimization
        // std::string outFeatureOptIni(output_folder + "IscanOptFeatureInit_" + to_string(mvTempScan[0].scanID) + ".txt");
        // exportIntegratedOptFeatures(outFeatureOptIni,0,1);
        std::string outFeatureOptRefined(output_folder_iscan +to_string(mnCount) + "_IscanOptFeatureRefine_" + to_string(mvTempScan[0].scanID) + ".txt");
        exportIntegratedOptFeatures(outFeatureOptRefined, 1, 1);

        // 3. All feature points
        std::string outFeatureIni(output_folder_iscan + to_string(mnCount) + "_IscanFeatureInit_" + to_string(mvTempScan[0].scanID) + ".txt");
        exportIntegratedFeatures(outFeatureIni, 0, 1);
        std::string outFeatureRefined(output_folder_iscan + to_string(mnCount) + "_IscanFeatureRefine_" + to_string(mvTempScan[0].scanID) + ".txt");
        exportIntegratedFeatures(outFeatureRefined, 1, 1);
    }

    //after optimization
    UpdateMergeTree();

    if (intermediate_result_flag)
    {
        std::string outFeatureOptRefinedMerged(output_folder_iscan + to_string(mnCount) + "_IscanOptFeatureRefineMerged_" + to_string(mvTempScan[0].scanID) + ".txt");
        exportIntegratedOptFeatures(outFeatureOptRefinedMerged, 1, 1);
    }

    // broadcast to map
    broadcastIScan();

    // export broadcast points
    if (intermediate_result_flag)
    {
        std::string outFeatureIscanDs(output_folder_iscan + to_string(mnCount) + "_IscanOptFeatureBroadcastDs_" + to_string(mvTempScan[0].scanID) + ".txt");
        exportIscanFeatureBroadcast(outFeatureIscanDs, 1);
        std::string outFeatureIscan(output_folder_iscan + to_string(mnCount) + "_IscanOptFeatureBroadcast_" + to_string(mvTempScan[0].scanID) + ".txt");
        exportIscanFeatureBroadcast(outFeatureIscan, 0);
    }

    // clear all used info for current integration
    resetIntegration();
    
    //system("read -p 'Press Enter to continue...' var");

}

/********************************************************
Function to derive pose information in the local frame
Input:
    mvTempScan
Ouput:
    vec_R/r_ini
    vec_R/r_ref
********************************************************/
/**
 * @brief Derive relative poses of each scan with respect to the local frame.
 * 
 * @details Sets pose of the first scan to origin and computes transformation for subsequent scans
 * using odometry or trajectory. Updates r/R_lu_local for each scan.
 */
void Mapping::derivePoseLocal()
{
     /*Definition of the local frame:
    The postion/orientation of the first scan are set to origin, remaining scans are integrated to the starting one.
    To ensure that point cloud in mapping is leveled, use nominal leveling rotation or the calibrated*/
    
    //@@@TEST
    //Eigen::Matrix3d initR = mpTraj->R_lu_b;mvTempScan[0].R_lu_lup;
    //Eigen::Matrix3d initR = computeAverageLevelR().toRotationMatrix();
    Eigen::Matrix3d init_R = mvTempScan[0].R_lu_lup;
    Eigen::Vector3d init_r = Eigen::Vector3d::Zero();

    //derive lu to local for each pose
    Eigen::Vector3d r_lu_local_t1 = init_r;
    Eigen::Matrix3d R_lu_local_t1 = init_R;
    mR_lu_local_t0 = R_lu_local_t1;
    mr_lu_local_t0 = r_lu_local_t1; // the pose at the end of previous integrated scan

    vector<Eigen::Vector3d> vec_r_ini;
    vector<Eigen::Matrix3d> vec_R_ini;
    vec_r_ini.push_back(r_lu_local_t1);
    vec_R_ini.push_back(R_lu_local_t1);

    // go through remaining scans 
    for (int nScan = 0; nScan < mNumIntegrationScan; nScan++)
    {
        // current scan to the first scan from odometry thread
        Eigen::Matrix3d R_lu_local_t2 = R_lu_local_t1 * mvTempScan[nScan].R_lut2_lut1;
        Eigen::Vector3d r_lu_local_t2 = r_lu_local_t1 + R_lu_local_t1 * mvTempScan[nScan].r_lut2_lut1;

        mvTempScan[nScan].r_lu_local_ini = r_lu_local_t2;
        mvTempScan[nScan].R_lu_local_ini = R_lu_local_t2;
        mvTempScan[nScan].r_lu_local_updated = r_lu_local_t2;
        mvTempScan[nScan].R_lu_local_updated = R_lu_local_t2;
        // used for next scan
        r_lu_local_t1 = r_lu_local_t2;
        R_lu_local_t1 = R_lu_local_t2;

        vec_r_ini.push_back(r_lu_local_t2);
        vec_R_ini.push_back(R_lu_local_t2);
    }
    mv_R_lu_local_ini = vec_R_ini;
    mv_r_lu_local_ini = vec_r_ini;

    //derive lu to local using reference trajectory 
    if (mpTraj)
    {
        //same initial for the reference traj
        Eigen::Vector3d r_lu_local_t1_ref = init_r;
        Eigen::Matrix3d R_lu_local_t1_ref = init_R;

        vector<Eigen::Vector3d> vec_r_ref;
        vector<Eigen::Matrix3d> vec_R_ref;
        vec_r_ref.push_back(r_lu_local_t1_ref);
        vec_R_ref.push_back(R_lu_local_t1_ref);

        Eigen::Vector3d r_lu_b_ref = mpTraj->r_lu_b;
        Eigen::Matrix3d R_lu_b_ref = mpTraj->R_lu_b;
        for (int nScan = 0; nScan < mNumIntegrationScan; nScan++)
        {

            Eigen::Vector3d r_b_m_t1, r_b_m_t2, r_lu_m_t1, r_lu_m_t2;
            Eigen::Matrix3d R_b_m_t1, R_b_m_t2, R_lu_m_t1, R_lu_m_t2;

            int idx;
            if (!(mpTraj->bopInterpolation((mvTempScan[nScan].mTimeEnd - mvTempScan[nScan].mTimeTracking) / 1000.0, -1, idx, r_b_m_t1, R_b_m_t1) && mpTraj->bopInterpolation(mvTempScan[nScan].mTimeEnd / 1000.0, -1, idx, r_b_m_t2, R_b_m_t2)))
                throw std::runtime_error("Wrong Time");

            r_lu_m_t1 = r_b_m_t1 + R_b_m_t1 * r_lu_b_ref;
            R_lu_m_t1 = R_b_m_t1 * R_lu_b_ref;

            r_lu_m_t2 = r_b_m_t2 + R_b_m_t2 * r_lu_b_ref;
            R_lu_m_t2 = R_b_m_t2 * R_lu_b_ref;

            Eigen::Matrix3d R_lut2_lut1_ref = R_lu_m_t1.inverse() * R_lu_m_t2;
            Eigen::Vector3d r_lut2_lut1_ref = R_lu_m_t1.inverse() * (r_lu_m_t2 - r_lu_m_t1);

            // current scan to the first scan from odometry thread
            Eigen::Matrix3d R_lu_local_t2_ref = R_lu_local_t1_ref * R_lut2_lut1_ref;
            Eigen::Vector3d r_lu_local_t2_ref = r_lu_local_t1_ref + R_lu_local_t1_ref * r_lut2_lut1_ref;

            mvTempScan[nScan].r_lu_local_ref = r_lu_local_t2_ref;
            mvTempScan[nScan].R_lu_local_ref = R_lu_local_t2_ref;
            // used for next scan
            r_lu_local_t1_ref = r_lu_local_t2_ref;
            R_lu_local_t1_ref = R_lu_local_t2_ref;

            vec_r_ref.push_back(r_lu_local_t2_ref);
            vec_R_ref.push_back(R_lu_local_t2_ref);
        }
        mv_R_lu_local_ref = vec_R_ref;
        mv_r_lu_local_ref = vec_r_ref;
    }
}


/*
Save tree points from each scan in the height range
*/
/**
 * @brief Backup tree points from each scan for later processing.
 * 
 * @details Stores 3D tree points per tree per scan while applying optional height filters.
 */
void Mapping::backupTreePoints()
{
    for (int nScan = 0; nScan < mNumIntegrationScan; nScan++)
    {
        int count = 0;
        vector<vector<PointType>> tree_pts_per_scan; 
        for (int nTree = 0; nTree < mvTempScan[nScan].vTreeRawPoints.size(); nTree++)
        {
            vector<PointType> pts_per_tree;
            vector<PointType> *pPointList = &(mvTempScan[nScan].vTreeRawPoints[nTree]);
            Eigen::Matrix3f R_temp = mvTempScan[nScan].R_lu_lup.cast<float>();
            int numPoint = pPointList->size();
            for (int nP = 0; nP < numPoint; nP++)
            {
                PointType p = pPointList->at(nP);
                Eigen::Vector3f r_I_lu(p.x, p.y, p.z);

                // Compute the coordinates in leveled lidar unit frame (the Z will be the height information)
                Eigen::Vector3f rIlup = R_temp * r_I_lu;
                // height range info
                //if (rIlup(2) < mvTempScan[nScan].gH + MaxHeightTreePoint && rIlup(2) > mvTempScan[nScan].gH + MinHeightTreePoint)
                //(TODO:keep all features because some area is not leveled)
                // if (rIlup(2) < mvTempScan[nScan].gH + MaxHeightTreePoint)
                {
                    pts_per_tree.push_back(p);
                    count++;
                }
            }
            tree_pts_per_scan.push_back(pts_per_tree);
        }
        raw_tree_points_per_iscan.push_back(tree_pts_per_scan);
        fDebug << mnCount << "\t" << "Scan index: " << mvTempScan[nScan].scanID  << ", tree points in range: " << count <<endl;
    }
}

/* function to integrate tree from the scans*/
/*majar steps
1. match tree from different scans
2. get the portion based on [MaxHeightTreePoint, MinHeightTreePoint], then compute the center and each tree and remove points that are 1 m away from center, 
    for tree with more than 10 points, keep it -> "tree with points in range: "
3. check number of points for each tree, if its larger than MinNumTreePoints * numScans, save it
*/
/**
 * @brief Integrate tree features from all scans.
 * 
 * @details Performs tree matching across scans, removes noisy trees, computes tree centers,
 * and updates valid trees based on point density.
 */
void Mapping::integrateTree()
{
#if 0
    // fill trees of first scan
    vector<int> treeCount(1000, 0); // how many scan has this tree
    for (int nT = 0; nT < mvTempScan[0].vTreeLoc.size(); nT++)
    {
        Eigen::Quaterniond q(mvTempScan[0].R_lu_local_ini);
        Eigen::Vector3d t = mvTempScan[0].r_lu_local_ini;

        mvIntegratedTree.push_back(q * mvTempScan[0].vTreeLoc[nT] + t);
        vector<pair<int, int>> treeInfo{make_pair(0, nT)};
        mvIntegratedTreeIds.push_back(treeInfo);
        treeCount[mvIntegratedTree.size() - 1] += 1;
    }

    // fMapLog << "Number of ini tree: " << mvIntegratedTree.size() << endl;
    int numIniTree = mvIntegratedTree.size();

    // go through remaining scans
    for (int nScan = 1; nScan < mNumIntegrationScan; nScan++)
    {
        int numTree = mvTempScan[nScan].vTreeLoc.size();
        // fMapLog << "---------------- Scan " << mvTempScan[nScan].scanID << " with " << numTree << " trees" << endl;

        Eigen::Quaterniond q(mvTempScan[nScan].R_lu_local_ini);
        Eigen::Vector3d t = mvTempScan[nScan].r_lu_local_ini;

        // predict the location of trees of scan
        vector<bool> vbMatch(numTree, false);
        vector<double> vecDistanceThreshold;
        vector<Eigen::Vector3d> vecEstimateCenter;
        for (int nT = 0; nT < numTree; nT++)
        {
            vecDistanceThreshold.push_back(max(2.0, 0.05 * mvTempScan[nScan].vTreeLoc[nT].norm()));
            vecEstimateCenter.push_back(q * mvTempScan[nScan].vTreeLoc[nT] + t);
        }

        // match tree
        std::vector<pair<int, int>> pairs;
        matchTreeScantoMap(&mvIntegratedTree, &vecEstimateCenter, &vecDistanceThreshold, pairs);
        int numPair = pairs.size();

        // iteratively map cur tree to map tree. Currently deactivate.
        vector<bool> validTreePair(numPair, true);

        // save mapped trees to the mvIntegratedTreeIds.
        int validPair = 0;
        for (int nPair = 0; nPair < numPair; nPair++)
        {
            int mapTreeId = pairs[nPair].first;
            int curTreeId = pairs[nPair].second;
            vbMatch[curTreeId] = true;
            mvIntegratedTreeIds[mapTreeId].push_back(make_pair(nScan, curTreeId));
            treeCount[mapTreeId] += 1;
            validPair++;
        }

        // save unmatched trees as new tree
        for (int nT = 0; nT < numTree; nT++)
        {
            if (vbMatch[nT])
                continue;

            mvIntegratedTree.push_back(q * mvTempScan[nScan].vTreeLoc[nT] + t);
            vector<pair<int, int>> treeInfo{make_pair(nScan, nT)};
            mvIntegratedTreeIds.push_back(treeInfo);
            treeCount[mvIntegratedTree.size() - 1] += 1;
        }

        // fMapLog << validPair << " pairs, MapTree: " << numIniTree << " -> " << mvIntegratedTree.size() << endl;

        numIniTree = mvIntegratedTree.size();
    } // end of go through scans

    // remove invalid trees
    //fMapLog << "--------------------------------------------------" << endl;
    fMapLog << "Number of ini tree: " << numIniTree;
#endif
	std::cout << "running_iscan_index: " << running_iscan_index << std::endl;
    std::cout << "Index: " << running_iscan_index
          << " — Size of mvIntegratedTreeIdsVector: " << mvIntegratedTreeIdsVector.size()
          << " — Size of mvIntegratedTreeIdsVector[running_iscan_index]: "
          << mvIntegratedTreeIdsVector[running_iscan_index].size() << std::endl;
	mvIntegratedTreeIds = mvIntegratedTreeIdsVector[running_iscan_index];
	mvIntegratedTree = mvIntegratedTreeVector[running_iscan_index];
    running_iscan_index++;
	int numIniTree = mvIntegratedTree.size();
    vector<Eigen::Vector3d> vIntegratedTree;           // locations of the tree (vector3d)
    vector<vector<pair<int, int>>> vIntegratedTreeIds; // scan id + tree id of each tree

    int validTree = 0;
    for (int nTree = 0; nTree < numIniTree; nTree++)
    {
        Eigen::Vector3f center_all = Eigen::Vector3f::Zero();
        vector<PointType> vPointRawPerTree;
        vector<Eigen::Vector3f> vPointLocalPerTree; // in local frame defined by the m  of first scan

        for (int nScan = 0; nScan < mvIntegratedTreeIds[nTree].size(); nScan++)
        {
            int scanId = mvIntegratedTreeIds[nTree][nScan].first;
            int treeId = mvIntegratedTreeIds[nTree][nScan].second;
            vector<PointType> *pPointList = &(mvTempScan[scanId].vTreeRawPoints[treeId]);
            int numPoint = pPointList->size();
            for (int nP = 0; nP < numPoint; nP++)
            {
                PointType p = pPointList->at(nP);
                Eigen::Vector3f r_I_lu(p.x, p.y, p.z);

                //compute the coordinates in the local frame (assuming no distortion)
                Eigen::Vector3f r_I_local = mvTempScan[scanId].R_lu_local_ini.cast<float>() * r_I_lu + mvTempScan[scanId].r_lu_local_ini.cast<float>();

                //Compute the coordinates in leveled lidar unit frame (the Z will be the height information)
                Eigen::Vector3f rilup = mvTempScan[scanId].R_lu_lup.cast<float>() * r_I_lu;

                // height range info
                //test
				//if (rilup(2) < mvTempScan[scanId].gH + MaxHeightTreePoint)
                //(TODO:keep all features because some area is not leveled)
                // if (rilup(2) < mvTempScan[scanId].gH + MaxHeightTreePoint && rilup(2) > mvTempScan[scanId].gH + MinHeightTreePoint)
                {
                    center_all += r_I_local;
                    p.intensity += float(100 * scanId); 
                    vPointRawPerTree.push_back(p);
                    vPointLocalPerTree.push_back(r_I_local);
                }
            }
        }

        int numPointAll = vPointRawPerTree.size();
        center_all = center_all / float(numPointAll);

        Eigen::Vector3f center = Eigen::Vector3f::Zero();
        vector<PointType> vPointRaw;
        //vector<Eigen::Vector3f> vecPointLocal;

        for (int i = 0; i < numPointAll; i++)
        {
            float dis_sqr = (vPointLocalPerTree[i](0) - center_all(0)) * (vPointLocalPerTree[i](0) - center_all(0)) + (vPointLocalPerTree[i](1) - center_all(1)) * (vPointLocalPerTree[i](1) - center_all(1));
            if (dis_sqr < 1.0 * 1.0)
            {
                vPointRaw.push_back(vPointRawPerTree[i]);
                center += vPointLocalPerTree[i];
            //    vecPointLocal.push_back(vPointLocalPerTree[i]);
            }
        }
        int numPoint = vPointRaw.size();
        center = center / float(numPoint);
		//(TODO:CHUNXI)
        //if (numPoint > 2) // remove the one with 0 observations
        {
            CylinderPara treePara;
            treePara.x = center.cast<double>();
            treePara.n << 0.0, 0.0, 1.0;
            treePara.r = 0.0;

            mvTreeCentersLocal.push_back(center);
            mvTreePointsRaw.push_back(vPointRaw);
            mvTreeParamsLocal.push_back(treePara);

            vIntegratedTree.push_back(mvIntegratedTree[nTree]);       // locations of the tree (vector3d)
            vIntegratedTreeIds.push_back(mvIntegratedTreeIds[nTree]); // scan id + tree id of each tree

#ifdef MAPTESTINVALID
            fMapLog << nTree << "\t" << numPoint << "\t" << treePara.x.transpose() << "\t" << treePara.n.transpose() << "\t" << treePara.r << endl;
#endif

            //mvValidTrees.push_back(true);

            // // enough points and inlier points
            // if (numPoint > MinNumTreePoints * double(mNumIntegrationScan))
            // //            if (float(numPoint) / float(numPointAll) > 0.8 && numPoint > 100)
            // {
            //     mvValidTrees.push_back(true);
            //     validTree++;
            // }
            // else
            //     mvValidTrees.push_back(false);
        }
    }
    
    //mNumValidTree = validTree;
#ifdef EXPORT_LOG
	fLog << "filtered vIntegratedTree size: " << vIntegratedTree.size() << std::endl;
	fLog << "mvIntegratedTree size: " << mvIntegratedTree.size() << std::endl;
#endif

    mvIntegratedTree = vIntegratedTree;
    mvIntegratedTreeIds = vIntegratedTreeIds;
    numIniTree = vIntegratedTree.size();
    // for (int nTree = 0; nTree < numIniTree; nTree++)
    // {
    //     fMapLog << mvTreeCenters[nTree].transpose() << "\t" << mvValidTrees[nTree] << "\t" << mvTreePoints[nTree].size() << endl;
    // }

    fMapLog << ", tree with points in the height range: " << numIniTree;

    mvValidTrees.clear();
    mNumValidTree = 0;
    //(TODO)
    for (int nTree = 0; nTree < mvIntegratedTree.size(); nTree++)
    {
        if (mvTreePointsRaw[nTree].size() > MinNumTreePointPerScan * double(mNumIntegrationScan))//25
        {
            mvValidTrees.push_back(true);
            mNumValidTree++;
        }
        else
        {
            mvValidTrees.push_back(false);
        }
    }

    fMapLog << ", tree with sufficient number of points: " << mNumValidTree << endl;
}

/* function to extract plane patches from the scans*/
/**
 * @brief Integrate planar ground features from all scans.
 * 
 * @details Transforms ground points to local frame, assigns them to a voxel grid, fits planar patches,
 * and filters out weak planes. Updates mvPlaneParamLocal.
 */
void Mapping::integrateGroundPoints()
{
    // compute ground points in lup and mapping
    double minX = 100.0, maxX = -100.0, minY = 100.0, maxY = -100.0;
    for (int nScan = 0; nScan < mNumIntegrationScan; nScan++)
    {
        int scanId = nScan;
        for (std::size_t nP = 0; nP < mvTempScan[nScan].vGroundPlanarRawPoints.size(); nP++)
        {
            PointType p = mvTempScan[nScan].vGroundPlanarRawPoints[nP];

            Eigen::Vector3f rilu(p.x, p.y, p.z);
            p.intensity += float(100 * scanId);
            mvGroundPointsRaw.push_back(p);

            Eigen::Vector3d r_I_local;
            r_I_local = mvTempScan[scanId].R_lu_local_ini * rilu.cast<double>() + mvTempScan[scanId].r_lu_local_ini;
            PointType pm = p;
            pm.x = float(r_I_local(0));
            pm.y = float(r_I_local(1));
            pm.z = float(r_I_local(2));
            mvGroundPointsLocal.push_back(pm);

            minX = min(minX, r_I_local(0));
            minY = min(minY, r_I_local(1));
            maxX = max(maxX, r_I_local(0));
            maxY = max(maxY, r_I_local(1));
        }
    }

    // partition to grid.. 2* (tree+grid) <= res
    double res = KGround_resolution;
    double tree_radius = 0.5; // 1m
    double grid_radius = 0.5;
    minX = res * floor(minX / res);
    minY = res * floor(minY / res);
    maxX = res * ceil(maxX / res);
    maxY = res * ceil(maxY / res);
    int col = round((maxX - minX) / res);
    int row = round((maxY - minY) / res);

    // compute center of each grid
    vector<pair<double, double>> vGridCenters;
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            vGridCenters.push_back(make_pair(double(j + 0.5) * res + minX, double(i + 0.5) * res + minY));
        }
    }

    // compare the tree center vs the grid center
    // the center: [x_min + grid_radius, x_max - grid_radius]
    // if x_tree > half. x_tree - tree_radius
    for (auto c : mvTreeCentersLocal)
    {
        double x_tree = double(c(0));
        double y_tree = double(c(1));
        int i = floor((y_tree - minY) / res);
        int j = floor((x_tree - minX) / res);
        if (i < 0 || i >= row || j < 0 || j >= col)
            continue;
        int index = i * col + j;

        // fMapLog << mvTreeCentersLocal.size() << "\t" << i << "\t" << j << "\t" << index << endl;

        // x coordinates. too close
        if (x_tree >= vGridCenters[index].first)
        {
            vGridCenters[index].first = max(x_tree - tree_radius - grid_radius, minX + res * double(j) + grid_radius);
        }
        else
        {
            vGridCenters[index].first = min(x_tree + tree_radius + grid_radius, minX + res * double(j + 1) - grid_radius);
        }

        // y coordinates. too close
        if (y_tree >= vGridCenters[index].second)
        {
            vGridCenters[index].second = max(y_tree - tree_radius - grid_radius, minY + res * double(i) + grid_radius);
        }
        else
        {
            vGridCenters[index].second = min(y_tree + tree_radius + grid_radius, minY + res * double(i + 1) - grid_radius);
        }
    }
    mvGridCenters = vGridCenters;

    // find points belong to each grid
    vector<int> temp;
    vector<vector<int>> vGridPointIndex(vGridCenters.size(), temp);
    for (int nP = 0; nP < mvGroundPointsLocal.size(); nP++)
    {
        double x_p = double(mvGroundPointsLocal[nP].x);
        double y_p = double(mvGroundPointsLocal[nP].y);
        int i = floor((y_p - minY) / res);
        int j = floor((x_p - minX) / res);
        int index = i * col + j;

        if (abs(x_p - vGridCenters[index].first) < grid_radius && abs(y_p - vGridCenters[index].second) < grid_radius)
        {
            vGridPointIndex[index].push_back(nP);
        }
    }
    mvGridPointIndex = vGridPointIndex;

    // find valid planar patches and get plane parameters in local
    vGridPointIndex.clear();
    vector<Eigen::Vector4d> vPlaneParam;
    int numPoint = 0;
    for (int nGrid = 0; nGrid < mvGridPointIndex.size(); nGrid++)
    {

        if (mvGridPointIndex[nGrid].size() < MinNumGroundPoints * double(mNumIntegrationScan)) // size limitation
            continue;
        Eigen::Vector4f params;
        Eigen::Vector3f center;

        std::vector<Eigen::Vector3f> allPoints;
        for (int nP = 0; nP < mvGridPointIndex[nGrid].size(); nP++)
        {
            PointType p = mvGroundPointsLocal[mvGridPointIndex[nGrid][nP]];
            allPoints.push_back(Eigen::Vector3f(p.x, p.y, p.z));
        }

        if (planeFitting(allPoints, params, center))
        {
            vGridPointIndex.push_back(mvGridPointIndex[nGrid]);
            vPlaneParam.push_back(params.cast<double>());
            numPoint += mvGridPointIndex[nGrid].size();
        }
    }
    mvPlaneParamLocal = vPlaneParam;
    mNumPlanarPatch = vGridPointIndex.size();
    mvGridPointIndex = vGridPointIndex;
    fMapLog << "Plane patch numbers: " << mNumPlanarPatch <<" "<< endl;

    // CONTROL
    // if no enough valid planar patch
    if (mNumPlanarPatch < 3)
    {
        fMapLog << "\tError! Num of planar patches are too few: " << mNumPlanarPatch << endl;
        f_mapping_debug << "Iscan " << mnCount << ": \tError! Num of planar patches are too few: " << mNumPlanarPatch << endl;
        Flag_Scan_Inter = false;
    }  
}

/********************************************************
 * Perform optimization for scan integration
 * Consider distortion in the computation.
 ********************************************************/
/**
 * @brief Optimize tree and ground features jointly using trajectory and cylinder model.
 * 
 * @param fix_traj_flag If true, fixes trajectory and only updates feature parameters
 * 
 * @details Uses Ceres solver to optimize scan poses and tree/ground features by minimizing residuals.
 * Removes noisy observations and filters poor tree fits based on inlier ratio.
 */
void Mapping::optimizeIntegratedScan(bool fix_traj_flag)
{
    double tree_std = mPara.iscanTreeStd;
    double ground_std = mPara.iscanGroundStd;
    int iter = 0;
    //(TODO/MODI): change max_iter from 2 to 1
    int max_iter = 2;

    //input from class
    vector<vector<PointType>> tree_raw_points_temp = mvTreePointsRaw;
    vector<vector<int>> ground_point_index_temp = mvGridPointIndex;
    vector<CylinderPara> tree_para_local = mvTreeParamsLocal;
    vector<Eigen::Vector4d> plane_para_local = mvPlaneParamLocal;
    vector<bool> valid_tree;
    //initial value for the pose
    int numScan = mNumIntegrationScan;
    int numEpoch = mNumIntegrationScan + 1; // one epoch corresponds to the starting

    vector<Eigen::Vector3d> vec_r_ini;
    vector<Eigen::Matrix3d> vec_R_ini;
    vector<Eigen::Vector3d> vec_r_refined;
    vector<Eigen::Matrix3d> vec_R_refined;
    for (int nEpoch = 0; nEpoch < numEpoch; nEpoch++)
    {
        if (nEpoch == 0)
        {
            vec_r_ini.push_back(mr_lu_local_t0);
            vec_R_ini.push_back(mR_lu_local_t0);
        }
        else
        {
            int scanId = nEpoch - 1;
            vec_r_ini.push_back(mvTempScan[scanId].r_lu_local_updated);
            vec_R_ini.push_back(mvTempScan[scanId].R_lu_local_updated);
        }
    }

    //iteratively conduct the optimization
    while (iter < max_iter)
    {
        fMapLog << "Iter " << ++iter;

        // if(iter==1){
        //     estimate_tree_normal_ = false; estimate_tree_radius_ = false; tree_outlier_multiplier_ = 5.0;
        // }else if(iter == 2){
        //     estimate_tree_normal_ = false; estimate_tree_radius_ = true; tree_outlier_multiplier_ = 1.0;
        // }
        //estimate_tree_normal_ = false; estimate_tree_radius_ = false; tree_outlier_multiplier_ = 3.0;

        //Number of tree with sufficient number of points
        mNumValidTree = 0;
        valid_tree.clear();
        for (int nTree = 0; nTree < mvIntegratedTree.size(); nTree++)
        {
            if(tree_raw_points_temp[nTree].size() > min(int(MinNumTreePointPerScan * double(mNumIntegrationScan)), mPara.minNumPtsPerIscanTree))
            {
                valid_tree.push_back(true);
                mNumValidTree++;
            }
            else
            {
                valid_tree.push_back(false);
            }
        }

        // initialization
        int numCylinder = mNumValidTree;
        int numPlane = mNumPlanarPatch;

        // unknowns
        int numPosUnknown = numEpoch * PosBlockSize;
        int numOriUnknown = numEpoch * OriBlockSize;
        int numCylinderUnknown = numCylinder * CylinderBlockSize;
        int numPlaneUnknown = numPlane * PlaneBlockSize;
        int numUnknown = numPosUnknown + numOriUnknown + numCylinderUnknown + numPlaneUnknown; // all unknowns
        double *poss = new double[numPosUnknown];
        double *poss_std = new double[numPosUnknown];
        double *oris = new double[numOriUnknown];
        double *oris_std = new double[numOriUnknown];
        double *cylinders = new double[numCylinderUnknown];
        double *cylinders_std = new double[numCylinderUnknown];
        double *planes = new double[numPlaneUnknown];
        double *planes_std = new double[numPlaneUnknown];

        // fill parameter
        // 1. trajectory information
        // 0: intial, 1-> scan id = 0;
        //currently there is no fixed pose (TODO: may need to check)
        //mr/R_lu_local_t0, mvTempScan[].r/R_lu_local_updated
        for (int nEpoch = 0; nEpoch < numEpoch; nEpoch++)
        {
            *(poss + nEpoch * PosBlockSize + 0) = vec_r_ini[nEpoch](0);
            *(poss + nEpoch * PosBlockSize + 1) = vec_r_ini[nEpoch](1);
            *(poss + nEpoch * PosBlockSize + 2) = vec_r_ini[nEpoch](2);

            Eigen::Quaterniond q(vec_R_ini[nEpoch]);
            *(oris + nEpoch * OriBlockSize + 0) = q.x();
            *(oris + nEpoch * OriBlockSize + 1) = q.y();
            *(oris + nEpoch * OriBlockSize + 2) = q.z();
            *(oris + nEpoch * OriBlockSize + 3) = q.w();

            if(fix_traj_flag)
            {
            *(poss_std + nEpoch * PosBlockSize + 0) = 1.0E-20;
            *(poss_std + nEpoch * PosBlockSize + 1) = 1.0E-20;
            *(poss_std + nEpoch * PosBlockSize + 2) = 1.0E-20;
            *(oris_std + nEpoch * OriBlockSize + 0) = 1.0E-20;
            *(oris_std + nEpoch * OriBlockSize + 1) = 1.0E-20;
            *(oris_std + nEpoch * OriBlockSize + 2) = 1.0E-20;
            *(oris_std + nEpoch * OriBlockSize + 3) = 1.0E-20;
            }
            else
            {
            *(poss_std + nEpoch * PosBlockSize + 0) = 1.0E+20;
            *(poss_std + nEpoch * PosBlockSize + 1) = 1.0E+20;
            *(poss_std + nEpoch * PosBlockSize + 2) = 1.0E+20;
            *(oris_std + nEpoch * OriBlockSize + 0) = 1.0E+20;
            *(oris_std + nEpoch * OriBlockSize + 1) = 1.0E+20;
            *(oris_std + nEpoch * OriBlockSize + 2) = 1.0E+20;
            *(oris_std + nEpoch * OriBlockSize + 3) = 1.0E+20;
            }
        }

        // 2. cylinder parameter
        //tree_para_local for true valid_tree, estimate X Y only
        int countCylinder = 0;
        for (int nTree = 0; nTree < valid_tree.size(); nTree++)
        {
            if (valid_tree[nTree])
            {
                *(cylinders + countCylinder * CylinderBlockSize + 0) = tree_para_local[nTree].x(0); // XYZ (Z is nominal)
                *(cylinders + countCylinder * CylinderBlockSize + 1) = tree_para_local[nTree].x(1);
                *(cylinders + countCylinder * CylinderBlockSize + 2) = tree_para_local[nTree].x(2);
                *(cylinders + countCylinder * CylinderBlockSize + 3) = tree_para_local[nTree].n(0);
                *(cylinders + countCylinder * CylinderBlockSize + 4) = tree_para_local[nTree].n(1);
                *(cylinders + countCylinder * CylinderBlockSize + 5) = tree_para_local[nTree].n(2);
                *(cylinders + countCylinder * CylinderBlockSize + 6) = tree_para_local[nTree].r;

                *(cylinders_std + countCylinder * CylinderBlockSize + 0) = 1.0E+20; // XYZ
                *(cylinders_std + countCylinder * CylinderBlockSize + 1) = 1.0E+20;
                *(cylinders_std + countCylinder * CylinderBlockSize + 2) = 1.0E-20;
                *(cylinders_std + countCylinder * CylinderBlockSize + 3) = estimate_tree_normal_?1.0E+20:1.0E-20; // ux, uy, uz
                *(cylinders_std + countCylinder * CylinderBlockSize + 4) = estimate_tree_normal_?1.0E+20:1.0E-20;
                *(cylinders_std + countCylinder * CylinderBlockSize + 5) = 1.0E-20;
                *(cylinders_std + countCylinder * CylinderBlockSize + 6) = estimate_tree_radius_?1.0E+20:1.0E-20; // r
                countCylinder++;
            }
        }

        // 3. plane parameters
        for (int nPlane = 0; nPlane < ground_point_index_temp.size(); nPlane++)
        {
            Eigen::Vector4d params = plane_para_local[nPlane];
            int fixIndex = planeParamTrans(params);

            // fMapLog << params.transpose() << "\t" << fixIndex << endl;
            *(planes + nPlane * 4 + 0) = params(0);
            *(planes + nPlane * 4 + 1) = params(1);
            *(planes + nPlane * 4 + 2) = params(2);
            *(planes + nPlane * 4 + 3) = params(3);

            *(planes_std + nPlane * 4 + 0) = 1.0E+20;
            *(planes_std + nPlane * 4 + 1) = 1.0E+20;
            *(planes_std + nPlane * 4 + 2) = 1.0E+20;
            *(planes_std + nPlane * 4 + 3) = 1.0E+20;

            // fixing the element equal to 1
            *(planes_std + nPlane * 4 + fixIndex) = 1.0E-20;
        }

        //------------- check fixed parameters
        int numFixedUnknown = 0;
        for (int i = 0; i < numPosUnknown; i++)
        {
            double temp;
            temp = *(poss_std + i);
            if (temp < SigmaFixed)
                numFixedUnknown += 1;
        }
        for (int i = 0; i < numOriUnknown; i++)
        {
            double temp;
            temp = *(oris_std + i);
            if (temp < SigmaFixed)
                numFixedUnknown += 1;
        }
        for (int i = 0; i < numCylinderUnknown; i++)
        {
            double temp;
            temp = *(cylinders_std + i);
            if (temp < SigmaFixed)
                numFixedUnknown += 1;
        }
        for (int i = 0; i < numPlaneUnknown; i++)
        {
            double temp;
            temp = *(planes_std + i);
            if (temp < SigmaFixed)
                numFixedUnknown += 1;
        }

        //--------------------------build ceres  ------------------------------
        ceres::LossFunction *loss_function = new ceres::HuberLoss(0.1);
        // ceres::LossFunction *loss_function = NULL; // new ceres::HuberLoss(0.1);
        ceres::LocalParameterization *q_parameterization = new ceres::EigenQuaternionParameterization(); // x, y, z, w
		ceres::Problem problem;
        ceres::CostFunction *costFunction;

        //--------------------------add parameter block ------------------------------
        for (int nEpoch = 0; nEpoch < numEpoch; nEpoch++)
        {
            problem.AddParameterBlock(poss + nEpoch * PosBlockSize, PosBlockSize);
        }
        for (int nEpoch = 0; nEpoch < numEpoch; nEpoch++)
        {
            problem.AddParameterBlock(oris + nEpoch * OriBlockSize, OriBlockSize, q_parameterization);
        }
        for (int nCylinder = 0; nCylinder < numCylinder; nCylinder++)
        {
            problem.AddParameterBlock(cylinders + nCylinder * CylinderBlockSize, CylinderBlockSize);
        }
        for (int nPlane = 0; nPlane < numPlane; nPlane++)
        {
            problem.AddParameterBlock(planes + nPlane * PlaneBlockSize, PlaneBlockSize);
        }

        //-------------------------add observation -------------------------------
        int cylinderId = 0;
        int numTreeObs = 0;
        for (int nTree = 0; nTree < valid_tree.size(); nTree++)
        {
            if (!valid_tree[nTree])
                continue;

            double *cylinder = cylinders + cylinderId * CylinderBlockSize; // corresponding cylinder block
            // for each tree point
            for (int nP = 0; nP < tree_raw_points_temp[nTree].size(); nP++)
            {
                int scanId = int(tree_raw_points_temp[nTree][nP].intensity / 100);
                double t_ratio = double(tree_raw_points_temp[nTree][nP].intensity - int(tree_raw_points_temp[nTree][nP].intensity));

                double *t_prev = poss + scanId * PosBlockSize;
                double *q_prev = oris + scanId * OriBlockSize;
                double *t_cur = poss + (scanId + 1) * PosBlockSize;
                double *q_cur = oris + (scanId + 1) * OriBlockSize;
                Eigen::Vector3d curr_point(tree_raw_points_temp[nTree][nP].x, tree_raw_points_temp[nTree][nP].y, tree_raw_points_temp[nTree][nP].z);
                // costFunction = CylinderFactorTwoEpoch::Create(curr_point, t_ratio, tree_std);
                // problem.AddResidualBlock(costFunction, loss_function, t_prev, q_prev, t_cur, q_cur, cylinder);

                costFunction = CylinderFactorOneEpoch::Create(curr_point, tree_std);
                problem.AddResidualBlock(costFunction, loss_function, t_prev, q_prev, cylinder);
                
                numTreeObs++;
            }
            cylinderId++;
        }

        int numGroundObs = 0;
        for (int nPlane = 0; nPlane < ground_point_index_temp.size(); nPlane++)
        {
            int planeId = nPlane;
            double *plane = planes + planeId * PlaneBlockSize; // corresponding plane block

            for (int nP = 0; nP < ground_point_index_temp[nPlane].size(); nP++)
            {
                int index = ground_point_index_temp[nPlane][nP];
                int scanId = int(mvGroundPointsRaw[index].intensity / 100);
                double t_ratio = double(mvGroundPointsRaw[index].intensity - int(mvGroundPointsRaw[index].intensity));
                double *t_prev = poss + scanId * PosBlockSize;
                double *q_prev = oris + scanId * OriBlockSize;
                double *t_cur = poss + (scanId + 1) * PosBlockSize;
                double *q_cur = oris + (scanId + 1) * OriBlockSize;

                Eigen::Vector3d curr_point(mvGroundPointsRaw[index].x, mvGroundPointsRaw[index].y, mvGroundPointsRaw[index].z);
                // costFunction = PlaneFactorTwoEpoch::Create(curr_point, t_ratio, ground_std);
                // problem.AddResidualBlock(costFunction, loss_function, t_prev, q_prev, t_cur, q_cur, plane);

                costFunction = PlaneFactorOneEpoch::Create(curr_point, ground_std);
                problem.AddResidualBlock(costFunction, loss_function, t_prev, q_prev, plane);
                numGroundObs++;
            }
        }

        // information related to odometry thread
        if (mbOdoIntegration)
        {
            for (int nScan = 0; nScan < numScan; nScan++)
            {
                // for each scan, epoches from [nScan, nScan+1] are used
                double *t_ini = poss + nScan * PosBlockSize;
                double *q_ini = oris + nScan * OriBlockSize;
                double *t_end = poss + (nScan + 1) * PosBlockSize;
                double *q_end = oris + (nScan + 1) * OriBlockSize;

                // transformation estimated by odometry
                Eigen::Quaterniond q_odo(mvTempScan[nScan].R_lut2_lut1);
                Eigen::Vector3d t_odo = mvTempScan[nScan].r_lut2_lut1;

                double acc_t;
                double acc_q;
                if (mPara.odo_from_trajectory_flag)
                {
                    acc_t = mPara.odoPositionStdTraj;
                    acc_q = deg2rad(mPara.odoOrientationStdTraj);
                }
                else
                {
                    acc_t = mPara.odoPositionStd;
                    acc_q = deg2rad(mPara.odoOrientationStd);
                }

                Eigen::Matrix<double, 6, 6> sqrt_information = Eigen::MatrixXd::Identity(6, 6);
                sqrt_information(0, 0) = 1.0 / acc_t;
                sqrt_information(1, 1) = 1.0 / acc_t;
                sqrt_information(2, 2) = 1.0 / acc_t;
                sqrt_information(3, 3) = 1.0 / acc_q;
                sqrt_information(4, 4) = 1.0 / acc_q;
                sqrt_information(5, 5) = 1.0 / acc_q;

                costFunction = PoseGraph3dErrorTerm::Create(t_odo, q_odo, sqrt_information);
				if (mPara.odo_from_trajectory_flag)
				{
					problem.AddResidualBlock(costFunction, nullptr, t_ini, q_ini, t_end, q_end);
				}
				else
				{
					problem.AddResidualBlock(costFunction, loss_function, t_ini, q_ini, t_end, q_end);
				}
                double dis_std = 0.5 ;
				costFunction = DistanceBetweenTwoEpochsConstraints::Create(t_odo.norm(),dis_std);
                problem.AddResidualBlock(costFunction, nullptr, t_ini, q_ini, t_end, q_end);

            }
        }

        //--------------------------Fix Parameters  ------------------------------
        // set constant unknowns
        Checkforconstantparams(&problem, poss, poss_std, PosBlockSize, numEpoch);
        Checkforconstantparams(&problem, oris, oris_std, OriBlockSize, numEpoch);
        Checkforconstantparams(&problem, cylinders, cylinders_std, CylinderBlockSize, numCylinder);
        Checkforconstantparams(&problem, planes, planes_std, PlaneBlockSize, numPlane);

        //-------------------------- Solve -------------------------------------------
        TicToc t_solver;
        ceres::Solver::Options solverOptions;
        solverOptions.max_num_iterations = 10;
        solverOptions.linear_solver_type = ceres::LinearSolverType::SPARSE_NORMAL_CHOLESKY; // ceres::LinearSolverType::SPARSE_NORMAL_CHOLESKY;
        solverOptions.minimizer_progress_to_stdout = false;
        ceres::Solver::Summary summary;
        ceres::Solve(solverOptions, &problem, &summary);
        fMapLog << ", Time: " << t_solver.toc() << " ms"<< endl;
        cout << "Time: " << t_solver.toc() << endl;

        //******************** Update parameters **********************
        // 1. update trajectory info
        vec_r_refined.clear();
        vec_R_refined.clear();
        for (int nEpoch = 0; nEpoch < numEpoch; nEpoch++)
        {
            Eigen::Vector3d r{*(poss + nEpoch * PosBlockSize + 0), *(poss + nEpoch * PosBlockSize + 1), *(poss + nEpoch * PosBlockSize + 2)};
            Eigen::Quaterniond q_refined(*(oris + nEpoch * OriBlockSize + 3), *(oris + nEpoch * OriBlockSize + 0), *(oris + nEpoch * OriBlockSize + 1), *(oris + nEpoch * OriBlockSize + 2));
            Eigen::Matrix3d R = q_refined.toRotationMatrix();
            vec_r_refined.push_back(r);
            vec_R_refined.push_back(R);
        }

        // 2. cylinder parameter. valid_tree -> false when unreasonable cylinder
        CylinderPara treePara;
        vector<CylinderPara> vTreeParamsUpdated = tree_para_local;
        countCylinder = 0;
        for (int nTree = 0; nTree < valid_tree.size(); nTree++)
        {
            if (valid_tree[nTree])
            {
                treePara.x(0) = *(cylinders + countCylinder * CylinderBlockSize + 0);
                treePara.x(1) = *(cylinders + countCylinder * CylinderBlockSize + 1);
                treePara.x(2) = *(cylinders + countCylinder * CylinderBlockSize + 2);
                treePara.n(0) = *(cylinders + countCylinder * CylinderBlockSize + 3);
                treePara.n(1) = *(cylinders + countCylinder * CylinderBlockSize + 4);
                treePara.n(2) = *(cylinders + countCylinder * CylinderBlockSize + 5);
                treePara.r = *(cylinders + countCylinder * CylinderBlockSize + 6);
                countCylinder++;
                // failure due to wrong normal vector, currently deactivated as only X and Y are estimated
                if (abs(treePara.n(0)) > 1.0 || abs(treePara.n(1)) > 1.0)
                    valid_tree[nTree] = false;
                else
                    vTreeParamsUpdated[nTree] = treePara;
            }
        }

        // 3. plane parameters
        vector<Eigen::Vector4d> vPlaneParamsUpdated;
        for (int nPlane = 0; nPlane < ground_point_index_temp.size(); nPlane++)
        {
            Eigen::Vector4d params{*(planes + nPlane * PlaneBlockSize + 0), *(planes + nPlane * PlaneBlockSize + 1), *(planes + nPlane * PlaneBlockSize + 2), *(planes + nPlane * PlaneBlockSize + 3)};
            vPlaneParamsUpdated.push_back(params);
        }

        //******************** Compute residual **********************
        double res_tree = 0.0, res_tree_ini = 0.0, res_ground = 0.0, res_ground_ini = 0.0;
        int numInlierTreePoints = 0;
        int num_tree_1 = 0, num_tree_2 = 0; //num_tree_2: final survived trees
        for (int nTree = 0; nTree < valid_tree.size(); nTree++)
        {
            if (!valid_tree[nTree])
                continue;
            
            num_tree_1++;

            vector<PointType> inlierTreePoints;
            
            // for each tree point
            for (int nP = 0; nP < tree_raw_points_temp[nTree].size(); nP++)
            {
                int scanId = int(tree_raw_points_temp[nTree][nP].intensity / 100);
                Eigen::Vector3d curr_point(tree_raw_points_temp[nTree][nP].x, tree_raw_points_temp[nTree][nP].y, tree_raw_points_temp[nTree][nP].z);
                double t_ratio = double(tree_raw_points_temp[nTree][nP].intensity - int(tree_raw_points_temp[nTree][nP].intensity));

                Eigen::Quaterniond q_prev_ini{vec_R_ini[scanId]};
                Eigen::Quaterniond q_cur_end_ini{vec_R_ini[scanId + 1]};
                Eigen::Quaterniond q_cur_ini = q_prev_ini.slerp(t_ratio, q_cur_end_ini);
                Eigen::Vector3d t_cur_ini = (1.0 - t_ratio) * vec_r_ini[scanId] + t_ratio * vec_r_ini[scanId + 1];
                // Eigen::Vector3d rIm_ini = q_cur_ini * curr_point + t_cur_ini;

                Eigen::Vector3d rIm_ini = vec_R_ini[scanId] * curr_point + vec_r_ini[scanId];
                double d_ini = computePoint2lineDistance(rIm_ini, tree_para_local[nTree].x, tree_para_local[nTree].x + tree_para_local[nTree].n);
                double res_ini = d_ini - tree_para_local[nTree].r;
                res_tree_ini = res_tree_ini + res_ini * res_ini;

                Eigen::Quaterniond q_prev{vec_R_refined[scanId]};
                Eigen::Quaterniond q_cur_end{vec_R_refined[scanId + 1]};
                Eigen::Quaterniond q_cur = q_prev.slerp(t_ratio, q_cur_end);
                Eigen::Vector3d t_cur = (1.0 - t_ratio) * vec_r_refined[scanId] + t_ratio * vec_r_refined[scanId + 1];
                // Eigen::Vector3d rIm = q_cur * curr_point + t_cur;
                Eigen::Vector3d rIm = vec_R_refined[scanId] * curr_point + vec_r_refined[scanId];

                double d_refine = computePoint2lineDistance(rIm, vTreeParamsUpdated[nTree].x, vTreeParamsUpdated[nTree].x + vTreeParamsUpdated[nTree].n);
                double res_refine = d_refine - vTreeParamsUpdated[nTree].r;
                res_tree = res_tree + res_refine * res_refine;
                //(TODO:CHUNXI)filter points using fitting cylinder 
                 //if (res_refine < tree_outlier_multiplier_ * tree_std)
                {
                    inlierTreePoints.push_back(tree_raw_points_temp[nTree][nP]);
                    numInlierTreePoints++;
                }
            }
            if(double(inlierTreePoints.size())/double(tree_raw_points_temp[nTree] .size())<ValidInlierRatio)
            {
                valid_tree[nTree] = false;
                tree_raw_points_temp[nTree].clear();
            }
            else
            {
                tree_raw_points_temp[nTree] = inlierTreePoints;
                num_tree_2++;
            }
        }
        fMapLog << "\tNumber of trees with sufficient number of points " << mNumValidTree << " (" << min(int(MinNumTreePointPerScan * double(mNumIntegrationScan)), mPara.minNumPtsPerIscanTree) <<")" << endl;
        fMapLog << "\tNumber of trees with reasonable cylinder parameters: " << num_tree_1 << endl;
        fMapLog << "\tNumber of trees with sufficient percent of inlier: " <<  num_tree_2<<endl;

        //check if the distribution of tree is good enough
        int num_isolated_tree = 0;
        for (int nTree = 0; nTree < valid_tree.size(); nTree++)
        {
            if (!valid_tree[nTree])
                continue;
            
            bool flag = true;
            for(int n_cand_tree = 0; n_cand_tree < valid_tree.size(); n_cand_tree++)
            {
                if (!valid_tree[n_cand_tree])
                    continue;

                if (n_cand_tree == nTree)
                    continue;

                if ((vTreeParamsUpdated[nTree].x(0) - vTreeParamsUpdated[n_cand_tree].x(0)) * (vTreeParamsUpdated[nTree].x(0) - vTreeParamsUpdated[n_cand_tree].x(0))
                + (vTreeParamsUpdated[nTree].x(1) - vTreeParamsUpdated[n_cand_tree].x(1)) * (vTreeParamsUpdated[nTree].x(1)- vTreeParamsUpdated[n_cand_tree].x(1)) < MinDistanceIsolatedTree *MinDistanceIsolatedTree)
                {
                    flag = false;
                    break;
                }
            }
            if(flag)            
                num_isolated_tree ++;
        }
        fMapLog << "\tNumber of isolated trees with sufficient percent of inlier: " << num_isolated_tree<<endl;

        int numInlierGroundPoints = 0;
        for (int nPlane = 0; nPlane < ground_point_index_temp.size(); nPlane++)
        {
            vector<int> inlierGridPointIndex;

            for (int nP = 0; nP < ground_point_index_temp[nPlane].size(); nP++)
            {
                int index = ground_point_index_temp[nPlane][nP];
                int scanId = int(mvGroundPointsRaw[index].intensity / 100);
                Eigen::Vector3d curr_point(mvGroundPointsRaw[index].x, mvGroundPointsRaw[index].y, mvGroundPointsRaw[index].z);
                double t_ratio = double(mvGroundPointsRaw[index].intensity - int(mvGroundPointsRaw[index].intensity));

                Eigen::Quaterniond q_prev_ini{vec_R_ini[scanId]};
                Eigen::Quaterniond q_cur_end_ini{vec_R_ini[scanId + 1]};
                Eigen::Quaterniond q_cur_ini = q_prev_ini.slerp(t_ratio, q_cur_end_ini);
                Eigen::Vector3d t_cur_ini = (1.0 - t_ratio) * vec_r_ini[scanId] + t_ratio * vec_r_ini[scanId + 1];
                // Eigen::Vector3d rIm_ini = q_cur_ini * curr_point + t_cur_ini;
                Eigen::Vector3d rIm_ini = vec_R_ini[scanId] * curr_point + vec_r_ini[scanId];

                Eigen::Vector3d normal_ini = plane_para_local[nPlane].head(3);
                double d_ini = (rIm_ini.transpose() * normal_ini + plane_para_local[nPlane](3)) / normal_ini.norm();
                res_ground_ini = res_ground_ini + d_ini * d_ini;

                Eigen::Quaterniond q_prev{vec_R_refined[scanId]};
                Eigen::Quaterniond q_cur_end{vec_R_refined[scanId + 1]};
                Eigen::Quaterniond q_cur = q_prev.slerp(t_ratio, q_cur_end);
                Eigen::Vector3d t_cur = (1.0 - t_ratio) * vec_r_refined[scanId] + t_ratio * vec_r_refined[scanId + 1];
                // Eigen::Vector3d rIm = q_cur * curr_point + t_cur;
                Eigen::Vector3d rIm = vec_R_refined[scanId] * curr_point + vec_r_refined[scanId];

                Eigen::Vector3d normal_updated = vPlaneParamsUpdated[nPlane].head(3);
                double d_refine = (rIm.transpose() * normal_updated + vPlaneParamsUpdated[nPlane](3)) / normal_updated.norm();
                res_ground = res_ground + d_refine * d_refine;

                if (d_refine < 1.0 * ground_std)
                {
                    inlierGridPointIndex.push_back(index);
                    numInlierGroundPoints++;
                }
            }
            ground_point_index_temp[nPlane] = inlierGridPointIndex;
        }

        cout << "\tError for tree/ground point: " << sqrt(res_tree_ini / double(numTreeObs)) << " -> " << sqrt(res_tree / double(numTreeObs)) << " m";
        cout << "\t" << sqrt(res_ground_ini / double(numGroundObs)) << " -> " << sqrt(res_ground / double(numGroundObs)) << " m" << endl;

        fMapLog << "\tError for tree/ground point: " << sqrt(res_tree_ini / double(numTreeObs)) << " -> " << sqrt(res_tree / double(numTreeObs)) << " m";
        fMapLog << "\t" << sqrt(res_ground_ini / double(numGroundObs)) << " -> " << sqrt(res_ground / double(numGroundObs)) << " m" << endl;
        fMapLog << "\tInlier tree point: " << numTreeObs << " -> " << numInlierTreePoints;
        fMapLog << "\tInlier ground point: " << numGroundObs << " -> " << numInlierGroundPoints << endl;



        delete[] poss;
        delete[] poss_std;
        delete[] oris;
        delete[] oris_std;
        delete[] cylinders;
        delete[] cylinders_std;
        delete[] planes;
        delete[] planes_std;

        //in case trajectory is fixed, only feature parameters are updated
        if(fix_traj_flag)
        {
            tree_para_local = vTreeParamsUpdated;
            plane_para_local = vPlaneParamsUpdated;
            break;
        }

        //CONTROL
        // if enough valid treesupdate the file in the class
        if (num_isolated_tree < Min_Num_Tree_)
        {
            fMapLog << "\tError! Num of trees in scan integration too few: " << num_isolated_tree  << endl;
            f_mapping_debug << "Iscan " << mnCount <<": \tNumber of trees in scan integration too few: " << num_isolated_tree  << endl;
            Flag_Scan_Inter = false;
            break;
        } 
        // if outlier ratio of tree points is too large
        else if (double(numInlierTreePoints) / double(numTreeObs) < 0.6)
        {
            fMapLog << "\tError! Too many outlier tree points: " << double(numInlierTreePoints) / double(numTreeObs)  << endl;
            f_mapping_debug << "Iscan " << mnCount <<": \tError! Too many outlier tree points: " << double(numInlierTreePoints) / double(numTreeObs)  << endl;
            Flag_Scan_Inter = false;
            break;
        }
        else
        {
            vec_R_ini = vec_R_refined;
            vec_r_ini = vec_r_refined;
            tree_para_local = vTreeParamsUpdated;
            plane_para_local = vPlaneParamsUpdated;
        }
    }


    // if success
    if(fix_traj_flag)
    {
        mvTreePointsRaw = tree_raw_points_temp;
        mvTreeParamsLocal = tree_para_local;
        mvGridPointIndex = ground_point_index_temp;
        mvPlaneParamLocal = plane_para_local;
        mvValidTrees= valid_tree; 
        mv_R_lu_local = mv_R_lu_local_ini;
        mv_r_lu_local = mv_r_lu_local_ini;
        return;
    }

    if (Flag_Scan_Inter)
    {
        mv_R_lu_local = vec_R_refined;
        mv_r_lu_local = vec_r_refined;
        for (int nScan = 0; nScan < numScan; nScan++)
        {
            int epochId = nScan + 1;
            mvTempScan[nScan].R_lu_local_updated = vec_R_refined[epochId];
            mvTempScan[nScan].r_lu_local_updated = vec_r_refined[epochId];
        }
        mvTreePointsRaw = tree_raw_points_temp;
        mvTreeParamsLocal = tree_para_local;
        mvGridPointIndex = ground_point_index_temp;
        mvPlaneParamLocal = plane_para_local;
        mvValidTrees= valid_tree; 
    }
    else
    {
        //dont update the pose results
        mv_R_lu_local = mv_R_lu_local_ini;
        mv_r_lu_local = mv_r_lu_local_ini;
    }


}


/********************************************************
 * For the invalid tree, 
 * Conduct the initialization of the T_lu_m for the first Iscan (change to meaningful location TODO)
 ********************************************************/
/**
 * @brief Update and merge invalid or outlier trees after scan integration.
 * 
 * @details Recomputes parameters for invalid trees and merges trees that are spatially close.
 * Ensures robustness against repeated observations and duplicates.
 */
void Mapping::UpdateMergeTree()
{
    //update the center of tree which is not valid
    for (int nTree = 0; nTree < mvValidTrees.size(); nTree++)
    {
        if (mvValidTrees[nTree])
            continue;

        if(mvTreePointsRaw[nTree].empty())
            continue;

        //for the remaining
        Eigen::Vector3d center = Eigen::Vector3d::Zero(); 
        // for each tree point
        for (int nP = 0; nP < mvTreePointsRaw[nTree].size(); nP++)
        {
            int scanId = int(mvTreePointsRaw[nTree][nP].intensity / 100);
            Eigen::Vector3d curr_point(mvTreePointsRaw[nTree][nP].x, mvTreePointsRaw[nTree][nP].y, mvTreePointsRaw[nTree][nP].z);
            double t_ratio = double(mvTreePointsRaw[nTree][nP].intensity - int(mvTreePointsRaw[nTree][nP].intensity));

            Eigen::Quaterniond q_prev{mv_R_lu_local[scanId]};
            Eigen::Quaterniond q_cur_end{mv_R_lu_local[scanId + 1]};
            Eigen::Quaterniond q_cur = q_prev.slerp(t_ratio, q_cur_end);
            Eigen::Vector3d t_cur = (1.0 - t_ratio) * mv_r_lu_local[scanId] + t_ratio * mv_r_lu_local[scanId + 1];
            // Eigen::Vector3d r_I_local = q_cur * curr_point + t_cur;

            Eigen::Vector3d r_I_local = mv_R_lu_local[scanId] * curr_point + mv_r_lu_local[scanId];

            center += r_I_local;
        }
        center = center / double(mvTreePointsRaw[nTree].size());
        CylinderPara treePara;
        treePara.x = center;
        treePara.n << 0.0, 0.0, 1.0;
        treePara.r = 0.0;
        mvTreeParamsLocal[nTree] = treePara;
    }

    // merge tree based on the x and y component
    vector<bool> visited_flags(mvValidTrees.size(), false);
    int num_merge_tree = 0;
    for (int i = 0; i < mvValidTrees.size(); i++)
    {
        if (visited_flags[i])
            continue;
        for (int j = i + 1; j < mvValidTrees.size(); j++)
        {
            double dis_sqr = (mvTreeParamsLocal[i].x(0) - mvTreeParamsLocal[j].x(0)) * (mvTreeParamsLocal[i].x(0) - mvTreeParamsLocal[j].x(0)) + (mvTreeParamsLocal[i].x(1) - mvTreeParamsLocal[j].x(1)) * (mvTreeParamsLocal[i].x(1) - mvTreeParamsLocal[j].x(1));
            if (dis_sqr > 0.2 * 0.2)
                continue;

            // if the distance is close, merge one with small number of points to the other
            int idx_ref, idx_can;
            if (mvTreePointsRaw[i].size() > mvTreePointsRaw[j].size())
            {
                idx_ref = i;
                idx_can = j;
            }
            else
            {
                idx_ref = j;
                idx_can = i;
            }

            mvTreePointsRaw[idx_ref].insert(mvTreePointsRaw[idx_ref].end(), mvTreePointsRaw[idx_can].begin(), mvTreePointsRaw[idx_can].end());
            mvTreePointsRaw[idx_can].clear();
            visited_flags[idx_can] = true;
            num_merge_tree++;
        }
    }
    fMapLog << "Merge tree after Iscan optimization: " << num_merge_tree << " trees merged"<<endl;
}


/********************************************************
 * Wrap up the info from the scan integration for local to mapping
 * Conduct the initialization of the T_lu_m for the first Iscan (change to meaningful location TODO)
 ********************************************************/
/**
 * @brief Broadcasts the result of scan integration to a new IntegratedScan instance.
 *
 * Sets initial pose information for the first Iscan based on whether the mapping is
 * being initialized or is referencing a global map. Updates tree and ground points,
 * performs downsampling, and computes transformations to local frames.
 */
void Mapping::broadcastIScan()
{
    fMapLog << "Save integrated scan, " << flush;

    pCurrentIScan = new IntegratedScan();
    pCurrentIScan->numScan = mNumIntegrationScan;
    pCurrentIScan->TimeInit = mvTempScan[0].mTimeEnd - mvTempScan[0].mTimeTracking; //time init corresponding to the first epoch
    pCurrentIScan->TimeEnd = mvTempScan.back().mTimeEnd; //time end corresponding to the last 
    pCurrentIScan->index = mnCount;

    // Pose for the timeInit
    if (!mbInit) //if current integration is initialization
    {
        // initialization of the mapping frame
        // Case 1. no global map-> initialize as the local frame
        if (!global_map_flag_)
        {
            pCurrentIScan->r_lu_m_ini = mv_r_lu_local[0];
            pCurrentIScan->R_lu_m_ini = mv_R_lu_local[0];
        }
        // Case 2. there is global map, localization is needed. Currently from traj (TODO)
        else
        {
            Eigen::Vector3d r_b_global_t1, r_lu_global_t1, r_lu_m_t1;
            Eigen::Matrix3d R_b_global_t1, R_lu_global_t1, R_lu_m_t1;

            int idx = 0;

            if (!mpTraj->bopInterpolation((mvTempScan[0].mTimeEnd - mvTempScan[0].mTimeTracking) / 1000.0, -1, idx, r_b_global_t1, R_b_global_t1))
                throw std::runtime_error("Wrong Time");

            //the initial value from trajectory, no necessary its accurate
            r_lu_m_t1 = r_b_global_t1 + R_b_global_t1 * mpTraj->r_lu_b - global_const_shift_;
            R_lu_m_t1 = R_b_global_t1 * mpTraj->R_lu_b;
            pCurrentIScan->r_lu_m_ini = r_lu_m_t1;
            pCurrentIScan->R_lu_m_ini = R_lu_m_t1;
        }
    }
    else //the pose related to t_end of previous integrated scan
    {
        // the r_lu_m for last individual scan from the previous integrated scan
        pCurrentIScan->r_lu_m_ini = mvpIScans.back()->v_r_mapping.back();
        pCurrentIScan->R_lu_m_ini = mvpIScans.back()->v_R_mapping.back();
    }

    // add new
    if (mPara.reoptimize_fetch_raw_points)
    {
        int num_points = 0;
        for (int nScan = 0; nScan < mvTempScan.size(); nScan++)
        {
            pCurrentIScan->raw_non_ground_points.push_back(mvTempScan[nScan].non_ground_raw_points_);
            num_points += mvTempScan[nScan].non_ground_raw_points_.size();
        }
        fDebug << num_points << " non-ground points" << endl;
    }

    for(int nScan = 0; nScan < mvTempScan.size(); nScan ++)
    {
        mvTempScan[nScan].cleanRedundantInfo();
    }
    pCurrentIScan->indScans = mvTempScan; //@@@@@

    // trajectory
    pCurrentIScan->v_r_local = mv_r_lu_local; //numScan + 1 epoch of refined poses to local
    pCurrentIScan->v_R_local = mv_R_lu_local;

    // for each tree
    vector<vector<PointType>> vTreePointsRaw;
    vector<CylinderPara> vTreeParamsLocal; //tree parameters
    vector<bool> vValidTrees; //whether included in the optimization

    for(int nT = 0; nT < mvTreePointsRaw.size(); nT++)
    {
        //threshold on number of points 
        if(mvTreePointsRaw[nT].size() > mPara.minNumPtsPerIscanTree)
        {
            vTreePointsRaw.push_back(mvTreePointsRaw[nT]);
            vTreeParamsLocal.push_back(mvTreeParamsLocal[nT]);
            vValidTrees.push_back(mvValidTrees[nT]);
        }
    }
    pCurrentIScan->vTreePointRaw = vTreePointsRaw; //  points in lu frame with scan info
    pCurrentIScan->vTreeParamLocal = vTreeParamsLocal; //parameters in local frame
    //pCurrentIScan->vValidTree = vValidTrees;   //flag showing whether it is used in optimization or not
    vector<int> treeStatus(vTreeParamsLocal.size(), 0);
    pCurrentIScan->vTreeStatus = treeStatus;
    pCurrentIScan->raw_tree_points = raw_tree_points_per_iscan;

    // ground point
    pCurrentIScan->vGroundPointRaw = mvGroundPointsRaw; // points in lu frame with scan info

    //compute tree/ground points in raw to local frame
    pCurrentIScan->transform2Start(mPara.ground_downsample_distance_scan_integration);

    
    //pCurrentIScan->updateTreeParamLocal();
    fMapLog << "Trees after removing insufficient ones: " <<  mvTreePointsRaw.size() << " -> " << vTreePointsRaw.size() << " (threshold: " << mPara.minNumPtsPerIscanTree <<")" <<endl;
    fMapLog << "Ground points after downsampling: " << pCurrentIScan->vGroundPointRaw.size() << " -> " << pCurrentIScan->pGroundPointLocalDs->size() << endl;

     //CONTROL
    if (mvTreePointsRaw.size() < Min_Num_Tree_)
    {
        f_mapping_debug << "Iscan " << mnCount <<": \tNumber of trees are too few: " << mvTreePointsRaw.size()  << endl;

        Flag_Scan_Inter = false;
    }
}

/**
 * @brief Localizes the current integrated scan using tree points with respect to map trees.
 *
 * Uses Iterative Closest Point (ICP) between tree centers in current scan and map trees
 * to update transformation from local frame to mapping frame.
 */
void Mapping::Localize()
{
    //current location in the mapping is pCurrentIScan->r_lu_m_ini 
    fMapLog << "Inital lu to mapping " << pCurrentIScan->r_lu_m_ini.transpose() << "\t" << rad2deg(Find_Rotation(pCurrentIScan->R_lu_m_ini)).transpose() << endl;

    pCurrentIScan->computeMapTreeParameter(IntegratedScan::INIT);

    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_tree_Iscan (new pcl::PointCloud<pcl::PointXYZ>);
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_tree_map (new pcl::PointCloud<pcl::PointXYZ>);

    pcl::PointXYZ p;
    vector<Eigen::Vector3d> vScanTreeCenterMapping;
    for (int nT = 0; nT < pCurrentIScan->vTreeParamLocal.size(); nT++)
    {

        vScanTreeCenterMapping.push_back(pCurrentIScan->vTreeParamMapping[nT].x);
        vScanTreeCenterMapping[nT](2) = 0.0;
        p.x = pCurrentIScan->vTreeParamMapping[nT].x(0);
        p.y = pCurrentIScan->vTreeParamMapping[nT].x(1);
        p.z = 0.0;
        cloud_tree_Iscan->push_back(p);
    }

    vector<Eigen::Vector3d> v_map_neighboring_tree;
    for (int nMapT = 0; nMapT < mpMapTree.size(); nMapT++)
    {
        Eigen::Vector3d tree_loc = mpMapTree[nMapT]->para.x;
        if((tree_loc.head(2) - pCurrentIScan->r_lu_m_ini.head(2)).norm() < 40.0)
        {
            v_map_neighboring_tree.push_back(tree_loc);
            p.x = tree_loc(0);
            p.y = tree_loc(1);
            p.z = 0.0;
            cloud_tree_map->push_back(p);
        }
    }


    std::string outLocalization(output_folder + "LocalizationResult" + ".txt");
    std::ofstream fLocalize(outLocalization, std::ifstream::out);
    fLocalize << fixed << std::setprecision(6);

    for (int i = 0; i < cloud_tree_Iscan->size(); i++)
    {
        fLocalize << "1\t" << cloud_tree_Iscan->points[i].x + global_const_shift_(0) << "\t" << cloud_tree_Iscan->points[i].y+ global_const_shift_(1) << "\t" 
        << cloud_tree_Iscan->points[i].z + global_const_shift_(2)+10.0<< endl;
    }
    for (int i = 0; i < cloud_tree_map->size(); i++)
    {
        fLocalize << "2\t" << cloud_tree_map->points[i].x + global_const_shift_(0)<< "\t" << cloud_tree_map->points[i].y + global_const_shift_(1)<< "\t" 
        << cloud_tree_map->points[i].z + global_const_shift_(2)+10.0 << endl;
    }


    pcl::IterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ> icp;
    icp.setInputSource(cloud_tree_Iscan);
    icp.setInputTarget(cloud_tree_map);
    icp.setMaxCorrespondenceDistance(5.0) ;//2.0

    pcl::PointCloud<pcl::PointXYZ> Final;
    icp.align(Final);
    
    fMapLog << "Number of Iscan tree / neighboring map tree: " << cloud_tree_Iscan->size() << " " << cloud_tree_map->size() <<endl;
    fMapLog << "has converged:" << icp.hasConverged() << " score: " << icp.getFitnessScore() << std::endl; 
    if(!icp.hasConverged() )
            throw std::runtime_error("Wrong initialization");

    fMapLog << icp.getFinalTransformation() << std::endl;

    for (int i = 0; i < Final.size(); i++)
    {
        fLocalize << "3\t" << Final.points[i].x + global_const_shift_(0) << "\t" << Final.points[i].y + global_const_shift_(1)<< "\t" 
        << Final.points[i].z + global_const_shift_(2)+10.0 << endl;
    }
    fLocalize.close();

    Eigen::Matrix3d R_update = icp.getFinalTransformation().block<3, 3>(0, 0).cast<double>();
    Eigen::Vector3d r_update = icp.getFinalTransformation().col(3).head(3).cast<double>();
    // fMapLog << R_update << std::endl;
    // fMapLog << r_update << std::endl;
    Eigen::Vector3d angles;
    angles = Find_Rotation(R_update);
    Eigen::Matrix3d orthoR;
    Compute_Rotation(0.0, 0.0, angles(2), orthoR);

    Eigen::Matrix3d R_lu_m_updated = orthoR * pCurrentIScan->R_lu_m_ini;
    Eigen::Vector3d r_lu_m_updated = orthoR * pCurrentIScan->r_lu_m_ini + r_update;
    pCurrentIScan->R_lu_m_ini = R_lu_m_updated;
    pCurrentIScan->r_lu_m_ini = r_lu_m_updated;
    //fMapLog << "Updated lu to mapping" << pCurrentIScan->r_lu_m_ini.transpose() << "\t" << rad2deg(Find_Rotation(pCurrentIScan->R_lu_m_ini)).transpose() << endl;

    //use ground points to estimate dZ
    pCurrentIScan->computeMapGroundPoints(IntegratedScan::INIT);
    // pcl::octree::OctreePointCloudSearch<PointType> octreeGroundPointsFromMap(0.2);
    // octreeGroundPointsFromMap.setInputCloud(mpMapGroundPoint);
    // octreeGroundPointsFromMap.addPointsFromInputCloud();
    std::vector<int> pointSearchInd;
    std::vector<float> pointSearchSqDis;
    float Z_dif_mean = 0.0;
    int count = 0;
    for (int nP = 0; nP < pCurrentIScan->pGroundPointMappingDs->size(); nP++)
    {
        PointType p = pCurrentIScan->pGroundPointMappingDs->points[nP];

        if (octreeGroundPointsFromMap.nearestKSearch(p, 1, pointSearchInd, pointSearchSqDis) > 0)
        {
            Z_dif_mean += mpMapGroundPoint->points[pointSearchInd[0]].z - pCurrentIScan->pGroundPointMappingDs->points[nP].z;
            count++;
        }
    }
    Z_dif_mean /= float(count);
    pCurrentIScan->r_lu_m_ini(2) = pCurrentIScan->r_lu_m_ini(2) + Z_dif_mean;
    fMapLog << "Updated lu to mapping" << pCurrentIScan->r_lu_m_ini.transpose() << "\t" << rad2deg(Find_Rotation(pCurrentIScan->R_lu_m_ini)).transpose() << endl;

    //system("read -p 'Press Enter to continue...' var");
}

/**
 * @brief Localizes the current integrated scan using only ground points (DTM-based).
 *
 * Similar to Localize(), but uses downsampled ground points and ICP for alignment.
 * Only z-axis adjustment is performed for simplicity.
 */
void Mapping::LocalizeDTMOnly()
{
    //current location in the mapping is pCurrentIScan->r_lu_m_ini 
    fMapLog << "Inital lu to mapping " << pCurrentIScan->r_lu_m_ini.transpose() << "\t" << rad2deg(Find_Rotation(pCurrentIScan->R_lu_m_ini)).transpose() << endl;

    //use ground points to estimate transformation
    double dis_threshold = 30*30;
    pCurrentIScan->computeMapGroundPoints(IntegratedScan::INIT);
    std::vector<int> pointSearchInd;
    std::vector<float> pointSearchSqDis;
    float Z_dif_mean = 0.0;
    int count = 0;
    pcl::PointXYZ pcl_p;
    pcl::PointCloud<pcl::PointXYZ>::Ptr map_ground (new pcl::PointCloud<pcl::PointXYZ>);
    pcl::PointCloud<pcl::PointXYZ>::Ptr current_scan_ground (new pcl::PointCloud<pcl::PointXYZ>);
    for (int nP = 0; nP < pCurrentIScan->pGroundPointMappingDs->size(); nP++)
    {
        PointType p = pCurrentIScan->pGroundPointMappingDs->points[nP];

        if (octreeGroundPointsFromMap.nearestKSearch(p, 1, pointSearchInd, pointSearchSqDis) > 0)
        {
            if(pointSearchSqDis[0] > dis_threshold)//if dis >30, then continue
            {
                continue;
            }
            pcl_p.x = mpMapGroundPoint->points[pointSearchInd[0]].x;
            pcl_p.y = mpMapGroundPoint->points[pointSearchInd[0]].y;
            pcl_p.z = mpMapGroundPoint->points[pointSearchInd[0]].z;
            map_ground->push_back(pcl_p);

            pcl_p.x = pCurrentIScan->pGroundPointMappingDs->points[nP].x;
            pcl_p.y = pCurrentIScan->pGroundPointMappingDs->points[nP].y;
            pcl_p.z = pCurrentIScan->pGroundPointMappingDs->points[nP].z;
            current_scan_ground->push_back(pcl_p);

            count++;
        }
    }

    pcl::IterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ> icp;
    icp.setInputSource(current_scan_ground);
    icp.setInputTarget(map_ground);
    icp.setMaxCorrespondenceDistance(5.0) ;//2.0

    pcl::PointCloud<pcl::PointXYZ> Final;
    icp.align(Final);
    
    fMapLog << "Number of current_scan_ground / map_ground: " << current_scan_ground->size() << " " << map_ground->size() <<endl;
    fMapLog << "has converged:" << icp.hasConverged() << " score: " << icp.getFitnessScore() << std::endl; 
    if(!icp.hasConverged() )
            throw std::runtime_error("Wrong initialization");

    fMapLog << icp.getFinalTransformation() << std::endl;

    Eigen::Matrix3d R_update = icp.getFinalTransformation().block<3, 3>(0, 0).cast<double>();
    Eigen::Vector3d r_update = icp.getFinalTransformation().col(3).head(3).cast<double>();

    std::string outLocalization(output_folder + "DTM_only_LocalizationResult" + ".txt");
    std::ofstream fLocalize(outLocalization, std::ifstream::out);
    fLocalize << fixed << std::setprecision(6);

    for (int i = 0; i < current_scan_ground->size(); i++)
    {
        fLocalize << "1\t" << map_ground->points[i].x + global_const_shift_(0) << "\t" << map_ground->points[i].y+ global_const_shift_(1) << "\t" 
        << map_ground->points[i].z + global_const_shift_(2)+10.0<< endl;

        Eigen::Vector3d transformed_p;
        transformed_p = R_update*current_scan_ground->points[i].getVector3fMap().cast<double>()+r_update;
        fLocalize << "2\t" << transformed_p(0) + global_const_shift_(0)<< "\t" << transformed_p(1) + global_const_shift_(1)<< "\t" 
        << transformed_p(2) + global_const_shift_(2)+10.0 << endl;
    }
    fLocalize.close();

    Eigen::Vector3d angles;
    angles = Find_Rotation(R_update);
    Eigen::Matrix3d orthoR;
    Compute_Rotation(0.0, 0.0, angles(2), orthoR);

    fMapLog << "icp results rotation"<<angles.transpose() << std::endl;
    fMapLog << "icp results translation"<< r_update.transpose() << std::endl;

    std::cout << "icp results rotation"<<angles.transpose() << std::endl;
    std::cout << "icp results translation"<< r_update.transpose() << std::endl;

    Eigen::Matrix3d R_lu_m_updated = orthoR * pCurrentIScan->R_lu_m_ini;
    Eigen::Vector3d r_lu_m_updated = orthoR * pCurrentIScan->r_lu_m_ini + r_update;
    pCurrentIScan->R_lu_m_ini = R_lu_m_updated;
    pCurrentIScan->r_lu_m_ini = r_lu_m_updated;

    fMapLog << "Updated lu to mapping" << pCurrentIScan->r_lu_m_ini.transpose() << "\t" << rad2deg(Find_Rotation(pCurrentIScan->R_lu_m_ini)).transpose() << endl;

    //system("read -p 'Press Enter to continue...' var");
}

/********************************************************
 * Clear information used in the computation.
 ********************************************************/
/**
 * @brief Clears all integration-related buffers and temporary scan data.
 */
void Mapping::resetIntegration()
{
    vector<ScanInfo>().swap(mvTempScan);
    vector<Eigen::Vector3d>().swap(mvIntegratedTree);
    vector<vector<pair<int, int>>>().swap(mvIntegratedTreeIds);
    vector<Eigen::Vector3f>().swap(mvTreeCentersLocal);
    vector<vector<PointType>>().swap(mvTreePointsRaw); // rIlu. Filtered tree points. Will be updated in the optimization

    vector<bool>().swap(mvValidTrees);              // whether included in the optimization
    vector<CylinderPara>().swap(mvTreeParamsLocal); // tree parameters

    std::vector<PointType>().swap(mvGroundPointsRaw);        // r_I_lu
    std::vector<PointType>().swap(mvGroundPointsLocal);      // r_I_local
    std::vector<pair<double, double>>().swap(mvGridCenters); // center of each grid in local frame
    vector<vector<int>>().swap(mvGridPointIndex);            // index of points in mvGroundPointsRaw
    vector<Eigen::Vector4d>().swap(mvPlaneParamLocal);

    vector<Eigen::Vector3d>().swap(mv_r_lu_local); // refined in the odometry
    vector<Eigen::Matrix3d>().swap(mv_R_lu_local);

    vector<vector<vector<PointType>>>().swap(raw_tree_points_per_iscan);
	tree_id_in_iscan_to_location_index_.clear();
	tree_id_in_iscan_to_vector_index_.clear();

    mNumPlanarPatch = 0;
	// running_iscan_index++;
}

/**
 * @brief Performs Iterative Closest Point (ICP) between map and current scan tree points.
 *
 * @param map_target List of target tree centers from map
 * @param current_scan_source List of tree centers from current scan
 * @param transformation Output rotation matrix
 * @param translation Output translation vector
 * @return true if ICP converged, false otherwise
 */
bool ICPforMapandCurrentScan(const vector<Eigen::Vector3d>& map_target,const vector<Eigen::Vector3d>& current_scan_source,
    Eigen::Matrix3d& transformation,Eigen::Vector3d& translation)
{
    // Convert point clouds to PCL point clouds
    pcl::PointCloud<pcl::PointXYZ> pcl_target;
    pcl::PointCloud<pcl::PointXYZ> pcl_source;
    for (const auto& point : map_target)
        pcl_target.push_back(pcl::PointXYZ(point(0),point(1),point(2)));
    for (const auto& point : current_scan_source)
        pcl_source.push_back(pcl::PointXYZ(point(0), point(1),point(2)));


    pcl::IterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ> icp;
	icp.setInputSource(pcl_source.makeShared());
	icp.setInputTarget(pcl_target.makeShared());
	icp.setMaxCorrespondenceDistance(1.0);
	icp.setMaximumIterations(50);

	pcl::PointCloud<pcl::PointXYZ> Final;
	icp.align(Final);

	if (!icp.hasConverged())
		return false;

    // Print transformation matrix
    transformation = icp.getFinalTransformation().block<3, 3>(0, 0).cast<double>();
    translation = icp.getFinalTransformation().block<3, 1>(0, 3).cast<double>();
    std::cout << "Loop closure, Transformation matrix:\n" << transformation << std::endl;
    std::cout << "Loop closure, Translation vector:\n" << translation << std::endl;
    return true;
}

/**
 * @brief Determines if loop closure is likely by checking spatial and temporal conditions.
 *
 * If sufficient time and spatial overlap occur, ICP is performed to confirm loop closure.
 * @return true if loop closure is detected, false otherwise
 */
bool Mapping::DetermineLoopClosureIfExist(Eigen::Vector3d scan_location,int iscan_index, const vector<Eigen::Vector3d>& current_scan_tree_locations,
    Eigen::Matrix3d& transformation,Eigen::Vector3d& translation)
{
    // it will be cnducted after 5 mins
	int number_scans_in_5_mins = (5 * 60 / 0.1);
    if(mvpIScans.size()*mNumInitialScan<number_scans_in_5_mins)
    {
        return false;
    }
	// less than 30 m , revist interval>5 mins;
	bool revisit =false;
    for(int i =0;i<mvpIScans.size();i++)
    {
		double dis = (mvpIScans[i]->r_lu_m_ini - scan_location).norm();
		if (dis < 30)
		{
			if (fabs(mvpIScans[i]->index - iscan_index) > 30)
			{
				revisit = true;
			}
		}
    }
	if (!revisit)
	{
		return false;
	}
    double x_min = DBL_MAX,x_max = -DBL_MAX,y_min = DBL_MAX,y_max = -DBL_MAX;
    // calculate boudning box of current scan locations
    for(int i = 0; i<current_scan_tree_locations.size();i++)
    {
        x_min = std::min(x_min,current_scan_tree_locations[i](0));
        y_min = std::min(x_min,current_scan_tree_locations[i](1));
        x_max = std::max(x_max,current_scan_tree_locations[i](0));
        y_max = std::max(y_max,current_scan_tree_locations[i](1));
    }

    std::vector<Eigen::Vector3d> map_trees_in_box;

    for (int nMapT = 0; nMapT < mpMapTree.size(); nMapT++) // prev scan
    {
        if(mpMapTree[nMapT]->status == MapTree::DEACTIVATE)
            continue;
        Eigen::Vector3d mapCenter = mpMapTree[nMapT]->para.x;

        if(mapCenter(0)>x_min && mapCenter(0)<x_max && mapCenter(1)>y_min && mapCenter(1)<y_max)
        {
            map_trees_in_box.push_back(mapCenter);
        }
    }
    // Eigen::Matrix3d transformation,Eigen::Vector3d translation;
    if(!ICPforMapandCurrentScan(map_trees_in_box,current_scan_tree_locations,transformation,translation))
    {
        return false;
    }
    return true;
}

/********************************************************
 * Main function for transform current integrated scan to map
 * Refine T_lu_m_ini for current scan
 * Output: R/r_local_m_updated
 ********************************************************/
/**
 * @brief Computes transformation from current scan to global map using tree alignment.
 *
 * Performs matching of trees in current scan and map trees, followed by optimization of
 * transformation using tree and ground points. Updates poses and cylinder parameters.
 */
void Mapping::computePoseToMap()
{
    Flag_Pose_to_Map = true;
    fMapLog << "--------------------------------------------------" << endl;
    fMapLog << "2. Compute Trans from Integrated Scan to Mapping" << endl;

    // ----------------------------------------------------------------
    // 1. compute the trans between local and mapping.
    // R_lu_m = R_l_m * R_lu_l.    r_lu_m = r_l_m + R_l_m * r_lu_l
    pCurrentIScan->R_local_m_ini = pCurrentIScan->R_lu_m_ini * pCurrentIScan->v_R_local[0].transpose();
    pCurrentIScan->r_local_m_ini = pCurrentIScan->r_lu_m_ini - pCurrentIScan->R_local_m_ini * pCurrentIScan->v_r_local[0];

    //in case the global map is available, T_local_m_ini is not reliable 
    if (!mbInit && global_map_flag_)
    {   
        //refine T_local_m_ini
        if(DTM_only_flag_)
        {
            LocalizeDTMOnly();
        }
        else
        {
            Localize();
        }
    }

    pCurrentIScan->R_local_m_ini = pCurrentIScan->R_lu_m_ini * pCurrentIScan->v_R_local[0].transpose();
    pCurrentIScan->r_local_m_ini = pCurrentIScan->r_lu_m_ini - pCurrentIScan->R_local_m_ini * pCurrentIScan->v_r_local[0];
    fMapLog << "Initial trans from local to mapping " << pCurrentIScan->r_local_m_ini.transpose() << "\t" << rad2deg(Find_Rotation(pCurrentIScan->R_local_m_ini)).transpose() << endl;
    // Initial value for the local to mapping to be refined
    pCurrentIScan->R_local_m_updated = pCurrentIScan->R_local_m_ini;
    pCurrentIScan->r_local_m_updated = pCurrentIScan->r_local_m_ini;

    // ----------------------------------------------------------------
    // If initialized or global map, conduct optimization
    if (mbInit || global_map_flag_)
    {
        int iter = 0;
        iter_count = 0;
        while(iter < 2)
        {
            // ----------------------------------------------------------------
            // 2. refine the T_l_m. update the r_lu_m_updated & R_lu_m_updated
            // 2.1 estimate location of tree using the initial r_l_m;
            pCurrentIScan->computeMapTreeParameter(IntegratedScan::REFINED);
            pCurrentIScan->computeMapTreePoints(IntegratedScan::REFINED);

            // 2.2 match trees from this individual scan to map tree
            int numTree = pCurrentIScan->vTreeParamLocal.size();
            vector<double> vDistanceThreshold;
            vector<Eigen::Vector3d> vScanTreeCenterMapping;
            for (int nT = 0; nT < numTree; nT++)
            {
                // similarly, threshold is a function of tree to lidar unit distance
                vDistanceThreshold.push_back(max(2.0/double(iter+1), 0.05 * pCurrentIScan->vTreeParamLocal[nT].x.norm()));
                //vDistanceThreshold.push_back(max(2.0, 0.05 * pCurrentIScan->vTreeParamLocal[nT].x.norm()));

                vScanTreeCenterMapping.push_back(pCurrentIScan->vTreeParamMapping[nT].x);
            }
			// //(TODO: real loop closure)
			// Eigen::Matrix3d transformation; Eigen::Vector3d translation;
			// bool revisit = DetermineLoopClosureIfExist(pCurrentIScan->r_local_m_updated, pCurrentIScan->index, vScanTreeCenterMapping, transformation, translation);
			// if (revisit)
			// {
			// 	pCurrentIScan->R_local_m_updated = transformation*pCurrentIScan->R_local_m_updated;
			// 	pCurrentIScan->r_local_m_updated = translation+ transformation*pCurrentIScan->r_local_m_updated;
			// 	pCurrentIScan->computeMapTreeParameter(IntegratedScan::REFINED);
			// 	pCurrentIScan->computeMapTreePoints(IntegratedScan::REFINED);
			// 	// 2.2 match trees from this individual scan to map tree
			// 	int numTree = pCurrentIScan->vTreeParamLocal.size();
			// 	vector<Eigen::Vector3d> vScanTreeCenterMapping_update;
			// 	for (int nT = 0; nT < numTree; nT++)
			// 	{
			// 		vScanTreeCenterMapping_update.push_back(pCurrentIScan->vTreeParamMapping[nT].x);
			// 	}
			// 	vScanTreeCenterMapping.swap(vScanTreeCenterMapping_update);
			// }

            std::vector<pair<int, int>> pairs; // map tree id & scan tree id
            matchTreeScantoMap(&vScanTreeCenterMapping, &vDistanceThreshold, pairs);
            fMapLog << "Number of tree matches between current Iscan to map: " << pairs.size() << endl;
            mvPairs = pairs;
            // some outlier removal or random search (Todo)

            // 2.3 conduct an optimization
            TicToc t_optimize;
            optimizeIScantoMap();
            fMapLog << "\tTime for this iteration: " << t_optimize.toc() << endl;
            iter++;
            if(!Flag_Pose_to_Map)
                break;
        }
    }
    
    //if not enough matches are found
    if(!Flag_Pose_to_Map)
    {
        fMapLog << "\tFailed in Iscan to map, reset transformation value to initial " <<endl;
        pCurrentIScan->R_local_m_updated = pCurrentIScan->R_local_m_ini;
        pCurrentIScan->r_local_m_updated = pCurrentIScan->r_local_m_ini;
        pCurrentIScan->flag_iscan_to_map =false;
    }

    // compute the T_lu_m individual scan
    pCurrentIScan->computeIndividualPoseMapping();

    // compute reference trajectory
    if (mbTraj)
        computeRefTraj();

    // compute the transformation between current Iscan and previous one
    if (mvpIScans.size() > 1)
    {
        //Eigen::Matrix3d R_local_t2_local_t1;
        //Eigen::Vector3d r_local_t2_local_t1;
        Eigen::Matrix3d R_local_m_t1 = mvpIScans[mvpIScans.size() - 2]->R_local_m_updated;
        Eigen::Matrix3d R_local_m_t2 = pCurrentIScan->R_local_m_updated;
        Eigen::Vector3d r_local_m_t1 = mvpIScans[mvpIScans.size() - 2]->r_local_m_updated;
        Eigen::Vector3d r_local_m_t2 = pCurrentIScan->r_local_m_updated;


        if (mPara.odo_from_trajectory_flag)
        {
            R_local_m_t1 = mvpIScans[mvpIScans.size() - 2]->v_R_mapping_ref[0] * mvpIScans[mvpIScans.size() - 2]->v_R_local[0].inverse();
            r_local_m_t1 = mvpIScans[mvpIScans.size() - 2]->v_r_mapping_ref[0] - R_local_m_t1 * mvpIScans[mvpIScans.size() - 2]->v_r_local[0];

            R_local_m_t2 = pCurrentIScan->v_R_mapping_ref[0] * pCurrentIScan->v_R_local[0].inverse();
            r_local_m_t2 = pCurrentIScan->v_r_mapping_ref[0] - R_local_m_t2 * pCurrentIScan->v_r_local[0];

        }

        //constraint in iscan level long term optimization
        pCurrentIScan->R_local_t2_local_t1 = R_local_m_t1.inverse() * R_local_m_t2;
        pCurrentIScan->r_local_t2_local_t1 = R_local_m_t1.inverse() * (r_local_m_t2 - r_local_m_t1);

    }
    // export trajectory in mapping
    int startIndex = mvpIScans.size() == 1 ? 0 : mvpIScans[mvpIScans.size() - 2]->indScans.back().scanID;
    for (int i = 0; i < pCurrentIScan->v_R_mapping.size(); i++)
    {
        int index = i == 0 ? startIndex : pCurrentIScan->indScans[i - 1].scanID;
        fMappingTraj << i << "\t" << pCurrentIScan->v_r_mapping[i].transpose() << "\t"
                     << "0\t" << index << "\t" << rad2deg(Find_Rotation(pCurrentIScan->v_R_mapping[i])).transpose() << endl;
        if (mbTraj)
        {
            fMappingTraj << i << "\t" << pCurrentIScan->v_r_mapping_ref[i].transpose() << "\t"
                         << "1\t" << index << "\t" << rad2deg(Find_Rotation(pCurrentIScan->v_R_mapping_ref[i])).transpose() << endl;
        }
    }
}

/********************************************************
 * Function to refine the pose of the integrated scan
 *  T_l_m_ini -> T_l_m_refined
 *  Updated: mvPairs,mpMapTree[].para
 ********************************************************/
/**
 * @brief Optimizes alignment between the current Iscan and map tree points.
 *
 * Uses Ceres optimization to refine local-to-map transformation and tree cylinder parameters
 * by minimizing residuals between observed and model tree/ground points.
 * @return true if optimization succeeds, false if insufficient tree matches are found
 */
bool Mapping::optimizeIScantoMap()
{
    double tree_std = mPara.mapTreeStd;
    double ground_std = mPara.mapGroundStd;
    double global_tree_loc_std = mPara.global_tree_loc_std;
    int max_num_map_point_per_tree = 1000;

    int numPair = mvPairs.size();
    cout << "Tree pairs: " << numPair << ", ground points: " << pCurrentIScan->pGroundPointLocalDs->size() << endl;

    // pcl::octree::OctreePointCloudSearch<PointType> octreeGroundPointsFromMap(0.2);
    if (!global_map_flag_)
    {
        // build octree for ground points
        TicToc t_build_tree;
        octreeGroundPointsFromMap.deleteTree();
        octreeGroundPointsFromMap.setInputCloud(mpMapGroundPoint);
        octreeGroundPointsFromMap.addPointsFromInputCloud();
        fMapLog << mpMapGroundPoint->size() << "\t" << octreeGroundPointsFromMap.getLeafCount() << endl;

        fMapLog << ", Build Tree time: " << t_build_tree.toc() << endl;
    }

    int iter = 0;
    int maxIter = 1;
    while (iter < maxIter)
    {
        fMapLog << "Iter " << ++iter <<flush;

        // R_local_m_updated is updated every iteration
        pCurrentIScan->computeMapGroundPoints(IntegratedScan::REFINED);

        // parameter to be estimated. before estimated
        Eigen::Quaterniond q_ini(pCurrentIScan->R_local_m_updated);
        Eigen::Vector3d t_ini = pCurrentIScan->r_local_m_updated;

        // initialize parameters
        double para_q[4] = {q_ini.x(), q_ini.y(), q_ini.z(), q_ini.w()};
        double para_t[3] = {t_ini[0], t_ini[1], t_ini[2]};

        Eigen::Map<Eigen::Quaterniond> q_last_curr(para_q); // requires array in [x, y, z, w]
        Eigen::Map<Eigen::Vector3d> t_last_curr(para_t);

        cout << "Initial trans from local to mapping: " << t_last_curr.transpose() << "\t" << rad2deg(Find_Rotation(q_last_curr.toRotationMatrix())).transpose() << endl;
        //fMapLog << "Initial trans from local to mapping: " << t_last_curr.transpose() << "\t" << rad2deg(Find_Rotation(q_last_curr.toRotationMatrix())).transpose() << endl;

        int numCylinder = mvPairs.size();
        int numCylinderUnknown = numCylinder * CylinderBlockSize;
        double *cylinders = new double[numCylinderUnknown];
        double *cylinders_std = new double[numCylinderUnknown];
        int numUnknown = 1 * (PosBlockSize + OriBlockSize) + numCylinderUnknown; // all unknowns

        // 2. cylinder parameter
        int countCylinder = 0;
        for (int nPair = 0; nPair < numPair; nPair++)
        {
            int mapTreeId = mvPairs[nPair].first;

            *(cylinders + countCylinder * CylinderBlockSize + 0) = mpMapTree[mapTreeId]->para.x(0); // XYZ (Z is nominal)
            *(cylinders + countCylinder * CylinderBlockSize + 1) = mpMapTree[mapTreeId]->para.x(1);
            *(cylinders + countCylinder * CylinderBlockSize + 2) = mpMapTree[mapTreeId]->para.x(2);
            *(cylinders + countCylinder * CylinderBlockSize + 3) = mpMapTree[mapTreeId]->para.n(0);
            *(cylinders + countCylinder * CylinderBlockSize + 4) = mpMapTree[mapTreeId]->para.n(1);
            *(cylinders + countCylinder * CylinderBlockSize + 5) = mpMapTree[mapTreeId]->para.n(2);
            *(cylinders + countCylinder * CylinderBlockSize + 6) = mpMapTree[mapTreeId]->para.r; 

            *(cylinders_std + countCylinder * CylinderBlockSize + 0) = 1.0E+20; // XYZ
            *(cylinders_std + countCylinder * CylinderBlockSize + 1) = 1.0E+20;
            *(cylinders_std + countCylinder * CylinderBlockSize + 2) = 1.0E-20;
            *(cylinders_std + countCylinder * CylinderBlockSize + 3) = 1.0E-20; // ux, uy, uz
            *(cylinders_std + countCylinder * CylinderBlockSize + 4) = 1.0E-20;
            *(cylinders_std + countCylinder * CylinderBlockSize + 5) = 1.0E-20;
            *(cylinders_std + countCylinder * CylinderBlockSize + 6) = 1.0E-20; // r

            if (mPara.bMapTreeDirec)
            {
                *(cylinders_std + countCylinder * CylinderBlockSize + 3) = 1.0E+20; // ux, uy
                *(cylinders_std + countCylinder * CylinderBlockSize + 4) = 1.0E+20;
            }
            if (mPara.bMapTreeRadius)
                *(cylinders_std + countCylinder * CylinderBlockSize + 6) = 1.0E+20; // r

            countCylinder++;
        }
        int numFixedUnknown = 0;
        for (int i = 0; i < numCylinderUnknown; i++)
        {
            double temp;
            temp = *(cylinders_std + i);
            if (temp < SigmaFixed)
                numFixedUnknown += 1;
        }

        // build ceres
        ceres::LossFunction *loss_function = new ceres::HuberLoss(0.1);
        ceres::LocalParameterization *q_parameterization = new ceres::EigenQuaternionParameterization(); // x, y, z, w
	    ceres::Problem problem;//= new ceres::Problem;
        ceres::CostFunction *costFunction;

        problem.AddParameterBlock(para_q, 4, q_parameterization);
        problem.AddParameterBlock(para_t, 3);
        for (int nCylinder = 0; nCylinder < numCylinder; nCylinder++)
        {
            problem.AddParameterBlock(cylinders + nCylinder * CylinderBlockSize, CylinderBlockSize);
        }

        //fMapLog << "\tNumber of unknowns: " << numUnknown << endl;
        //fMapLog << "\tNumber of fixed unknowns: " << numFixedUnknown << endl;

        //-------------------------add observation -------------------------------
        int numMapTreeObs = 0, numScanTreeObs = 0;
        for (int nPair = 0; nPair < numPair; nPair++)
        {
            int cylinderId = nPair;
            int mapTreeId = mvPairs[nPair].first;
            int scanTreeId = mvPairs[nPair].second;
            double *cylinder = cylinders + cylinderId * CylinderBlockSize;

            // two types of observation
            // 1. Tree points from map
            int increment = 1;
            if(mpMapTree[mapTreeId]->numPoint != 0)
            {
                increment = ceil(mpMapTree[mapTreeId]->numPoint/max_num_map_point_per_tree);
                increment = max(1, increment);
                //fDebug << increment << " " << mpMapTree[mapTreeId]->numPoint << " " << max_num_map_point_per_tree <<endl;
            }
            for (int nT = 0; nT < mpMapTree[mapTreeId]->visibleInfo.size(); nT++)
            {
                int iscanId = mpMapTree[mapTreeId]->visibleInfo[nT].first;
                int treeId = mpMapTree[mapTreeId]->visibleInfo[nT].second;
                for (int nP = 0; nP < mvpIScans[iscanId]->vTreePointMapping[treeId].size(); nP = nP + increment)
                {
                    // points in mapping
                    Eigen::Vector3d curr_point(mvpIScans[iscanId]->vTreePointMapping[treeId][nP].x,
                                               mvpIScans[iscanId]->vTreePointMapping[treeId][nP].y,
                                               mvpIScans[iscanId]->vTreePointMapping[treeId][nP].z);
                    
                    //apply the scale to std
                    costFunction = CylinderFactorObjectPoint::Create(curr_point, tree_std / sqrt(increment));
                    problem.AddResidualBlock(costFunction, loss_function, cylinder);
                    numMapTreeObs++;
                }
            }

            // 2. Tree points from current integrated scan
            for (int nP = 0; nP < pCurrentIScan->vTreePointMapping[scanTreeId].size(); nP++)
            {
                // points in local ***
                Eigen::Vector3d curr_point(pCurrentIScan->vTreePointLocal[scanTreeId][nP].x,
                                           pCurrentIScan->vTreePointLocal[scanTreeId][nP].y,
                                           pCurrentIScan->vTreePointLocal[scanTreeId][nP].z);

                costFunction = CylinderFactorOneEpoch::Create(curr_point, tree_std);
                problem.AddResidualBlock(costFunction, loss_function, para_t, para_q, cylinder);
                numScanTreeObs++;
            }
        }

        int numGroundObs = 0;
        std::vector<int> pointSearchInd;
        std::vector<float> pointSearchSqDis;
        for (int nP = 0; nP < pCurrentIScan->pGroundPointMappingDs->size(); nP++)
        {
            PointType p = pCurrentIScan->pGroundPointMappingDs->points[nP];
            if (octreeGroundPointsFromMap.nearestKSearch(p, 5, pointSearchInd, pointSearchSqDis) > 0)
            {
                Eigen::Matrix<double, 5, 3> matA0;
                Eigen::Matrix<double, 5, 1> matB0 = -1 * Eigen::Matrix<double, 5, 1>::Ones();
                if (pointSearchSqDis[4] < 1.0)
                {
                    for (int j = 0; j < 5; j++)
                    {
                        matA0(j, 0) = mpMapGroundPoint->points[pointSearchInd[j]].x;
                        matA0(j, 1) = mpMapGroundPoint->points[pointSearchInd[j]].y;
                        matA0(j, 2) = mpMapGroundPoint->points[pointSearchInd[j]].z;
                    }
                    // find the norm of plane
                    Eigen::Vector3d norm = matA0.colPivHouseholderQr().solve(matB0);
                    double negative_OA_dot_norm = 1 / norm.norm();
                    norm.normalize();

                    // Here n(pa, pb, pc) is unit norm of plane
                    bool planeValid = true;
                    for (int j = 0; j < 5; j++)
                    {
                        // if OX * n > 0.2, then plane is not fit well
                        if (fabs(norm(0) * mpMapGroundPoint->points[pointSearchInd[j]].x +
                                 norm(1) * mpMapGroundPoint->points[pointSearchInd[j]].y +
                                 norm(2) * mpMapGroundPoint->points[pointSearchInd[j]].z + negative_OA_dot_norm) > 0.2)
                        {
                            planeValid = false;
                            break;
                        }
                    }
                    if (planeValid)
                    {
                        Eigen::Vector4d params;
                        params << norm, negative_OA_dot_norm;
                        Eigen::Vector3d curr_point(pCurrentIScan->pGroundPointLocalDs->points[nP].x,
                                                   pCurrentIScan->pGroundPointLocalDs->points[nP].y,
                                                   pCurrentIScan->pGroundPointLocalDs->points[nP].z);

                        costFunction = PlaneFactorOneEpochFixPlane::Create(curr_point, params, ground_std);
                        problem.AddResidualBlock(costFunction, loss_function, para_t, para_q);
                        numGroundObs++;
                    }
                }
            }
        }

        //fDebug << "REference control" <<endl;
        // for tree from global map, add constraint in the X and Y coordinates
        int num_global_tree = 0;
        for (int nPair = 0; nPair < numPair; nPair++)
        {
            int mapTreeId = mvPairs[nPair].first;
            int cylinderId = nPair;

            if (mpMapTree[mapTreeId]->map_tree_flag)
            {
                double *cylinder = cylinders + cylinderId * CylinderBlockSize;
                //costFunction = CylinderPriorConstraints::Create(mpMapTree[mapTreeId]->para_ref_.x, global_tree_loc_std);
                costFunction = CylinderPriorConstraintsPointToLine::Create(mpMapTree[mapTreeId]->para_ref_.x, global_tree_loc_std);
                //costFunction = CylinderFactorObjectPoint::Create(mpMapTree[mapTreeId]->para_ref_.x, global_tree_loc_std);
                //Eigen::Vector3d temp_check{mpMapTree[mapTreeId]->para_ref_.x(0) + 0.1, mpMapTree[mapTreeId]->para_ref_.x(1) + 0.1, mpMapTree[mapTreeId]->para_ref_.x(2) + 0.1};
                //costFunction = CylinderPriorConstraintsPointToLine::Create(temp_check, global_tree_loc_std);
 
                problem.AddResidualBlock(costFunction, nullptr, cylinder);                
                num_global_tree++;

                // fDebug << num_global_tree << "\t" << mpMapTree[mapTreeId]->para_ref_.x(0) << "\t" << mpMapTree[mapTreeId]->para_ref_.x(1) << "\t" << mpMapTree[mapTreeId]->para_ref_.x(2)
                // << "\t" << global_tree_loc_std <<endl;
            }
        }

        //--------------------------Fix Parameters  ------------------------------
        Checkforconstantparams(&problem, cylinders, cylinders_std, CylinderBlockSize, numCylinder);

        // solve
        // ceres::Problem::EvaluateOptions evalop;
        // double initialCost;
        // std::vector<double> resBlocks_ini, resBlocks_out;
        // problem->Evaluate(evalop, &initialCost, &resBlocks_ini, NULL, 0);
        TicToc t_solver;
        ceres::Solver::Options solverOptions;
        solverOptions.max_num_iterations = 10;
        solverOptions.linear_solver_type = ceres::LinearSolverType::SPARSE_NORMAL_CHOLESKY; // ceres::LinearSolverType::SPARSE_NORMAL_CHOLESKY;
        solverOptions.minimizer_progress_to_stdout = false;
        ceres::Solver::Summary summary;
        ceres::Solve(solverOptions, &problem, &summary);
        fMapLog << ", Time: " << t_solver.toc() << endl;
        cout << "Time: " << t_solver.toc() << endl;
        // double finalCost;
        // problem->Evaluate(evalop, &finalCost, &resBlocks_out, NULL, 0);

        fMapLog << "\tTree obs from current scan/map: " << numScanTreeObs << " " << numMapTreeObs << "\tGround obs: " << numGroundObs << endl;
        fMapLog << "\tNumber of pair with global tree: " << num_global_tree<< endl;

        cout << "Tree obs from current scan/map: " << numScanTreeObs << " " << numMapTreeObs << "\tGround obs: " << numGroundObs << endl;

        //fMapLog << "\tRefined trans from local to mapping: " << t_last_curr.transpose() << "\t" << rad2deg(Find_Rotation(q_last_curr.toRotationMatrix())).transpose() << endl;
        //fMapLog << "\tNorm of the quaternion: " << q_ini.norm() << "\t" << q_last_curr.norm() <<endl;

        //******************** Update parameters **********************
        // 1. trajectory info
        pCurrentIScan->r_local_m_updated = t_last_curr;

        //condition the rotation matrix to avoid numerical issue
        Eigen::Matrix3d tempR, orthoR;
        tempR = q_last_curr.toRotationMatrix();
        Eigen::Vector3d angles;
        angles =Find_Rotation(tempR);
        Compute_Rotation(angles(0), angles(1), angles(2), orthoR);
        pCurrentIScan->R_local_m_updated = orthoR;
        fMapLog << "\tRefined trans from local to mapping: " << pCurrentIScan->r_local_m_updated.transpose() << "\t" << rad2deg(Find_Rotation(pCurrentIScan->R_local_m_updated)).transpose() << endl;

        // 2. cylinder parameter
        CylinderPara treePara;
        vector<CylinderPara> vTreeParamsUpdated;
        countCylinder = 0;
        for (int nPair = 0; nPair < numPair; nPair++)
        {
            int mapTreeId = mvPairs[nPair].first;
            // some criteria of map tree to define what parameter to estimate (TBD)
            treePara.x(0) = *(cylinders + countCylinder * CylinderBlockSize + 0);
            treePara.x(1) = *(cylinders + countCylinder * CylinderBlockSize + 1);
            treePara.x(2) = *(cylinders + countCylinder * CylinderBlockSize + 2);
            treePara.n(0) = *(cylinders + countCylinder * CylinderBlockSize + 3);
            treePara.n(1) = *(cylinders + countCylinder * CylinderBlockSize + 4);
            treePara.n(2) = *(cylinders + countCylinder * CylinderBlockSize + 5);
            treePara.r = *(cylinders + countCylinder * CylinderBlockSize + 6);
            countCylinder++;
            vTreeParamsUpdated.push_back(treePara);
        }

        //******************** Compute residual **********************
        double res_scantree = 0.0, res_scantree_ini = 0.0, res_maptree = 0.0, res_maptree_ini = 0.0, res_ground = 0.0, res_ground_ini = 0.0;
        double res_scantree_ceres = 0.0, res_scantree_ini_ceres = 0.0, res_maptree_ceres = 0.0, res_maptree_ini_ceres = 0.0, res_ground_ceres = 0.0, res_ground_ini_ceres = 0.0;
        int obsCount = 0;
        std::ofstream fFeatureIscanToMapBefore;
        std::ofstream fFeatureIscanToMapAfter;
        
        if(intermediate_result_flag)
        {
            std::string outFeatureIni(output_folder_iscan_map + to_string(pCurrentIScan->index) + "_IscanToMap_" + to_string(pCurrentIScan->indScans[0].scanID) + "_iter" + to_string(iter_count) + "_ini.txt");
            std::string outFeatureRef(output_folder_iscan_map + to_string(pCurrentIScan->index) + "_IscanToMap_" + to_string(pCurrentIScan->indScans[0].scanID) + "_iter" + to_string(iter_count) + "_ref.txt");
            fFeatureIscanToMapBefore.open(outFeatureIni);
            fFeatureIscanToMapAfter.open(outFeatureRef);
        }

        for (int nPair = 0; nPair < numPair; nPair++)
        {
            int cylinderId = nPair;
            int mapTreeId = mvPairs[nPair].first;
            int scanTreeId = mvPairs[nPair].second;
            // double *cylinder = cylinders + cylinderId * CylinderBlockSize;
            
            double res_per_Map_tree = 0.0;
            int num_point_map_tree = 0;
            // two types of observation
            // 1. Tree points from map
            for (int nT = 0; nT < mpMapTree[mapTreeId]->visibleInfo.size(); nT++)
            {
                int iscanId = mpMapTree[mapTreeId]->visibleInfo[nT].first;
                int treeId = mpMapTree[mapTreeId]->visibleInfo[nT].second;
                for (int nP = 0; nP < mvpIScans[iscanId]->vTreePointMapping[treeId].size(); nP++)
                {
                    // points in mapping
                    Eigen::Vector3d curr_point(mvpIScans[iscanId]->vTreePointMapping[treeId][nP].x,
                                               mvpIScans[iscanId]->vTreePointMapping[treeId][nP].y,
                                               mvpIScans[iscanId]->vTreePointMapping[treeId][nP].z);

                    double d_ini = computePoint2lineDistance(curr_point, mpMapTree[mapTreeId]->para.x, mpMapTree[mapTreeId]->para.x + mpMapTree[mapTreeId]->para.n);
                    double res_ini = d_ini - mpMapTree[mapTreeId]->para.r;
                    res_maptree_ini = res_maptree_ini + res_ini * res_ini;

                    double d_ref = computePoint2lineDistance(curr_point, vTreeParamsUpdated[cylinderId].x, vTreeParamsUpdated[cylinderId].x + vTreeParamsUpdated[cylinderId].n);
                    double res_ref = d_ref - vTreeParamsUpdated[cylinderId].r;
                    res_maptree = res_maptree + res_ref * res_ref;
                    
                    res_per_Map_tree = res_per_Map_tree + res_ref * res_ref;
                    num_point_map_tree++;

                    // res_maptree_ini_ceres = res_maptree_ini_ceres + resBlocks_ini[obsCount] * resBlocks_ini[obsCount];
                    // res_maptree_ceres = res_maptree_ceres + resBlocks_out[obsCount] * resBlocks_out[obsCount];
                    // obsCount++;

                    if(intermediate_result_flag)
                    {
                        fFeatureIscanToMapBefore<< nPair << "\t" << curr_point.transpose() <<"\t" << -1 <<endl;
                        fFeatureIscanToMapAfter<< nPair << "\t" << curr_point.transpose() << "\t" << -1 <<endl;
                    }

                }
            }
            if(mpMapTree[mapTreeId]->visibleInfo.size() == 0)
                mpMapTree[mapTreeId]->rmse = 0;
            else
            mpMapTree[mapTreeId]->rmse = sqrt(res_per_Map_tree / double(num_point_map_tree));

            // 2. Tree points from current integrated scan
            double res_per_Iscan_tree = 0.0;
            int num_point_Iscan_tree = 0;
            for (int nP = 0; nP < pCurrentIScan->vTreePointMapping[scanTreeId].size(); nP++)
            {
                // points in local
                Eigen::Vector3d curr_point(pCurrentIScan->vTreePointLocal[scanTreeId][nP].x,
                                           pCurrentIScan->vTreePointLocal[scanTreeId][nP].y,
                                           pCurrentIScan->vTreePointLocal[scanTreeId][nP].z);

                double d_ini = computePoint2lineDistance(q_ini * curr_point + t_ini, mpMapTree[mapTreeId]->para.x, mpMapTree[mapTreeId]->para.x + mpMapTree[mapTreeId]->para.n);
                double res_ini = d_ini - mpMapTree[mapTreeId]->para.r;
                res_scantree_ini = res_scantree_ini + res_ini * res_ini;

                double d_ref = computePoint2lineDistance(q_last_curr * curr_point + t_last_curr, vTreeParamsUpdated[cylinderId].x, vTreeParamsUpdated[cylinderId].x + vTreeParamsUpdated[cylinderId].n);
                double res_ref = d_ref - vTreeParamsUpdated[cylinderId].r;
                res_scantree = res_scantree + res_ref * res_ref;

                // res_scantree_ini_ceres = res_scantree_ini_ceres + resBlocks_ini[obsCount] * resBlocks_ini[obsCount];
                // res_scantree_ceres = res_scantree_ceres + resBlocks_out[obsCount] * resBlocks_out[obsCount];
                // obsCount++;

                res_per_Iscan_tree = res_per_Iscan_tree + res_ref * res_ref;
                num_point_Iscan_tree++;

                if (intermediate_result_flag)
                {
                    fFeatureIscanToMapBefore << nPair << "\t" << (q_ini * curr_point + t_ini).transpose() << "\t" << 0  << endl;
                    fFeatureIscanToMapAfter << nPair << "\t" << (q_last_curr * curr_point + t_last_curr).transpose() << "\t" << 0  << endl;
                }
            }
            mpMapTree[mapTreeId]->candRmse = sqrt(res_per_Iscan_tree / double(num_point_Iscan_tree));

        }

        for (int nP = 0; nP < pCurrentIScan->pGroundPointMappingDs->size(); nP++)
        {
            PointType p = pCurrentIScan->pGroundPointMappingDs->points[nP];
            if (octreeGroundPointsFromMap.nearestKSearch(p, 5, pointSearchInd, pointSearchSqDis) > 0)
            {
                Eigen::Matrix<double, 5, 3> matA0;
                Eigen::Matrix<double, 5, 1> matB0 = -1 * Eigen::Matrix<double, 5, 1>::Ones();
                if (pointSearchSqDis[4] < 1.0)
                {
                    for (int j = 0; j < 5; j++)
                    {
                        matA0(j, 0) = mpMapGroundPoint->points[pointSearchInd[j]].x;
                        matA0(j, 1) = mpMapGroundPoint->points[pointSearchInd[j]].y;
                        matA0(j, 2) = mpMapGroundPoint->points[pointSearchInd[j]].z;
                    }
                    // find the norm of plane
                    Eigen::Vector3d norm = matA0.colPivHouseholderQr().solve(matB0);
                    double negative_OA_dot_norm = 1 / norm.norm();
                    norm.normalize();

                    // Here n(pa, pb, pc) is unit norm of plane
                    bool planeValid = true;
                    for (int j = 0; j < 5; j++)
                    {
                        // if OX * n > 0.2, then plane is not fit well
                        if (fabs(norm(0) * mpMapGroundPoint->points[pointSearchInd[j]].x +
                                 norm(1) * mpMapGroundPoint->points[pointSearchInd[j]].y +
                                 norm(2) * mpMapGroundPoint->points[pointSearchInd[j]].z + negative_OA_dot_norm) > 0.2)
                        {
                            planeValid = false;
                            break;
                        }
                    }
                    if (planeValid)
                    {
                        Eigen::Vector4d params;
                        params << norm, negative_OA_dot_norm;
                        Eigen::Vector3d curr_point(pCurrentIScan->pGroundPointLocalDs->points[nP].x,
                                                   pCurrentIScan->pGroundPointLocalDs->points[nP].y,
                                                   pCurrentIScan->pGroundPointLocalDs->points[nP].z);

                        double d_ini = ((q_ini * curr_point + t_ini).transpose() * norm + negative_OA_dot_norm);
                        res_ground_ini = res_ground_ini + d_ini * d_ini;

                        double d_refine = ((q_last_curr * curr_point + t_last_curr).transpose() * norm + negative_OA_dot_norm);
                        res_ground = res_ground + d_refine * d_refine;

                        if (intermediate_result_flag)
                        {
                            for (int j = 0; j < 5; j++)
                            {
                                fFeatureIscanToMapBefore << -1 << "\t" << mpMapGroundPoint->points[pointSearchInd[j]].x << "\t"
                                                         << mpMapGroundPoint->points[pointSearchInd[j]].y << "\t"
                                                         << mpMapGroundPoint->points[pointSearchInd[j]].z << "\t" << -1 << endl;
                                fFeatureIscanToMapAfter << -1 << "\t" << mpMapGroundPoint->points[pointSearchInd[j]].x << "\t"
                                                        << mpMapGroundPoint->points[pointSearchInd[j]].y << "\t"
                                                        << mpMapGroundPoint->points[pointSearchInd[j]].z << "\t" << -1 << endl;
                            }
                            fFeatureIscanToMapBefore << -1 << "\t" << (q_ini * curr_point + t_ini).transpose() << "\t" << 0 << endl;
                            fFeatureIscanToMapAfter << -1 << "\t" << (q_last_curr * curr_point + t_last_curr).transpose() << "\t" << 0 << endl;
                        }
                        // res_ground_ini_ceres = res_ground_ini_ceres + resBlocks_ini[obsCount] * resBlocks_ini[obsCount];
                        // res_ground_ceres = res_ground_ceres + resBlocks_out[obsCount] * resBlocks_out[obsCount];
                        // obsCount++;
                    }
                }
            }
        }

        cout << "\tError for scan/map tree & ground point: " << sqrt(res_scantree_ini / double(numScanTreeObs)) << " -> " << sqrt(res_scantree / double(numScanTreeObs)) << " m";
        cout << "\t" << sqrt(res_maptree_ini / double(numMapTreeObs)) << " -> " << sqrt(res_maptree / double(numMapTreeObs)) << " m" << endl;
        cout << "\t" << sqrt(res_ground_ini / double(numGroundObs)) << " -> " << sqrt(res_ground / double(numGroundObs)) << " m" << endl;

        fMapLog << "\tError for scan/map tree & ground point: " << sqrt(res_scantree_ini / double(numScanTreeObs)) << " -> " << sqrt(res_scantree / double(numScanTreeObs)) << " m";
        fMapLog << "\t" << sqrt(res_maptree_ini / double(numMapTreeObs)) << " -> " << sqrt(res_maptree / double(numMapTreeObs)) << " m";
        fMapLog << "\t" << sqrt(res_ground_ini / double(numGroundObs)) << " -> " << sqrt(res_ground / double(numGroundObs)) << " m" << endl;

        // cout << "\tError for scan/map tree & ground point: " << sqrt(res_scantree_ini_ceres / double(numScanTreeObs)) * tree_std << " -> " << sqrt(res_scantree_ceres / double(numScanTreeObs)) * tree_std << " m";
        // cout << "\t" << sqrt(res_maptree_ini_ceres / double(numMapTreeObs)) * tree_std << " -> " << sqrt(res_maptree_ceres / double(numMapTreeObs)) * tree_std << " m" << endl;
        // cout << "\t" << sqrt(res_ground_ini_ceres / double(numGroundObs)) * ground_std << " -> " << sqrt(res_ground_ceres / double(numGroundObs)) * ground_std << " m" << endl;

        // fMapLog << "\tError for scan/map tree & ground point: " << sqrt(res_scantree_ini_ceres / double(numScanTreeObs)) * tree_std << " -> " << sqrt(res_scantree_ceres / double(numScanTreeObs)) * tree_std << " m";
        // fMapLog << "\t" << sqrt(res_maptree_ini_ceres / double(numMapTreeObs)) * tree_std << " -> " << sqrt(res_maptree_ceres / double(numMapTreeObs)) * tree_std << " m" << endl;
        // fMapLog << "\t" << sqrt(res_ground_ini_ceres / double(numGroundObs)) * ground_std << " -> " << sqrt(res_ground_ceres / double(numGroundObs)) * ground_std << " m" << endl;

        // update tree parameter, check if it's valid
        vector<pair<int, int>> updatedPairs; // map tree id & scan tree id
        vector<CylinderPara> vValidTreeParamsUpdated;

        for (int nPair = 0; nPair < numPair; nPair++)
        {
            //check if tree parameter is valid 
            int mapTreeId = mvPairs[nPair].first;
            int cylinderId = nPair;
            
            //if the parameter is valid, save this pair
            if(mpMapTree[mapTreeId]->is_valid_para(vTreeParamsUpdated[cylinderId]))
            {
                updatedPairs.push_back(mvPairs[nPair]);
                vValidTreeParamsUpdated.push_back(vTreeParamsUpdated[cylinderId]);
            }
        }

        mvPairs = updatedPairs;
        numPair = mvPairs.size();
        mvValidTreeParamsUpdated = vValidTreeParamsUpdated;
        fMapLog << "\tValid pairs " << numPair << endl;

        // system("read -p 'Press Enter to continue...' var");
        delete []cylinders;
        delete []cylinders_std;
        //delete problem;
        //delete costFunction;

        iter_count++;
        if (intermediate_result_flag)
        {
            fFeatureIscanToMapBefore.close();
            fFeatureIscanToMapAfter.close();
        }

        // CONTROL
        if (numPair < Min_Num_Tree_)
        {
            f_mapping_debug << "Iscan " << mnCount << ": \tNumber of tree matches are too few: " << numPair << endl;
            Flag_Pose_to_Map = false;
            return false;
        }
    }
    return true;
}

/********************************************************
 * Function to check if the matched trees are valid, if yes, add them to mapTree
 * Input: mvPairs, mvValidTreeParamsUpdated
 * Output: mvPairs (updated)
 ********************************************************/
/**
 * @brief Validate and add matched integrated scan trees to the existing map trees.
 * 
 * This function goes through each matched pair in `mvPairs`, validates the tree match using residuals, 
 * and adds the corresponding tree from the current integrated scan to the map tree if valid.
 * 
 * @details If the tree is added successfully, the tree status is updated in `vTreeStatus`. The valid pairs
 * are retained in `mvPairs`.
 */
void Mapping::addIscanTreetoMapTree()
{
    // go through all survived pairs
    int numValidPair = 0;
    vector<pair<int, int>> updatedPairs; // map tree id & scan tree id
    for (int nPair = 0; nPair < mvPairs.size(); nPair++)
    {
        int mapTreeId = mvPairs[nPair].first;
        int IscanTreeId = mvPairs[nPair].second;
        pair<int, int> cand_treeinfo = make_pair(mnCount, IscanTreeId);

        //if this match is valid
        if (mpMapTree[mapTreeId]->add_iscan_tree(mvValidTreeParamsUpdated[nPair], cand_treeinfo))
        {
            numValidPair++;
            pCurrentIScan->vTreeStatus[IscanTreeId] = 1;
            updatedPairs.push_back(mvPairs[nPair]);
        }
    }
    fMapLog << "\tValid pairs after checking residuals " << mvPairs.size() << " -> " << numValidPair << endl;

    mvPairs = updatedPairs;
}
/********************************************************
 * Function to update all points/features using T_l_m_updated,
 * and then add features to map:
 *  1. tree features that are not corresponding to available trees
 *  2. all ground points
 ********************************************************/
/**
 * @brief Add features (trees and ground points) from the current integrated scan to the global map.
 * 
 * @details This function performs:
 * 1. Matching of trees between integrated scan and map tree.
 * 2. Adding unmatched trees as new map trees.
 * 3. Downsampling and adding ground points to the global map.
 * 
 * Requires transformation T_l_m_updated to be computed prior to this step.
 */

void Mapping::addFeaturetoMap()
{
    fMapLog << "--------------------------------------------------" << endl;
    fMapLog << "3. Add Feature to Map" << endl;

    // tree feautures
    // compute tree parameter/points using the updated T_l_m
    pCurrentIScan->computeMapTreePoints(IntegratedScan::REFINED);
    pCurrentIScan->computeMapTreeParameter(IntegratedScan::REFINED);

    // 1. Some trees are already matched to current tree
    if( (mbInit || global_map_flag_ ) && Flag_Pose_to_Map )//&& (!DTM_only_flag_)
        addIscanTreetoMapTree();

    // 2. Trees that are not matched, initialize a new maptree
    for (int nTree = 0; nTree < pCurrentIScan->vTreeStatus.size(); nTree++)
    {    
        if (pCurrentIScan->vTreeStatus[nTree] == 0) // not added to map yet
        {
            MapTree *pTempTree = new MapTree(mpMapTree.size(), pCurrentIScan->vTreeParamMapping[nTree],make_pair(mnCount, nTree), this);
            
#ifdef FIT_CYLINDER_TREE
			if (pTempTree->poor_state_tree_ || pTempTree->fit_tree_failed_)
			{
				delete pTempTree;
				continue;
			}
#endif
            mpMapTree.push_back(pTempTree);
        }
    }

#ifdef MAPTESTINVALID
    for (int nT = 0; nT < mpMapTree.size(); nT++)
        fMapLog << nT << "\t" << mpMapTree[nT]->para.x.transpose() << "\t" << mpMapTree[nT]->para.n.transpose() << "\t" << mpMapTree[nT]->para.r << endl;
#endif

    // ground points
    pCurrentIScan->computeMapGroundPoints(IntegratedScan::REFINED);
    if (!global_map_flag_)
    {
        TicToc t_update_ground_map;

        *mpMapGroundPoint += *pCurrentIScan->pGroundPointMappingDs;

        // pcl::VoxelGrid<PointType> downSizeFilter;
        // downSizeFilter.setInputCloud(mpMapGroundPoint);
        // downSizeFilter.setLeafSize(0.1, 0.1, 0.1);
        // downSizeFilter.filter(*mpMapGroundPoint);
        //(TODO:pcl have integrater flow problem)
        //(TODO:filter size from 0.1 to 0.2)
        *mpMapGroundPoint=DownSamplePointCloudBasedOnDistance(mpMapGroundPoint,0.2);

        fMapLog << "Update ground map time: " << t_update_ground_map.toc() << endl;
    }

    fMapLog << "\tGround points: " << mpMapGroundPoint->size();
    fMapLog << "\tMap trees: " << mpMapTree.size() << endl;
}

/**********************************************
match tree of a scan to tree map:
input:
    mapTree: map trees
    preScanTree: estimated trees of a scan (predicted by R_lup_m and r_lup_m)
    threshold: matching threshold
output:
    pairs (map id <-> scan tree id)
***********************************************/
/**
 * @brief Match trees of the current scan to a list of map tree centers using a distance threshold.
 * 
 * @param mapTree Vector of 3D centers for map trees.
 * @param preScanTree Vector of 3D estimated centers from the scan.
 * @param threshold Matching threshold (per tree).
 * @param pairs Output matched pairs of (map tree index, scan tree index).
 */

void Mapping::matchTreeScantoMap(const vector<Eigen::Vector3d> *mapTree, const vector<Eigen::Vector3d> *preScanTree, const vector<double> *threshold, std::vector<pair<int, int>> &pairs)
{
    // std::vector<int> pairMapTreeId;
    int numTree = preScanTree->size();

    vector<pair<int, int>> tempPairs;
    for (int nT = 0; nT < numTree; nT++)
    {
        Eigen::Vector3d estimateCenter = preScanTree->at(nT);
        double MaxDistance = threshold->at(nT);
        // fMapLog<< estimateCenter.transpose()<<" " << MaxDistance <<endl;
        double score_best = 100.0;
        int id_best = -1;
        for (int nMapT = 0; nMapT < mapTree->size(); nMapT++) // prev scan
        {
            Eigen::Vector3d mapCenter = mapTree->at(nMapT);
            double d_2d = (estimateCenter(0) - mapCenter(0)) * (estimateCenter(0) - mapCenter(0)) + (estimateCenter(1) - mapCenter(1)) * (estimateCenter(1) - mapCenter(1));

            if (d_2d < MaxDistance * MaxDistance)
            {
                double score = d_2d;
                if (score < score_best)
                {
                    score_best = score;
                    id_best = nMapT;
                }
            }
        }
        if (id_best >= 0)
            tempPairs.push_back(make_pair(id_best, nT));
    }

    // backward matching
    bool bMatchFb = false;
    if (bMatchFb)
    {
        for (int nP = 0; nP < tempPairs.size(); nP++)
        {
            int mapTreeId = tempPairs[nP].first;
            Eigen::Vector3d mapCenter = mapTree->at(mapTreeId);
            double score_best = 100.0;
            int id_best = -1;

            for (int nT = 0; nT < numTree; nT++)
            {
                Eigen::Vector3d estimateCenter = preScanTree->at(nT);
                double MaxDistance = threshold->at(nT);

                double d_2d = (estimateCenter(0) - mapCenter(0)) * (estimateCenter(0) - mapCenter(0)) + (estimateCenter(1) - mapCenter(1)) * (estimateCenter(1) - mapCenter(1));

                if (d_2d < MaxDistance * MaxDistance)
                {
                    double score = d_2d;
                    if (score < score_best)
                    {
                        score_best = score;
                        id_best = nT;
                    }
                }
            }

            if (id_best == tempPairs[nP].second)
                pairs.push_back(tempPairs[nP]);
        }
    }
    else
    {
        pairs = tempPairs;
    }
}

/**********************************************
match tree of a scan to tree map (mpMapTree). Use the tree parameter for computation
input:
    preScanTree: estimated trees of a scan (predicted by R_lup_m and r_lup_m)
    threshold: matching threshold
output:
    pairs (map id <-> scan tree id)
***********************************************/
/**
 * @brief Match trees of the current scan to a list of map tree centers using a distance threshold.
 * 
 * @param mapTree Vector of 3D centers for map trees.
 * @param preScanTree Vector of 3D estimated centers from the scan.
 * @param threshold Matching threshold (per tree).
 * @param pairs Output matched pairs of (map tree index, scan tree index).
 */

void Mapping::matchTreeScantoMap(const vector<Eigen::Vector3d> *ScanTree, const vector<double> *threshold, std::vector<pair<int, int>> &pairs)
{
    // std::vector<int> pairMapTreeId;
    int numTree = ScanTree->size();

    vector<pair<int, int>> tempPairs;
    for (int nT = 0; nT < numTree; nT++)
    {
        Eigen::Vector3d estimateCenter = ScanTree->at(nT);
        double MaxDistance = threshold->at(nT);
        // fMapLog<< estimateCenter.transpose()<<" " << MaxDistance <<endl;
        double score_best = 100.0;
        int id_best = -1;
        for (int nMapT = 0; nMapT < mpMapTree.size(); nMapT++) // prev scan
        {
            if(mpMapTree[nMapT]->status == MapTree::DEACTIVATE)
                continue;
            // Eigen::Vector3d mapCenter = mpMapTree[nMapT]->para.x;
            // double d_2d = (estimateCenter(0) - mapCenter(0)) * (estimateCenter(0) - mapCenter(0)) + (estimateCenter(1) - mapCenter(1)) * (estimateCenter(1) - mapCenter(1));
            double dis = mpMapTree[nMapT]->compute_dis(estimateCenter);
            if (dis < MaxDistance )
            {
                double score = dis;
                if (score < score_best)
                {
                    score_best = score;
                    id_best = nMapT;
                }
            }
        }
        if (id_best >= 0)
            tempPairs.push_back(make_pair(id_best, nT));
    }

    // backward matching
    bool bMatchFb = true;
    if (bMatchFb)
    {
        for (int nP = 0; nP < tempPairs.size(); nP++)
        {
            int mapTreeId = tempPairs[nP].first;
            Eigen::Vector3d mapCenter = mpMapTree[mapTreeId]->para.x;
            double score_best = 100.0;
            int id_best = -1;

            for (int nT = 0; nT < numTree; nT++)
            {
                Eigen::Vector3d estimateCenter = ScanTree->at(nT);
                double MaxDistance = threshold->at(nT);
                // double d_2d = (estimateCenter(0) - mapCenter(0)) * (estimateCenter(0) - mapCenter(0)) + (estimateCenter(1) - mapCenter(1)) * (estimateCenter(1) - mapCenter(1));

                double dis = mpMapTree[mapTreeId]->compute_dis(estimateCenter);

                if (dis < MaxDistance)
                {
                    double score = dis;
                    if (score < score_best)
                    {
                        score_best = score;
                        id_best = nT;
                    }
                }
            }

            if (id_best == tempPairs[nP].second)
                pairs.push_back(tempPairs[nP]);
        }
    }
    else
    {
        pairs = tempPairs;
    }
}

/**********************************************
match tree of a scan to tree map:
input:
    mapTree: map trees
    scanTree: trees locations of a scan (R_tree_lup and r_tree_lup)
    pairs: matching pairs
output:
    validTreePair (map id <-> scan tree id)
    R,t: updated 2d similarity
***********************************************/
/**
 * @brief Estimate a 2D similarity transformation (rotation + translation) between matched scan and map tree centers.
 * 
 * @param mapTree Map tree center locations.
 * @param scanTree Scan tree center locations.
 * @param pairs Matched pairs of trees between map and scan.
 * @param validTreePair Output flags indicating valid pairs used in estimation.
 * @param R Output rotation matrix (3x3).
 * @param t Output translation vector (3x1).
 * @param initResult Optional initial 2D transform (cos, sin, tx, ty).
 */
void Mapping::compute2dTransScantoMap(const vector<Eigen::Vector3d> *mapTree, const vector<Eigen::Vector3d> *scanTree, const vector<pair<int, int>> *pairs, std::vector<bool> &validTreePair,
                                      Eigen::Matrix3d &R, Eigen::Vector3d &t, Eigen::Vector4d initResult = Eigen::Vector4d(1, 0, 0, 0))
{
    // compute
    int numPair = pairs->size();
    int numValidPair = numPair;
    Eigen::Vector4d tempResult = initResult;
    double DistanceThreshold = 0.5;
    int iter = 0;
    int maxIter = 4;
    double scale_, theta_, tx_, ty_;

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

            Eigen::Vector3d center_prev = mapTree->at(pairs->at(i).first);
            Eigen::Vector3d center_cur = scanTree->at(pairs->at(i).second);

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
        fMapLog << "\t2d iter " << iter << "\t" << s1 << " -> " << s2 << endl;
        fMapLog << "\tPair: " << numValidPair << " -> ";
#endif
        count = 0;
        numValidPair = 0;
        bool bCont = false;

        for (int i = 0; i < numPair; ++i)
        {
            if (!validTreePair[i])
                continue;

            double d_sqr = res2(count * 2) * res2(count * 2) + res2(count * 2 + 1) * res2(count * 2 + 1);
            // cout << res2(count * 2) << ", " << res2(count * 2 + 1) << ", " << sqrt(d_sqr) << endl;

            count++;

            if (d_sqr > DistanceThreshold * DistanceThreshold)
            {
                validTreePair[i] = false;
                bCont = true;
            }
            else
                numValidPair++;
        }
        scale_ = sqrt(result_(0) * result_(0) + result_(1) * result_(1));
        theta_ = atan2(result_(1), result_(0));
        tx_ = result_(2);
        ty_ = result_(3);
#ifdef EXPORT_LOG
        fMapLog << numValidPair << endl;
        fMapLog << "Result: " << scale_ << "\t" << rad2deg(theta_) << "\t" << tx_ << "\t" << ty_ << endl;
#endif
        DistanceThreshold = 0.9 * DistanceThreshold;

        if (!bCont)
            break;
    }

    Eigen::AngleAxisd rollAngle_(theta_, Eigen::Vector3d::UnitZ());
    Eigen::Matrix3d R_kap_ = rollAngle_.matrix();

    R = R_kap_;
    t = R * Eigen::Vector3d(tx_, ty_, 0);
}

/**********************************************
compute info related to the map tree
input:
    mapTree: map trees
    mvpIScans
output: center of each tree, numPoint of each tree
***********************************************/
/**
 * @brief Compute the center and number of points for each tree in `mpMapTree`.
 * 
 * @details It averages all visible tree points to calculate the center and updates the number of points.
 * If the tree is in ESTABLISHED status, it computes the surface ratio as well.
 */

void Mapping::computeMapTreeInfo()
{
    for (int nMapT = 0; nMapT < mpMapTree.size(); nMapT++)
    {
        if(mpMapTree[nMapT]->status< 0)
        {
            continue;
        }
        int mapTreeId = nMapT;
        Eigen::Vector3d center = Eigen::Vector3d::Zero();
        int numPoint = 0;
    
        for (int nT = 0; nT < mpMapTree[mapTreeId]->visibleInfo.size(); nT++)
        {
            int iscanId = mpMapTree[mapTreeId]->visibleInfo[nT].first;
            int treeId = mpMapTree[mapTreeId]->visibleInfo[nT].second;
            for (int nP = 0; nP < mvpIScans[iscanId]->vTreePointMapping[treeId].size(); nP++)
            {
                // points in mapping
                Eigen::Vector3d curr_point(mvpIScans[iscanId]->vTreePointMapping[treeId][nP].x,
                                        mvpIScans[iscanId]->vTreePointMapping[treeId][nP].y,
                                        mvpIScans[iscanId]->vTreePointMapping[treeId][nP].z);

               center += curr_point;
               numPoint++;
            }
        }
        center = center / double(numPoint);

        mpMapTree[nMapT]->center = center;
        mpMapTree[nMapT]->numPoint = numPoint;
        if(mpMapTree[nMapT]->status == MapTree::ESTABLISHED)
            mpMapTree[nMapT]->compute_surface_ratio();
    }
}

/**
 * @brief Merge map trees that are spatially and structurally similar.
 * 
 * @return int The number of merged trees.
 * 
 * @details This function first checks proximity and residuals between tree point sets. If compatible,
 * trees are merged via `AddMapTree()`. Two rounds are used — one for well-established trees and one
 * for INIT or TBD trees.
 */

int Mapping::mergeMapTree()
{
    // reset visit flag for the tree
    // for (int nMapT = 0; nMapT < mpMapTree.size(); nMapT++)
    // {
    //     mpMapTree[nMapT]->bVisit = false;
    // }

    //int mergedTreeId = 0;
    int num_merged_trees = 0;
    //for each well established tree, check whether other tree is neighboring
    for(int nMapT_Ref = 0; nMapT_Ref < mpMapTree.size(); nMapT_Ref++)
    {
        //(TODO/MODI)（disbale modification now）
        //if it's a well established tree and not visited yet mpMapTree[nMapT_Ref]->bVisit == true ||
        if( mpMapTree[nMapT_Ref]->status == MapTree::INIT || mpMapTree[nMapT_Ref]->status == MapTree::TBD
        || mpMapTree[nMapT_Ref]->status == MapTree::GLOBAL_INIT || mpMapTree[nMapT_Ref]->status == MapTree::DEACTIVATE ) 
        {
            continue;
        }
        
        //mpMapTree[nMapT_Ref]->groupId = mergedTreeId;
        //mpMapTree[nMapT_Ref]->bVisit = true;
        //num_merged_trees++;
        for(int nMapT_Cand = 0; nMapT_Cand < mpMapTree.size(); nMapT_Cand++)
        {
            // //if the tree is not visited yet
            // if(mpMapTree[nMapT_Cand]->bVisit == true)
            // {
            //     continue; 
            // }

            if(nMapT_Cand == nMapT_Ref)
                continue;

            //if the tree a global tree
            if(mpMapTree[nMapT_Cand]->status == MapTree::GLOBAL_INIT || mpMapTree[nMapT_Cand]->status == MapTree::DEACTIVATE )
            {
                continue; 
            }

            //first compute the center point to tree distance, if too large, skip
            if( mpMapTree[nMapT_Ref]->compute_dis(mpMapTree[nMapT_Cand]->center)> 1.0)
            {
                continue;
            }

            //then compute each point to its distance 
            double dis_rmse = computeTreeToTreeResidual(nMapT_Ref,nMapT_Cand);
            //fMapLog << "Pair " << nMapT_Ref << " - " << nMapT_Cand << " :" << dis_rmse<<endl;
            
            if( (dis_rmse < 0.2  && mpMapTree[nMapT_Ref]->status == MapTree::ESTABLISHED)
            || (dis_rmse < 0.1  && mpMapTree[nMapT_Ref]->status == MapTree::SOLID)
            ||(dis_rmse < 0.3  && mpMapTree[nMapT_Ref]->status == MapTree::FITTED))
            {
                //mpMapTree[nMapT_Cand]->groupId = mergedTreeId;
                //mpMapTree[nMapT_Cand]->bVisit = true; 
                num_merged_trees++;

                mpMapTree[nMapT_Ref]->AddMapTree(mpMapTree[nMapT_Cand]);
            }
            //(TODO/MODI)
            // if(dis_rmse < 0.3  &&(mpMapTree[nMapT_Ref]->status == MapTree::INIT||mpMapTree[nMapT_Ref]->status == MapTree::TBD))
            // {
            //     num_merged_trees++;
            //     mpMapTree[nMapT_Ref]->AddMapTree(mpMapTree[nMapT_Cand]);
            // }
        } //end of go through each tree
        //mergedTreeId++;
    }

    //another round to merge TBD and INIT
    for(int nMapT_Ref = 0; nMapT_Ref < mpMapTree.size(); nMapT_Ref++)
    {
        if(!(mpMapTree[nMapT_Ref]->status == MapTree::INIT || mpMapTree[nMapT_Ref]->status == MapTree::TBD))
        {
            continue;
        }
        for(int nMapT_Cand = 0; nMapT_Cand < mpMapTree.size(); nMapT_Cand++)
        {
            if(nMapT_Cand == nMapT_Ref)
                continue;

            if(!(mpMapTree[nMapT_Ref]->status == MapTree::INIT || mpMapTree[nMapT_Ref]->status == MapTree::TBD))
            {
                continue;
            }

            //first compute the center point to tree distance, if too large, skip
            if( mpMapTree[nMapT_Ref]->compute_dis(mpMapTree[nMapT_Cand]->center)> 1.0)
            {
                continue;
            }

            //then compute each point to its distance 
            double dis_rmse = computeTreeToTreeResidual(nMapT_Ref,nMapT_Cand);
            //fMapLog << "Pair " << nMapT_Ref << " - " << nMapT_Cand << " :" << dis_rmse<<endl;
            
            if(dis_rmse < 0.2)
            {
                num_merged_trees++;
                mpMapTree[nMapT_Ref]->AddMapTree(mpMapTree[nMapT_Cand]);
            }
        } //end of go through each tree
    }

    fMapLog << "Total # of trees and trees being merged: " <<  mpMapTree.size() << " / " << num_merged_trees <<endl;
    return num_merged_trees;
}

/**
 * @brief Compute RMSE between tree points of a candidate and a reference map tree.
 * 
 * @param refId Reference tree ID in `mpMapTree`.
 * @param candId Candidate tree ID in `mpMapTree`.
 * @return double RMSE of point-to-tree surface distance.
 */

double Mapping::computeTreeToTreeResidual(int refId, int candId)
{
    double dis_rmse = 0.0;
    //(TODO) move candidate tree center to ref tree center
	Eigen::Vector3d move_vector = mpMapTree[candId]->para.x - mpMapTree[refId]->para.x;
    for (int nT = 0; nT < mpMapTree[candId]->visibleInfo.size(); nT++)
    {
        int iscanId = mpMapTree[candId]->visibleInfo[nT].first;
        int treeId = mpMapTree[candId]->visibleInfo[nT].second;
        for (int nP = 0; nP < mvpIScans[iscanId]->vTreePointMapping[treeId].size(); nP++)
        {
            // points in mapping
            Eigen::Vector3d curr_point(mvpIScans[iscanId]->vTreePointMapping[treeId][nP].x,
                                       mvpIScans[iscanId]->vTreePointMapping[treeId][nP].y,
                                       mvpIScans[iscanId]->vTreePointMapping[treeId][nP].z);
			//(TODO)
			curr_point -= move_vector;
            double dis = mpMapTree[refId]->compute_dis(curr_point);
            dis_rmse = dis_rmse + dis * dis;
        }
    }
    
    dis_rmse = sqrt(dis_rmse / double(mpMapTree[candId]->numPoint));
    return(dis_rmse);
}

//the outlier ratio is calculated
/**
 * @brief Compute RMSE and outlier ratio between a candidate and reference tree.
 * 
 * @param refId Reference map tree index.
 * @param candId Candidate map tree index.
 * @param residual_distribution Output: vector containing ratio of [negative outliers, inliers, positive outliers].
 * @return double RMSE of residual distances.
 */

double Mapping::computeTreeToTreeResidual(int refId, int candId, vector<double> &residual_distribution )
{
    double threshold = 0.2;
    int num_negative_outlier = 0, num_positive_outlier = 0, num_inlier = 0;
    double dis_rmse = 0.0;
    for (int nT = 0; nT < mpMapTree[candId]->visibleInfo.size(); nT++)
    {
        int iscanId = mpMapTree[candId]->visibleInfo[nT].first;
        int treeId = mpMapTree[candId]->visibleInfo[nT].second;
        for (int nP = 0; nP < mvpIScans[iscanId]->vTreePointMapping[treeId].size(); nP++)
        {
            // points in mapping
            Eigen::Vector3d curr_point(mvpIScans[iscanId]->vTreePointMapping[treeId][nP].x,
                                       mvpIScans[iscanId]->vTreePointMapping[treeId][nP].y,
                                       mvpIScans[iscanId]->vTreePointMapping[treeId][nP].z);
            double dis = mpMapTree[refId]->compute_dis(curr_point);
            dis_rmse = dis_rmse + dis * dis;
            if(abs(dis) < threshold)
            {
                num_inlier++;
            }
            else
            {
                dis < 0 ? num_negative_outlier++ :  num_positive_outlier++;
            }
        }
    }
    
    dis_rmse = sqrt(dis_rmse / double(mpMapTree[candId]->numPoint));
    residual_distribution.push_back(double(num_negative_outlier)/double(mpMapTree[candId]->numPoint));
    residual_distribution.push_back(double(num_inlier)/double(mpMapTree[candId]->numPoint));
    residual_distribution.push_back(double(num_positive_outlier)/double(mpMapTree[candId]->numPoint));
    return(dis_rmse);
}

/**
 * @brief Perform long-term loop closure optimization using matched map and scan features.
 * 
 * @details Includes:
 * - Feature extraction
 * - Raw-level or structural optimization
 * - Pose and map tree update
 * - Export of trajectory and feature data
 */

void Mapping::LoopClosure()
{
    fMapLog << "--------------------------------------------------" << endl;
    fMapLog << "Long term optimization.." <<endl;

    //iscan before iscan_index_start will not change
    DeriveIndexLC();

    //derive planar patch 
    DerivePlanarPatchLC();

    //export the feature
    std::string outFeatureBeforeLC(output_folder_loopclouse + prefix +  to_string(count_lc) + "_FeatureBeforeLC_" + to_string(mvpIScans[iscan_index_end_]->indScans.back().scanID) + ".las");
    // ExportFeatureLC(outFeatureBeforeLC);
    std::string outTrajBeforeLC(output_folder_loopclouse + prefix + to_string(count_lc) + "_TrajBeforeLC_" + to_string(mvpIScans[iscan_index_end_]->indScans.back().scanID) + ".txt");
    ExportTrajectoryLC(outTrajBeforeLC);

    // conduct optimization
	if (raw_point_flag_)
	{
		OptimizeRawLevelLC();
		// OptimizeRawLevelLCUsingSurfaceElements();
	}
	else
	{
		OptimizeLC();
	}

    //update the pose and points for the Iscan and mapTree
    UpdatePosePointsLC();
    
    //export the feature
    std::string outFeatureAfterLC(output_folder_loopclouse + prefix + to_string(count_lc) + "_FeatureAfterLC_" + to_string(mvpIScans[iscan_index_end_]->indScans.back().scanID) + ".las");
    // ExportFeatureLC(outFeatureAfterLC);
    std::string outTrajAfterLC(output_folder_loopclouse + prefix+ to_string(count_lc) + "_TrajAfterLC_" + to_string(mvpIScans[iscan_index_end_]->indScans.back().scanID) + ".txt");
    ExportTrajectoryLC(outTrajAfterLC);

    std::string outTrajRefLC(output_folder_loopclouse + prefix + to_string(count_lc) + "_TrajRefLC_" + to_string(mvpIScans[iscan_index_end_]->indScans.back().scanID) + ".txt");
    ExportTrajectoryRefLC(outTrajRefLC);
    
    //system("read -p 'Press Enter to continue...' var");

    //update trees 
    UpdateMapTreeLC();

    //update info for trees included/or not included in optimization
    //update status as well
    fMapLog<< "Update tree status" <<endl;
    std::vector<int> pointSearchInd;
    std::vector<float> pointSearchSqDis;
    for (int nTree = 0; nTree < mpMapTree.size(); nTree++)
    {
        mpMapTree[nTree]->ComputeAttributes();

        //derive the projection on ground
        PointType p;
        p.x = mpMapTree[nTree]->center(0);
        p.y = mpMapTree[nTree]->center(1);
        p.z = mpMapTree[nTree]->center(2);
        double z_temp;
        if (octreeGroundPointsFromMap.nearestKSearch(p, 1, pointSearchInd, pointSearchSqDis) > 0)
        {
            z_temp = mpMapGroundPoint->points[pointSearchInd[0]].z;
        }
        else
        {
            z_temp = p.z;
        }

        Eigen::Vector3d ground_point(p.x, p.y, z_temp);
        mpMapTree[nTree]->ground_point = ground_point;
    }
    CountMapTreeType();

    count_lc++;
    iscan_index_prev_loop_closure_ = iscan_index_end_;
    //system("read -p 'Press Enter to continue...' var");
}

/**
 * @brief Identify the range of Iscans to be refined during loop closure optimization.
 * 
 * @details It updates `iscan_index_start_`, `iscan_index_end_`, and `iscan_index_feature_start_`
 * based on usage in active map trees.
 */

void Mapping::DeriveIndexLC()
{
    // the end of the Iscan_end from previous LC
    iscan_index_start_ = 0;
    iscan_index_end_ = mnCount - 1;

    //test 
    iscan_index_prev_loop_closure_ = -1;

    // go through the trees that will be used in LC to find iscan start
    iscan_index_feature_start_ =  iscan_index_prev_loop_closure_ + 1;

    //for a tree that includes newly added Iscan
    int num_used_tree = 0;
    vector<int> used_trees; //# vector of unknown tree -> mapTreeID

    for (int nTree = 0; nTree < mpMapTree.size(); nTree++)
    {
        if (mpMapTree[nTree]->end_iscan_index_ >iscan_index_prev_loop_closure_  && (mpMapTree[nTree]->status == MapTree::ESTABLISHED ||
                                              mpMapTree[nTree]->status == MapTree::SOLID || mpMapTree[nTree]->status == MapTree::FITTED))
        {
            used_trees.push_back(nTree);
            num_used_tree++;
            iscan_index_feature_start_ = min(iscan_index_feature_start_, mpMapTree[nTree]->start_iscan_index_);
        }
    }
    used_trees_ = used_trees;

    fMapLog << "Start index for including feature, refining pose, End index: " << iscan_index_feature_start_ << " / " << iscan_index_start_
            << " / " << iscan_index_end_ << endl;
    //system("read -p 'Press Enter to continue...' var");
}

/**
 * @brief Compute normal vector and covariance matrix from a set of 3D points.
 * 
 * @param data Input point cloud as a vector of 3D points.
 * @param parameters Output: Plane parameters (nx, ny, nz, cx, cy, cz).
 * @param cov Output: Covariance matrix of the input points.
 * @return true If normal computation succeeds.
 * @return false If the input point set is too small (<3).
 */

bool ComputePointsNormal(const std::vector< Eigen::Vector3d>& data, std::vector<double> &parameters, Eigen::Matrix3d& cov)
{
	if (data.size() < 3)
	{
		return false;
	}
	Eigen::Matrix3d covariance_matrix = Eigen::Matrix3d::Zero();
	Eigen::Matrix<double, 1, 9, Eigen::RowMajor> inter_matrix = Eigen::Matrix<double, 1, 9, Eigen::RowMajor>::Zero();
	size_t point_count = data.size();
	for (int i = 0; i < point_count; ++i)
	{
		inter_matrix[0] += data[i](0) * data[i](0);
		inter_matrix[1] += data[i](0) * data[i](1);
		inter_matrix[2] += data[i](0) * data[i](2);
		inter_matrix[3] += data[i](1) * data[i](1);
		inter_matrix[4] += data[i](1) * data[i](2);
		inter_matrix[5] += data[i](2) * data[i](2);
		inter_matrix[6] += data[i](0);
		inter_matrix[7] += data[i](1);
		inter_matrix[8] += data[i](2);
	}

	inter_matrix /= static_cast<double> (point_count);
	covariance_matrix.coeffRef(0) = inter_matrix[0] - inter_matrix[6] * inter_matrix[6];
	covariance_matrix.coeffRef(1) = inter_matrix[1] - inter_matrix[6] * inter_matrix[7];
	covariance_matrix.coeffRef(2) = inter_matrix[2] - inter_matrix[6] * inter_matrix[8];
	covariance_matrix.coeffRef(4) = inter_matrix[3] - inter_matrix[7] * inter_matrix[7];
	covariance_matrix.coeffRef(5) = inter_matrix[4] - inter_matrix[7] * inter_matrix[8];
	covariance_matrix.coeffRef(8) = inter_matrix[5] - inter_matrix[8] * inter_matrix[8];
	covariance_matrix.coeffRef(3) = covariance_matrix.coeff(1);
	covariance_matrix.coeffRef(6) = covariance_matrix.coeff(2);
	covariance_matrix.coeffRef(7) = covariance_matrix.coeff(5);

	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::Matrix3d> solver(covariance_matrix, Eigen::Matrix3d::Identity(), Eigen::ComputeEigenvectors | Eigen::Ax_lBx);

	parameters.clear();
	parameters.resize(6);
	parameters[0] = solver.eigenvectors().col(0)[0];
	parameters[1] = solver.eigenvectors().col(0)[1];
	parameters[2] = solver.eigenvectors().col(0)[2];
	parameters[3] = inter_matrix(6);
	parameters[4] = inter_matrix(7);
	parameters[5] = inter_matrix(8);

	Eigen::VectorXd eigenvalues = solver.eigenvalues();
	std::sort(eigenvalues.data(), eigenvalues.data() + eigenvalues.size(), std::greater<double>());

	double lamda1 = eigenvalues(0);
	double lamda2 = eigenvalues(1);
	double lamda3 = eigenvalues(2);
	double planarity = (std::sqrt(lamda2) - std::sqrt(lamda3)) / std::sqrt(lamda1);
	cov = covariance_matrix;
	return true;
}

/**
 * @brief Derives planar patches from ground points for loop closure.
 *
 * This function extracts reference ground points either from the global map or
 * previous IScans and partitions them into grids. For each grid cell, it checks 
 * if enough ground points from new IScans are present. If yes, it computes the 
 * local planar parameters for later optimization during loop closure.
 */

void Mapping::DerivePlanarPatchLC()
{
    //reference ground point cloud, from global map and/or Iscans before min(iscan_index_feature_start_, iscan_index_start_)
    pcl::PointCloud<PointType>::Ptr reference_ground_points (new pcl::PointCloud<PointType>);
    if(global_map_flag_)
    {
        pcl::copyPointCloud(*mpMapGroundPoint, *reference_ground_points);
    }
    else
    {
        for (int iscan_id = 0; iscan_id < min(iscan_index_feature_start_, iscan_index_start_); iscan_id++)
        {
            for (int nP = 0; nP < mvpIScans[iscan_id]->pGroundPointMappingDs->size(); nP++)
            {
                reference_ground_points->push_back( mvpIScans[iscan_id]->pGroundPointMappingDs->points[nP]);
            }
        }

        if (reference_ground_points->size() > 0)
        {
            //(TODO:filter size from 0.1 to 0.2)
            double ds_size = 0.2;

            // pcl::VoxelGrid<PointType> downSizeFilter;
            // downSizeFilter.setInputCloud(reference_ground_points);
            // downSizeFilter.setLeafSize(ds_size, ds_size, ds_size);
            // downSizeFilter.filter(*reference_ground_points);
            //(TODO:pcl have integrater overflow problem)
            *reference_ground_points=DownSamplePointCloudBasedOnDistance(reference_ground_points,ds_size);
        }
    }
    fMapLog << "Number of points in the reference ground points: " << reference_ground_points->size() <<endl;
    
    //define the area of conducting loop closure
    double minX = 100000.0, maxX = -100000.0, minY = 1000000.0, maxY = -1000000.0;
    for (int iscan_id = iscan_index_prev_loop_closure_+1; iscan_id <= iscan_index_end_; iscan_id ++)
    {
        minX = min(minX, mvpIScans[iscan_id]->r_local_m_updated(0)-50.0);
        maxX = max(maxX, mvpIScans[iscan_id]->r_local_m_updated(0)+50.0);
        minY = min(minY, mvpIScans[iscan_id]->r_local_m_updated(1)-50.0);
        maxY = max(maxY, mvpIScans[iscan_id]->r_local_m_updated(1)+50.0);
    }
    fMapLog << "Area for optimization: " << minX << "\t" << maxX  <<"\t" << minY << "\t" << maxY <<endl;

    // partition to grid.. 2* (tree+grid) <= res
    double res = KGround_resolution;
    //double tree_radius = 0.5; // 1m
    //double grid_radius = 0.29;
    double grid_radius = mPara.patch_size_lc * 0.5;
    minX = res * floor(minX / res);
    minY = res * floor(minY / res);
    maxX = res * ceil(maxX / res);
    maxY = res * ceil(maxY / res);
    int col = round((maxX - minX) / res);
    int row = round((maxY - minY) / res);

    // compute center of each grid
    vector<pair<double, double>> grid_centers;
    int num_grid = grid_centers.size();
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            grid_centers.push_back(make_pair(double(j + 0.5) * res + minX, double(i + 0.5) * res + minY));
        }
    }

    // find points belong to each grid
    vector<vector<pair<int, int>>> point_index_per_grid(grid_centers.size(), vector<pair<int, int>>(0));
    vector<vector<Eigen::Vector3d>> point_reference_per_grid(grid_centers.size(), vector<Eigen::Vector3d>(0));
    vector<int> num_point_from_new_iscan_per_grid(grid_centers.size(),0);
    // ground points from Iscan that will be included
    for (int iscan_id = iscan_index_feature_start_; iscan_id <= iscan_index_end_; iscan_id++)
    {
        if (raw_point_flag_)
        {
            //compute raw points in mapping
            mvpIScans[iscan_id]->computeRawGroundPointToMapping();

            for (int nP = 0; nP < mvpIScans[iscan_id]->vGroundPointMapping.size(); nP++)
            {
                double x_p = double(mvpIScans[iscan_id]->vGroundPointMapping[nP].x);
                double y_p = double(mvpIScans[iscan_id]->vGroundPointMapping[nP].y);
                int i = floor((y_p - minY) / res);
                int j = floor((x_p - minX) / res);
                if (i < 0 || i >= row || j < 0 || j >= col)
                    continue;

                int index = i * col + j;

                if (abs(x_p - grid_centers[index].first) < grid_radius && abs(y_p - grid_centers[index].second) < grid_radius)
                {
                    point_index_per_grid[index].push_back(make_pair(iscan_id, nP));
                    if (iscan_id >= iscan_index_prev_loop_closure_ + 1)
                    {
                        num_point_from_new_iscan_per_grid[index]++;
                    }
                }
            }
        }
        else
        {
            //if iscan level, the iscan is failed
            if(!mvpIScans[iscan_id]->flag_iscan_to_map)
                continue;

            for (int nP = 0; nP < mvpIScans[iscan_id]->pGroundPointMappingDs->size(); nP++)
            {
                double x_p = double(mvpIScans[iscan_id]->pGroundPointMappingDs->points[nP].x);
                double y_p = double(mvpIScans[iscan_id]->pGroundPointMappingDs->points[nP].y);
                int i = floor((y_p - minY) / res);
                int j = floor((x_p - minX) / res);
                if (i < 0 || i >= row || j < 0 || j >= col)
                    continue;

                int index = i * col + j;

                if (abs(x_p - grid_centers[index].first) < grid_radius && abs(y_p - grid_centers[index].second) < grid_radius)
                {
                    point_index_per_grid[index].push_back(make_pair(iscan_id, nP));
                    if (iscan_id >= iscan_index_prev_loop_closure_ + 1)
                    {
                        num_point_from_new_iscan_per_grid[index]++;
                    }
                }
            }
        }
    }
    int reference_point_number = 0;
    // ground points from reference point cloud
    for (int nP = 0; nP < reference_ground_points->size(); nP++)
    {
        double x_p = double(reference_ground_points->points[nP].x);
        double y_p = double(reference_ground_points->points[nP].y);
        double z_p = double(reference_ground_points->points[nP].z);
        int i = floor((y_p - minY) / res);
        int j = floor((x_p - minX) / res);

        if (i < 0 || i >= row || j < 0 || j >= col)
            continue;

        int index = i * col + j;

        if (abs(x_p - grid_centers[index].first) < grid_radius && abs(y_p - grid_centers[index].second) < grid_radius)
        {
            point_reference_per_grid[index].push_back(Eigen::Vector3d(x_p, y_p, z_p));
            reference_point_number++;
        }
    }
    fMapLog << "\treference_point_number: " << reference_point_number <<endl;
    // find valid planar patches and get plane parameters in local
    vector<vector<pair<int,int>>> point_index_per_grid_valid;
    vector<vector<Eigen::Vector3d>> point_ref_per_grid_valid;
    vector<Eigen::Vector4d> plane_params;
    std::vector<Eigen::Vector3d> patch_points;
    PointType p;

    for (int nGrid = 0; nGrid < grid_centers.size(); nGrid++)
    {
        patch_points.clear();
        if (num_point_from_new_iscan_per_grid[nGrid] < 100) // size limitation for points from new scan
            continue;

        Eigen::Vector4f params;
        Eigen::Vector3f center;
        float Z_mean = 0.0;
        for (int nP = 0; nP < point_index_per_grid[nGrid].size(); nP++)
        {
            pair<int, int> info = point_index_per_grid[nGrid][nP];
            if (raw_point_flag_)
            {
                p = mvpIScans[info.first]->vGroundPointMapping[info.second];
            }
            else
            {
                p = mvpIScans[info.first]->pGroundPointMappingDs->points[info.second];
            }
            Z_mean += p.z;
            patch_points.push_back(p.getVector3fMap().cast<double>());
        }
        Z_mean /= float(point_index_per_grid[nGrid].size());
        std::vector<double> parameters;
		Eigen::Matrix3d cov;
		// if(ComputePointsNormal(patch_points, parameters,cov))
        // {
        //     Eigen::Vector3d plane_normal(parameters[0], parameters[1], parameters[2]);
        //     plane_normal = plane_normal.normalized();
        //     params << plane_normal(0), plane_normal(1), plane_normal(2), -Z_mean;
        // }
        // else
        {
            params << 0.0, 0.0, 1.0, -Z_mean;
        }   
        plane_params.push_back(params.cast<double>());
        point_index_per_grid_valid.push_back(point_index_per_grid[nGrid]);
        point_ref_per_grid_valid.push_back(point_reference_per_grid[nGrid]);
    }
    point_ref_per_grid_ = point_ref_per_grid_valid;
    point_index_per_grid_ = point_index_per_grid_valid;
    plane_params_ = plane_params;

    fMapLog << "\tValid planar patch: " << point_index_per_grid.size() << " -> " << point_index_per_grid_.size() <<endl;
}

/**
 * @brief Performs loop closure optimization using trees (cylinders) and planar patches.
 *
 * This method formulates a nonlinear optimization problem using the Ceres Solver,
 * where the cost functions incorporate:
 *  - Tree observations (cylindrical feature constraints)
 *  - Ground plane observations
 *  - Global tree priors
 *  - Odometry constraints between IScans
 *
 * It updates the trajectory and tree parameters based on minimized residuals
 * and records the before-after errors for analysis.
 */

void Mapping::OptimizeLC()
{
    double tree_std = mPara.mapTreeStd;
    double ground_std = mPara.mapGroundStd;
    double ref_ground_std = ref_ground_std_; // 0.02 0.005
    double global_tree_loc_std = mPara.global_tree_loc_std;
    int max_num_map_point_per_tree = Max_Point_Per_Tree_Iscan_LC;

    int iter = 0;
    int maxIter = 1;
    while (iter < maxIter)
    {
        fMapLog << "Iter " << ++iter;

        fMapLog << "\tNumber of included trees " << used_trees_.size();

        // initialization
        int numEpoch = iscan_index_end_ - iscan_index_start_ + 1;
        int numCylinder = used_trees_.size();
        int numPlane = plane_params_.size();

        // unknowns
        int numPosUnknown = numEpoch * PosBlockSize;
        int numOriUnknown = numEpoch * OriBlockSize;
        int numCylinderUnknown = numCylinder * CylinderBlockSize;
        int numPlaneUnknown = numPlane * PlaneBlockSize;
        int numUnknown = numPosUnknown + numOriUnknown + numCylinderUnknown + numPlaneUnknown; // all unknowns
        double *poss = new double[numPosUnknown];
        double *poss_std = new double[numPosUnknown];
        double *oris = new double[numOriUnknown];
        double *oris_std = new double[numOriUnknown];
        double *cylinders = new double[numCylinderUnknown];
        double *cylinders_std = new double[numCylinderUnknown];
        double *planes = new double[numPlaneUnknown];
        double *planes_std = new double[numPlaneUnknown];

        // fill parameter
        // 1. trajectory information
        for (int nEpoch = 0; nEpoch < numEpoch; nEpoch++)
        {
            int iscan_id = nEpoch + iscan_index_start_;
            *(poss + nEpoch * PosBlockSize + 0) = mvpIScans[iscan_id]->r_local_m_updated(0);
            *(poss + nEpoch * PosBlockSize + 1) = mvpIScans[iscan_id]->r_local_m_updated(1);
            *(poss + nEpoch * PosBlockSize + 2) = mvpIScans[iscan_id]->r_local_m_updated(2);

            Eigen::Quaterniond q(mvpIScans[iscan_id]->R_local_m_updated);
            *(oris + nEpoch * OriBlockSize + 0) = q.x();
            *(oris + nEpoch * OriBlockSize + 1) = q.y();
            *(oris + nEpoch * OriBlockSize + 2) = q.z();
            *(oris + nEpoch * OriBlockSize + 3) = q.w();

            *(poss_std + nEpoch * PosBlockSize + 0) = 1.0E+20;
            *(poss_std + nEpoch * PosBlockSize + 1) = 1.0E+20;
            *(poss_std + nEpoch * PosBlockSize + 2) = 1.0E+20;
            *(oris_std + nEpoch * OriBlockSize + 0) = 1.0E+20;
            *(oris_std + nEpoch * OriBlockSize + 1) = 1.0E+20;
            *(oris_std + nEpoch * OriBlockSize + 2) = 1.0E+20;
            *(oris_std + nEpoch * OriBlockSize + 3) = 1.0E+20;
        }

        // 2. cylinder parameter
        int countCylinder = 0;
        for (int nTree = 0; nTree < numCylinder; nTree++)
        {
            int map_tree_id = used_trees_[nTree];
            *(cylinders + nTree * CylinderBlockSize + 0) = mpMapTree[map_tree_id]->para.x(0); // XYZ (Z is nominal)
            *(cylinders + nTree * CylinderBlockSize + 1) = mpMapTree[map_tree_id]->para.x(1);
            *(cylinders + nTree * CylinderBlockSize + 2) = mpMapTree[map_tree_id]->para.x(2);
            *(cylinders + nTree * CylinderBlockSize + 3) = mpMapTree[map_tree_id]->para.n(0);
            *(cylinders + nTree * CylinderBlockSize + 4) = mpMapTree[map_tree_id]->para.n(1);
            *(cylinders + nTree * CylinderBlockSize + 5) = mpMapTree[map_tree_id]->para.n(2);
            *(cylinders + nTree * CylinderBlockSize + 6) = mpMapTree[map_tree_id]->para.r;

            *(cylinders_std + nTree * CylinderBlockSize + 0) = 1.0E+20; // XYZ
            *(cylinders_std + nTree * CylinderBlockSize + 1) = 1.0E+20;
            *(cylinders_std + nTree * CylinderBlockSize + 2) = 1.0E-20;
            *(cylinders_std + nTree * CylinderBlockSize + 3) = 1.0E+20; // ux, uy, uz
            *(cylinders_std + nTree * CylinderBlockSize + 4) = 1.0E+20;
            *(cylinders_std + nTree * CylinderBlockSize + 5) = 1.0E-20;
            *(cylinders_std + nTree * CylinderBlockSize + 6) = 1.0E+20; // r
        }

        // 3. plane parameters
        for (int nPlane = 0; nPlane < numPlane; nPlane++)
        {
            Eigen::Vector4d params = plane_params_[nPlane];
            int fixIndex = planeParamTrans(params);

            *(planes + nPlane * 4 + 0) = params(0);
            *(planes + nPlane * 4 + 1) = params(1);
            *(planes + nPlane * 4 + 2) = params(2);
            *(planes + nPlane * 4 + 3) = params(3);

            *(planes_std + nPlane * 4 + 0) = 1.0E+20;
            *(planes_std + nPlane * 4 + 1) = 1.0E+20;
            *(planes_std + nPlane * 4 + 2) = 1.0E+20;
            *(planes_std + nPlane * 4 + 3) = 1.0E+20;

            // fixing the element equal to 1
            *(planes_std + nPlane * 4 + fixIndex) = 1.0E-20;
        }

        //------------- check fixed parameters
        int numFixedUnknown = 0;
        for (int i = 0; i < numPosUnknown; i++)
        {
            double temp;
            temp = *(poss_std + i);
            if (temp < SigmaFixed)
                numFixedUnknown += 1;
        }
        for (int i = 0; i < numOriUnknown; i++)
        {
            double temp;
            temp = *(oris_std + i);
            if (temp < SigmaFixed)
                numFixedUnknown += 1;
        }
        for (int i = 0; i < numCylinderUnknown; i++)
        {
            double temp;
            temp = *(cylinders_std + i);
            if (temp < SigmaFixed)
                numFixedUnknown += 1;
        }
        for (int i = 0; i < numPlaneUnknown; i++)
        {
            double temp;
            temp = *(planes_std + i);
            if (temp < SigmaFixed)
                numFixedUnknown += 1;
        }

        //--------------------------build ceres  ------------------------------
        ceres::LossFunction *loss_function = new ceres::HuberLoss(0.1);
        ceres::LocalParameterization *q_parameterization = new ceres::EigenQuaternionParameterization(); // x, y, z, w
        ceres::Problem problem ;//= new ceres::Problem;
        ceres::CostFunction *costFunction;
        //--------------------------add parameter block ------------------------------
        for (int nEpoch = 0; nEpoch < numEpoch; nEpoch++)
        {
            problem.AddParameterBlock(poss + nEpoch * PosBlockSize, PosBlockSize);
        }
        for (int nEpoch = 0; nEpoch < numEpoch; nEpoch++)
        {
            problem.AddParameterBlock(oris + nEpoch * OriBlockSize, OriBlockSize, q_parameterization);
        }
        for (int nCylinder = 0; nCylinder < numCylinder; nCylinder++)
        {
            problem.AddParameterBlock(cylinders + nCylinder * CylinderBlockSize, CylinderBlockSize);
        }
        for (int nPlane = 0; nPlane < numPlane; nPlane++)
        {
            problem.AddParameterBlock(planes + nPlane * PlaneBlockSize, PlaneBlockSize);
        }

        //-------------------------add observation -------------------------------
        // 1. Tree observation
        int numFixTreeObs = 0;
        int numUnknownTreeObs = 0;
        std::set<uint64_t> downsample_index_set;
        int original_tree_points_number = 0;
        int downsample_tree_points_number =0;
        for (int nCylinder = 0; nCylinder < numCylinder; nCylinder++)
        {
            int cylinderId = nCylinder;
            int map_tree_id = used_trees_[nCylinder];
            double *cylinder = cylinders + cylinderId * CylinderBlockSize; // corresponding cylinder block
            
            int increment = 1;
            if(mpMapTree[map_tree_id]->numPoint != 0)
            {
                increment = ceil(mpMapTree[map_tree_id]->numPoint/max_num_map_point_per_tree);
                increment = max(1, increment);
            }
            for (int nT = 0; nT < mpMapTree[map_tree_id]->visibleInfo.size(); nT++)
            {
                int iscanId = mpMapTree[map_tree_id]->visibleInfo[nT].first;
                int treeId = mpMapTree[map_tree_id]->visibleInfo[nT].second;
                int epoch_id = iscanId - iscan_index_start_;

                // if iscan level, the iscan is failed
                if (!mvpIScans[iscanId]->flag_iscan_to_map)
                    continue;

                //modification: dont apply uniform scale to each iscan
                int num_p = mvpIScans[iscanId]->vTreePointMapping[treeId].size();
                //int incre_individual = num_p < Min_Tree_Points_Per_Iscan_LC? 1 : increment;
                int incre_individual = increment;
                
                for (int nP = 0; nP < mvpIScans[iscanId]->vTreePointMapping[treeId].size(); nP = nP + incre_individual)
                {
                    original_tree_points_number++;
                    //if from the iscan that is fix
                    if(iscanId < iscan_index_start_)
                    {
                        // points in mapping
                        Eigen::Vector3d r_I_m(mvpIScans[iscanId]->vTreePointMapping[treeId][nP].x,
                                                mvpIScans[iscanId]->vTreePointMapping[treeId][nP].y,
                                                mvpIScans[iscanId]->vTreePointMapping[treeId][nP].z);
                        uint64_t key=GetVoxelkeyIndex(r_I_m,kDownsamplingSize);
                        if(downsample_index_set.find(key)==downsample_index_set.end())
                        {
                            downsample_index_set.insert(key);
                        }
                        else
                        {
                            continue;
                        }
                        //apply the scale to std
                        downsample_tree_points_number++;
                        costFunction = CylinderFactorObjectPoint::Create(r_I_m, tree_std / sqrt(increment));
                        //costFunction = CylinderFactorObjectPoint::Create(r_I_m, tree_std);
                        problem.AddResidualBlock(costFunction, loss_function, cylinder);
                        numFixTreeObs++;
                    }
                    else
                    {
                        double *t_cur = poss + epoch_id * PosBlockSize;
                        double *q_cur = oris + epoch_id * OriBlockSize;
                        Eigen::Vector3d r_I_local(mvpIScans[iscanId]->vTreePointLocal[treeId][nP].x,
                                                mvpIScans[iscanId]->vTreePointLocal[treeId][nP].y,
                                                mvpIScans[iscanId]->vTreePointLocal[treeId][nP].z);
                        Eigen::Vector3d point_in_mapping = mvpIScans[iscanId]->R_local_m_updated * r_I_local
                        + mvpIScans[iscanId]->r_local_m_updated;
                        uint64_t key=GetVoxelkeyIndex(point_in_mapping,kDownsamplingSize);
                        if(downsample_index_set.find(key)==downsample_index_set.end())
                        {
                            downsample_index_set.insert(key);
                        }
                        else
                        {
                            continue;
                        }
                        downsample_tree_points_number++;
                        costFunction = CylinderFactorOneEpoch::Create(r_I_local, tree_std / sqrt(increment));
                        //costFunction = CylinderFactorOneEpoch::Create(r_I_local, tree_std);
                        problem.AddResidualBlock(costFunction, loss_function, t_cur, q_cur, cylinder);
                        numUnknownTreeObs++;
                    }
                }
            }
        }

        std::cout<<"LC original_tree_points_number: "<<original_tree_points_number<<std::endl;
        std::cout<<"LC downsample_tree_points_number: "<<downsample_tree_points_number<<std::endl;


        // 2. ground observation
        int numRefGroundObs = 0;
        int numFixGroundObs = 0;
        int numUnknownGroundObs = 0;
        int original_ground_points_number = 0;
        int downsample_ground_points_number =0;
        std::set<uint64_t> downsample_index_set_ground;
        for (int nPlane = 0; nPlane < numPlane; nPlane++)
        {
            int planeId = nPlane;
            double *plane = planes + planeId * PlaneBlockSize; // corresponding plane block
            for (int nP = 0; nP < point_index_per_grid_[nPlane].size(); nP++)
            {
                int iscanId = point_index_per_grid_[nPlane][nP].first;
                int ptId = point_index_per_grid_[nPlane][nP].second;
                int epoch_id = iscanId - iscan_index_start_;
                original_ground_points_number++;
                // if from the iscan that is fix
                if (iscanId < iscan_index_start_)
                {
                    // points in mapping
                    Eigen::Vector3d r_I_m(mvpIScans[iscanId]->pGroundPointMappingDs->points[ptId].x,
                                          mvpIScans[iscanId]->pGroundPointMappingDs->points[ptId].y,
                                          mvpIScans[iscanId]->pGroundPointMappingDs->points[ptId].z);

                    uint64_t key=GetVoxelkeyIndex(r_I_m,KDownsampleLCGround);
                    if(downsample_index_set_ground.find(key)==downsample_index_set_ground.end())
                    {
                        downsample_index_set_ground.insert(key);
                    }
                    else
                    {
                        continue;
                    }
                    downsample_ground_points_number++;
                    costFunction = PlaneFactorObjectPoint::Create(r_I_m, ground_std);
                    problem.AddResidualBlock(costFunction, loss_function, plane);
                    numFixGroundObs++;
                }
                else
                {
                    double *t_cur = poss + epoch_id * PosBlockSize;
                    double *q_cur = oris + epoch_id * OriBlockSize;

                    Eigen::Vector3d r_I_local(mvpIScans[iscanId]->pGroundPointLocalDs->points[ptId].x,
                                              mvpIScans[iscanId]->pGroundPointLocalDs->points[ptId].y,
                                              mvpIScans[iscanId]->pGroundPointLocalDs->points[ptId].z);
                    Eigen::Vector3d point_in_mapping = mvpIScans[iscanId]->R_local_m_updated * r_I_local
                    + mvpIScans[iscanId]->r_local_m_updated;
                    uint64_t key=GetVoxelkeyIndex(point_in_mapping,KDownsampleLCGround);
                    if(downsample_index_set_ground.find(key)==downsample_index_set_ground.end())
                    {
                        downsample_index_set_ground.insert(key);
                    }
                    else
                    {
                        continue;
                    }
                    downsample_ground_points_number++;
                    costFunction = PlaneFactorOneEpoch::Create(r_I_local, ground_std);
                    problem.AddResidualBlock(costFunction, loss_function, t_cur, q_cur, plane);
                    numUnknownGroundObs++;
                }
            }

            for (int nP = 0; nP < point_ref_per_grid_[nPlane].size(); nP++)
            {

                costFunction = PlaneFactorObjectPoint::Create(point_ref_per_grid_[nPlane][nP], ref_ground_std);
                problem.AddResidualBlock(costFunction, loss_function, plane);
                numRefGroundObs++;
            }
        }

        std::cout<<"LC original_ground_points_number: "<<original_ground_points_number<<std::endl;
        std::cout<<"LC downsample_ground_points_number: "<<downsample_ground_points_number<<std::endl;

        // 3. Global tree
        int num_global_tree = 0;
        for (int nCylinder = 0; nCylinder < numCylinder; nCylinder++)
        {
            int cylinderId = nCylinder;
            int map_tree_id = used_trees_[nCylinder];
            double *cylinder = cylinders + cylinderId * CylinderBlockSize; // corresponding cylinder block
           
            if (mpMapTree[map_tree_id]->map_tree_flag)
            {
                double *cylinder = cylinders + cylinderId * CylinderBlockSize;
                //costFunction = CylinderPriorConstraints::Create(mpMapTree[mapTreeId]->para_ref_.x, global_tree_loc_std);
                costFunction = CylinderPriorConstraintsPointToLine::Create(mpMapTree[map_tree_id]->para_ref_.x, global_tree_loc_std);
                problem.AddResidualBlock(costFunction, nullptr, cylinder);                
                num_global_tree++;
            }
        }

        // 4. Edge constraint
        if (1)
        {
            for (int nEpoch = 0; nEpoch < numEpoch - 1; nEpoch++)
            {
                // for each scan, epoches from [nScan, nScan+1] are used
                double *t_1 = poss + nEpoch * PosBlockSize;
                double *q_1 = oris + nEpoch * OriBlockSize;
                double *t_2 = poss + (nEpoch + 1) * PosBlockSize;
                double *q_2 = oris + (nEpoch + 1) * OriBlockSize;

                int iscan_id_1 = nEpoch + iscan_index_start_;
                int iscan_id_2 = nEpoch + iscan_index_start_ + 1;

                // transformation estimated by odometry
                Eigen::Quaterniond q_prior(mvpIScans[iscan_id_2]->R_local_t2_local_t1);
                Eigen::Vector3d t_prior = mvpIScans[iscan_id_2]->r_local_t2_local_t1;

                // derive sqrt of information matrix
                // Sample: from CERES const Eigen::Matrix<double, 6, 6> sqrt_information =constraint.information.llt().matrixL();
                                
                double acc_t, acc_q; /// 0.1 deg
                if (mPara.odo_from_trajectory_flag)
                {
                    acc_t = mPara.odoPositionStdTraj;
                    acc_q = deg2rad(mPara.odoOrientationStdTraj);
                }
                else
                {
                    acc_t = mPara.odoPositionStd;
                    acc_q = deg2rad(mPara.odoOrientationStd);
                }

                Eigen::Matrix<double, 6, 6> sqrt_information = Eigen::MatrixXd::Identity(6, 6);
                sqrt_information(0, 0) = 1.0 / acc_t;
                sqrt_information(1, 1) = 1.0 / acc_t;
                sqrt_information(2, 2) = 1.0 / acc_t;
                sqrt_information(3, 3) = 1.0 / acc_q;
                sqrt_information(4, 4) = 1.0 / acc_q;
                sqrt_information(5, 5) = 1.0 / acc_q;

                costFunction = PoseGraph3dErrorTerm::Create(t_prior, q_prior, sqrt_information);
                
				if (mPara.odo_from_trajectory_flag)
				{
					problem.AddResidualBlock(costFunction, nullptr, t_1, q_1, t_2, q_2);
				}
				else
				{
					problem.AddResidualBlock(costFunction, loss_function, t_1, q_1, t_2, q_2);
				}

                double dis_std = 0.5 ;
				costFunction = DistanceBetweenTwoEpochsConstraints::Create(t_prior.norm(),dis_std);
                problem.AddResidualBlock(costFunction, nullptr, t_1, q_1, t_2, q_2);


                //        Eigen::Matrix3d R_lu_local_t2 = R_lu_local_t1 * mvTempScan[0].R_lut2_lut1;
                //    Eigen::Vector3d r_lu_local_t2 = r_lu_local_t1 + R_lu_local_t1 * mvTempScan[0].r_lut2_lut1;
            }
        }
        //--------------------------Fix Parameters  ------------------------------
        // set constant unknowns
        Checkforconstantparams(&problem, poss, poss_std, PosBlockSize, numEpoch);
        Checkforconstantparams(&problem, oris, oris_std, OriBlockSize, numEpoch);
        Checkforconstantparams(&problem, cylinders, cylinders_std, CylinderBlockSize, numCylinder);
        Checkforconstantparams(&problem, planes, planes_std, PlaneBlockSize, numPlane);

        TicToc t_solver;
        ceres::Solver::Options solverOptions;
        solverOptions.max_num_iterations = 10;
        solverOptions.linear_solver_type = ceres::LinearSolverType::SPARSE_NORMAL_CHOLESKY; // ceres::LinearSolverType::SPARSE_NORMAL_CHOLESKY;
        solverOptions.minimizer_progress_to_stdout = false;
        ceres::Solver::Summary summary;
        ceres::Solve(solverOptions, &problem, &summary);
        fMapLog << ", Time: " << t_solver.toc() << " ms" << endl;
        fMapLog << "\tFix/unknown tree observations: " << numFixTreeObs << " " << numUnknownTreeObs << endl;
        fMapLog << "\tRef/fix/unknown ground observations: " << numRefGroundObs << " " << numFixGroundObs << " " << numUnknownGroundObs<< endl;
        fMapLog << "\tNumber of pair with global tree: " << num_global_tree<< endl;

        cout << "Time: " << t_solver.toc() << endl;

        //******************** Update parameters **********************
        // 1. trajectory info
        for (int nEpoch = 0; nEpoch < numEpoch; nEpoch++)
        {
            Eigen::Vector3d r{*(poss + nEpoch * PosBlockSize + 0), *(poss + nEpoch * PosBlockSize + 1), *(poss + nEpoch * PosBlockSize + 2)};
            Eigen::Quaterniond q_refined(*(oris + nEpoch * OriBlockSize + 3), *(oris + nEpoch * OriBlockSize + 0), *(oris + nEpoch * OriBlockSize + 1), *(oris + nEpoch * OriBlockSize + 2));
            Eigen::Matrix3d tempR, orthoR;

            tempR = q_refined.toRotationMatrix();
            Eigen::Vector3d angles;
            angles = Find_Rotation(tempR);
            Compute_Rotation(angles(0), angles(1), angles(2), orthoR);

            int iscan_id = nEpoch + iscan_index_start_;
            mvpIScans[iscan_id]->r_local_m_temp = mvpIScans[iscan_id]->r_local_m_updated;
            mvpIScans[iscan_id]->R_local_m_temp = mvpIScans[iscan_id]->R_local_m_updated;

            mvpIScans[iscan_id]->r_local_m_updated = r;
            mvpIScans[iscan_id]->R_local_m_updated = orthoR;
        }

        //******************** Residuals  **********************
        // 1. Tree observation
        double res_tree_map = 0.0, res_tree_unknown = 0.0, res_tree_map_ini = 0.0, res_tree_unknown_ini = 0.0;
        CylinderPara treePara;
        int num_valid_tree_para = 0, num_invalid_tree_para = 0;
        for (int nCylinder = 0; nCylinder < numCylinder; nCylinder++)
        {
            int cylinderId = nCylinder;
            int map_tree_id = used_trees_[nCylinder];
            double *cylinder = cylinders + cylinderId * CylinderBlockSize; // corresponding cylinder block
            treePara.x(0) = *(cylinders + cylinderId * CylinderBlockSize + 0);
            treePara.x(1) = *(cylinders + cylinderId * CylinderBlockSize + 1);
            treePara.x(2) = *(cylinders + cylinderId * CylinderBlockSize + 2);
            treePara.n(0) = *(cylinders + cylinderId * CylinderBlockSize + 3);
            treePara.n(1) = *(cylinders + cylinderId * CylinderBlockSize + 4);
            treePara.n(2) = *(cylinders + cylinderId * CylinderBlockSize + 5);
            treePara.r = *(cylinders + cylinderId * CylinderBlockSize + 6);

            //compute residual
            int increment = 1;
            if(mpMapTree[map_tree_id]->numPoint != 0)
            {
                increment = ceil(mpMapTree[map_tree_id]->numPoint/max_num_map_point_per_tree);
                increment = max(1, increment);
            }

            for (int nT = 0; nT < mpMapTree[map_tree_id]->visibleInfo.size(); nT++)
            {
                int iscan_id = mpMapTree[map_tree_id]->visibleInfo[nT].first;
                int treeId = mpMapTree[map_tree_id]->visibleInfo[nT].second;
                int epoch_id = iscan_id - iscan_index_start_;

                // if iscan level, the iscan is failed
                if (!mvpIScans[iscan_id]->flag_iscan_to_map)
                    continue;

                int num_p = mvpIScans[iscan_id]->vTreePointMapping[treeId].size();
                //int incre_individual = num_p < Min_Tree_Points_Per_Iscan_LC? 1 : increment;
                int incre_individual = increment;
                for (int nP = 0; nP < mvpIScans[iscan_id]->vTreePointMapping[treeId].size(); nP = nP + incre_individual)
                {
                    //if from the iscan that is fix
                    if(iscan_id < iscan_index_start_)
                    {
                        // points in mapping
                        Eigen::Vector3d r_I_m(mvpIScans[iscan_id]->vTreePointMapping[treeId][nP].x,
                                                mvpIScans[iscan_id]->vTreePointMapping[treeId][nP].y,
                                                mvpIScans[iscan_id]->vTreePointMapping[treeId][nP].z);
                        
                        //apply the scale to std
                        double d_ini = computePoint2CylinderDistance(r_I_m,mpMapTree[map_tree_id]->para);
                        double d_ref = computePoint2CylinderDistance(r_I_m,treePara);
                        res_tree_map_ini = res_tree_map_ini + d_ini * d_ini;
                        res_tree_map = res_tree_map + d_ref * d_ref;
                    }
                    else
                    {
                        Eigen::Vector3d r_I_local(mvpIScans[iscan_id]->vTreePointLocal[treeId][nP].x,
                                                  mvpIScans[iscan_id]->vTreePointLocal[treeId][nP].y,
                                                  mvpIScans[iscan_id]->vTreePointLocal[treeId][nP].z);

                        double d_ini = computePoint2CylinderDistance(mvpIScans[iscan_id]->R_local_m_temp * r_I_local + mvpIScans[iscan_id]->r_local_m_temp,
                                                                     mpMapTree[map_tree_id]->para);
                        double d_ref = computePoint2CylinderDistance(mvpIScans[iscan_id]->R_local_m_updated * r_I_local + mvpIScans[iscan_id]->r_local_m_updated,
                                                                     treePara);
                        res_tree_unknown_ini = res_tree_unknown_ini + d_ini * d_ini;
                        res_tree_unknown = res_tree_unknown + d_ref * d_ref;
                    }
                }
            }

            //update tree parameter
            mpMapTree[map_tree_id]->UpdateTreeParam(treePara)?num_valid_tree_para++: num_invalid_tree_para++;
        }
        fMapLog << "\tNumber of trees with valid and invalid parameters: " << num_valid_tree_para << " / " << num_invalid_tree_para <<endl;
        fMapLog << "\tError for tree point from fix Iscans: " << sqrt(res_tree_map_ini / double(numFixTreeObs)) << " -> " << sqrt(res_tree_map / double(numFixTreeObs)) << " m" <<endl;
        fMapLog << "\tError for tree point from estimated Iscans: " << sqrt(res_tree_unknown_ini / double(numUnknownTreeObs)) << " -> " << sqrt(res_tree_unknown / double(numUnknownTreeObs)) << " m" <<endl;

        // 2. ground observation
        double res_ground_ref = 0.0, res_ground_fix = 0.0, res_ground_unknown = 0.0, res_ground_ref_ini = 0.0, res_ground_fix_ini = 0.0, res_ground_unknown_ini = 0.0;
        for (int nPlane = 0; nPlane < numPlane; nPlane++)
        {
            int planeId = nPlane;
            Eigen::Vector4d updated_para{*(planes + planeId * PlaneBlockSize + 0), *(planes + planeId * PlaneBlockSize + 1),
                                         *(planes + planeId * PlaneBlockSize + 2), *(planes + planeId * PlaneBlockSize + 3)};

            for (int nP = 0; nP < point_index_per_grid_[nPlane].size(); nP++)
            {
                int iscan_id = point_index_per_grid_[nPlane][nP].first;
                int ptId = point_index_per_grid_[nPlane][nP].second;
                int epoch_id = iscan_id - iscan_index_start_;

                // if from the iscan that is fix
                if (iscan_id < iscan_index_start_)
                {
                    // points in mapping
                    Eigen::Vector3d r_I_m(mvpIScans[iscan_id]->pGroundPointLocalDs->points[ptId].x,
                                          mvpIScans[iscan_id]->pGroundPointLocalDs->points[ptId].y,
                                          mvpIScans[iscan_id]->pGroundPointLocalDs->points[ptId].z);
                    double d_ini = (r_I_m.transpose() * plane_params_[nPlane].head(3) + plane_params_[nPlane](3)) / plane_params_[nPlane].head(3).norm();
                    double d_ref = (r_I_m.transpose() * updated_para.head(3) + updated_para(3)) / updated_para.head(3).norm();
                    res_ground_fix_ini = res_ground_fix_ini + d_ini * d_ini;
                    res_ground_fix = res_ground_fix + d_ref * d_ref;
                }
                else
                {
                    Eigen::Vector3d r_I_local(mvpIScans[iscan_id]->pGroundPointLocalDs->points[ptId].x,
                                              mvpIScans[iscan_id]->pGroundPointLocalDs->points[ptId].y,
                                              mvpIScans[iscan_id]->pGroundPointLocalDs->points[ptId].z);
                    double d_ini = ((mvpIScans[iscan_id]->R_local_m_temp * r_I_local + mvpIScans[iscan_id]->r_local_m_temp).transpose() * plane_params_[nPlane].head(3) + plane_params_[nPlane](3)) / plane_params_[nPlane].head(3).norm();
                    double d_ref = ((mvpIScans[iscan_id]->R_local_m_updated * r_I_local + mvpIScans[iscan_id]->r_local_m_updated).transpose() * updated_para.head(3) + updated_para(3)) / updated_para.head(3).norm();
                    res_ground_unknown_ini = res_ground_unknown_ini + d_ini * d_ini;
                    res_ground_unknown = res_ground_unknown + d_ref * d_ref;
                }
            }
            for (int nP = 0; nP < point_ref_per_grid_[nPlane].size(); nP++)
            {
                Eigen::Vector3d r_I_m = point_ref_per_grid_[nPlane][nP];
                double d_ini = (r_I_m.transpose() * plane_params_[nPlane].head(3) + plane_params_[nPlane](3)) / plane_params_[nPlane].head(3).norm();
                double d_ref = (r_I_m.transpose() * updated_para.head(3) + updated_para(3)) / updated_para.head(3).norm();
                res_ground_ref_ini = res_ground_ref_ini + d_ini * d_ini;
                res_ground_ref = res_ground_ref + d_ref * d_ref;
            }
        }
        fMapLog << "\tError for ground point for ref planar patches: " << sqrt(res_ground_ref_ini / double(numRefGroundObs)) << " -> " << sqrt(res_ground_ref / double(numRefGroundObs)) << " m" << endl;
        fMapLog << "\tError for ground point for fix planar patches: " << sqrt(res_ground_fix_ini / double(numFixGroundObs)) << " -> " << sqrt(res_ground_fix / double(numFixGroundObs)) << " m" << endl;
        fMapLog << "\tError for ground point for unknown planar patches: " << sqrt(res_ground_unknown_ini / double(numUnknownGroundObs)) << " -> " << sqrt(res_ground_unknown / double(numUnknownGroundObs)) << " m" << endl;


    }//end of iteration
}

// uint64_t GetVoxelkeyIndex(const Eigen::Vector3d& point,double voxel_size)
// {
// 	double cell_size = voxel_size;
// 	uint64_t length = 65535;
// 	uint64_t x_index = static_cast<uint64_t>(point(0) / cell_size);
// 	uint64_t y_index = static_cast<uint64_t>(point(1) / cell_size);
// 	uint64_t z_index = static_cast<uint64_t>(point(2) / cell_size);
// 	uint64_t x_index64 = static_cast<uint64_t>(x_index%length);
// 	uint64_t y_index64 = static_cast<uint64_t>(y_index%length);
// 	uint64_t z_index64 = static_cast<uint64_t>(z_index%length);
// 	uint64_t key_index = (x_index64 << 32) + (y_index64 << 16) + z_index64;
// 	return key_index;
// }

/**
 * @brief Computes the normal vector of a set of 3D points and evaluates planarity.
 * 
 * @param[in] data Vector of 3D points.
 * @param[out] parameters Vector of 6 elements: normal vector (3) and centroid (3).
 * @return Planarity score as per eigenvalue-based descriptor.
 *
 * This function calculates the covariance matrix and derives its eigenvectors 
 * to estimate the normal direction and center. It returns the planarity measure
 * using the square roots of the eigenvalues.
 */
double ComputePointsNormal(const std::vector< Eigen::Vector3d>& data, std::vector<double> &parameters)
{
	if (data.size() < 3)
	{
		return false;
	}
	Eigen::Matrix3d covariance_matrix = Eigen::Matrix3d::Zero();
	Eigen::Matrix<double, 1, 9, Eigen::RowMajor> inter_matrix = Eigen::Matrix<double, 1, 9, Eigen::RowMajor>::Zero();
	size_t point_count = data.size();
	for (int i = 0; i < point_count; ++i)
	{
		inter_matrix[0] += data[i](0) * data[i](0);
		inter_matrix[1] += data[i](0) * data[i](1);
		inter_matrix[2] += data[i](0) * data[i](2);
		inter_matrix[3] += data[i](1) * data[i](1);
		inter_matrix[4] += data[i](1) * data[i](2);
		inter_matrix[5] += data[i](2) * data[i](2);
		inter_matrix[6] += data[i](0);
		inter_matrix[7] += data[i](1);
		inter_matrix[8] += data[i](2);
	}

	inter_matrix /= static_cast<double> (point_count);
	covariance_matrix.coeffRef(0) = inter_matrix[0] - inter_matrix[6] * inter_matrix[6];
	covariance_matrix.coeffRef(1) = inter_matrix[1] - inter_matrix[6] * inter_matrix[7];
	covariance_matrix.coeffRef(2) = inter_matrix[2] - inter_matrix[6] * inter_matrix[8];
	covariance_matrix.coeffRef(4) = inter_matrix[3] - inter_matrix[7] * inter_matrix[7];
	covariance_matrix.coeffRef(5) = inter_matrix[4] - inter_matrix[7] * inter_matrix[8];
	covariance_matrix.coeffRef(8) = inter_matrix[5] - inter_matrix[8] * inter_matrix[8];
	covariance_matrix.coeffRef(3) = covariance_matrix.coeff(1);
	covariance_matrix.coeffRef(6) = covariance_matrix.coeff(2);
	covariance_matrix.coeffRef(7) = covariance_matrix.coeff(5);

	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::Matrix3d> solver(covariance_matrix, Eigen::Matrix3d::Identity(), Eigen::ComputeEigenvectors | Eigen::Ax_lBx);

	parameters.clear();
	parameters.resize(6);
	parameters[0] = solver.eigenvectors().col(0)[0];
	parameters[1] = solver.eigenvectors().col(0)[1];
	parameters[2] = solver.eigenvectors().col(0)[2];
	parameters[3] = inter_matrix(6);
	parameters[4] = inter_matrix(7);
	parameters[5] = inter_matrix(8);

	Eigen::VectorXd eigenvalues = solver.eigenvalues();
	std::sort(eigenvalues.data(), eigenvalues.data() + eigenvalues.size(), std::greater<double>());

	double lamda1 = eigenvalues(0);
	double lamda2 = eigenvalues(1);
	double lamda3 = eigenvalues(2);
	double planarity = (std::sqrt(lamda2) - std::sqrt(lamda3)) / std::sqrt(lamda1);
	return planarity;
}

/**
 * @brief Performs raw-level loop closure optimization using planar surface elements.
 *
 * This method voxelizes the environment using raw tree points, extracts dominant 
 * planar patches (based on eigenvalue decomposition), and builds a Ceres optimization 
 * problem that refines poses and plane parameters. It uses odometry constraints 
 * and surface observations to optimize alignment across IScans.
 */

void Mapping::OptimizeRawLevelLCUsingSurfaceElements()
{
    //map saving the start pose index of a Iscan
    unordered_map<int, int> iscan_pose_map;
    int pose_count = 0;
    for (int nIscan = iscan_index_start_; nIscan <=  iscan_index_end_; nIscan++)
    {
        iscan_pose_map[nIscan] = pose_count;
        pose_count += mvpIScans[nIscan]->numScan;
    }
    pose_count ++;
    struct SurfacePoint
    {
        Eigen::Vector3d point;
        int epoch_id;
	};
    std::map<uint64_t, std::vector<SurfacePoint>> voxel_key_to_inner_points;
	std::map<uint64_t, std::vector<Eigen::Vector3d>> voxel_key_to_inner_points_mapping;
    double voxel_size = 0.10; //0.2 0.1
    int voxel_points_number_threshold = 10; //30 10
    for(int iscanId =0 ;iscanId< mvpIScans.size();iscanId++)
    {
        for(const auto& tree : mvpIScans[iscanId]->vTreePointRaw)
        {
            for(const auto& tree_p : tree)
            {
                int scanId = int(tree_p.intensity / 100);
                int epoch_id = iscan_pose_map[iscanId] + scanId;
                SurfacePoint s_p;
                s_p.point=tree_p.getVector3fMap().cast<double>();
                s_p.epoch_id = epoch_id;

				Eigen::Vector3d point_in_mapping = mvpIScans[iscanId]->v_R_mapping[scanId] * s_p.point 
					+ mvpIScans[iscanId]->v_r_mapping[scanId];


                uint64_t key = GetVoxelkeyIndex(point_in_mapping, voxel_size);
		        voxel_key_to_inner_points[key].push_back(s_p);
				voxel_key_to_inner_points_mapping[key].push_back(point_in_mapping);
            }
        }

   //     for(const auto& ground_p : mvpIScans[iscanId]->vGroundPointRaw)
   //     {
   //         int scanId = int(ground_p.intensity / 100);
   //         int epoch_id = iscan_pose_map[iscanId] + scanId;
   //         SurfacePoint s_p;
   //         s_p.point=ground_p.getVector3fMap().cast<double>();
   //         s_p.epoch_id = epoch_id;

			//Eigen::Vector3d point_in_mapping = mvpIScans[iscanId]->v_R_mapping[scanId] * s_p.point
			//	+ mvpIScans[iscanId]->v_r_mapping[scanId];

   //         uint64_t key = GetVoxelkeyIndex(point_in_mapping, voxel_size);
   //         voxel_key_to_inner_points[key].push_back(s_p);
			//voxel_key_to_inner_points_mapping[key].push_back(point_in_mapping);
   //     }
    }
    std::vector<Eigen::Vector3d> voxel_points;
    std::vector<Eigen::Vector4d> surface_paramenters_vector;
    std::vector<uint64_t> surface_keys;
    int surface_element_number = 0;
    for (const auto& item : voxel_key_to_inner_points_mapping)
	{
		if (item.second.size() < voxel_points_number_threshold)
		{
			continue;
		}
		voxel_points.clear();
		for (int i = 0; i < item.second.size(); i++)
		{
			voxel_points.push_back(item.second[i]);
		}
		std::vector<double> parameters;
		double planarity = ComputePointsNormal(voxel_points, parameters);
		// if (planarity < planarity_threshold)
		// {
		// 	//continue;
		// }
		Eigen::Vector3d plane_center(parameters[3], parameters[4], parameters[5]);
		Eigen::Vector3d plane_normal(parameters[0], parameters[1], parameters[2]);
        double d = -plane_normal.dot(plane_center); // 计算平面常数项
        Eigen::Vector4d surface_ele(plane_normal(0),plane_normal(1),plane_normal(2),d);
		surface_paramenters_vector.push_back(surface_ele);
        surface_keys.push_back(item.first);
        surface_element_number++;
	}
	std::map<uint64_t, std::vector<Eigen::Vector3d>> empty_map;
	voxel_key_to_inner_points_mapping.swap(empty_map);

    double surface_std = mPara.mapGroundStd;

    int iter = 0;
    int maxIter = 1;
    while (iter < maxIter)
    {
        fMapLog << "Iter " << ++iter;
        fMapLog << "\tNumber of included trees " << used_trees_.size();
        
        // initialization
        int numEpoch = pose_count;
        int numSurface = surface_paramenters_vector.size();

        // unknowns
        int numPosUnknown = numEpoch * PosBlockSize;
        int numOriUnknown = numEpoch * OriBlockSize;
        int numPlaneUnknown = numSurface * PlaneBlockSize;
        int numUnknown = numPosUnknown + numOriUnknown +  numPlaneUnknown; // all unknowns
        double *poss = new double[numPosUnknown];
        double *poss_std = new double[numPosUnknown];
        double *oris = new double[numOriUnknown];
        double *oris_std = new double[numOriUnknown];
        double *planes = new double[numPlaneUnknown];
        double *planes_std = new double[numPlaneUnknown];

        // fill parameter
        // 1. trajectory information
        vector<int> check_(numEpoch, 0);
        for (int nIscan = iscan_index_start_; nIscan <= iscan_index_end_; nIscan++)
        {
            int iscan_id = nIscan;
            int pose_start = iscan_pose_map[iscan_id];

            // for the last Iscan, use all pose from v_R/r_mapping
            int num_pose_per_iscan = (nIscan == iscan_index_end_) ? mvpIScans[nIscan]->numScan + 1 : mvpIScans[nIscan]->numScan;
            for (int nPose = 0; nPose < num_pose_per_iscan; nPose++)
            {
                int pose_id = pose_start + nPose;
                *(poss + pose_id * PosBlockSize + 0) = mvpIScans[iscan_id]->v_r_mapping[nPose](0);
                *(poss + pose_id * PosBlockSize + 1) = mvpIScans[iscan_id]->v_r_mapping[nPose](1);
                *(poss + pose_id * PosBlockSize + 2) = mvpIScans[iscan_id]->v_r_mapping[nPose](2);

                Eigen::Quaterniond q(mvpIScans[iscan_id]->v_R_mapping[nPose]);
                *(oris + pose_id * OriBlockSize + 0) = q.x();
                *(oris + pose_id * OriBlockSize + 1) = q.y();
                *(oris + pose_id * OriBlockSize + 2) = q.z();
                *(oris + pose_id * OriBlockSize + 3) = q.w();

                *(poss_std + pose_id * PosBlockSize + 0) = 1.0E+20;
                *(poss_std + pose_id * PosBlockSize + 1) = 1.0E+20;
                *(poss_std + pose_id * PosBlockSize + 2) = 1.0E+20;
                *(oris_std + pose_id * OriBlockSize + 0) = 1.0E+20;
                *(oris_std + pose_id * OriBlockSize + 1) = 1.0E+20;
                *(oris_std + pose_id * OriBlockSize + 2) = 1.0E+20;
                *(oris_std + pose_id * OriBlockSize + 3) = 1.0E+20;

                check_[pose_id]++;
            }
        }
        for (int i = 0; i < numEpoch; i++)
        {
            if (check_[i] != 1)
            {
                throw std::runtime_error("Prolem in filling pose");
            }
        }

        // 2. surface parameters
        for (int nPlane = 0; nPlane < numSurface; nPlane++)
        {
            Eigen::Vector4d params = surface_paramenters_vector[nPlane];
            // int fixIndex = planeParamTrans(params);

            *(planes + nPlane * 4 + 0) = params(0);
            *(planes + nPlane * 4 + 1) = params(1);
            *(planes + nPlane * 4 + 2) = params(2);
            *(planes + nPlane * 4 + 3) = params(3);

            *(planes_std + nPlane * 4 + 0) = 1.0E+20;
            *(planes_std + nPlane * 4 + 1) = 1.0E+20;
            *(planes_std + nPlane * 4 + 2) = 1.0E+20;
            *(planes_std + nPlane * 4 + 3) = 1.0E+20;

            // fixing the element equal to 1
            // *(planes_std + nPlane * 4 + fixIndex) = 1.0E-20;
        }

        //--------------------------build ceres  ------------------------------
        ceres::LossFunction *loss_function = new ceres::HuberLoss(0.1);
        ceres::LocalParameterization *q_parameterization = new ceres::EigenQuaternionParameterization(); // x, y, z, w
        ceres::Problem problem ;//= new ceres::Problem;
        ceres::CostFunction *costFunction;
        //--------------------------add parameter block ------------------------------
        for (int nEpoch = 0; nEpoch < numEpoch; nEpoch++)
        {
            problem.AddParameterBlock(poss + nEpoch * PosBlockSize, PosBlockSize);
        }
        for (int nEpoch = 0; nEpoch < numEpoch; nEpoch++)
        {
            problem.AddParameterBlock(oris + nEpoch * OriBlockSize, OriBlockSize, q_parameterization);
        }
        for (int nPlane = 0; nPlane < numSurface; nPlane++)
        {
            problem.AddParameterBlock(planes + nPlane * PlaneBlockSize, PlaneBlockSize);
        }

        // 1. surface observation
        int numUnknownGroundObs = 0;
        int numPlane = numSurface;
        for (int nPlane = 0; nPlane < numPlane; nPlane++)
        {
            int planeId = nPlane;
            double *plane = planes + planeId * PlaneBlockSize; // corresponding plane block
            uint64_t surface_key = surface_keys[nPlane];
            const auto & inner_points = voxel_key_to_inner_points[surface_key];
            for (int nP = 0; nP < inner_points.size(); nP++)
            {
                int epoch_id = inner_points[nP].epoch_id;
                Eigen::Vector3d r_I_lu=inner_points[nP].point;

                double *t_start = poss + epoch_id * PosBlockSize;
                double *q_start = oris + epoch_id * OriBlockSize;

                costFunction = PlaneFactorOneEpoch::Create(r_I_lu, surface_std);
                problem.AddResidualBlock(costFunction, loss_function, t_start, q_start, plane);
                numUnknownGroundObs++;
            }
        }
#if 1
		for (int nIscan = iscan_index_start_; nIscan <= iscan_index_end_; nIscan++)
		{
			int iscan_id = nIscan;
			int pose_start = iscan_pose_map[iscan_id];

			//for the first iscan, start from the second one.
			//int start_scan = (nIscan == iscan_index_start_) ? 1 : 0;
			for (int nScan = 0; nScan < mvpIScans[nIscan]->numScan; nScan++)
			{
				int pose_id_start = pose_start + nScan;
				int pose_id_end = pose_id_start + 1;

				double *t_1 = poss + pose_id_start * PosBlockSize;
				double *q_1 = oris + pose_id_start * OriBlockSize;
				double *t_2 = poss + pose_id_end * PosBlockSize;
				double *q_2 = oris + pose_id_end * OriBlockSize;

				// transformation estimated by odometry
				Eigen::Quaterniond q_prior(mvpIScans[iscan_id]->indScans[nScan].R_lut2_lut1);
				Eigen::Vector3d t_prior = mvpIScans[iscan_id]->indScans[nScan].r_lut2_lut1;

				// derive sqrt of information matrix
				// Sample: from CERES const Eigen::Matrix<double, 6, 6> sqrt_information =constraint.information.llt().matrixL();
				double acc_t;
				double acc_q;
				if (mPara.odo_from_trajectory_flag)
				{
					acc_t = mPara.odoPositionStdTraj;
					acc_q = deg2rad(mPara.odoOrientationStdTraj);
				}
				else
				{
					acc_t = mPara.odoPositionStd;
					acc_q = deg2rad(mPara.odoOrientationStd);
				}

				Eigen::Matrix<double, 6, 6> sqrt_information = Eigen::MatrixXd::Identity(6, 6);
				sqrt_information(0, 0) = 1.0 / acc_t;
				sqrt_information(1, 1) = 1.0 / acc_t;
				sqrt_information(2, 2) = 1.0 / acc_t;
				sqrt_information(3, 3) = 1.0 / acc_q;
				sqrt_information(4, 4) = 1.0 / acc_q;
				sqrt_information(5, 5) = 1.0 / acc_q;

				costFunction = PoseGraph3dErrorTerm::Create(t_prior, q_prior, sqrt_information);
				problem.AddResidualBlock(costFunction, nullptr, t_1, q_1, t_2, q_2);

				double dis_std = 0.05; // 0.5
				costFunction = DistanceBetweenTwoEpochsConstraints::Create(t_prior.norm(), dis_std);
				problem.AddResidualBlock(costFunction, nullptr, t_1, q_1, t_2, q_2);
			}
		}
#endif
        
        //--------------------------Fix Parameters  ------------------------------
        // set constant unknowns
        Checkforconstantparams(&problem, poss, poss_std, PosBlockSize, numEpoch);
        Checkforconstantparams(&problem, oris, oris_std, OriBlockSize, numEpoch);
        Checkforconstantparams(&problem, planes, planes_std, PlaneBlockSize, numPlane);

        TicToc t_solver;
        ceres::Solver::Options solverOptions;
        solverOptions.max_num_iterations = 10;
        solverOptions.minimizer_progress_to_stdout = true;
		// if (numFixTreeObs + numUnknownTreeObs > 1000000)
		{
			cout << "Using iterative schur" << endl;
			solverOptions.trust_region_strategy_type = ceres::TrustRegionStrategyType::LEVENBERG_MARQUARDT;
			solverOptions.preconditioner_type = ceres::PreconditionerType::SCHUR_JACOBI; //schur jacob  with iterative schur only
			solverOptions.linear_solver_type = ceres::LinearSolverType::ITERATIVE_SCHUR;//ITERATIVE_SCHUR; //SPARSE_NORMAL_CHOLESKY;
			solverOptions.max_num_consecutive_invalid_steps = 20;
			solverOptions.use_nonmonotonic_steps = true;
		}
		// else
		// {
		// 	cout << "Using SPARSE_NORMAL_CHOLESKY" << endl;
		// 	solverOptions.linear_solver_type = ceres::LinearSolverType::SPARSE_NORMAL_CHOLESKY; // ceres::LinearSolverType::SPARSE_NORMAL_CHOLESKY;
		// }
        ceres::Solver::Summary summary;
        ceres::Solve(solverOptions, &problem, &summary);
        fMapLog << ", Time: " << t_solver.toc() << " ms" << endl;
        fMapLog << "\t" << pose_count << " pose from " << iscan_index_end_ - iscan_index_start_ + 1 << " Iscans" << endl;

        cout << "Time: " << t_solver.toc() << endl;

        ///////////////////////////////////////////////////////////////////
        //******************** Update parameters **********************
        // 1. trajectory info
        for (int nIscan = iscan_index_start_; nIscan <= iscan_index_end_; nIscan++)
        {
            int iscan_id = nIscan;
            int pose_start = iscan_pose_map[iscan_id];

            // for the last Iscan, use all pose from v_R/r_mapping
            int num_pose_per_iscan = mvpIScans[nIscan]->numScan + 1;
            for (int nPose = 0; nPose < num_pose_per_iscan; nPose++)
            {
                int pose_id = pose_start + nPose;

                Eigen::Vector3d r{*(poss + pose_id * PosBlockSize + 0), *(poss + pose_id * PosBlockSize + 1), *(poss + pose_id * PosBlockSize + 2)};
                Eigen::Quaterniond q_refined(*(oris + pose_id * OriBlockSize + 3), *(oris + pose_id * OriBlockSize + 0),
                                            *(oris + pose_id * OriBlockSize + 1), *(oris + pose_id * OriBlockSize + 2));
                Eigen::Matrix3d tempR, orthoR;
                tempR = q_refined.toRotationMatrix();
                Eigen::Vector3d angles;
                angles = Find_Rotation(tempR);
                Compute_Rotation(angles(0), angles(1), angles(2), orthoR);
                mvpIScans[iscan_id]->v_r_mapping[nPose] = r;
                mvpIScans[iscan_id]->v_R_mapping[nPose] = orthoR;
            }
        }
    }

}

/**
 * @brief Writes a point cloud with optional metadata to a LAS file using PDAL.
 *
 * This function constructs a PDAL point table, fills in the point and attribute data
 * (if provided), and writes the result to a `.las` file.
 *
 * @param ofname_cloud Output LAS file path.
 * @param cloud Vector of 3D point positions.
 * @param cloud_info Optional vector of metadata (e.g., intensity, classification, tree ID).
 * @param x_min X offset for writing point cloud.
 * @param y_min Y offset for writing point cloud.
 * @param z_min Z offset for writing point cloud.
 */

void outputPointCloudByPDAL(std::string &ofname_cloud, const std::vector<Eigen::Vector3d> &cloud, const std::vector<CloudInfo> &cloud_info,
	double x_min, double y_min, double z_min)
{
	Options xjOptions;
	xjOptions.add("filename", ofname_cloud);
	xjOptions.add("extra_dims", "all");
	//pdal::Option las_opt_extra("use_eb_vlr", true);
	//xjOptions.add(las_opt_extra);
	xjOptions.add("offset_x", (double)x_min);
	xjOptions.add("offset_y", (double)y_min);
	xjOptions.add("offset_z", (double)z_min);
	xjOptions.add("scale_x", (double)0.00001);   //JD debug
	xjOptions.add("scale_y", (double)0.00001);   //JD debug
	xjOptions.add("scale_z", (double)0.00001);   //JD debug

	PointTable table;
	table.layout()->registerDim(Dimension::Id::X);
	table.layout()->registerDim(Dimension::Id::Y);
	table.layout()->registerDim(Dimension::Id::Z);
	table.layout()->registerDim(Dimension::Id::PointId);
	table.layout()->registerDim(Dimension::Id::PointSourceId);
	table.layout()->registerDim(Dimension::Id::OriginId);
	table.layout()->registerDim(Dimension::Id::Intensity);
	table.layout()->registerDim(Dimension::Id::GpsTime);
	table.layout()->registerDim(Dimension::Id::EchoRange);
	table.layout()->registerDim(Dimension::Id::Classification);
	table.layout()->registerDim(Dimension::Id::ScanDirectionFlag);

	table.layout()->registerDim(Dimension::Id::Red);
	table.layout()->registerDim(Dimension::Id::Green);
	table.layout()->registerDim(Dimension::Id::Blue);
	table.layout()->registerDim(Dimension::Id::ClusterID);
	table.layout()->registerDim(Dimension::Id::InternalTime);
	bool no_attibute = cloud_info.empty();
	PointViewPtr view(new PointView(table));
	int id = 0;
	for (int i = 0; i < cloud.size(); i++)
	{
		double x = cloud[i].x() + x_min;
		double y = cloud[i].y() + y_min;
		double z = cloud[i].z() + z_min;


		view->setField(pdal::Dimension::Id::X, id, x);
		view->setField(pdal::Dimension::Id::Y, id, y);
		view->setField(pdal::Dimension::Id::Z, id, z);
		if (no_attibute)
		{
			id++;
			continue;
		}
		view->setField(pdal::Dimension::Id::Intensity, id, cloud_info[i].intensity);
		view->setField(pdal::Dimension::Id::PointId, id, cloud_info[i].point_id);
		view->setField(pdal::Dimension::Id::OriginId, id, cloud_info[i].origin_id);
		view->setField(pdal::Dimension::Id::GpsTime, id, cloud_info[i].gps_time);
		view->setField(pdal::Dimension::Id::PointSourceId, id, cloud_info[i].point_source_id);
		view->setField(pdal::Dimension::Id::Intensity, id, cloud_info[i].intensity);
		view->setField(pdal::Dimension::Id::GpsTime, id, cloud_info[i].gps_time);
		view->setField(pdal::Dimension::Id::EchoRange, id, cloud_info[i].echo_range);
		view->setField(pdal::Dimension::Id::Classification, id, cloud_info[i].classification);
		view->setField(pdal::Dimension::Id::ScanDirectionFlag, id, cloud_info[i].scan_direction_flag);

		view->setField(pdal::Dimension::Id::Red, id, cloud_info[i].red);
		view->setField(pdal::Dimension::Id::Green, id, cloud_info[i].green);
		view->setField(pdal::Dimension::Id::Blue, id, cloud_info[i].blue);

		view->setField(pdal::Dimension::Id::ClusterID, id, cloud_info[i].tree_ID);
		view->setField(pdal::Dimension::Id::InternalTime, id, cloud_info[i].time_ratio);
		id++;

	}

	BufferReader xjBufferReader;
	xjBufferReader.addView(view);
	//xjBufferReader.addView(view);

	StageFactory factory;
	Stage *writer = factory.createStage("writers.las");
	writer->setInput(xjBufferReader);
	writer->setOptions(xjOptions);
	writer->prepare(table);
	writer->execute(table);
}

struct VectorSizeComparator {
    bool operator()(const std::pair<int, std::vector<Eigen::Vector3d>>& vec1, const std::pair<int, std::vector<Eigen::Vector3d>>& vec2) const {
        return vec1.second.size() > vec2.second.size();
    }
};

/**
 * @brief Calculates the most dominant cylindrical portion (cluster) of each tree based on scan contributions.
 *
 * This function clusters tree points by scan ID and identifies the dominant cluster for each tree. Then,
 * it fits a cylinder to that cluster to estimate radius and fitting error.
 * The points are also exported as a LAS file for visualization or debugging.
 *
 * @return A vector of RadiusInfo structs, each containing radius, error, and success flag for each tree.
 */

std::vector<RadiusInfo> Mapping::CalculateMostPortionRadius()
{
	tree_portion_points_.clear();
	tree_portion_points_info_.clear();
    std::vector<RadiusInfo> radius_info_results(used_trees_.size());
    for (int nCylinder = 0; nCylinder < used_trees_.size(); nCylinder++)
    {
        int map_tree_id = used_trees_[nCylinder];


        std::map<int,std::vector<Eigen::Vector3d>> scan_id_to_points;
        std::map<int,int> scan_id_to_cluster_id;
        int tree_points_number = 0;
        int initial_cluster_id = -1;
        for (int nT = 0; nT < mpMapTree[map_tree_id]->visibleInfo.size(); nT++)
        {
            int iscanId = mpMapTree[map_tree_id]->visibleInfo[nT].first;
            int treeId = mpMapTree[map_tree_id]->visibleInfo[nT].second;

            int num_p = mvpIScans[iscanId]->vTreePointMapping[treeId].size();
            for (int nP = 0; nP < mvpIScans[iscanId]->vTreePointMapping[treeId].size(); nP ++)
            {
                // points in mapping
                Eigen::Vector3d r_I_m(mvpIScans[iscanId]->vTreePointMapping[treeId][nP].x,
                                    mvpIScans[iscanId]->vTreePointMapping[treeId][nP].y,
                                    mvpIScans[iscanId]->vTreePointMapping[treeId][nP].z);
                int local_scanId = int(mvpIScans[iscanId]->vTreePointRaw[treeId][nP].intensity / 100);
				int scanId = mvpIScans[iscanId]->indScans[local_scanId].scanID;
                scan_id_to_points[scanId].push_back(r_I_m);
                scan_id_to_cluster_id[scanId]=initial_cluster_id;
                tree_points_number++;
            }
        }
        int cluster_id = 0;
        std::map<int,int> cluster_id_to_points_number;
        int interval_threshold = 25;
        std::queue<int> seed_scan_ids;
		std::vector<std::pair<int, std::vector<Eigen::Vector3d>>> sorted_scan_id_to_points(scan_id_to_points.begin(), scan_id_to_points.end());
		std::sort(sorted_scan_id_to_points.begin(), sorted_scan_id_to_points.end(), VectorSizeComparator());
        for(const auto& item: sorted_scan_id_to_points)
        {
            if(scan_id_to_cluster_id[item.first]==-1)
            {
				if(item.second.size()<5) continue;
                cluster_id_to_points_number[cluster_id] = item.second.size();
                scan_id_to_cluster_id[item.first]=cluster_id;
                seed_scan_ids.push(item.first);
                do
                {
                    int current_scan_id =seed_scan_ids.front();
		            seed_scan_ids.pop();
                    for(int id = current_scan_id - interval_threshold; id<=current_scan_id + interval_threshold;id++)
                    {
                        if(scan_id_to_points.find(id)==scan_id_to_points.end() || scan_id_to_cluster_id[id]!=-1)
                        {
                            continue; 
                        }
						if(scan_id_to_points[id].size()<5) continue;
                        cluster_id_to_points_number[cluster_id] += scan_id_to_points[id].size();
                        scan_id_to_cluster_id[id]=cluster_id;
                        seed_scan_ids.push(id);
                    }
                
                } while (!seed_scan_ids.empty());
                if(cluster_id_to_points_number[cluster_id]>(tree_points_number*0.5))
                {
                    break;
                }
                cluster_id++;
            }
            // else
            // {
            //     continue;
            // }
        }
        int max_points_number = INT_MIN;
        int final_cluster_id;
        for(const auto & item:cluster_id_to_points_number) 
        {
            if(item.second>max_points_number)
            {
                max_points_number = item.second;
                final_cluster_id = item.first;
            }
        }

        std::vector<Eigen::Vector3d> points_to_calculate_radius;
		CloudInfo p_info;
        for(const auto& item: scan_id_to_cluster_id)
        {
            if(item.second == final_cluster_id)
            {
                points_to_calculate_radius.insert(points_to_calculate_radius.end(), scan_id_to_points[item.first].begin(),
                  scan_id_to_points[item.first].end());

                for (const auto& p : scan_id_to_points[item.first])
                {
                    tree_portion_points_.push_back(p);
                    p_info.origin_id = item.first;
                    tree_portion_points_info_.push_back(p_info);
                }
            }

        }
        std::vector<double> cylinder_parameters;
        double final_error;
        RadiusInfo radius_info;
		std::vector<PointType> points_convert;
		PointType point_t;
		Eigen::Vector3d offset = points_to_calculate_radius.front();
		for (const auto& p : points_to_calculate_radius)
		{
			point_t.getVector3fMap() = (p - offset).cast<float>();
			points_convert.push_back(point_t);
		}
		//tree_portion_points_.insert(tree_portion_points_.end(), points_to_calculate_radius.begin(), points_to_calculate_radius.end());
        if(Forestry_SLAM::FittingCylinderModel(points_convert,cylinder_parameters,final_error)==1)
        {
            radius_info.success =true;
            radius_info.radius = cylinder_parameters[6];
            radius_info.error = final_error;
        }
        radius_info_results[nCylinder]=radius_info;
    }
	clock_t now_time = clock();
	std::string file_name_test = "D:\\DATA_extend\\CPT_VLP_16HR_20230302\\test_tree_points_"+to_string(now_time)+".las";
	outputPointCloudByPDAL(file_name_test, tree_portion_points_, tree_portion_points_info_,0,0,0);
    return radius_info_results;
}
//pose for each scan is used for estimation
/**
 * @brief Performs loop closure optimization using raw scan-level data.
 *
 * This function builds a Ceres optimization problem with tree-cylinder and ground-plane constraints,
 * as well as odometry constraints between consecutive scans. It updates trajectory poses and tree parameters
 * after optimization, supporting fine-grained alignment in raw mapping mode.
 *
 * Optimization uses downsampling for computational efficiency, and optionally incorporates radius
 * refinements from `CalculateMostPortionRadius()`.
 */

void Mapping::OptimizeRawLevelLC()
{
    double tree_std = mPara.mapTreeStd;
    double ground_std = mPara.mapGroundStd;
    double ref_ground_std = ref_ground_std_;
    double global_tree_loc_std = mPara.global_tree_loc_std;
    int max_num_map_point_per_tree = Max_Point_Per_Tree_Raw_LC;

    int iter = 0;
    int maxIter = 1;
    while (iter < maxIter)
    {
        fMapLog << "Iter " << ++iter;
        fMapLog << "\tNumber of included trees " << used_trees_.size();

        //map saving the start pose index of a Iscan
        unordered_map<int, int> iscan_pose_map;
        int pose_count = 0;
        for (int nIscan = iscan_index_start_; nIscan <=  iscan_index_end_; nIscan++)
        {
            iscan_pose_map[nIscan] = pose_count;
            pose_count += mvpIScans[nIscan]->numScan;
        }
        pose_count ++;
        
        // initialization
        int numEpoch = pose_count;
        int numCylinder = used_trees_.size();
        int numPlane = plane_params_.size();

        // unknowns
        int numPosUnknown = numEpoch * PosBlockSize;
        int numOriUnknown = numEpoch * OriBlockSize;
        int numCylinderUnknown = numCylinder * CylinderBlockSize;
        int numPlaneUnknown = numPlane * PlaneBlockSize;
        int numUnknown = numPosUnknown + numOriUnknown + numCylinderUnknown + numPlaneUnknown; // all unknowns
        double *poss = new double[numPosUnknown];
        double *poss_std = new double[numPosUnknown];
        double *oris = new double[numOriUnknown];
        double *oris_std = new double[numOriUnknown];
        double *cylinders = new double[numCylinderUnknown];
        double *cylinders_std = new double[numCylinderUnknown];
        double *planes = new double[numPlaneUnknown];
        double *planes_std = new double[numPlaneUnknown];

        // fill parameter
        // 1. trajectory information
        vector<int> check_(numEpoch, 0);
        for (int nIscan = iscan_index_start_; nIscan <= iscan_index_end_; nIscan++)
        {
            int iscan_id = nIscan;
            int pose_start = iscan_pose_map[iscan_id];

            // for the last Iscan, use all pose from v_R/r_mapping
            int num_pose_per_iscan = (nIscan == iscan_index_end_) ? mvpIScans[nIscan]->numScan + 1 : mvpIScans[nIscan]->numScan;
            for (int nPose = 0; nPose < num_pose_per_iscan; nPose++)
            {
                int pose_id = pose_start + nPose;
                *(poss + pose_id * PosBlockSize + 0) = mvpIScans[iscan_id]->v_r_mapping[nPose](0);
                *(poss + pose_id * PosBlockSize + 1) = mvpIScans[iscan_id]->v_r_mapping[nPose](1);
                *(poss + pose_id * PosBlockSize + 2) = mvpIScans[iscan_id]->v_r_mapping[nPose](2);

                Eigen::Quaterniond q(mvpIScans[iscan_id]->v_R_mapping[nPose]);
                *(oris + pose_id * OriBlockSize + 0) = q.x();
                *(oris + pose_id * OriBlockSize + 1) = q.y();
                *(oris + pose_id * OriBlockSize + 2) = q.z();
                *(oris + pose_id * OriBlockSize + 3) = q.w();

                *(poss_std + pose_id * PosBlockSize + 0) = 1.0E+20;
                *(poss_std + pose_id * PosBlockSize + 1) = 1.0E+20;
                *(poss_std + pose_id * PosBlockSize + 2) = 1.0E+20;
                *(oris_std + pose_id * OriBlockSize + 0) = 1.0E+20;
                *(oris_std + pose_id * OriBlockSize + 1) = 1.0E+20;
                *(oris_std + pose_id * OriBlockSize + 2) = 1.0E+20;
                *(oris_std + pose_id * OriBlockSize + 3) = 1.0E+20;

                check_[pose_id]++;
            }
        }
        for (int i = 0; i < numEpoch; i++)
        {
            if (check_[i] != 1)
            {
                throw std::runtime_error("Prolem in filling pose");
            }
        }
		std::vector<RadiusInfo> radius_info_results;
		if (radius_stategy_)
		{
			radius_info_results = CalculateMostPortionRadius();
		}
        // 2. cylinder parameter
        int countCylinder = 0;
        int updated_radius_count = 0;
        for (int nTree = 0; nTree < numCylinder; nTree++)
        {
            int map_tree_id = used_trees_[nTree];
            *(cylinders + nTree * CylinderBlockSize + 0) = mpMapTree[map_tree_id]->para.x(0); // XYZ (Z is nominal)
            *(cylinders + nTree * CylinderBlockSize + 1) = mpMapTree[map_tree_id]->para.x(1);
            *(cylinders + nTree * CylinderBlockSize + 2) = mpMapTree[map_tree_id]->para.x(2);
            *(cylinders + nTree * CylinderBlockSize + 3) = mpMapTree[map_tree_id]->para.n(0);
            *(cylinders + nTree * CylinderBlockSize + 4) = mpMapTree[map_tree_id]->para.n(1);
            *(cylinders + nTree * CylinderBlockSize + 5) = mpMapTree[map_tree_id]->para.n(2);
            *(cylinders + nTree * CylinderBlockSize + 6) = mpMapTree[map_tree_id]->para.r;

            *(cylinders_std + nTree * CylinderBlockSize + 0) = 1.0E+20; // XYZ
            *(cylinders_std + nTree * CylinderBlockSize + 1) = 1.0E+20;
            *(cylinders_std + nTree * CylinderBlockSize + 2) = 1.0E-20;
            *(cylinders_std + nTree * CylinderBlockSize + 3) = 1.0E+20; // ux, uy, uz
            *(cylinders_std + nTree * CylinderBlockSize + 4) = 1.0E+20;
            *(cylinders_std + nTree * CylinderBlockSize + 5) = 1.0E-20;
            *(cylinders_std + nTree * CylinderBlockSize + 6) = 1.0E+20; // r

			if (radius_stategy_)
			{
				if (radius_info_results[nTree].success && radius_info_results[nTree].error < 0.3&& radius_info_results[nTree].radius < 1)
				{
					*(cylinders + nTree * CylinderBlockSize + 6) = radius_info_results[nTree].radius;
					fMapLog << "updated tree radius, error: " << radius_info_results[nTree].radius << "," << radius_info_results[nTree].error << std::endl;
					*(cylinders_std + nTree * CylinderBlockSize + 6) = 1.0E-20; // r
					updated_radius_count++;
				}
			}
        }
		radius_stategy_ = false;//only apply in second iteration.
        std::cout<<"MODIFICATION: tree number/updated_radius_count: "<<numCylinder<<"/"<<updated_radius_count<<std::endl;

        // 3. plane parameters
        for (int nPlane = 0; nPlane < numPlane; nPlane++)
        {
            Eigen::Vector4d params = plane_params_[nPlane];
            int fixIndex = planeParamTrans(params);

            *(planes + nPlane * 4 + 0) = params(0);
            *(planes + nPlane * 4 + 1) = params(1);
            *(planes + nPlane * 4 + 2) = params(2);
            *(planes + nPlane * 4 + 3) = params(3);

            *(planes_std + nPlane * 4 + 0) = 1.0E+20;
            *(planes_std + nPlane * 4 + 1) = 1.0E+20;
            *(planes_std + nPlane * 4 + 2) = 1.0E+20;
            *(planes_std + nPlane * 4 + 3) = 1.0E+20;

            // fixing the element equal to 1
            *(planes_std + nPlane * 4 + fixIndex) = 1.0E-20;
        }

        //------------- check fixed parameters
        int numFixedUnknown = 0;
        for (int i = 0; i < numPosUnknown; i++)
        {
            double temp;
            temp = *(poss_std + i);
            if (temp < SigmaFixed)
                numFixedUnknown += 1;
        }
        for (int i = 0; i < numOriUnknown; i++)
        {
            double temp;
            temp = *(oris_std + i);
            if (temp < SigmaFixed)
                numFixedUnknown += 1;
        }
        for (int i = 0; i < numCylinderUnknown; i++)
        {
            double temp;
            temp = *(cylinders_std + i);
            if (temp < SigmaFixed)
                numFixedUnknown += 1;
        }
        for (int i = 0; i < numPlaneUnknown; i++)
        {
            double temp;
            temp = *(planes_std + i);
            if (temp < SigmaFixed)
                numFixedUnknown += 1;
        }

        //--------------------------build ceres  ------------------------------
        ceres::LossFunction *loss_function = new ceres::HuberLoss(0.1);
        ceres::LocalParameterization *q_parameterization = new ceres::EigenQuaternionParameterization(); // x, y, z, w
        ceres::Problem problem ;//= new ceres::Problem;
        ceres::CostFunction *costFunction;
        //--------------------------add parameter block ------------------------------
        for (int nEpoch = 0; nEpoch < numEpoch; nEpoch++)
        {
            problem.AddParameterBlock(poss + nEpoch * PosBlockSize, PosBlockSize);
        }
        for (int nEpoch = 0; nEpoch < numEpoch; nEpoch++)
        {
            problem.AddParameterBlock(oris + nEpoch * OriBlockSize, OriBlockSize, q_parameterization);
        }
        for (int nCylinder = 0; nCylinder < numCylinder; nCylinder++)
        {
            problem.AddParameterBlock(cylinders + nCylinder * CylinderBlockSize, CylinderBlockSize);
        }
        for (int nPlane = 0; nPlane < numPlane; nPlane++)
        {
            problem.AddParameterBlock(planes + nPlane * PlaneBlockSize, PlaneBlockSize);
        }

        //-------------------------add observation -------------------------------
        // 1. Tree observation
        int numFixTreeObs = 0;
        int numUnknownTreeObs = 0;
        std::set<uint64_t> downsample_index_set_tree;
        int original_tree_points_number = 0;
        int downsample_points_number = 0;
        for (int nCylinder = 0; nCylinder < numCylinder; nCylinder++)
        {
            int cylinderId = nCylinder;
            int map_tree_id = used_trees_[nCylinder];
            double *cylinder = cylinders + cylinderId * CylinderBlockSize; // corresponding cylinder block
            
            int increment = 1;
            if(mpMapTree[map_tree_id]->numPoint != 0)
            {
                increment = ceil(mpMapTree[map_tree_id]->numPoint/max_num_map_point_per_tree);
                increment = max(1, increment);
            }
            for (int nT = 0; nT < mpMapTree[map_tree_id]->visibleInfo.size(); nT++)
            {
                int iscanId = mpMapTree[map_tree_id]->visibleInfo[nT].first;
                int treeId = mpMapTree[map_tree_id]->visibleInfo[nT].second;

                int num_p = mvpIScans[iscanId]->vTreePointMapping[treeId].size();
                int incre_individual = num_p < Min_Tree_Points_Per_Iscan_LC? 1 : increment;
                for (int nP = 0; nP < mvpIScans[iscanId]->vTreePointMapping[treeId].size(); nP = nP + incre_individual)
                {
                    original_tree_points_number++;
                    //if from the iscan that is fix
                    if(iscanId < iscan_index_start_)
                    {
                        // points in mapping
                        Eigen::Vector3d r_I_m(mvpIScans[iscanId]->vTreePointMapping[treeId][nP].x,
                                              mvpIScans[iscanId]->vTreePointMapping[treeId][nP].y,
                                              mvpIScans[iscanId]->vTreePointMapping[treeId][nP].z);
                        uint64_t key=GetVoxelkeyIndex(r_I_m,KDownsampleSizeRawLevelTree);
                        if(downsample_index_set_tree.find(key)==downsample_index_set_tree.end())
                        {
                            downsample_index_set_tree.insert(key);
                        }
                        else
                        {
                            continue;
                        }
                        downsample_points_number++;
                        //apply the scale to std
                        costFunction = CylinderFactorObjectPoint::Create(r_I_m, tree_std / sqrt(increment));
                        //costFunction = CylinderFactorObjectPoint::Create(r_I_m, tree_std);
                        problem.AddResidualBlock(costFunction, loss_function, cylinder);
                        numFixTreeObs++;
                    }
                    else
                    {
                        int scanId = int(mvpIScans[iscanId]->vTreePointRaw[treeId][nP].intensity / 100);
                        double t_ratio = double(mvpIScans[iscanId]->vTreePointRaw[treeId][nP].intensity -
                                                int(mvpIScans[iscanId]->vTreePointRaw[treeId][nP].intensity));
                        Eigen::Vector3d r_I_lu(mvpIScans[iscanId]->vTreePointRaw[treeId][nP].x,
                                               mvpIScans[iscanId]->vTreePointRaw[treeId][nP].y,
                                               mvpIScans[iscanId]->vTreePointRaw[treeId][nP].z);

                        Eigen::Vector3d point_in_mapping = mvpIScans[iscanId]->v_R_mapping[scanId] * r_I_lu 
                            + mvpIScans[iscanId]->v_r_mapping[scanId];
                        uint64_t key=GetVoxelkeyIndex(point_in_mapping,KDownsampleSizeRawLevelTree);
                        if(downsample_index_set_tree.find(key)==downsample_index_set_tree.end())
                        {
                            downsample_index_set_tree.insert(key);
                        }
                        else
                        {
                            continue;
                        }
                        downsample_points_number++;
                        int epoch_id = iscan_pose_map[iscanId] + scanId;
                        double *t_start = poss + epoch_id * PosBlockSize;
                        double *q_start = oris + epoch_id * OriBlockSize;
                        double *t_end = poss + (epoch_id + 1) * PosBlockSize;
                        double *q_end = oris + (epoch_id + 1) * OriBlockSize;

                        // costFunction = CylinderFactorTwoEpoch::Create(r_I_lu, t_ratio, tree_std / sqrt(increment));                        
                        // //costFunction = CylinderFactorTwoEpoch::Create(r_I_lu, t_ratio, tree_std);
                        // problem.AddResidualBlock(costFunction, loss_function, t_start, q_start, t_end, q_end, cylinder);

                        costFunction = CylinderFactorOneEpoch::Create(r_I_lu, tree_std / sqrt(increment));
                        problem.AddResidualBlock(costFunction, loss_function, t_start, q_start, cylinder);
                        numUnknownTreeObs++;
                    }
                }
            }
        }
        std::cout<<"RawLC original_tree_points_number: "<<original_tree_points_number<<std::endl;
        std::cout<<"RawLC downsample_tree_points_number: "<<downsample_points_number<<std::endl;

        // 2. ground observation
        int numRefGroundObs = 0;
        int numFixGroundObs = 0;
        int numUnknownGroundObs = 0;
        std::set<uint64_t> downsample_index_set_ground;
        for (int nPlane = 0; nPlane < numPlane; nPlane++)
        {
            int planeId = nPlane;
            double *plane = planes + planeId * PlaneBlockSize; // corresponding plane block
            for (int nP = 0; nP < point_index_per_grid_[nPlane].size(); nP++)
            {
                int iscanId = point_index_per_grid_[nPlane][nP].first;
                int ptId = point_index_per_grid_[nPlane][nP].second;
                int epoch_id = iscanId - iscan_index_start_;

                // if from the iscan that is fix
                if (iscanId < iscan_index_start_)
                {
                    // points in mapping
                    Eigen::Vector3d r_I_m(mvpIScans[iscanId]->vGroundPointMapping[ptId].x,
                                          mvpIScans[iscanId]->vGroundPointMapping[ptId].y,
                                          mvpIScans[iscanId]->vGroundPointMapping[ptId].z);

                    uint64_t key=GetVoxelkeyIndex(r_I_m,KDownsampleSizeRawLevelGround);
                    if(downsample_index_set_ground.find(key)==downsample_index_set_ground.end())
                    {
                        downsample_index_set_ground.insert(key);
                    }
                    else
                    {
                        continue;
                    }
                    downsample_points_number++;

                    costFunction = PlaneFactorObjectPoint::Create(r_I_m, ground_std);
                    problem.AddResidualBlock(costFunction, loss_function, plane);
                    numFixGroundObs++;
                }
                else
                {
                    int scanId = int(mvpIScans[iscanId]->vGroundPointRaw[ptId].intensity / 100);
                    double t_ratio = double(mvpIScans[iscanId]->vGroundPointRaw[ptId].intensity -
                                            int(mvpIScans[iscanId]->vGroundPointRaw[ptId].intensity));
                    Eigen::Vector3d r_I_lu(mvpIScans[iscanId]->vGroundPointRaw[ptId].x,
                                           mvpIScans[iscanId]->vGroundPointRaw[ptId].y,
                                           mvpIScans[iscanId]->vGroundPointRaw[ptId].z);

                    Eigen::Vector3d point_in_mapping = mvpIScans[iscanId]->v_R_mapping[scanId] * r_I_lu 
                        + mvpIScans[iscanId]->v_r_mapping[scanId];
                    uint64_t key=GetVoxelkeyIndex(point_in_mapping,KDownsampleSizeRawLevelGround);
                    if(downsample_index_set_ground.find(key)==downsample_index_set_ground.end())
                    {
                        downsample_index_set_ground.insert(key);
                    }
                    else
                    {
                        continue;
                    }
                    downsample_points_number++;

                    int epoch_id = iscan_pose_map[iscanId] + scanId;
                    double *t_start = poss + epoch_id * PosBlockSize;
                    double *q_start = oris + epoch_id * OriBlockSize;
                    double *t_end = poss + (epoch_id + 1) * PosBlockSize;
                    double *q_end = oris + (epoch_id + 1) * OriBlockSize;

                    // costFunction = PlaneFactorTwoEpoch::Create(r_I_lu, t_ratio, ground_std);
                    // problem.AddResidualBlock(costFunction, loss_function, t_start, q_start, t_end, q_end, plane);

                    costFunction = PlaneFactorOneEpoch::Create(r_I_lu, ground_std);
                    problem.AddResidualBlock(costFunction, loss_function, t_start, q_start, plane);
                    numUnknownGroundObs++;
                }
            }

            for (int nP = 0; nP < point_ref_per_grid_[nPlane].size(); nP++)
            {

                costFunction = PlaneFactorObjectPoint::Create(point_ref_per_grid_[nPlane][nP], ref_ground_std);
                problem.AddResidualBlock(costFunction, loss_function, plane);
                numRefGroundObs++;
            }
        }
      


        // 3. Global tree
        int num_global_tree = 0;
        for (int nCylinder = 0; nCylinder < numCylinder; nCylinder++)
        {
            int cylinderId = nCylinder;
            int map_tree_id = used_trees_[nCylinder];
            double *cylinder = cylinders + cylinderId * CylinderBlockSize; // corresponding cylinder block
           
            if (mpMapTree[map_tree_id]->map_tree_flag)
            {
                double *cylinder = cylinders + cylinderId * CylinderBlockSize;
                //costFunction = CylinderPriorConstraints::Create(mpMapTree[mapTreeId]->para_ref_.x, global_tree_loc_std);
                costFunction = CylinderPriorConstraintsPointToLine::Create(mpMapTree[map_tree_id]->para_ref_.x, global_tree_loc_std);
                problem.AddResidualBlock(costFunction, nullptr, cylinder);                
                num_global_tree++;
            }
        }

        // 4. odometry constraints
        //based on test, if odometry is not from traj
        int odometry_constraint_count = 0;
        if (mPara.odo_from_trajectory_flag)
        {
            for (int nIscan = iscan_index_start_; nIscan <= iscan_index_end_; nIscan++)
            {
                int iscan_id = nIscan;
                int pose_start = iscan_pose_map[iscan_id];

                //for the first iscan, start from the second one.
                //int start_scan = (nIscan == iscan_index_start_) ? 1 : 0;
                for (int nScan = 0; nScan < mvpIScans[nIscan]->numScan; nScan++)
                {
                    int pose_id_start = pose_start + nScan;
                    int pose_id_end = pose_id_start + 1 ;

                    double *t_1 = poss + pose_id_start * PosBlockSize;
                    double *q_1 = oris + pose_id_start * OriBlockSize;
                    double *t_2 = poss + pose_id_end * PosBlockSize;
                    double *q_2 = oris + pose_id_end * OriBlockSize;

                    // transformation estimated by odometry
                    Eigen::Quaterniond q_prior(mvpIScans[iscan_id]->indScans[nScan].R_lut2_lut1);
                    Eigen::Vector3d t_prior = mvpIScans[iscan_id]->indScans[nScan].r_lut2_lut1;

                    // derive sqrt of information matrix
                    // Sample: from CERES const Eigen::Matrix<double, 6, 6> sqrt_information =constraint.information.llt().matrixL();
                    double acc_t; 
                    double acc_q; 
                    if(mPara.odo_from_trajectory_flag)
                    {
                        acc_t = mPara.odoPositionStdTraj;
                        acc_q = deg2rad(mPara.odoOrientationStdTraj); 
                    }
                    else
                    {
                        acc_t = mPara.odoPositionStd;
                        acc_q = deg2rad(mPara.odoOrientationStd); 
                    }

                    Eigen::Matrix<double, 6, 6> sqrt_information = Eigen::MatrixXd::Identity(6, 6);
                    sqrt_information(0, 0) = 1.0 / acc_t;
                    sqrt_information(1, 1) = 1.0 / acc_t;
                    sqrt_information(2, 2) = 1.0 / acc_t;
                    sqrt_information(3, 3) = 1.0 / acc_q;
                    sqrt_information(4, 4) = 1.0 / acc_q;
                    sqrt_information(5, 5) = 1.0 / acc_q;

                    costFunction = PoseGraph3dErrorTerm::Create(t_prior, q_prior, sqrt_information);
                    problem.AddResidualBlock(costFunction, nullptr, t_1, q_1, t_2, q_2);

                double dis_std = 0.5 ;
				costFunction = DistanceBetweenTwoEpochsConstraints::Create(t_prior.norm(),dis_std);
                problem.AddResidualBlock(costFunction, nullptr, t_1, q_1, t_2, q_2);
                    
                    odometry_constraint_count++;
                }
            }

        }

        //--------------------------Fix Parameters  ------------------------------
        // set constant unknowns
        Checkforconstantparams(&problem, poss, poss_std, PosBlockSize, numEpoch);
        Checkforconstantparams(&problem, oris, oris_std, OriBlockSize, numEpoch);
        Checkforconstantparams(&problem, cylinders, cylinders_std, CylinderBlockSize, numCylinder);
        Checkforconstantparams(&problem, planes, planes_std, PlaneBlockSize, numPlane);

        TicToc t_solver;
        ceres::Solver::Options solverOptions;
        solverOptions.max_num_iterations = 12;
        solverOptions.minimizer_progress_to_stdout = true;
		if (numFixTreeObs + numUnknownTreeObs > 1000000)
		{
			cout << "Using iterative schur" << endl;
			solverOptions.trust_region_strategy_type = ceres::TrustRegionStrategyType::LEVENBERG_MARQUARDT;
			solverOptions.preconditioner_type = ceres::PreconditionerType::SCHUR_JACOBI; //schur jacob  with iterative schur only
			solverOptions.linear_solver_type = ceres::LinearSolverType::ITERATIVE_SCHUR;//ITERATIVE_SCHUR; //SPARSE_NORMAL_CHOLESKY;
			solverOptions.max_num_consecutive_invalid_steps = 20;
			solverOptions.use_nonmonotonic_steps = true;
		}
		else
		{
			cout << "Using SPARSE_NORMAL_CHOLESKY" << endl;
			solverOptions.linear_solver_type = ceres::LinearSolverType::SPARSE_NORMAL_CHOLESKY; // ceres::LinearSolverType::SPARSE_NORMAL_CHOLESKY;
		}
        ceres::Solver::Summary summary;
        ceres::Solve(solverOptions, &problem, &summary);
        fMapLog << ", Time: " << t_solver.toc() << " ms" << endl;
        fMapLog << "\t" << pose_count << " pose from " << iscan_index_end_ - iscan_index_start_ + 1 << " Iscans" << endl;

        fMapLog << "\tFix/unknown tree observations: " << numFixTreeObs << " " << numUnknownTreeObs << endl;
        fMapLog << "\tRef/fix/unknown ground observations: " << numRefGroundObs << " " << numFixGroundObs << " " << numUnknownGroundObs << endl;
        fMapLog << "\tNumber of pair with global tree: " << num_global_tree << endl;
        if (mPara.odo_from_trajectory_flag)
            fMapLog << "\t" << odometry_constraint_count << " constraints from odometry" << endl;

        cout << "Time: " << t_solver.toc() << endl;

        ///////////////////////////////////////////////////////////////////


        //******************** Update parameters **********************
        // 1. trajectory info
        for (int nIscan = iscan_index_start_; nIscan <= iscan_index_end_; nIscan++)
        {
            int iscan_id = nIscan;
            int pose_start = iscan_pose_map[iscan_id];

            // for the last Iscan, use all pose from v_R/r_mapping
            int num_pose_per_iscan = mvpIScans[nIscan]->numScan + 1;
            for (int nPose = 0; nPose < num_pose_per_iscan; nPose++)
            {
                int pose_id = pose_start + nPose;

                Eigen::Vector3d r{*(poss + pose_id * PosBlockSize + 0), *(poss + pose_id * PosBlockSize + 1), *(poss + pose_id * PosBlockSize + 2)};
                Eigen::Quaterniond q_refined(*(oris + pose_id * OriBlockSize + 3), *(oris + pose_id * OriBlockSize + 0),
                                            *(oris + pose_id * OriBlockSize + 1), *(oris + pose_id * OriBlockSize + 2));
                Eigen::Matrix3d tempR, orthoR;
                tempR = q_refined.toRotationMatrix();
                Eigen::Vector3d angles;
                angles = Find_Rotation(tempR);
                Compute_Rotation(angles(0), angles(1), angles(2), orthoR);
                mvpIScans[iscan_id]->v_r_mapping[nPose] = r;
                mvpIScans[iscan_id]->v_R_mapping[nPose] = orthoR;
            }
        }

        //******************** Residuals  **********************
        // 1. Tree observation
        double res_tree_map = 0.0, res_tree_unknown = 0.0, res_tree_map_ini = 0.0, res_tree_unknown_ini = 0.0;
        CylinderPara treePara;
        int num_valid_tree_para = 0, num_invalid_tree_para = 0;
        for (int nCylinder = 0; nCylinder < numCylinder; nCylinder++)
        {
            int cylinderId = nCylinder;
            int map_tree_id = used_trees_[nCylinder];
            double *cylinder = cylinders + cylinderId * CylinderBlockSize; // corresponding cylinder block
            treePara.x(0) = *(cylinders + cylinderId * CylinderBlockSize + 0);
            treePara.x(1) = *(cylinders + cylinderId * CylinderBlockSize + 1);
            treePara.x(2) = *(cylinders + cylinderId * CylinderBlockSize + 2);
            treePara.n(0) = *(cylinders + cylinderId * CylinderBlockSize + 3);
            treePara.n(1) = *(cylinders + cylinderId * CylinderBlockSize + 4);
            treePara.n(2) = *(cylinders + cylinderId * CylinderBlockSize + 5);
            treePara.r = *(cylinders + cylinderId * CylinderBlockSize + 6);

            //compute residual
            int increment = 1;
            if(mpMapTree[map_tree_id]->numPoint != 0)
            {
                increment = ceil(mpMapTree[map_tree_id]->numPoint/max_num_map_point_per_tree);
                increment = max(1, increment);
            }

            for (int nT = 0; nT < mpMapTree[map_tree_id]->visibleInfo.size(); nT++)
            {
                int iscan_id = mpMapTree[map_tree_id]->visibleInfo[nT].first;
                int treeId = mpMapTree[map_tree_id]->visibleInfo[nT].second;

                int num_p = mvpIScans[iscan_id]->vTreePointMapping[treeId].size();
                int incre_individual = num_p < Min_Tree_Points_Per_Iscan_LC? 1 : increment;
                for (int nP = 0; nP < mvpIScans[iscan_id]->vTreePointMapping[treeId].size(); nP = nP + incre_individual)
                {
                    //if from the iscan that is fix
                    if(iscan_id < iscan_index_start_)
                    {
                        // points in mapping
                        Eigen::Vector3d r_I_m(mvpIScans[iscan_id]->vTreePointMapping[treeId][nP].x,
                                                mvpIScans[iscan_id]->vTreePointMapping[treeId][nP].y,
                                                mvpIScans[iscan_id]->vTreePointMapping[treeId][nP].z);
                        
                        //apply the scale to std
                        double d_ini = computePoint2CylinderDistance(r_I_m,mpMapTree[map_tree_id]->para);
                        double d_ref = computePoint2CylinderDistance(r_I_m,treePara);
                        res_tree_map_ini = res_tree_map_ini + d_ini * d_ini;
                        res_tree_map = res_tree_map + d_ref * d_ref;
                    }
                    else
                    {
                        int scanId = int(mvpIScans[iscan_id]->vTreePointRaw[treeId][nP].intensity / 100);
                        double t_ratio = double(mvpIScans[iscan_id]->vTreePointRaw[treeId][nP].intensity -
                                                int(mvpIScans[iscan_id]->vTreePointRaw[treeId][nP].intensity));
                        Eigen::Vector3d r_I_lu(mvpIScans[iscan_id]->vTreePointRaw[treeId][nP].x,
                                               mvpIScans[iscan_id]->vTreePointRaw[treeId][nP].y,
                                               mvpIScans[iscan_id]->vTreePointRaw[treeId][nP].z);

                        Eigen::Vector3d r_I_m_ini(mvpIScans[iscan_id]->vTreePointMapping[treeId][nP].x,
                                                  mvpIScans[iscan_id]->vTreePointMapping[treeId][nP].y,
                                                  mvpIScans[iscan_id]->vTreePointMapping[treeId][nP].z);

                        Eigen::Quaterniond q_prev{mvpIScans[iscan_id]->v_R_mapping[scanId]};
                        Eigen::Quaterniond q_cur_end{mvpIScans[iscan_id]->v_R_mapping[scanId + 1]};
                        Eigen::Quaterniond q_cur = q_prev.slerp(t_ratio, q_cur_end);
                        Eigen::Vector3d t_cur = (1.0 - t_ratio) * mvpIScans[iscan_id]->v_r_mapping[scanId] + t_ratio * mvpIScans[iscan_id]->v_r_mapping[scanId + 1];


                        double d_ini = computePoint2CylinderDistance(r_I_m_ini,  mpMapTree[map_tree_id]->para);
                        // double d_ref = computePoint2CylinderDistance(q_cur * r_I_lu + t_cur, treePara);
                        double d_ref = computePoint2CylinderDistance(mvpIScans[iscan_id]->v_R_mapping[scanId] * r_I_lu +
                         mvpIScans[iscan_id]->v_r_mapping[scanId], treePara);
                        res_tree_unknown_ini = res_tree_unknown_ini + d_ini * d_ini;
                        res_tree_unknown = res_tree_unknown + d_ref * d_ref;
                    }
                }
            }

            //update tree parameter
            mpMapTree[map_tree_id]->UpdateTreeParam(treePara)?num_valid_tree_para++: num_invalid_tree_para++;
        }
        fMapLog << "\tNumber of trees with valid and invalid parameters: " << num_valid_tree_para << " / " << num_invalid_tree_para<<endl;
        fMapLog << "\tError for tree point from fix scans: " << sqrt(res_tree_map_ini / double(numFixTreeObs)) << " -> " << sqrt(res_tree_map / double(numFixTreeObs)) << " m" <<endl;
        fMapLog << "\tError for tree point from estimated scans: " << sqrt(res_tree_unknown_ini / double(numUnknownTreeObs)) << " -> " << sqrt(res_tree_unknown / double(numUnknownTreeObs)) << " m" <<endl;

        // 2. ground observation
        double res_ground_ref = 0.0, res_ground_fix = 0.0, res_ground_unknown = 0.0, res_ground_ref_ini = 0.0, res_ground_fix_ini = 0.0, res_ground_unknown_ini = 0.0;
        for (int nPlane = 0; nPlane < numPlane; nPlane++)
        {
            int planeId = nPlane;
            Eigen::Vector4d updated_para{*(planes + planeId * PlaneBlockSize + 0), *(planes + planeId * PlaneBlockSize + 1),
                                         *(planes + planeId * PlaneBlockSize + 2), *(planes + planeId * PlaneBlockSize + 3)};

            for (int nP = 0; nP < point_index_per_grid_[nPlane].size(); nP++)
            {
                int iscan_id = point_index_per_grid_[nPlane][nP].first;
                int ptId = point_index_per_grid_[nPlane][nP].second;
                int epoch_id = iscan_id - iscan_index_start_;

                // if from the iscan that is fix
                if (iscan_id < iscan_index_start_)
                {
                    Eigen::Vector3d r_I_m(mvpIScans[iscan_id]->vGroundPointMapping[ptId].x,
                                          mvpIScans[iscan_id]->vGroundPointMapping[ptId].y,
                                          mvpIScans[iscan_id]->vGroundPointMapping[ptId].z);

                    double d_ini = (r_I_m.transpose() * plane_params_[nPlane].head(3) + plane_params_[nPlane](3)) / plane_params_[nPlane].head(3).norm();
                    double d_ref = (r_I_m.transpose() * updated_para.head(3) + updated_para(3)) / updated_para.head(3).norm();
                    res_ground_fix_ini = res_ground_fix_ini + d_ini * d_ini;
                    res_ground_fix = res_ground_fix + d_ref * d_ref;
                }
                else
                {
                    int scanId = int(mvpIScans[iscan_id]->vGroundPointRaw[ptId].intensity / 100);
                    double t_ratio = double(mvpIScans[iscan_id]->vGroundPointRaw[ptId].intensity -
                                            int(mvpIScans[iscan_id]->vGroundPointRaw[ptId].intensity));
                    Eigen::Vector3d r_I_lu(mvpIScans[iscan_id]->vGroundPointRaw[ptId].x,
                                           mvpIScans[iscan_id]->vGroundPointRaw[ptId].y,
                                           mvpIScans[iscan_id]->vGroundPointRaw[ptId].z);

                    Eigen::Vector3d r_I_m_ini(mvpIScans[iscan_id]->vGroundPointMapping[ptId].x,
                                              mvpIScans[iscan_id]->vGroundPointMapping[ptId].y,
                                              mvpIScans[iscan_id]->vGroundPointMapping[ptId].z);

                    Eigen::Quaterniond q_prev{mvpIScans[iscan_id]->v_R_mapping[scanId]};
                    Eigen::Quaterniond q_cur_end{mvpIScans[iscan_id]->v_R_mapping[scanId + 1]};
                    Eigen::Quaterniond q_cur = q_prev.slerp(t_ratio, q_cur_end);
                    Eigen::Vector3d t_cur = (1.0 - t_ratio) * mvpIScans[iscan_id]->v_r_mapping[scanId] + t_ratio * mvpIScans[iscan_id]->v_r_mapping[scanId + 1];

                    double d_ini = (r_I_m_ini.transpose() * plane_params_[nPlane].head(3) + plane_params_[nPlane](3)) / plane_params_[nPlane].head(3).norm();
                    // double d_ref = ((q_cur * r_I_lu + t_cur).transpose() * updated_para.head(3) + updated_para(3)) / updated_para.head(3).norm();
                    double d_ref = ((mvpIScans[iscan_id]->v_R_mapping[scanId] * r_I_lu + mvpIScans[iscan_id]->v_r_mapping[scanId]).transpose() * 
                    updated_para.head(3) + updated_para(3)) / updated_para.head(3).norm();
                    res_ground_unknown_ini = res_ground_unknown_ini + d_ini * d_ini;
                    res_ground_unknown = res_ground_unknown + d_ref * d_ref;
                }
            }
            for (int nP = 0; nP < point_ref_per_grid_[nPlane].size(); nP++)
            {
                Eigen::Vector3d r_I_m = point_ref_per_grid_[nPlane][nP];
                double d_ini = (r_I_m.transpose() * plane_params_[nPlane].head(3) + plane_params_[nPlane](3)) / plane_params_[nPlane].head(3).norm();
                double d_ref = (r_I_m.transpose() * updated_para.head(3) + updated_para(3)) / updated_para.head(3).norm();
                res_ground_ref_ini = res_ground_ref_ini + d_ini * d_ini;
                res_ground_ref = res_ground_ref + d_ref * d_ref;
            }
        }
        fMapLog << "\tError for ground point for ref planar patches: " << sqrt(res_ground_ref_ini / double(numRefGroundObs)) << " -> " << sqrt(res_ground_ref / double(numRefGroundObs)) << " m" << endl;
        fMapLog << "\tError for ground point for fix planar patches: " << sqrt(res_ground_fix_ini / double(numFixGroundObs)) << " -> " << sqrt(res_ground_fix / double(numFixGroundObs)) << " m" << endl;
        fMapLog << "\tError for ground point for unknown planar patches: " << sqrt(res_ground_unknown_ini / double(numUnknownGroundObs)) << " -> " << sqrt(res_ground_unknown / double(numUnknownGroundObs)) << " m" << endl;


    }//end of iteration
}

/**
 * @brief Updates pose, ground, and tree point mappings for IScans after loop closure optimization.
 *
 * Depending on whether raw mapping is used, it either updates per-scan poses and remaps raw data
 * or recomputes mapping-level features for each scan.
 */

void Mapping::UpdatePosePointsLC()
{
    // TODO incase we are doing two LC together
    if (raw_point_flag_)
    {
        //traj info have updated

        // update the pose, tree points, and ground points (ds) for each updated Iscan
        for (int iscan_id = iscan_index_start_; iscan_id <= iscan_index_end_; iscan_id++)
        {
            mvpIScans[iscan_id]->computeRawGroundPointToMapping();
            mvpIScans[iscan_id]->computeRawTreePointToMapping();
        }

        //TODO Check update mpMapGroundPoint
    }
    else
    {
        // update the relative transformation between iscans (this is not used)
        for (int iscan_id = iscan_index_start_; iscan_id < iscan_index_end_; iscan_id++)
        {
            Eigen::Matrix3d R_local_m_t1 = mvpIScans[iscan_id]->R_local_m_updated;
            Eigen::Matrix3d R_local_m_t2 = mvpIScans[iscan_id + 1]->R_local_m_updated;
            Eigen::Vector3d r_local_m_t1 = mvpIScans[iscan_id]->r_local_m_updated;
            Eigen::Vector3d r_local_m_t2 = mvpIScans[iscan_id + 1]->r_local_m_updated;

            mvpIScans[iscan_id + 1]->R_local_t2_local_t1_updated = R_local_m_t1.inverse() * R_local_m_t2;
            mvpIScans[iscan_id + 1]->r_local_t2_local_t1_updated = R_local_m_t1.inverse() * (r_local_m_t2 - r_local_m_t1);
        }

        // update the pose, tree points, and ground points (ds) for each updated Iscan
        for (int iscan_id = iscan_index_start_; iscan_id <= iscan_index_end_; iscan_id++)
        {
            mvpIScans[iscan_id]->computeIndividualPoseMapping();
            mvpIScans[iscan_id]->computeMapTreePoints(IntegratedScan::REFINED);
            mvpIScans[iscan_id]->computeMapGroundPoints(IntegratedScan::REFINED);
        }    
        
        //update ground point map if no global map
        if(!global_map_flag_)
        {
            mpMapGroundPoint.reset(new pcl::PointCloud<PointType>);

            // TODO in case iscan_index_start_ != 0
            for (int iscan_id = 0; iscan_id <= iscan_index_end_; iscan_id++)
            {
                *mpMapGroundPoint += *mvpIScans[iscan_id]->pGroundPointMappingDs;
            }
            // pcl::VoxelGrid<PointType> downSizeFilter;
            // downSizeFilter.setInputCloud(mpMapGroundPoint);
            // downSizeFilter.setLeafSize(0.1, 0.1, 0.1);
            // downSizeFilter.filter(*mpMapGroundPoint);
            //(TODO:pcl have integrater overflow problem)
            //(TODO:filter size from 0.1 to 0.2)
            *mpMapGroundPoint=DownSamplePointCloudBasedOnDistance(mpMapGroundPoint,0.2);
            fMapLog << "Regenerate ground point map: " << mpMapGroundPoint->size() << endl;
        }
    }
    


}

/**
 * @brief Updates tree parameters for trees not included in the loop closure optimization.
 *
 * After optimization, this function updates all trees that were not part of the loop closure
 * process by re-running `OptimizeLC()` on the remaining set.
 * It also clears temporary planar patch and pose tracking variables.
 */

void Mapping::UpdateMapTreeLC()
{
    //for trees that are included, the parameters are updated. But pose for all Iscans are updated, this will affect the trees that are not included
    fMapLog << "Update remaining trees " <<endl;
    //deactivate planar patch
    vector<vector<Eigen::Vector3d>>().swap(point_ref_per_grid_);
    vector<vector<pair<int, int>>>().swap(point_index_per_grid_);
    vector<Eigen::Vector4d>().swap(plane_params_);

    //deactivate pose
    iscan_index_start_ = iscan_index_end_ + 1; // fix all pose

    //find trees that are not included
    unordered_set<int> used_tree_set;
    for (int nTree = 0; nTree < used_trees_.size(); nTree++)
    {
        used_tree_set.insert(used_trees_[nTree]);
    }

    int num_used_tree = 0;
    vector<int> used_trees; //# vector of unknown tree -> mapTreeID
    for (int nTree = 0; nTree < mpMapTree.size(); nTree++)
    {
        if (mpMapTree[nTree]->status == MapTree::ESTABLISHED || mpMapTree[nTree]->status == MapTree::SOLID ||
            mpMapTree[nTree]->status == MapTree::FITTED || mpMapTree[nTree]->status == MapTree::TBD)
        {
            // not in the used list
            if (used_tree_set.find(nTree) == used_tree_set.end())
            {
                used_trees.push_back(nTree);

                num_used_tree++;
            }
        }
    }
    used_trees_ = used_trees;
    
    //conduct lsa
    OptimizeLC();



}

/**
 * @brief Exports optimized tree and ground points as a LAS file with metadata.
 *
 * This function outputs the tree features and ground planar patches (from both reference
 * and current IScans) into a LAS file. Each point is annotated with scan origin,
 * tree ID (or 0 for ground), and classification flag.
 *
 * @param outPass Output LAS file path.
 */

void Mapping::ExportFeatureLC(const std::string outPass)
{
    PointTable table;
    table.layout()->registerDim(Dimension::Id::X);
    table.layout()->registerDim(Dimension::Id::Y);
    table.layout()->registerDim(Dimension::Id::Z);
    table.layout()->registerDim(Dimension::Id::OriginId);  //iScan index + 1 (0 for reference point cloud)
    table.layout()->registerDim(Dimension::Id::PointId);    //feature id ( 1-> num of tree), 0 for ground
    table.layout()->registerDim(Dimension::Id::ClassFlags); //tree status, 10->ground

    double offsetX = 0.0;
    double offsetY = 0.0;
    double offsetZ = 0.0;
    Options options;
    options.add("filename", outPass);
    options.add("scale_x", (double)1e-3);
    options.add("scale_y", (double)1e-3);
    options.add("scale_z", (double)1e-3);
    options.add("offset_x", offsetX);
    options.add("offset_y", offsetY);
    options.add("offset_z", offsetZ);
    options.add("minor_version", 2);
    options.add("extra_dims", "all");

    PointViewPtr view(new PointView(table));
    int num_point = 0;

    //std::ofstream fFeatureLC(outPass, std::ifstream::out);
    //fFeatureLC << fixed << std::setprecision(8);
    // update info for trees included in optimization
    for (auto i : used_trees_)
    {
        int mapTreeId = i;
        for (int nT = 0; nT < mpMapTree[mapTreeId]->visibleInfo.size(); nT++)
        {
            int iscan_id = mpMapTree[mapTreeId]->visibleInfo[nT].first;
            int treeId = mpMapTree[mapTreeId]->visibleInfo[nT].second;
            for (int nP = 0; nP < mvpIScans[iscan_id]->vTreePointMapping[treeId].size(); nP++)
            {
                // points in mapping
                Eigen::Vector3d r_I_m(mvpIScans[iscan_id]->vTreePointMapping[treeId][nP].x,
                                      mvpIScans[iscan_id]->vTreePointMapping[treeId][nP].y,
                                      mvpIScans[iscan_id]->vTreePointMapping[treeId][nP].z);

                int flag = (iscan_id < iscan_index_start_) ? 0 : 1;

                view->setField(pdal::Dimension::Id::X, num_point, r_I_m(0));
                view->setField(pdal::Dimension::Id::Y, num_point, r_I_m(1));
                view->setField(pdal::Dimension::Id::Z, num_point, r_I_m(2));
                view->setField(pdal::Dimension::Id::PointId, num_point, mapTreeId + 1);
                view->setField(pdal::Dimension::Id::OriginId, num_point, iscan_id + 1);
                view->setField(pdal::Dimension::Id::ClassFlags, num_point, static_cast<int>(mpMapTree[mapTreeId]->status));
                num_point++;


                //fFeatureLC << mapTreeId << "\t" << r_I_m.transpose() << "\t" << flag << "\t" << iscan_id << "\t" << mpMapTree[mapTreeId]->status
                //           << endl;
            }
        }
    }

    for (int nPlane = 0; nPlane < point_index_per_grid_.size(); nPlane++)
    {
        for (int nP = 0; nP < point_index_per_grid_[nPlane].size(); nP++)
        {
            int iscan_id = point_index_per_grid_[nPlane][nP].first;
            int ptId = point_index_per_grid_[nPlane][nP].second;
            Eigen::Vector3d r_I_m;
            if (raw_point_flag_)
            {
                r_I_m << mvpIScans[iscan_id]->vGroundPointMapping[ptId].x,
                    mvpIScans[iscan_id]->vGroundPointMapping[ptId].y,
                    mvpIScans[iscan_id]->vGroundPointMapping[ptId].z;
            }
            else
            {
                r_I_m << mvpIScans[iscan_id]->pGroundPointMappingDs->points[ptId].x,
                    mvpIScans[iscan_id]->pGroundPointMappingDs->points[ptId].y,
                    mvpIScans[iscan_id]->pGroundPointMappingDs->points[ptId].z;
            }
            view->setField(pdal::Dimension::Id::X, num_point, r_I_m(0));
            view->setField(pdal::Dimension::Id::Y, num_point, r_I_m(1));
            view->setField(pdal::Dimension::Id::Z, num_point, r_I_m(2));
            view->setField(pdal::Dimension::Id::PointId, num_point, 0);
            view->setField(pdal::Dimension::Id::OriginId, num_point, iscan_id+1);
            view->setField(pdal::Dimension::Id::ClassFlags, num_point, 10);
            num_point++;

            //fFeatureLC << -1 << "\t" << r_I_m.transpose() << "\t" << 0 << "\t" << iscan_id << "\t" << 0 << endl;
        }

        for (int nP = 0; nP < point_ref_per_grid_[nPlane].size(); nP++)
        {
           Eigen::Vector3d r_I_m =  point_ref_per_grid_[nPlane][nP];
            view->setField(pdal::Dimension::Id::X, num_point, r_I_m(0));
            view->setField(pdal::Dimension::Id::Y, num_point, r_I_m(1));
            view->setField(pdal::Dimension::Id::Z, num_point, r_I_m(2));
            view->setField(pdal::Dimension::Id::PointId, num_point, 0);
            view->setField(pdal::Dimension::Id::OriginId, num_point, 0);
            view->setField(pdal::Dimension::Id::ClassFlags, num_point, 10);
            num_point++;
            //fFeatureLC << -1 << "\t" << point_ref_per_grid_[nPlane][nP].transpose() << "\t" << 0 << "\t" << -1 << "\t" << 0 << endl;
        }
    }
    //fFeatureLC.close();

    BufferReader reader;
    reader.addView(view);
    StageFactory factory;

    // Set second argument to 'true' to let factory take ownership of
    // stage and facilitate clean up.
    Stage *writer = factory.createStage("writers.las");

    writer->setInput(reader);
    writer->setOptions(options);
    writer->prepare(table);
    writer->execute(table);
}

/**
 * @brief Exports updated poses (positions and angles) for each scan after loop closure.
 *
 * This function logs the refined trajectory for all scans included in loop closure.
 * The format includes iscan ID, position, scan index, mapping flag, and rotation angles.
 *
 * @param outPass Output file path (TXT).
 */

void Mapping::ExportTrajectoryLC(const std::string outPass)
{

    std::ofstream fTrajectoryLC(outPass, std::ifstream::out);
    fTrajectoryLC << fixed << std::setprecision(8);
    for (int nIscan = iscan_index_start_; nIscan <= iscan_index_end_; nIscan++)
    {
        // index at t_ini(0), i.e., t_end(last) from previous
        int startIndex = nIscan == 0 ? 0 : mvpIScans[nIscan - 1]->indScans.back().scanID;

        for (int nPose = 0; nPose < mvpIScans[nIscan]->v_R_mapping.size(); nPose++)
        {
            int index = nPose == 0 ? startIndex : mvpIScans[nIscan]->indScans[nPose - 1].scanID;
            fTrajectoryLC << nIscan<<"\t" << mvpIScans[nIscan]->v_r_mapping[nPose].transpose() << "\t"  << index << "\t" 
                         << mvpIScans[nIscan]->flag_iscan_to_map << "\t" << rad2deg(Find_Rotation(mvpIScans[nIscan]->v_R_mapping[nPose])).transpose() << endl;
        }
    }
    fTrajectoryLC.close();
}

/**
 * @brief Exports reference (pre-optimized) poses for trajectory visualization.
 *
 * Similar to `ExportTrajectoryLC`, but outputs reference trajectory stored before optimization.
 * Useful for comparison and debugging.
 *
 * @param outPass Output file path (TXT).
 */

void Mapping::ExportTrajectoryRefLC(const std::string outPass)
{

    std::ofstream fTrajectoryLC(outPass, std::ifstream::out);
    fTrajectoryLC << fixed << std::setprecision(8);
    for (int nIscan = iscan_index_start_; nIscan <= iscan_index_end_; nIscan++)
    {
        // index at t_ini(0), i.e., t_end(last) from previous
        int startIndex = nIscan == 0 ? 0 : mvpIScans[nIscan - 1]->indScans.back().scanID;

        for (int nPose = 0; nPose < mvpIScans[nIscan]->v_R_mapping.size(); nPose++)
        {
            int index = nPose == 0 ? startIndex : mvpIScans[nIscan]->indScans[nPose - 1].scanID;
            fTrajectoryLC << nIscan<<"\t" << mvpIScans[nIscan]->v_r_mapping_ref[nPose].transpose() << "\t"  << index << "\t" 
                         << mvpIScans[nIscan]->flag_iscan_to_map << "\t" << rad2deg(Find_Rotation(mvpIScans[nIscan]->v_R_mapping_ref[nPose])).transpose() << endl;
        }
    }
    fTrajectoryLC.close();
}


//load DTM, and form the mpMapGroundPoint and octreeGroundPointsFromMap
/**
 * @brief Load the global ground map from a file.
 *
 * This function reads the global ground map points from a specified file path,
 * removes any global constant shifts, stores the points in `mpMapGroundPoint`,
 * and constructs a PCL octree for spatial queries.
 *
 * @return true if the map was loaded successfully, false otherwise.
 */
bool Mapping::LoadGlobalGroundMap()
{
    mpMapGroundPoint.reset(new pcl::PointCloud<PointType>);

    if(mPara.global_ground_map_pass.empty())
        return false;

    ifstream f_global_ground_map;
    f_global_ground_map.open(mPara.global_ground_map_pass);

    if (!f_global_ground_map)
    {
        fMapLog << "failed to open" << mPara.global_ground_map_pass << endl;
        throw std::runtime_error("Wrong global ground map file");
        return false;
    }

    PointType p;
    string line;
    while (getline(f_global_ground_map, line))
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

            if (words.size() < 3)
            {
                fMapLog << line;
                fMapLog << "wrong global ground" << words.size() << endl;
                throw std::runtime_error("Wrong global ground map format");
                return false;
            }

            p.x = float(stod(words[0]) - global_const_shift_(0));
            p.y = float(stod(words[1]) - global_const_shift_(1));
            p.z = float(stod(words[2]) - global_const_shift_(2));
            mpMapGroundPoint->push_back(p);
        }
    }

    fMapLog << "Number of points in global map: " << mpMapGroundPoint->size() <<endl;
    
    octreeGroundPointsFromMap.setInputCloud(mpMapGroundPoint);
    octreeGroundPointsFromMap.addPointsFromInputCloud();   
    return true;
}


//load map tree
/**
 * @brief Load the global tree map and add valid trees into the map.
 *
 * Uses the EOP octree to check spatial validity and maps tree locations to heights
 * using the nearest ground point plus a constant height offset.
 *
 * @return true if the global tree map was loaded successfully, false otherwise.
 */
bool Mapping::LoadGlobalTreeMap()
{

    if(mPara.global_tree_map_pass.empty())
        return false;

    ifstream f_global_tree_map;
    f_global_tree_map.open(mPara.global_tree_map_pass);

    if (!f_global_tree_map)
    {
        fMapLog << "failed to open" << mPara.global_tree_map_pass << endl;
        throw std::runtime_error("Wrong global tree map file");
        return false;
    }
            
    Eigen::Vector3d loc;

    PointType p;
    pcl::octree::OctreePointCloudSearch<PointType> OctreeEopSearch(0.2);
    pcl::PointCloud<PointType>::Ptr eop_points (new pcl::PointCloud<PointType>);

    // build eop map
    if (global_map_flag_)
    {
        for (int i = 0; i < mpTraj->bopList.size(); i = i + 10)
        {
            p.x = float(mpTraj->bopList[i].pos.XO - global_const_shift_(0));
            p.y = float(mpTraj->bopList[i].pos.YO - global_const_shift_(1));
            p.z = 0.0;
            eop_points->push_back(p);
        }

        OctreeEopSearch.setInputCloud(eop_points);
        OctreeEopSearch.addPointsFromInputCloud();
    }

    std::vector<int> pointSearchInd;
    std::vector<float> pointSearchSqDis;

    string line;
    while (getline(f_global_tree_map, line))
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

            if (words.size() < 2)
            {
                fMapLog << line;
                fMapLog << "wrong global tree" << words.size() << endl;
                throw std::runtime_error("Wrong global tree map format");
                return false;
            }

            loc(0) = stod(words[0]) - global_const_shift_(0);
            loc(1) = stod(words[1]) - global_const_shift_(1);
            loc(2) = 0.0;

            p.x = loc(0); p.y =loc(1); p.z = loc(2);

            //check if the tree is too far, continue if its too far
            if (global_map_flag_)
            {

                if (OctreeEopSearch.nearestKSearch(p, 1, pointSearchInd, pointSearchSqDis) > 0)
                {
                    if ((eop_points->points[pointSearchInd[0]].x - p.x) * (eop_points->points[pointSearchInd[0]].x - p.x) 
                    + (eop_points->points[pointSearchInd[0]].y - p.y) *(eop_points->points[pointSearchInd[0]].y - p.y) > 80.0 *80.0)
                    {
                        continue;
                    }
                }
            }

            //
           
            if (octreeGroundPointsFromMap.nearestKSearch(p, 1, pointSearchInd, pointSearchSqDis) > 0)
            {
                loc(2) = mpMapGroundPoint->points[pointSearchInd[0]].z + 3.5; // current tree location is derived based on the 2-5 meter range, ground height + 3.5
            }
            

            CylinderPara treePara;
            treePara.x = loc;
            treePara.n << 0.0, 0.0, 1.0;
            treePara.r = 0.0;

            MapTree *pTempTree = new MapTree(mpMapTree.size(), treePara, this);
            mpMapTree.push_back(pTempTree);

        }
    }

    fMapLog << "Number of trees in global tree map: " << mpMapTree.size() <<endl;
    OctreeEopSearch.deleteTree();
    eop_points.reset();
    return true;
}

/***************************************************************
export all feature points from slam tracking thread: mvIntegratedTree
(vTreeRawPoints), pGroundPoints (ground planar points)
flag:
    1. init trajectory from odometry
    2. refined trajectory after integration
    3. reference trajectory
bpoint:
    true:export ground point
***************************************************************/
/**
 * @brief Export all integrated features (trees and optionally ground points) into a TXT file.
 *
 * @param outPass Output file path
 * @param flag Specifies which trajectory version to use (0 = init, 1 = refined, 2 = reference)
 * @param bGroundPoint Flag indicating whether to export ground points
 */
void Mapping::exportIntegratedFeatures(const std::string outPass, int flag, bool bGroundPoint)
{
    std::ofstream fPointFile(outPass, std::ifstream::out);
    fPointFile << fixed << std::setprecision(8);

    vector<Eigen::Vector3d> *vr=nullptr;
    vector<Eigen::Matrix3d> *vR=nullptr;
    if (flag == 0)
    {
        vr = &(mv_r_lu_local_ini);
        vR = &(mv_R_lu_local_ini);
    }
    else if (flag == 1)
    {
        vr = &(mv_r_lu_local);
        vR = &(mv_R_lu_local);
    }
    else if (flag == 2)
    {
        vr = &(mv_r_lu_local_ref);
        vR = &(mv_R_lu_local_ref);
    }

    // ground points
    if (bGroundPoint)
    {
        for (int nScan = 0; nScan < mNumIntegrationScan; nScan++)
        {
            int numPoint = mvTempScan[nScan].vGroundPlanarRawPoints.size();
            int scanId = nScan;
            fDebug<< numPoint <<endl;
            for (int nP = 0; nP < numPoint; nP++)
            {
                Eigen::Vector3d r_I_lu(mvTempScan[nScan].vGroundPlanarRawPoints[nP].x, mvTempScan[nScan].vGroundPlanarRawPoints[nP].y, mvTempScan[nScan].vGroundPlanarRawPoints[nP].z);

                double t_ratio = double(mvTempScan[nScan].vGroundPlanarRawPoints[nP].intensity - int(mvTempScan[nScan].vGroundPlanarRawPoints[nP].intensity));

                Eigen::Quaterniond q_cur_end(vR->at(scanId + 1));
                Eigen::Quaterniond q_prev(vR->at(scanId)); // = Eigen::Quaterniond::Identity(); // intepolated
                Eigen::Quaterniond q_cur = q_prev.slerp(t_ratio, q_cur_end);
                Eigen::Vector3d t_cur = (1.0 - t_ratio) * vr->at(scanId) + t_ratio * vr->at(scanId + 1);

                Eigen::Vector3d r_I_local;
                // r_I_local = q_cur * r_I_lu + t_cur;
                r_I_local = vR->at(scanId) * r_I_lu + vr->at(scanId);

                fPointFile << -1 << "\t" << r_I_local.transpose() << "\t" << 0 << "\t" << float(scanId) << "\t" << endl;
            }
        }
    }


    // tree feautres
    for (int nTree = 0; nTree < mvIntegratedTree.size(); nTree++)
    {
        for (int nScan = 0; nScan < mvIntegratedTreeIds[nTree].size(); nScan++)
        {
            int scanId = mvIntegratedTreeIds[nTree][nScan].first;
            int treeId = mvIntegratedTreeIds[nTree][nScan].second;
            vector<PointType> *pPointList = &(mvTempScan[scanId].vTreeRawPoints[treeId]);
            int numPoint = pPointList->size();
            for (int nP = 0; nP < numPoint; nP++)
            {
                Eigen::Vector3d r_I_lu(pPointList->at(nP).x, pPointList->at(nP).y, pPointList->at(nP).z);
                double t_ratio = double(pPointList->at(nP).intensity - int(pPointList->at(nP).intensity));

                Eigen::Quaterniond q_cur_end(vR->at(scanId + 1));
                Eigen::Quaterniond q_prev(vR->at(scanId)); // = Eigen::Quaterniond::Identity(); // intepolated
                Eigen::Quaterniond q_cur = q_prev.slerp(t_ratio, q_cur_end);
                Eigen::Vector3d t_cur = (1.0 - t_ratio) * vr->at(scanId) + t_ratio * vr->at(scanId + 1);

                Eigen::Vector3d r_I_local;
                // r_I_local = q_cur * r_I_lu + t_cur;
                r_I_local = vR->at(scanId) * r_I_lu + vr->at(scanId);

                if (r_I_local(2) < 10.0)
                {
                    fPointFile << nTree << "\t" << r_I_local.transpose() << "\t" << mvIntegratedTreeIds[nTree].size() << "\t" << float(scanId) << "\t" << endl;
                }
            }
        }
    }

    fPointFile.close();
}

/***************************************************************
Export all points used in optimization: valid mvTreePoints &
mvGridPointIndex of vGroundPlanarRawPoints while considering distortion
flag:
    1. init trajectory from odometry
    2. refined trajectory after integration
    3. reference trajectory
bpoint:
    true:export ground point
***************************************************************/
/**
 * @brief Export features (tree and ground points) used in optimization.
 *
 * @param outPass Output file path
 * @param flag Specifies which trajectory version to use (0 = init, 1 = refined, 2 = reference)
 * @param bGroundPoint Flag indicating whether to export ground points
 */

void Mapping::exportIntegratedOptFeatures(const std::string outPass, int flag, bool bGroundPoint)
{
    std::ofstream fPointFile(outPass, std::ifstream::out);
    fPointFile << fixed << std::setprecision(8);
    
    vector<Eigen::Vector3d> *vr=nullptr; 
    vector<Eigen::Matrix3d> *vR=nullptr;
    if (flag == 0)
    {
        vr = &(mv_r_lu_local_ini);
        vR = &(mv_R_lu_local_ini);
    }
    else if (flag == 1)
    {
        vr = &(mv_r_lu_local);
        vR = &(mv_R_lu_local);
    }
    else if (flag == 2)
    {
        vr = &(mv_r_lu_local_ref);
        vR = &(mv_R_lu_local_ref);
    }

    for (int nTree = 0; nTree < mvIntegratedTree.size(); nTree++)
    {
        if (!mvValidTrees[nTree])
            continue;

        for (int nP = 0; nP < mvTreePointsRaw[nTree].size(); nP++)
        {
            int scanId = int(mvTreePointsRaw[nTree][nP].intensity / 100);
            Eigen::Vector3d r_I_lu(mvTreePointsRaw[nTree][nP].x, mvTreePointsRaw[nTree][nP].y, mvTreePointsRaw[nTree][nP].z);

            double t_ratio = double(mvTreePointsRaw[nTree][nP].intensity - int(mvTreePointsRaw[nTree][nP].intensity));

            Eigen::Quaterniond q_cur_end(vR->at(scanId+1));
            Eigen::Quaterniond q_prev(vR->at(scanId)); // = Eigen::Quaterniond::Identity(); // intepolated
            Eigen::Quaterniond q_cur = q_prev.slerp(t_ratio, q_cur_end);
            Eigen::Vector3d t_cur = (1.0 - t_ratio) * vr->at(scanId) + t_ratio * vr->at(scanId+1);

            Eigen::Vector3d r_I_local;
            // r_I_local = q_cur * r_I_lu + t_cur;
            r_I_local = vR->at(scanId) * r_I_lu + vr->at(scanId);

            fPointFile << nTree << "\t" << r_I_local.transpose() << "\t" << mvIntegratedTreeIds[nTree].size() << "\t" << float(scanId) 
            << "\t" <<mvTreePointsRaw[nTree].size() << endl;
        }
    }

    // ground points
    if (bGroundPoint)
    {
        for (int nPlane = 0; nPlane < mvGridPointIndex.size(); nPlane++)
        {
            for (int nP = 0; nP < mvGridPointIndex[nPlane].size(); nP++)
            {
                int index = mvGridPointIndex[nPlane][nP];
                int scanId = int(mvGroundPointsRaw[index].intensity / 100);
                Eigen::Vector3d r_I_lu(mvGroundPointsRaw[index].x, mvGroundPointsRaw[index].y, mvGroundPointsRaw[index].z);

                double t_ratio = double(mvGroundPointsRaw[index].intensity - int(mvGroundPointsRaw[index].intensity));

                Eigen::Quaterniond q_cur_end(vR->at(scanId + 1));
                Eigen::Quaterniond q_prev(vR->at(scanId)); // = Eigen::Quaterniond::Identity(); // intepolated
                Eigen::Quaterniond q_cur = q_prev.slerp(t_ratio, q_cur_end);
                Eigen::Vector3d t_cur = (1.0 - t_ratio) * vr->at(scanId) + t_ratio * vr->at(scanId + 1);

                Eigen::Vector3d r_I_local;
                // r_I_local = q_cur * r_I_lu + t_cur;
                r_I_local = vR->at(scanId) * r_I_lu + vr->at(scanId);

                fPointFile << -1 << "\t" << r_I_local.transpose() << "\t" << -1 << "\t" << float(scanId) << "\t" << mvGridPointIndex[nPlane].size() << endl;
            }
        }
    }
    fPointFile.close();
}
/***************************************************************
Export all points broadcast to map: 

bpoint:
    true:export ground point
***************************************************************/
/**
 * @brief Export features from the current IScan for broadcasting.
 *
 * @param outPass Output file path
 * @param ds_flag Use downsampled ground points if true, else use raw points
 */
void Mapping::exportIscanFeatureBroadcast(const std::string outPass, int ds_flag)
{
    std::ofstream fPointFile(outPass, std::ifstream::out);
    fPointFile << fixed << std::setprecision(8);
    
    for (int nTree = 0; nTree < pCurrentIScan->vTreePointLocal.size(); nTree++)
    {
        // for each tree point
        for (int nP = 0; nP < pCurrentIScan->vTreePointLocal[nTree].size(); nP++)
        {
            fPointFile << nTree << "\t" << pCurrentIScan->vTreePointLocal[nTree][nP].x << "\t" << pCurrentIScan->vTreePointLocal[nTree][nP].y
                       << "\t" << pCurrentIScan->vTreePointLocal[nTree][nP].z << endl;
        }
    }

    if (ds_flag)
    {
        for (int nP = 0; nP < pCurrentIScan->pGroundPointLocalDs->size(); nP++)
        {
            fPointFile << -1 << "\t" << pCurrentIScan->pGroundPointLocalDs->points[nP].x << "\t" << pCurrentIScan->pGroundPointLocalDs->points[nP].y
                       << "\t" << pCurrentIScan->pGroundPointLocalDs->points[nP].z << endl;
        }
    }
    else
    {
        for (int nP = 0; nP < pCurrentIScan->vGroundPointLocal.size(); nP++)
        {
            fPointFile << -1 << "\t" << pCurrentIScan->vGroundPointLocal[nP].x << "\t" << pCurrentIScan->vGroundPointLocal[nP].y
                       << "\t" << pCurrentIScan->vGroundPointLocal[nP].z << endl;
        }
    }
    fPointFile.close();
}

/***************************************************************
 * Export tree parameters
 **************************************************************/
void Mapping::exportTreeModel(const std::string outPass)
{
    std::ofstream fTreeFile(outPass, std::ifstream::out);
    fTreeFile << fixed << std::setprecision(8);
    for (int nTree = 0; nTree < mvIntegratedTree.size(); nTree++)
    {
        // height range: 0 - 5
        for (int i = 0; i < 50; i++)
        {
            double lambda = (double(i) * 0.1 - mvTreeParamsLocal[nTree].x(2)) / mvTreeParamsLocal[nTree].n(2);
            Eigen::Vector3d p = mvTreeParamsLocal[nTree].x + lambda * mvTreeParamsLocal[nTree].n;
            fTreeFile << nTree << "\t" << p.transpose() << "\t" << mvValidTrees[nTree] << endl;
        }
    }
    fTreeFile.close();
}

/***************************************************************
 * Export the line representing axis of the map tree in the mapping frame
 * Option:
 *  Mapping frame
 *  Global frame (only valid if trajectory is available)
 **************************************************************/
/**
 * @brief Export tree axis lines for map trees in either the local mapping frame or global frame.
 *
 * @param outPass Output file path
 * @param exportOption MAPPING_FRAME or GLOBAL_FRAME
 */
void Mapping::exportMapTreeModel(const std::string outPass, eExportOption exportOption )
{
    std::ofstream fTreeFile(outPass, std::ifstream::out);
    fTreeFile << fixed << std::setprecision(8);
    fTreeFile << "TreeID\t" << "X\t"<< "Y\t"<< "Z\t" << "r\t" <<
     "Status\t" << "SurfaceRatio\t"<< "NumberIscan\t" << "NumerPts" <<endl;
    for (int nT = 0; nT < mpMapTree.size(); nT++)
    {
        if(mpMapTree[nT]->status>= 0)
        {        
            //  height range: 0 - 5
            for (int i = -25; i < 25; i++)
            {
                double lambda = (double(i) * 0.1 + mpMapTree[nT]->center(2) - mpMapTree[nT]->para.x(2)) / mpMapTree[nT]->para.n(2);
                Eigen::Vector3d p = mpMapTree[nT]->para.x + lambda * mpMapTree[nT]->para.n;
                Eigen::Vector3d p_;
                if (exportOption == MAPPING_FRAME)
                    p_ = p;
                else if (exportOption == GLOBAL_FRAME)
                    p_ = r_m_global + R_m_global* p;

                fTreeFile << nT + 1 << "\t" << p_.transpose() << "\t" << mpMapTree[nT]->para.r << "\t" << mpMapTree[nT]->status
                        << "\t" << mpMapTree[nT]->surfaceRatio << "\t" << mpMapTree[nT]->visibleInfo.size() 
                        << "\t" << mpMapTree[nT]->numPoint << endl; //"\t" << mpMapTree[nT]->groupId  << endl;
            }
        }
    }
    fTreeFile.close();
}


/***************************************************************
 * Export map tree parameters
 **************************************************************/
/**
 * @brief Export cylinder parameters of all map trees.
 *
 * @param outPass Output file path
 * @param exportOption MAPPING_FRAME or GLOBAL_FRAME
 */
void Mapping::exportMapTreeParam(const std::string outPass, eExportOption exportOption)
{
    std::ofstream fTreeFile(outPass, std::ifstream::out);
    fTreeFile << fixed << std::setprecision(6);

    fTreeFile << "TreeID\t" << "X\t"<< "Y\t"<< "Z\t" << "Nx\t"<< "Ny\t"<< "Nz\t" << "r\t" << "Cx\t" << "Cy\t" << "Cz\t" 
    << "Status\t" << "SurfaceRatio\t"<< "NumberIscan\t" << "NumerPts\t" << "Residual\t" <<  "Ref_Tree_Flag\t" 
    << "Outlier_negative\t" << "Inlier\t" <<  "Outlier_Positive" <<endl;
    for (int nT = 0; nT < mpMapTree.size(); nT++)
    {
        if (mpMapTree[nT]->status >= 0)
        {
            Eigen::Vector3d x_, x = mpMapTree[nT]->para.x;
            Eigen::Vector3d n_, n = mpMapTree[nT]->para.n;
            Eigen::Vector3d c_, c = mpMapTree[nT]->center;
            if (exportOption == MAPPING_FRAME)
            {
                x_ = x;
                n_ = n;
                c_ = c;
            }
            else if (exportOption == GLOBAL_FRAME)
            {
                x_ = r_m_global + R_m_global * x;
                n_ = R_m_global * n;
                n_ = n_ / n_(2);
                c_ = r_m_global + R_m_global * c;
            }
            fTreeFile << nT + 1 << "\t" << x_(0) << "\t" << x_(1) << "\t" << x_(2) << "\t"
                      << n_(0) << "\t" << n_(1) << "\t" << n_(2) << "\t"
                      << mpMapTree[nT]->para.r << "\t" << c_(0) << "\t" << c_(1) << "\t" << c_(2) << "\t"
                      << mpMapTree[nT]->status << "\t" << mpMapTree[nT]->surfaceRatio << "\t" << mpMapTree[nT]->visibleInfo.size() << "\t" 
                      << mpMapTree[nT]->numPoint << "\t" << mpMapTree[nT]->rmse <<  "\t" << mpMapTree[nT]->map_tree_flag <<  "\t" 
                      << mpMapTree[nT]->residual_distribution[0] << "\t" << mpMapTree[nT]->residual_distribution[1] <<  "\t" 
                      << mpMapTree[nT]->residual_distribution[2]   <<endl;
        }
    }
    fTreeFile.close();
}

/***************************************************************
 * Export all tree points and ground points in mapping 
 * (tree points in map tree, ground points after downsampling in the integrated scan)
 * Option:
 *  Mapping frame
 *  Global frame (only valid if trajectory is available)
 * Input:
 *  mpMapTree, mvpIScans[iscanId]->vTreePointMapping, mvpIScans[iscanId]->pGroundPointMappingDs
 **************************************************************/
/**
 * @brief Export all tree and ground feature points used in mapping.
 *
 * @param outPass Output file path
 * @param exportOption MAPPING_FRAME or GLOBAL_FRAME*/

void Mapping::exportFeatureMapping(const std::string outPass,  eExportOption exportOption)
{
	//feature_las_points_.clear();
    std::ofstream fPointFile(outPass, std::ifstream::out);
    fPointFile << fixed << std::setprecision(8);
    fPointFile << "TreeID\t" << "X\t"<< "Y\t"<< "Z\t" << "IscanID\t"<< "status\t" <<"angle\t"<< "dis\t" <<endl;

    for (int nMapT = 0; nMapT < mpMapTree.size(); nMapT++)
    {
        int mapTreeId = nMapT;

        //for global tree or deactivated tree
        if(mpMapTree[mapTreeId]->status < 0)
            continue;

        // 1. Tree points from map
        for (int nT = 0; nT < mpMapTree[mapTreeId]->visibleInfo.size(); nT++)
        {   
            int iscanId = mpMapTree[mapTreeId]->visibleInfo[nT].first;
            int treeId = mpMapTree[mapTreeId]->visibleInfo[nT].second;
            for (int nP = 0; nP < mvpIScans[iscanId]->vTreePointMapping[treeId].size(); nP++)
            {
                // points in mapping
                Eigen::Vector3d curr_point(mvpIScans[iscanId]->vTreePointMapping[treeId][nP].x,
                                           mvpIScans[iscanId]->vTreePointMapping[treeId][nP].y,
                                           mvpIScans[iscanId]->vTreePointMapping[treeId][nP].z);

                Eigen::Vector3d p_;
				if (exportOption == MAPPING_FRAME)
				{
					p_ = curr_point;
				}
				else if (exportOption == GLOBAL_FRAME)
				{
					p_ = r_m_global + R_m_global * curr_point;
					//feature_las_points_.push_back(p_);
				}

                fPointFile << mapTreeId << "\t" << p_.transpose() << "\t" << float(iscanId) << "\t" << mpMapTree[mapTreeId]->status << "\t"
                    << mpMapTree[mapTreeId]->compute_angle(curr_point) << "\t" << abs(mpMapTree[mapTreeId]->compute_dis(curr_point)) 
                    <<  endl;
            }
        }
    }

    for (int nIscan = 0; nIscan < mvpIScans.size(); nIscan++)
    {
        int iscanId = nIscan;
        if(raw_point_flag_)
        {
            for (int nP = 0; nP < mvpIScans[iscanId]->vGroundPointMapping.size(); nP++)
            {

                Eigen::Vector3d curr_point(mvpIScans[iscanId]->vGroundPointMapping[nP].x,
                                           mvpIScans[iscanId]->vGroundPointMapping[nP].y,
                                           mvpIScans[iscanId]->vGroundPointMapping[nP].z);
                Eigen::Vector3d p_;
                if (exportOption == MAPPING_FRAME)
                    p_ = curr_point;
                else if (exportOption == GLOBAL_FRAME)
                    p_ = r_m_global + R_m_global * curr_point;

                fPointFile << -1 << "\t" << p_.transpose() << "\t" << float(iscanId) << "\t" << 0 << "\t" << 0 << "\t" << 0 << endl; //"\t" << 0 << "\t" << 0 << endl;
            }
        }
        else
        {
            for (int nP = 0; nP < mvpIScans[iscanId]->pGroundPointMappingDs->size(); nP++)
            {

                Eigen::Vector3d curr_point(mvpIScans[iscanId]->pGroundPointMappingDs->points[nP].x,
                                           mvpIScans[iscanId]->pGroundPointMappingDs->points[nP].y,
                                           mvpIScans[iscanId]->pGroundPointMappingDs->points[nP].z);
                Eigen::Vector3d p_;
                if (exportOption == MAPPING_FRAME)
                    p_ = curr_point;
                else if (exportOption == GLOBAL_FRAME)
                    p_ = r_m_global + R_m_global * curr_point;

                fPointFile << -1 << "\t" << p_.transpose() << "\t" << float(iscanId) << "\t" << 0 << "\t" << 0 << "\t" << 0 << endl; //"\t" << 0 << "\t" << 0 << endl;
            }
        }
    }
    fPointFile.close();
}

/**
 * @brief Export the mapping features (trees and ground) to a LAS file.
 *
 * Includes extra PDAL fields like GpsTime, EchoRange, OriginId, etc.
 *
 * @param outPass Output LAS file path
 * @param exportOption MAPPING_FRAME or GLOBAL_FRAME
 */
void Mapping::ExportFeatureMappingLas(const std::string outPass, eExportOption exportOption)
{
    PointTable table;
    table.layout()->registerDim(Dimension::Id::X);
    table.layout()->registerDim(Dimension::Id::Y);
    table.layout()->registerDim(Dimension::Id::Z);
    table.layout()->registerDim(Dimension::Id::PointSourceId);
    table.layout()->registerDim(Dimension::Id::Intensity);

    table.layout()->registerDim(Dimension::Id::OriginId);
    table.layout()->registerDim(Dimension::Id::GpsTime);
    table.layout()->registerDim(Dimension::Id::EchoRange);
    table.layout()->registerDim(Dimension::Id::PointId);
    table.layout()->registerDim(Dimension::Id::ClassFlags);
    table.layout()->registerDim(Dimension::Id::NNDistance);

    double offsetX = 0.0;
    double offsetY = 0.0;
    double offsetZ = 0.0;

    if (exportOption == GLOBAL_FRAME)
    {
       offsetX = trunc(r_m_global(0));
       offsetY = trunc(r_m_global(1));
       offsetZ = trunc(r_m_global(2));  
    }

    Options options;
    options.add("filename", outPass);
    options.add("scale_x", (double)1e-3);
    options.add("scale_y", (double)1e-3);
    options.add("scale_z", (double)1e-3);
    options.add("offset_x", offsetX);
    options.add("offset_y", offsetY);
    options.add("offset_z", offsetZ);
    options.add("minor_version", 4);
    options.add("extra_dims", "all");

    PointViewPtr view(new PointView(table));
    int num_point = 0;
    for (int nMapT = 0; nMapT < mpMapTree.size(); nMapT++)
    {
        int mapTreeId = nMapT;

        //for global tree or deactivated tree
        if(mpMapTree[mapTreeId]->status < 0)
            continue;

        // 1. Tree points from map
        for (int nT = 0; nT < mpMapTree[mapTreeId]->visibleInfo.size(); nT++)
        {   
            int iscanId = mpMapTree[mapTreeId]->visibleInfo[nT].first;
            int treeId = mpMapTree[mapTreeId]->visibleInfo[nT].second;
            for (int nP = 0; nP < mvpIScans[iscanId]->vTreePointMapping[treeId].size(); nP++)
            {
                int scan_index_in_iscan = int(mvpIScans[iscanId]->vTreePointRaw[treeId][nP].intensity / 100);
                int scan_id = mvpIScans[iscanId]->indScans[scan_index_in_iscan].scanID;
                int laser_beam_id = int(mvpIScans[iscanId]->vTreePointRaw[treeId][nP].intensity)%100;
                Eigen::Vector3d r_I_lu(mvpIScans[iscanId]->vTreePointRaw[treeId][nP].x,
                                       mvpIScans[iscanId]->vTreePointRaw[treeId][nP].y,
                                       mvpIScans[iscanId]->vTreePointRaw[treeId][nP].z);

                double t_ratio = double(mvpIScans[iscanId]->vTreePointRaw[treeId][nP].intensity -
                                        int(mvpIScans[iscanId]->vTreePointRaw[treeId][nP].intensity));

                double t = mvpIScans[iscanId]->indScans[scan_index_in_iscan].mTimeInit + t_ratio * mvpIScans[iscanId]->indScans[scan_index_in_iscan].mTimeTracking;


                // points in mapping
                Eigen::Vector3d r_I_m(mvpIScans[iscanId]->vTreePointMapping[treeId][nP].x,
                                           mvpIScans[iscanId]->vTreePointMapping[treeId][nP].y,
                                           mvpIScans[iscanId]->vTreePointMapping[treeId][nP].z);

                Eigen::Vector3d p_;
				if (exportOption == MAPPING_FRAME)
				{
					p_ = r_I_m;
				}
				else if (exportOption == GLOBAL_FRAME)
				{
					p_ = r_m_global + R_m_global * r_I_m;
					//feature_las_points_.push_back(p_);
				}
                view->setField(pdal::Dimension::Id::X, num_point, p_(0));
                view->setField(pdal::Dimension::Id::Y, num_point, p_(1));
                view->setField(pdal::Dimension::Id::Z, num_point, p_(2));
                view->setField(pdal::Dimension::Id::GpsTime, num_point, t);
                view->setField(pdal::Dimension::Id::EchoRange, num_point, r_I_lu.norm());
                view->setField(pdal::Dimension::Id::PointId, num_point, mapTreeId + 1);
                view->setField(pdal::Dimension::Id::OriginId, num_point, scan_id);
                view->setField(pdal::Dimension::Id::PointSourceId, num_point, laser_beam_id);
                view->setField(pdal::Dimension::Id::ClassFlags, num_point, static_cast<int>(mpMapTree[mapTreeId]->status));
                view->setField(pdal::Dimension::Id::NNDistance, num_point, mpMapTree[mapTreeId]->compute_dis(r_I_m));

                num_point++;
                // fPointFile << mapTreeId << "\t" << p_.transpose() << "\t" << float(iscanId) << "\t" << mpMapTree[mapTreeId]->status << "\t"
                //     << mpMapTree[mapTreeId]->compute_angle(r_I_m) << "\t" << abs(mpMapTree[mapTreeId]->compute_dis(curr_point)) 
                //     <<  endl;
            }
        }
    }

    for (int nIscan = 0; nIscan < mvpIScans.size(); nIscan++)
    {
        int iscanId = nIscan;
        if(raw_point_flag_)
        {
            for (int nP = 0; nP < mvpIScans[iscanId]->vGroundPointMapping.size(); nP++)
            {
                int scan_index_in_iscan = int(mvpIScans[iscanId]->vGroundPointRaw[nP].intensity / 100);
                int scan_id = mvpIScans[iscanId]->indScans[scan_index_in_iscan].scanID;
                int laser_beam_id = int(mvpIScans[iscanId]->vGroundPointRaw[nP].intensity) % 100;
                Eigen::Vector3d r_I_lu(mvpIScans[iscanId]->vGroundPointRaw[nP].x,
                                      mvpIScans[iscanId]->vGroundPointRaw[nP].y,
                                      mvpIScans[iscanId]->vGroundPointRaw[nP].z);
                double t_ratio = double(mvpIScans[iscanId]->vGroundPointRaw[nP].intensity -
                                        int(mvpIScans[iscanId]->vGroundPointRaw[nP].intensity));

                Eigen::Vector3d r_I_m(mvpIScans[iscanId]->vGroundPointMapping[nP].x,
                                      mvpIScans[iscanId]->vGroundPointMapping[nP].y,
                                      mvpIScans[iscanId]->vGroundPointMapping[nP].z);
                double t = mvpIScans[iscanId]->indScans[scan_index_in_iscan].mTimeInit + t_ratio * mvpIScans[iscanId]->indScans[scan_index_in_iscan].mTimeTracking;

                Eigen::Vector3d p_;
                if (exportOption == MAPPING_FRAME)
                {
                    p_ = r_I_m;
                }
                else if (exportOption == GLOBAL_FRAME)
                {
                    p_ = r_m_global + R_m_global * r_I_m;
                }
                view->setField(pdal::Dimension::Id::X, num_point, p_(0));
                view->setField(pdal::Dimension::Id::Y, num_point, p_(1));
                view->setField(pdal::Dimension::Id::Z, num_point, p_(2));
                view->setField(pdal::Dimension::Id::GpsTime, num_point, t);
                view->setField(pdal::Dimension::Id::EchoRange, num_point, r_I_lu.norm());
                view->setField(pdal::Dimension::Id::PointId, num_point, 0);
                view->setField(pdal::Dimension::Id::OriginId, num_point, scan_id);
                view->setField(pdal::Dimension::Id::PointSourceId, num_point, laser_beam_id);
                view->setField(pdal::Dimension::Id::ClassFlags, num_point, 10);
                view->setField(pdal::Dimension::Id::NNDistance, num_point, 0.0);
                num_point++;
                //fPointFile << -1 << "\t" << p_.transpose() << "\t" << float(iscanId) << "\t" << 0 << "\t" << 0 << "\t" << 0 << endl; //"\t" << 0 << "\t" << 0 << endl;
            }
        }
        else
        {
            for (int nP = 0; nP < mvpIScans[iscanId]->pGroundPointMappingDs->size(); nP++)
            {

                Eigen::Vector3d curr_point(mvpIScans[iscanId]->pGroundPointMappingDs->points[nP].x,
                                           mvpIScans[iscanId]->pGroundPointMappingDs->points[nP].y,
                                           mvpIScans[iscanId]->pGroundPointMappingDs->points[nP].z);
                Eigen::Vector3d p_;
                if (exportOption == MAPPING_FRAME)
                    p_ = curr_point;
                else if (exportOption == GLOBAL_FRAME)
                    p_ = r_m_global + R_m_global * curr_point;
                view->setField(pdal::Dimension::Id::X, num_point, p_(0));
                view->setField(pdal::Dimension::Id::Y, num_point, p_(1));
                view->setField(pdal::Dimension::Id::Z, num_point, p_(2));
                view->setField(pdal::Dimension::Id::GpsTime, num_point, 0.0);
                view->setField(pdal::Dimension::Id::EchoRange, num_point, 0.0);
                view->setField(pdal::Dimension::Id::PointId, num_point, 0);
                view->setField(pdal::Dimension::Id::OriginId, num_point, iscanId);
                view->setField(pdal::Dimension::Id::PointSourceId, num_point, 0);
                view->setField(pdal::Dimension::Id::ClassFlags, num_point, 10);
                view->setField(pdal::Dimension::Id::NNDistance, num_point, 0.0);
                //fPointFile << -1 << "\t" << p_.transpose() << "\t" << float(iscanId) << "\t" << 0 << "\t" << 0 << "\t" << 0 << endl; //"\t" << 0 << "\t" << 0 << endl;
            }
        }
    }

    BufferReader reader;
    reader.addView(view);
    StageFactory factory;

    // Set second argument to 'true' to let factory take ownership of
    // stage and facilitate clean up.
    Stage *writer = factory.createStage("writers.las");

    writer->setInput(reader);
    writer->setOptions(options);
    writer->prepare(table);
    writer->execute(table);
}

/**
 * @brief Exports backup tree point cloud features to a LAS file.
 *
 * This method exports raw tree points (used as backups) for each IScan, including metadata like
 * time, echo range, scan IDs, and tree mapping information. Outputs are transformed into either
 * mapping or global frame.
 *
 * @param outPass Path to output LAS file.
 * @param exportOption Choose between MAPPING_FRAME or GLOBAL_FRAME for coordinate frame.
 */

void Mapping::ExportBackupTreeFeatureMappingLas(const std::string outPass, eExportOption exportOption)
{
    PointTable table;
    table.layout()->registerDim(Dimension::Id::X);
    table.layout()->registerDim(Dimension::Id::Y);
    table.layout()->registerDim(Dimension::Id::Z);
    table.layout()->registerDim(Dimension::Id::PointSourceId);
    table.layout()->registerDim(Dimension::Id::Intensity);

    table.layout()->registerDim(Dimension::Id::OriginId);
    table.layout()->registerDim(Dimension::Id::GpsTime);
    table.layout()->registerDim(Dimension::Id::EchoRange);
    table.layout()->registerDim(Dimension::Id::PointId);
    table.layout()->registerDim(Dimension::Id::ClassFlags);
    table.layout()->registerDim(Dimension::Id::NNDistance);

    double offsetX = 0.0;
    double offsetY = 0.0;
    double offsetZ = 0.0;

    if (exportOption == GLOBAL_FRAME)
    {
       offsetX = trunc(r_m_global(0));
       offsetY = trunc(r_m_global(1));
       offsetZ = trunc(r_m_global(2));  
    }

    Options options;
    options.add("filename", outPass);
    options.add("scale_x", (double)1e-3);
    options.add("scale_y", (double)1e-3);
    options.add("scale_z", (double)1e-3);
    options.add("offset_x", offsetX);
    options.add("offset_y", offsetY);
    options.add("offset_z", offsetZ);
    options.add("minor_version", 4);
    options.add("extra_dims", "all");

    PointViewPtr view(new PointView(table));
    int num_point = 0;
    for (int nIscan = 0; nIscan < mvpIScans.size(); nIscan++)
    {
        int iscanId = nIscan;
        for (int nScan = 0; nScan < mvpIScans[iscanId]->raw_tree_points_mapping.size(); nScan++) //individual scan
        {
            for (int nTree = 0; nTree < mvpIScans[iscanId]->raw_tree_points_mapping[nScan].size(); nTree++) //individual tree
            {
                int treeId = nTree;
                // for each tree point
                for (int nP = 0; nP < mvpIScans[iscanId]->raw_tree_points_mapping[nScan][treeId].size(); nP++)
                {
                    // points in mapping
                    Eigen::Vector3d r_I_m(mvpIScans[iscanId]->raw_tree_points_mapping[nScan][treeId][nP].x,
                                          mvpIScans[iscanId]->raw_tree_points_mapping[nScan][treeId][nP].y,
                                          mvpIScans[iscanId]->raw_tree_points_mapping[nScan][treeId][nP].z);

                    Eigen::Vector3d p_;
                    if (exportOption == MAPPING_FRAME)
                    {
                        p_ = r_I_m;
                    }
                    else if (exportOption == GLOBAL_FRAME)
                    {
                        p_ = r_m_global + R_m_global * r_I_m;
                    }
                    Eigen::Vector3d r_I_lu(mvpIScans[iscanId]->raw_tree_points[nScan][treeId][nP].x,
                                           mvpIScans[iscanId]->raw_tree_points[nScan][treeId][nP].y,
                                           mvpIScans[iscanId]->raw_tree_points[nScan][treeId][nP].z);
                    double t_ratio = double(mvpIScans[iscanId]->raw_tree_points_mapping[nScan][treeId][nP].intensity -
                                        int(mvpIScans[iscanId]->raw_tree_points_mapping[nScan][treeId][nP].intensity));
                    int laser_beam_id = int(mvpIScans[iscanId]->raw_tree_points_mapping[nScan][treeId][nP].intensity)%100;
                    int scan_index_in_iscan = nScan;
                    int scan_id = mvpIScans[iscanId]->indScans[scan_index_in_iscan].scanID;

                    double t = mvpIScans[iscanId]->indScans[scan_index_in_iscan].mTimeInit + t_ratio * mvpIScans[iscanId]->indScans[scan_index_in_iscan].mTimeTracking;

                    view->setField(pdal::Dimension::Id::X, num_point, p_(0));
                    view->setField(pdal::Dimension::Id::Y, num_point, p_(1));
                    view->setField(pdal::Dimension::Id::Z, num_point, p_(2));
                    view->setField(pdal::Dimension::Id::GpsTime, num_point, t);
                    view->setField(pdal::Dimension::Id::EchoRange, num_point, r_I_lu.norm());
                    view->setField(pdal::Dimension::Id::PointId, num_point, mvpIScans[iscanId]->matched_map_tree_id[nScan][treeId]+1);
                    view->setField(pdal::Dimension::Id::OriginId, num_point, scan_id);
                    view->setField(pdal::Dimension::Id::PointSourceId, num_point, laser_beam_id);
                    num_point++;
                }
            }
        }
    }

    BufferReader reader;
    reader.addView(view);
    StageFactory factory;

    // Set second argument to 'true' to let factory take ownership of
    // stage and facilitate clean up.
    Stage *writer = factory.createStage("writers.las");

    writer->setInput(reader);
    writer->setOptions(options);
    writer->prepare(table);
    writer->execute(table);
}


/***************************************************************
 * Export all tree points and ground points in mapping using the reference trajectory
 * Input: mpMapTree, mvpIScans[iscanId]->vGroundPointRaw/vTreePointRaw
 * Option:
 *  Mapping frame
 *  Global frame (only valid if trajectory is available)
 * (tree points in map tree, all ground points in the integrated scan)
 **************************************************************/
void Mapping::exportFeatureMappingRef(const std::string outPass, eExportOption exportOption)
{
    std::ofstream fPointFile(outPass, std::ifstream::out);
    fPointFile << fixed << std::setprecision(8);
    for (int nIscan = 0; nIscan < mvpIScans.size(); nIscan++)
    {
        int iscanId = nIscan;
        for (int nP = 0; nP < mvpIScans[iscanId]->vGroundPointRaw.size(); nP++)
        {

            int scanId = int(mvpIScans[iscanId]->vGroundPointRaw[nP].intensity / 100);
            Eigen::Vector3d curr_point(mvpIScans[iscanId]->vGroundPointRaw[nP].x, mvpIScans[iscanId]->vGroundPointRaw[nP].y, mvpIScans[iscanId]->vGroundPointRaw[nP].z);
            double t_ratio = double(mvpIScans[iscanId]->vGroundPointRaw[nP].intensity - int(mvpIScans[iscanId]->vGroundPointRaw[nP].intensity));

            Eigen::Quaterniond q_prev{mvpIScans[iscanId]->v_R_mapping_ref[scanId]};
            Eigen::Quaterniond q_cur_end{mvpIScans[iscanId]->v_R_mapping_ref[scanId + 1]};
            Eigen::Quaterniond q_cur = q_prev.slerp(t_ratio, q_cur_end);
            Eigen::Vector3d t_cur = (1.0 - t_ratio) * mvpIScans[iscanId]->v_r_mapping_ref[scanId] + t_ratio * mvpIScans[iscanId]->v_r_mapping_ref[scanId + 1];
            // Eigen::Vector3d rIm = q_cur * curr_point + t_cur;
            Eigen::Vector3d rIm = mvpIScans[iscanId]->v_R_mapping_ref[scanId] * curr_point + mvpIScans[iscanId]->v_r_mapping_ref[scanId];
            
            Eigen::Vector3d p_;
            if (exportOption == MAPPING_FRAME)
                p_ = rIm;
            else if (exportOption == GLOBAL_FRAME)
                p_ = r_m_global_ref_traj + R_m_global_ref_traj * rIm;

            fPointFile << -1 << "\t" << p_.transpose() << "\t" << float(iscanId) << endl;
        }
    }
    for (int nMapT = 0; nMapT < mpMapTree.size(); nMapT++)
    {
        int mapTreeId = nMapT;
        
        //for global tree or deactivated tree
        if(mpMapTree[mapTreeId]->status < 0)
            continue;

        for (int nT = 0; nT < mpMapTree[mapTreeId]->visibleInfo.size(); nT++)
        {
            int iscanId = mpMapTree[mapTreeId]->visibleInfo[nT].first;
            int treeId = mpMapTree[mapTreeId]->visibleInfo[nT].second;
            for (int nP = 0; nP < mvpIScans[iscanId]->vTreePointRaw[treeId].size(); nP++)
            {
                int scanId = int(mvpIScans[iscanId]->vTreePointRaw[treeId][nP].intensity / 100);

                Eigen::Vector3d curr_point(mvpIScans[iscanId]->vTreePointRaw[treeId][nP].x, mvpIScans[iscanId]->vTreePointRaw[treeId][nP].y, 
                mvpIScans[iscanId]->vTreePointRaw[treeId][nP].z);
                double t_ratio = double(mvpIScans[iscanId]->vTreePointRaw[treeId][nP].intensity - int(mvpIScans[iscanId]->vTreePointRaw[treeId][nP].intensity));

                Eigen::Quaterniond q_prev{mvpIScans[iscanId]->v_R_mapping_ref[scanId]};
                Eigen::Quaterniond q_cur_end{mvpIScans[iscanId]->v_R_mapping_ref[scanId + 1]};
                Eigen::Quaterniond q_cur = q_prev.slerp(t_ratio, q_cur_end);
                Eigen::Vector3d t_cur = (1.0 - t_ratio) * mvpIScans[iscanId]->v_r_mapping_ref[scanId] + t_ratio * mvpIScans[iscanId]->v_r_mapping_ref[scanId + 1];
                // Eigen::Vector3d rIm = q_cur * curr_point + t_cur;
                Eigen::Vector3d rIm = mvpIScans[iscanId]->v_R_mapping_ref[scanId] * curr_point + mvpIScans[iscanId]->v_r_mapping_ref[scanId];

                Eigen::Vector3d p_;
                if (exportOption == MAPPING_FRAME)
                    p_ = rIm;
                else if (exportOption == GLOBAL_FRAME)
                    p_ = r_m_global + R_m_global * rIm;
                fPointFile << mapTreeId << "\t" << p_.transpose() << "\t" << float(iscanId) << endl;
            }
        }
    }
    fPointFile.close();
}

/***************************************************************
 * Export ground map
 * Input: mpMapGroundPoint
 * Option:
 *  Mapping frame
 *  Global frame (only valid if trajectory is available)
 **************************************************************/
void Mapping::exportGroundMap(const std::string outPass, eExportOption exportOption)
{
    std::ofstream fPointFile(outPass, std::ifstream::out);
    fPointFile << fixed << std::setprecision(8);

    for (int nP = 0; nP < mpMapGroundPoint->size(); nP++)
    {
        Eigen::Vector3d curr_point(mpMapGroundPoint->points[nP].x,
                                   mpMapGroundPoint->points[nP].y,
                                   mpMapGroundPoint->points[nP].z);
        Eigen::Vector3d p_;
        if (exportOption == MAPPING_FRAME)
            p_ = curr_point;
        else if (exportOption == GLOBAL_FRAME)
            p_ = r_m_global + R_m_global * curr_point;
        fPointFile << p_.transpose() << "\t" << endl;
    }
    fPointFile.close();
}

/**
 * @brief Compares estimated and reference global trajectories and exports both.
 * 
 * @param outPass Output file path.
 * 
 * Useful for debugging and evaluating alignment between SLAM-estimated
 * trajectory and reference GNSS/IMU trajectory.
 */
void Mapping::exportTrajGlobalCheck(const std::string outPass)
{
    std::ofstream fTrajGlobal(outPass, std::ifstream::out);
    fTrajGlobal << fixed << std::setprecision(6);
    fTrajGlobal << "Flag\t" << "X\t"<< "Y\t"<< "Z\t" << "Scan_Index\t"<< "ome\t" <<"phi\t"<< "kap\t" << "Num_Tree"<<endl;

    for (int nIscan = 0; nIscan < mvpIScans.size(); nIscan++)
    {
        // index at t_ini(0), i.e., t_end(last) from previous
        int startIndex = nIscan == 0 ? mvpIScans[0]->indScans[0].scanID -1 : mvpIScans[nIscan - 1]->indScans.back().scanID;
        int num_tree_start = nIscan == 0 ? 0 : mvpIScans[nIscan - 1]->indScans.back().num_tree;
        for (int nPose = 0; nPose < mvpIScans[nIscan]->v_R_mapping.size(); nPose++)
        {
            int index = nPose == 0 ? startIndex : mvpIScans[nIscan]->indScans[nPose - 1].scanID;
            int num_tree = nPose == 0 ? num_tree_start : mvpIScans[nIscan]->indScans[nPose - 1].num_tree;
            Eigen::Vector3d r_lu_g, r_lu_g_ref;
            Eigen::Matrix3d R_lu_g, R_lu_g_ref;

            r_lu_g = r_m_global + R_m_global * mvpIScans[nIscan]->v_r_mapping[nPose];
            R_lu_g = R_m_global * mvpIScans[nIscan]->v_R_mapping[nPose];

            r_lu_g_ref = r_m_global_ref_traj + R_m_global_ref_traj * mvpIScans[nIscan]->v_r_mapping_ref[nPose];
            R_lu_g_ref = R_m_global_ref_traj * mvpIScans[nIscan]->v_R_mapping_ref[nPose];

            fTrajGlobal << "0\t" << r_lu_g.transpose() << "\t"
                         << index << "\t" << rad2deg(Find_Rotation(R_lu_g)).transpose() << "\t" << num_tree << endl;

            fTrajGlobal << "1\t" << r_lu_g_ref.transpose() << "\t"
                         << index << "\t" << rad2deg(Find_Rotation(R_lu_g_ref)).transpose() <<"\t" << num_tree << endl;

            // fTrajGlobal << mnCount << "\t" << mvpIScans[nIscan]->v_r_global_ref[nPose].transpose() << "\t"
            //             << "2\t" << index << "\t" << rad2deg(Find_Rotation(mvpIScans[nIscan]->v_R_global_ref[nPose])).transpose() << endl;
        }
    }
    fTrajGlobal.close();
}
/**
 * @brief Exports SLAM and GNSS trajectories to file and optionally accumulates into point clouds.
 * 
 * @param outPass Output file path.
 * @param slam_trajectory Optional pointer to SLAM trajectory point cloud.
 * @param gnss_trajectory Optional pointer to GNSS trajectory point cloud.
 * 
 * This function also applies offset alignment and stores the transformation matrix
 * for further correction steps.
 */

void Mapping::exportTrajGlobal(const std::string outPass,pcl::PointCloud<pcl::PointXYZ>* slam_trajectory,pcl::PointCloud<pcl::PointXYZ>* gnss_trajectory)
{
    std::ofstream fTrajGlobal(outPass, std::ifstream::out);
    fTrajGlobal << fixed << std::setprecision(6);
    fTrajGlobal << "Flag\t" << "X\t"<< "Y\t"<< "Z\t" << "Scan_Index\t"<< "ome\t" <<"phi\t"<< "kap"<<endl;

	pcl::PointXYZ point;
	Eigen::Vector3d offset;
	int count = 0;
    Eigen::Vector3d r_lu_g_prev;
    for (int nIscan = 0; nIscan < mvpIScans.size(); nIscan++)
    {
        // index at t_ini(0), i.e., t_end(last) from previous
        int startIndex = nIscan == 0 ? mvpIScans[0]->indScans[0].scanID -1 : mvpIScans[nIscan - 1]->indScans.back().scanID;
        int start_pose = nIscan == 0 ? 0:1;
        for (int nPose = start_pose; nPose < mvpIScans[nIscan]->v_R_mapping.size(); nPose++)
        {
            int index = nPose == 0 ? startIndex : mvpIScans[nIscan]->indScans[nPose - 1].scanID;
            Eigen::Vector3d r_lu_g, r_lu_g_ref;
            Eigen::Matrix3d R_lu_g, R_lu_g_ref;

            r_lu_g = r_m_global + R_m_global * mvpIScans[nIscan]->v_r_mapping[nPose];
			if (slam_trajectory)
			{
				if (count == 0)
				{
					offset = r_lu_g;
					offset_trans_traj_ = offset;
				}
				point.getVector3fMap() = (r_lu_g - offset).cast<float>();
				slam_trajectory->push_back(point);
			}

            R_lu_g = R_m_global * mvpIScans[nIscan]->v_R_mapping[nPose];

            r_lu_g_ref = r_m_global_ref_traj + R_m_global_ref_traj * mvpIScans[nIscan]->v_r_mapping_ref[nPose];
            R_lu_g_ref = R_m_global_ref_traj * mvpIScans[nIscan]->v_R_mapping_ref[nPose];
			if (gnss_trajectory)
			{
				point.getVector3fMap() = (r_lu_g_ref - offset).cast<float>();
				gnss_trajectory->push_back(point);
			}

            fTrajGlobal << "0\t" << r_lu_g.transpose() << "\t"
                         << index << "\t" << rad2deg(Find_Rotation(R_lu_g)).transpose() << endl;

            fTrajGlobal << "1\t" << r_lu_g_ref.transpose() << "\t"
                         << index << "\t" << rad2deg(Find_Rotation(R_lu_g_ref)).transpose() << endl;

			count++;

            if(count>1)
            {
                double successive_scan_dis = (r_lu_g_prev - r_lu_g).norm();
                if(successive_scan_dis>0.2)
                {
                    fMapLog<< " large distance in "<<mvpIScans[nIscan]->index<<",distance: "<<successive_scan_dis<<std::endl;
                }
            }
            r_lu_g_prev = r_lu_g;

        }
    }
    fTrajGlobal.close();


}

/**
 * @brief Apply ICP between SLAM and GNSS trajectories and update the global transformation.
 *
 * @param slam_trajectory SLAM trajectory point cloud
 * @param gnss_trajectory GNSS trajectory point cloud
 */
void Mapping::UpdateMappingToGlobalTransformation(const pcl::PointCloud<pcl::PointXYZ>& slam_trajectory, const pcl::PointCloud<pcl::PointXYZ>& gnss_trajectory)
{
	Eigen::Matrix4f icp_result = Trajectory_ICP(slam_trajectory, gnss_trajectory);
	fMapLog << " icp_result to transform traj: " << std::endl;
	fMapLog << icp_result << endl;
	Eigen::Matrix3d R_icp = icp_result.block<3, 3>(0, 0).cast<double>();
	Eigen::Vector3d r_icp = icp_result.block<3, 1>(0, 3).cast<double>();
	fMapLog << "R_icp: " << R_icp << endl;
	fMapLog << "r_icp: " << r_icp << endl;
	R_m_global = R_icp*R_m_global;
	r_m_global = R_icp*r_m_global - R_icp*offset_trans_traj_ + r_icp+ offset_trans_traj_;
	fMapLog << "after transform traj" << endl;
	fMapLog << "R_m_global: " << R_m_global << endl;
	fMapLog << "r_m_global: " << r_m_global << endl;
}




/***************************************************************
 * Compute the average of R_lu_lup from all mvTempScan
 * Q is an Mx4 matrix of quaternions. Qavg is the average quaternion Based on
 * Markley, F. Landis, Yang Cheng, John Lucas Crassidis, and Yaakov Oshman.
 * "Averaging quaternions." Journal of Guidance, Control, and Dynamics 30,
 * no. 4 (2007): 1193-1197.
 * https://www.mathworks.com/matlabcentral/fileexchange/40098-tolgabirdal-averaging_quaternions
 **************************************************************/
/**
 * @brief Compute and return the average orientation across scans using quaternion averaging.
 *
 * @return Average orientation as a quaternion
 */
Eigen::Quaterniond Mapping::computeAverageLevelR()
{
    // currently, equal weights
    double w = 1.0;
    double wSum = 0.0;
    Eigen::Matrix4d A = Eigen::Matrix4d::Zero();
    for (int nScan = 0; nScan < mvTempScan.size(); nScan++)
    {
        Eigen::Quaterniond q{mvTempScan[nScan].R_lu_lup};
        Eigen::Vector4d q_vec = q.coeffs();
        if (q_vec(0) < 0)
            q_vec *= -1;
        A = A + w * q_vec * q_vec.transpose();
        wSum += w;
        //  fMapLog << rad2deg(Find_Rotation(q.toRotationMatrix())).transpose() <<endl;
    }
    A = (1.0 / wSum) * A;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> saes(A);
    Eigen::Vector4d q_vec_avg = saes.eigenvectors().col(3);
    Eigen::Quaterniond q_avg{q_vec_avg(3), q_vec_avg(0), q_vec_avg(1), q_vec_avg(2)};
    //    fMapLog << "Average: " << rad2deg(Find_Rotation(q_avg.toRotationMatrix())).transpose() <<endl;

    return q_avg;
}

/***************************************************************
 * Compute the ref pose for each individual scan in the ISCAN
 **************************************************************/
/**
 * @brief Compute the reference trajectory for the current scan.
 *
 * Interpolates the GNSS trajectory (if available) to derive poses for each scan.
 */
void Mapping::computeRefTraj()
{
    // int startIndex = mbInit ? 0 : 1;

    // first pose from the previous
    if (!mbInit)
    {
        pCurrentIScan->v_R_mapping_ref.push_back(pCurrentIScan->v_R_mapping[0]);
        pCurrentIScan->v_r_mapping_ref.push_back(pCurrentIScan->v_r_mapping[0]);
    }
    else
    {
        //since we add current scan first
        //pCurrentIScan->v_R_mapping_ref.push_back(mvpIScans.back()->v_R_mapping_ref.back());
        //pCurrentIScan->v_r_mapping_ref.push_back(mvpIScans.back()->v_r_mapping_ref.back());
        pCurrentIScan->v_R_mapping_ref.push_back(mvpIScans[mvpIScans.size()-2]->v_R_mapping_ref.back());
        pCurrentIScan->v_r_mapping_ref.push_back(mvpIScans[mvpIScans.size()-2]->v_r_mapping_ref.back());
    }

    //
    Eigen::Vector3d r_lu_b_ref = mpTraj->r_lu_b;
    Eigen::Matrix3d R_lu_b_ref = mpTraj->R_lu_b;
    for (int nScan = 0; nScan < mNumIntegrationScan; nScan++)
    {

        Eigen::Vector3d r_b_global_t1, r_b_global_t2, r_lu_global_t1, r_lu_global_t2;
        Eigen::Matrix3d R_b_global_t1, R_b_global_t2, R_lu_global_t1, R_lu_global_t2;

        int idx;
        if (!(mpTraj->bopInterpolation((pCurrentIScan->indScans[nScan].mTimeEnd - pCurrentIScan->indScans[nScan].mTimeTracking) / 1000.0, -1, idx, r_b_global_t1, R_b_global_t1) && mpTraj->bopInterpolation(pCurrentIScan->indScans[nScan].mTimeEnd / 1000.0, -1, idx, r_b_global_t2, R_b_global_t2)))
            throw std::runtime_error("Wrong Time");

        r_lu_global_t1 = r_b_global_t1 + R_b_global_t1 * r_lu_b_ref;
        R_lu_global_t1 = R_b_global_t1 * R_lu_b_ref;

        r_lu_global_t2 = r_b_global_t2 + R_b_global_t2 * r_lu_b_ref;
        R_lu_global_t2 = R_b_global_t2 * R_lu_b_ref;

        Eigen::Matrix3d R_lut2_lut1_ref = R_lu_global_t1.inverse() * R_lu_global_t2;
        Eigen::Vector3d r_lut2_lut1_ref = R_lu_global_t1.inverse() * (r_lu_global_t2 - r_lu_global_t1);

        // fMapLog << mvTempScan[nScan].scanID << "\t" << mvTempScan[nScan].r_lup_m_ini.transpose() << "\t"<< rad2deg(Find_Rotation(mvTempScan[nScan].R_lup_m_ini)).transpose() <<endl;
        // fMapLog << mvTempScan[nScan].scanID << "\t" << r_lupt2_lutp1_ref.transpose() << "\t"<< rad2deg(Find_Rotation(R_lupt2_lutp1_ref)).transpose() <<endl;
        // current scan to the first scan from odometry thread
        Eigen::Matrix3d R_lu_m_t2_ref = pCurrentIScan->v_R_mapping_ref.back() * R_lut2_lut1_ref;
        Eigen::Vector3d r_lu_m_t2_ref = pCurrentIScan->v_r_mapping_ref.back() + pCurrentIScan->v_R_mapping_ref.back() * r_lut2_lut1_ref;

        pCurrentIScan->v_R_mapping_ref.push_back(R_lu_m_t2_ref);
        pCurrentIScan->v_r_mapping_ref.push_back(r_lu_m_t2_ref);

        //check only
        // if (nScan == 0)
        // {
        //     pCurrentIScan->v_R_global_ref.push_back(R_lu_global_t1);
        //     pCurrentIScan->v_r_global_ref.push_back(r_lu_global_t1);
        // }

        // pCurrentIScan->v_R_global_ref.push_back(R_lu_global_t2);
        // pCurrentIScan->v_r_global_ref.push_back(r_lu_global_t2);
    }
}

/***************************************************************
export estimated trajectory for the integrated scans
flag:
init trajectory from odometry
refined trajectory after integration
reference trajectory
***************************************************************/
/*void Mapping::exportTrajectory(const std::string outPass)
{
    std::ofstream fTrajectoryFile(outPass, std::ifstream::out);
    fTrajectoryFile << fixed << std::setprecision(8);
    // fTrajectoryFile << numScan <<" " << mvTempScan.size() <<endl;

    for (int nScan = 0; nScan < mNumIntegrationScan; nScan++)
    {
        fTrajectoryFile << nScan <<"\t" << mvTempScan[nScan].r_lu_local_ini.transpose() << "\t" << 0 << "\t" << nScan << "\t" << rad2deg(Find_Rotation(mvTempScan[nScan].R_lu_local_ini)).transpose() << endl;
        fTrajectoryFile << nScan <<"\t" <<mvTempScan[nScan].r_lu_local_updated.transpose() << "\t" << 1 << "\t" << nScan << "\t" << rad2deg(Find_Rotation(mvTempScan[nScan].R_lu_local_updated)).transpose() << endl;
        if (mpTraj)
        {
            fTrajectoryFile << nScan <<"\t" << mvTempScan[nScan].r_lu_local_ref.transpose() << "\t" << 2 << "\t" << nScan << "\t" << rad2deg(Find_Rotation(mvTempScan[nScan].R_lu_local_ref)).transpose() << endl;
        }
    }
    fTrajectoryFile.close();
}*/

void Mapping::exportIntegratedTrajectory(const std::string outPass)
{
    std::ofstream fTrajectoryFile(outPass, std::ifstream::out);
    fTrajectoryFile << fixed << std::setprecision(8);
    // fTrajectoryFile << numScan <<" " << mvTempScan.size() <<endl;

    for (int nEpoch = 0; nEpoch < mNumIntegrationScan+1; nEpoch++)
    {
        fTrajectoryFile << 0 <<"\t" << mv_r_lu_local_ini[nEpoch].transpose() << "\t" << nEpoch  << "\t" << rad2deg(Find_Rotation(mv_R_lu_local_ini[nEpoch])).transpose() << endl;
        fTrajectoryFile << 1 <<"\t" <<mv_r_lu_local[nEpoch].transpose()<< "\t" << nEpoch << "\t" << rad2deg(Find_Rotation(mv_R_lu_local[nEpoch])).transpose() << endl;
        if (mpTraj)
        {
        fTrajectoryFile << 2 <<"\t" <<mv_r_lu_local_ref[nEpoch].transpose()<< "\t" << nEpoch << "\t" << rad2deg(Find_Rotation(mv_R_lu_local_ref[nEpoch])).transpose() << endl;
        }
    }
    fTrajectoryFile.close();
}


// void Mapping::exportIntegratedTree(const std::string outPass)
// {
//     std::ofstream fPointFile(outPass, std::ifstream::out);
//     fPointFile << fixed << std::setprecision(8);
//     for (int nTree = 0; nTree < mvIntegratedTree.size(); nTree++)
//     {
//         for (int nScan = 0; nScan < mvIntegratedTreeIds[nTree].size(); nScan++)
//         {
//             int scanId = mvIntegratedTreeIds[nTree][nScan].first;
//             int treeId = mvIntegratedTreeIds[nTree][nScan].second;

//             vector<PointType> *pPointList = &(mvTempScan[scanId].vTreeRawPoints[treeId]);
//             int numPoint = pPointList->size();
//             for (int nP = 0; nP < numPoint; nP++)
//             {
//                 Eigen::Vector3f rilu(pPointList->at(nP).x, pPointList->at(nP).y, pPointList->at(nP).z);
//                 Eigen::Vector3d rib2;
//                 // rib2 = mvTempScan[scanId].R_lup_m_updated * (mvTempScan[scanId].R_lu_lup * rilu.cast<double>()) + mvTempScan[scanId].r_lup_m_updated;
//                 rib2 = mvTempScan[scanId].R_lu_local_ini * rilu.cast<double>() + mvTempScan[scanId].r_lu_local_ini;

//                 if (rib2(2) < 10.0)
//                 {
//                     fPointFile << rib2(0) << "\t" << rib2(1) << "\t" << rib2(2) << "\t" << mvIntegratedTreeIds[nTree].size() << "\t" << nTree << "\t" << float(scanId) << "\t" << endl;
//                 }
//             }
//         }
//     }
//     for (int nScan = 0; nScan < mvTempScan.size(); nScan++)
//     {
//         vector<PointType> *pPointList = &(mvTempScan[nScan].vGroundRawPoints);

//         int numPoint = pPointList->size();
//         int scanId = nScan;
//         for (int nP = 0; nP < numPoint; nP++)
//         {
//             Eigen::Vector3f rilu(pPointList->at(nP).x, pPointList->at(nP).y, pPointList->at(nP).z);
//             Eigen::Vector3d rib2;
//             //            rib2 = mvTempScan[scanId].R_lup_m_updated * (mvTempScan[scanId].R_lu_lup * rilu.cast<double>()) + mvTempScan[scanId].r_lup_m_updated;
//             rib2 = mvTempScan[scanId].R_lu_local_ini * rilu.cast<double>() + mvTempScan[scanId].r_lu_local_ini;

//             fPointFile << rib2(0) << "\t" << rib2(1) << "\t" << rib2(2) << "\t" << 0 << "\t" << 0 << "\t" << float(scanId) << "\t" << endl;
//         }
//     }
//     fPointFile.close();
// }

// void Mapping::exportFinalTree(const std::string outPass)
// {
//     std::ofstream fPointFile(outPass, std::ifstream::out);
//     fPointFile << fixed << std::setprecision(8);
//     for (int nTree = 0; nTree < mvTreePoints.size(); nTree++)
//     {
//         for (int nP = 0; nP < mvTreePoints[nTree].size(); nP++)
//         {
//             Eigen::Vector3f rilup(mvTreePoints[nTree][nP].x, mvTreePoints[nTree][nP].y, mvTreePoints[nTree][nP].z);
//             int scanId = int(mvTreePoints[nTree][nP].intensity / 100);
//             Eigen::Vector3d rib;
//             rib = mvTempScan[scanId].R_lu_local_ini * rilup.cast<double>() + mvTempScan[scanId].r_lu_local_ini;

//             fPointFile << rib(0) << "\t" << rib(1) << "\t" << rib(2) << "\t" << nTree << "\t" << mvValidTrees[nTree] << "\t" << float(scanId) << "\t" << endl;
//         }
//     }
//     fPointFile.close();
// }

// void Mapping::exportIntegratedTreeCenter(const std::string outPass)
// {
//     std::ofstream fPointFile(outPass, std::ifstream::out);
//     fPointFile << fixed << std::setprecision(8);
//     for (int nTree = 0; nTree < mvIntegratedTree.size(); nTree++)
//     {
//         for (int nScan = 0; nScan < mvIntegratedTreeIds[nTree].size(); nScan++)
//         {
//             int scanId = mvIntegratedTreeIds[nTree][nScan].first;
//             int treeId = mvIntegratedTreeIds[nTree][nScan].second;
//             Eigen::Vector3d rib2;
//             Eigen::Vector3d rr = mvTempScan[scanId].vTreeLoc[treeId];
//             rib2 = mvTempScan[scanId].R_lu_local_ini * (mvTempScan[scanId].vTreeLoc[treeId]) + mvTempScan[scanId].r_lu_local_ini;
//             fPointFile << 1 << "\t" << rib2(0) << "\t" << rib2(1) << "\t" << rib2(2) << "\t" << mvIntegratedTreeIds[nTree].size() << "\t" << nTree << "\t" << float(scanId) << "\t" << endl;
//             fPointFile << 2 << "\t" << rr(0) << "\t" << rr(1) << "\t" << rr(2) << "\t" << mvIntegratedTreeIds[nTree].size() << "\t" << nTree << "\t" << float(scanId) << "\t" << endl;
//         }
//         fPointFile << 3 << "\t" << mvTreeCentersLocal[nTree](0) << "\t" << mvTreeCentersLocal[nTree](1) << "\t" << mvTreeCentersLocal[nTree](2) << "\t" << mvValidTrees[nTree] << "\t" << mvIntegratedTreeIds[nTree].size() << "\t" << nTree << "\t" << endl;
//     }

//     fPointFile.close();
// }
// void Mapping::exportTree(bool validTree)
// {
//     int treeCount = 0;
//     for (std::size_t nTree = 0; nTree < mVecTrees.size(); nTree++)
//     {
//         if (validTree)
//         {
//             if (mVecTrees[nTree].numCluster <= 3)
//                 continue;
//         }
//         treeCount++;
//         for (std::size_t nCluster = 0; nCluster < mVecTrees[nTree].numCluster; nCluster++)
//         {

//             int scanID = mVecTrees[nTree].vecScanId[nCluster];
//             for (std::size_t nP = 0; nP < mVecTrees[nTree].mvPoints[nCluster].size(); nP++)
//             {
//                 Eigen::Vector3f rilu(mVecTrees[nTree].mvPoints[nCluster][nP].x, mVecTrees[nTree].mvPoints[nCluster][nP].y, mVecTrees[nTree].mvPoints[nCluster][nP].z);
//                 Eigen::Vector3d rib = mVecScan[scanID].r_lu_m + mVecScan[scanID].R_lu_m * rilu.cast<double>();
//                 fMappedPoints << nTree << "\t" << rib(0) << "\t" << rib(1) << "\t" << rib(2) << "\t" << scanID << "\t" << mVecScan[scanID].partID << endl;
//             }
//         }
//     }
//     cout << "Numer of tree " << treeCount << endl;
// }

/**
 * @brief Exports integrated tree and ground points using reference trajectory.
 * 
 * @param outPass Output file path to save exported points.
 * 
 * This function uses the reference trajectory from GNSS/IMU to convert all local scan
 * points to the mapping frame and writes them into a file with metadata indicating whether 
 * the point is a tree or ground point.
 */
void Mapping::exportIntegratedTreeRef(const std::string outPass)
{
    fMapLog << outPass << endl;

    Eigen::Vector3d r_b_m_, r_lu_m_;
    Eigen::Matrix3d R_b_m_, R_lu_m_;
    int startIndex = -1;
    double startTime = mvTempScan[0].mTimeInit;
    if (!mpTraj->bopInterpolation(startTime / 1000.0, -1, startIndex, r_b_m_, R_b_m_))
        return;

    std::ofstream fPointFile(outPass, std::ifstream::out);
    fPointFile << fixed << std::setprecision(8);
    for (int nTree = 0; nTree < mvIntegratedTree.size(); nTree++)
    {
        for (int nScan = 0; nScan < mvIntegratedTreeIds[nTree].size(); nScan++)
        {
            int scanId = mvIntegratedTreeIds[nTree][nScan].first;
            int treeId = mvIntegratedTreeIds[nTree][nScan].second;

            vector<PointType> *pPointList = &(mvTempScan[scanId].vTreeRawPoints[treeId]);
            int numPoint = pPointList->size();
            for (int nP = 0; nP < numPoint; nP++)
            {
                Eigen::Vector3f rilup(pPointList->at(nP).x, pPointList->at(nP).y, pPointList->at(nP).z);
                Eigen::Vector3d rilu = mvTempScan[scanId].R_lu_lup.inverse() * rilup.cast<double>();

                // time = t_end_prev + ratio * time_Tracking
                double time = (mvTempScan[scanId].mTimeEnd - mvTempScan[scanId].mTimeTracking) / 1000.0 + mvTempScan[scanId].mTimeTracking / 1000.0 * double(pPointList->at(nP).intensity - int(pPointList->at(nP).intensity));
                Eigen::Vector3d r_b_m_, r_lu_m_;
                Eigen::Matrix3d R_b_m_, R_lu_m_;
                int idx;
                if (mpTraj->bopInterpolation(time, startIndex, idx, r_b_m_, R_b_m_))
                {
                    r_lu_m_ = r_b_m_ + R_b_m_ * mpTraj->r_lu_b;
                    R_lu_m_ = R_b_m_ * mpTraj->R_lu_b;
                    Eigen::Vector3d rib2;
                    rib2 = R_lu_m_ * rilu.cast<double>() + r_lu_m_;
                    fPointFile << 0 << "\t" << rib2(0) << "\t" << rib2(1) << "\t" << rib2(2) << "\t" << mvIntegratedTreeIds[nTree].size() << "\t" << nTree << "\t" << float(scanId) << "\t" << endl;
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

    for (int nScan = 0; nScan < mNumIntegrationScan; nScan++)
    {
        vector<PointType> *pPointList = &(mvTempScan[nScan].vGroundRawPoints);

        int numPoint = pPointList->size();
        int scanId = nScan;
        for (int nP = 0; nP < numPoint; nP++)
        {
            Eigen::Vector3f rilup(pPointList->at(nP).x, pPointList->at(nP).y, pPointList->at(nP).z);
            Eigen::Vector3d rilu = mvTempScan[scanId].R_lu_lup.inverse() * rilup.cast<double>();
            double time = (mvTempScan[scanId].mTimeEnd - mvTempScan[scanId].mTimeTracking) / 1000.0 + mvTempScan[scanId].mTimeTracking / 1000.0 * double(pPointList->at(nP).intensity - int(pPointList->at(nP).intensity));
            Eigen::Vector3d r_b_m_, r_lu_m_;
            Eigen::Matrix3d R_b_m_, R_lu_m_;
            int idx;
            if (mpTraj->bopInterpolation(time, startIndex, idx, r_b_m_, R_b_m_))
            {
                r_lu_m_ = r_b_m_ + R_b_m_ * mpTraj->r_lu_b;
                R_lu_m_ = R_b_m_ * mpTraj->R_lu_b;
                Eigen::Vector3d rib2;
                rib2 = R_lu_m_ * rilu.cast<double>() + r_lu_m_;
                fPointFile << 1 << "\t" << rib2(0) << "\t" << rib2(1) << "\t" << rib2(2) << "\t" << 0 << "\t" << 0 << "\t" << float(scanId) << "\t" << endl;
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

/**
 * @brief Exports downsampled ground points for each grid to a file.
 * 
 * @param outPass Output path to save exported ground points.
 * 
 * Ground points are selected using grid-based filtering. Each line of the output
 * contains point coordinates and the grid index to which it belongs.
 */
void Mapping::exportGroundPoints(const std::string outPass)
{
    std::ofstream fGroundPoints(outPass, std::ifstream::out);
    fGroundPoints << fixed << std::setprecision(4);

    // for (int nP = 0; nP <mvGroundPointsLocal.size(); nP++)
    // {
    //         fGroundPoints << -1 << "\t" << mvGroundPointsLocal[nP].x << "\t" << mvGroundPointsLocal[nP].y
    //         << "\t" << mvGroundPointsLocal[nP].z << "\t" << mvGroundPointsLocal[nP].intensity<< endl;
    // }
    for (int nGrid = 0; nGrid < mvGridPointIndex.size(); nGrid++)
    {
        for (int nP = 0; nP < mvGridPointIndex[nGrid].size(); nP++)
        {
            int index = mvGridPointIndex[nGrid][nP];
            fGroundPoints << nGrid << "\t" << mvGroundPointsLocal[index].x << "\t" << mvGroundPointsLocal[index].y
                          << "\t" << mvGroundPointsLocal[index].z << "\t" << mvGroundPointsLocal[index].intensity << endl;
        }
    }
    /*for (int nScan = 0; nScan < mvTempScan.size(); nScan++)
    {
        int scanId = nScan;
        for (std::size_t nP = 0; nP < mvTempScan[nScan].pGroundPoints->size(); nP++)
        {
            Eigen::Vector3f rilu(mvTempScan[nScan].pGroundPoints->points[nP].x, mvTempScan[nScan].pGroundPoints->points[nP].y, mvTempScan[nScan].pGroundPoints->points[nP].z);
             Eigen::Vector3d rib2;
            //            rib2 = mvTempScan[scanId].R_lup_m_updated * (mvTempScan[scanId].R_lu_lup * rilu.cast<double>()) + mvTempScan[scanId].r_lup_m_updated;
            rib2 = mvTempScan[scanId].R_lup_m_ini * rilu.cast<double>() + mvTempScan[scanId].r_lup_m_ini;

            fGroundPoints << rib2(0) << "\t" << rib2(1) << "\t" << rib2(2) << "\t" << 0 << "\t" << 0 << "\t" << float(scanId) << "\t" << endl;

        }
    }
    */
    fGroundPoints.close();
}

// /*deprecated*/
// void Mapping::exportIntegratedTreeIni(const std::string outPass)
// {
//     fMapLog << outPass << endl;

//     std::ofstream fPointFile(outPass, std::ifstream::out);
//     fPointFile << fixed << std::setprecision(8);
//     for (int nTree = 0; nTree < mvIntegratedTree.size(); nTree++)
//     {
//         for (int nScan = 0; nScan < mvIntegratedTreeIds[nTree].size(); nScan++)
//         {
//             int scanId = mvIntegratedTreeIds[nTree][nScan].first;
//             int treeId = mvIntegratedTreeIds[nTree][nScan].second;

//             vector<PointType> *pPointList = &(mvTempScan[scanId].vTreeRawPoints[treeId]);
//             int numPoint = pPointList->size();
//             for (int nP = 0; nP < numPoint; nP++)
//             {
//                 Eigen::Vector3f rilu(pPointList->at(nP).x, pPointList->at(nP).y, pPointList->at(nP).z);
//                 Eigen::Vector3d rib2;
//                 rib2 = mvTempScan[scanId].R_lu_m * rilu.cast<double>() + mvTempScan[scanId].r_lu_m;
//                 fPointFile << rib2(0) << "\t" << rib2(1) << "\t" << rib2(2) << "\t" << mvIntegratedTreeIds[nTree].size() << "\t" << nTree << "\t" << float(scanId) << "\t" << endl;
//             }
//         }
//     }

//     fPointFile.close();
// }

/**
 * @brief Count different types of trees in the current map (e.g., ESTABLISHED, SOLID).
 */
void Mapping::CountMapTreeType()
{
    unordered_map<int, int> tree_type;
    for (int nT = 0; nT < mpMapTree.size(); nT++)
    {
        tree_type[mpMapTree[nT]->status]++;
    }

    fMapLog << "\tGlobal - " << tree_type[-1] << "  Init - " << tree_type[0] << "  TBD - " << tree_type[1] << "  Fitted - " << tree_type[2]
            << "  Established - " << tree_type[3] << "  Solid - " << tree_type[4] << "  Deactivated - " << tree_type[-2] << endl;
}

/**
 * @brief Compute surface coverage ratio of a tree by binning angles.
 *
 * @param paraVer Whether to use current tree or candidate tree
 */
void MapTree::compute_surface_ratio(eParaVer paraVer)
{

    vector<int> countBin(Num_Bin, 0);
    for (int nT = 0; nT < visibleInfo.size(); nT++)
    {
        int iscanId = visibleInfo[nT].first;
        int treeId = visibleInfo[nT].second;
        for (int nP = 0; nP < pMap->mvpIScans[iscanId]->vTreePointMapping[treeId].size(); nP++)
        {
            // points in mapping
            Eigen::Vector3d curr_point(pMap->mvpIScans[iscanId]->vTreePointMapping[treeId][nP].x,
                                       pMap->mvpIScans[iscanId]->vTreePointMapping[treeId][nP].y,
                                       pMap->mvpIScans[iscanId]->vTreePointMapping[treeId][nP].z);
            double kap = compute_angle(curr_point,paraVer);
            countBin[int((kap - -PI) / (2 * PI / double(Num_Bin)))]++;
        }
    }

    if (paraVer == CANDIDATE)
    {
        int iscanId = candTreeInfo.first;
        int treeId = candTreeInfo.second;
        for (int nP = 0; nP < pMap->mvpIScans[iscanId]->vTreePointMapping[treeId].size(); nP++)
        {
            // points in mapping
            Eigen::Vector3d curr_point(pMap->mvpIScans[iscanId]->vTreePointMapping[treeId][nP].x,
                                       pMap->mvpIScans[iscanId]->vTreePointMapping[treeId][nP].y,
                                       pMap->mvpIScans[iscanId]->vTreePointMapping[treeId][nP].z);
            double kap = compute_angle(curr_point, paraVer);
            countBin[int((kap - -PI) / (2 * PI / double(Num_Bin)))]++;
        }
    }
	if (paraVer == FIT_CYLINDER)
	{
		int iscanId = candTreeInfo.first;
		int treeId = candTreeInfo.second;
		for (int nP = 0; nP < pMap->mvpIScans[iscanId]->vTreePointMapping[treeId].size(); nP++)
		{
			// points in mapping
			Eigen::Vector3d curr_point(pMap->mvpIScans[iscanId]->vTreePointMapping[treeId][nP].x,
				pMap->mvpIScans[iscanId]->vTreePointMapping[treeId][nP].y,
				pMap->mvpIScans[iscanId]->vTreePointMapping[treeId][nP].z);
			double kap = compute_angle(curr_point, paraVer);
			int index = static_cast<int>((kap - -PI) / (2 * PI / double(Num_Bin)));
			countBin[index]++;
		}
	}

    int countValidBin = 0;
    for (int nBin = 0; nBin < Num_Bin; nBin++)
    {
        if (countBin[nBin] > Min_Number_Pts_Per_Bin)
            countValidBin++;
    }
    surfaceRatio = double(countValidBin) / double(Num_Bin);
}
/**
 * @brief Count the number of points currently assigned to this map tree.
 */
void MapTree::count_point()
{
    numPoint = 0;
    for (int nT = 0; nT < visibleInfo.size(); nT++)
    {
        int iscanId = visibleInfo[nT].first;
        int treeId = visibleInfo[nT].second;
        numPoint += pMap->mvpIScans[iscanId]->vTreePointMapping[treeId].size();
    }
}
/**
 * @brief Compute and update the centroid of the tree from all its points.
 */
void MapTree::ComputeCentroid()
{
    Eigen::Vector3d centeroid = Eigen::Vector3d::Zero();

    for (int nT = 0; nT < visibleInfo.size(); nT++)
    {
        int iscanId = visibleInfo[nT].first;
        int treeId = visibleInfo[nT].second;
        for (int nP = 0; nP < pMap->mvpIScans[iscanId]->vTreePointMapping[treeId].size(); nP++)
        {
            // points in mapping
            Eigen::Vector3d curr_point(pMap->mvpIScans[iscanId]->vTreePointMapping[treeId][nP].x,
                                       pMap->mvpIScans[iscanId]->vTreePointMapping[treeId][nP].y,
                                       pMap->mvpIScans[iscanId]->vTreePointMapping[treeId][nP].z);

            centeroid += curr_point;
        }
    }
    center = centeroid / double(numPoint);
}

/*****************************************************************
 * Add Iscan tree to current map tree, update status of a tree
 * ***************************************************************/
/**
 * @brief Add a new candidate iscan tree into the current map tree.
 *
 * Updates tree parameters and status if the residual is within threshold.
 *
 * @param updatedPara New candidate cylinder parameters
 * @param tree_info Scan and tree ID pair
 * @return true if the addition is successful
 */
bool MapTree::add_iscan_tree(CylinderPara updatedPara, pair<int, int> tree_info)
{
    // receive the new param, derive R_trans_candidate
    cand_para = updatedPara;
    candTreeInfo = tree_info;
    CylinderPara tree_para = cand_para;
    MapTree::eParaVer para_model = CANDIDATE;
#ifdef FIT_CYLINDER_TREE
    para_model = FIT_CYLINDER;
    if (visibleInfo.size() != 0)
	{
		CalculateAllPointsFitCylinderParameters();
		para = fit_para_;
	}
	std::vector<double> tree_paras = CalculateIncomingPointsFitCylinderParameters();
	if (!fit_tree_failed_)
	{
		if (IsPoorStateTree())
		{
			poor_state_tree_ = true;
			return false;
		}
	}
	else
	{
		return false;
	}
    tree_para = fit_para_;
#endif
	derive_trans_matrix(para_model);

    double res_can = compute_iscan_tree_residual(tree_info, para_model);

    //in case no points from map tree (GLOBAL_INIT), residual will be 0.0
    double res_map = compute_map_tree_residual(para_model);

    //fDebug << Id << " - " << tree_info.first << " - " << tree_info.second << " : cand rmse - " << res_can << "/" << candRmse << ", map rmse - " << res_map << "/" << rmse << endl;

    // Case 1. If this map tree is well established
    if (status == ESTABLISHED)
    {
        //fDebug << "ESTABLISHED - " << flush;
        // compute the rmse of candidate tree to the new parameter
        if (res_can > 0.2)
        {
            //fDebug << "Fail " << endl;
            return false;
        }
        // if sucess
        compute_surface_ratio(para_model);
                //fDebug << "Compute ratio - " << surfaceRatio << " - " << flush;

        update();
        is_established();
        //fDebug << "Status - " << status << " - " << flush;
    }
    // Case 2. If this map tree is a solid cylinder
    else if (status == SOLID)
    {
        //fDebug << "SOLID - " << flush;

        if (tree_para.r > 0.04 || res_can > 0.2)
        {
            //fDebug << "Fail " << endl;
            return false;
        }

        // if sucess
        compute_surface_ratio(para_model);
                //fDebug << "Compute ratio - " << surfaceRatio << " - " << flush;

        update();
        is_established();
        //fDebug << "Status - " << status << " - " << flush;
    }
    // Case 3.
    else
    {
        //fDebug << "Not decided - " << flush;

        if (res_can > 0.4 || res_map > 0.4)
        {
            //fDebug << "Fail " << endl;
            return false;
        }

        //
        compute_surface_ratio(para_model);
        //fDebug << "Compute ratio - " << surfaceRatio << " - " << flush;
        // if radius is large but small ratio, failed matches
        if (surfaceRatio < 0.4 && tree_para.r > 0.4)
        {
            //fDebug << "Fail " << endl;
            return false;
        }
        // if sucess
        update();
        is_established();
        //fDebug << "Status - " << status << " - " << flush;
    }
    
    //fDebug << "Success" << endl;

    return true;
}


/*Conduct fitting for the candidate iscan tree */
std::vector<double> MapTree::CalculateIncomingPointsFitCylinderParameters()
{
	// update parameter
	int iscanId = candTreeInfo.first;
	int treeId = candTreeInfo.second;

	int num_points_cand_tree = pMap->mvpIScans[iscanId]->vTreePointMapping[treeId].size();

	Eigen::Vector3d cur_center = Eigen::Vector3d::Zero();
	std::vector<PointType> tree_points;
	pcl::PointCloud<PointType> out_point_cloud;
	for (int nP = 0; nP < pMap->mvpIScans[iscanId]->vTreePointMapping[treeId].size(); nP++)
	{
		// points in mapping
		Eigen::Vector3d curr_point(pMap->mvpIScans[iscanId]->vTreePointMapping[treeId][nP].x,
			pMap->mvpIScans[iscanId]->vTreePointMapping[treeId][nP].y,
			pMap->mvpIScans[iscanId]->vTreePointMapping[treeId][nP].z);

		cur_center += curr_point;
		PointType tree_point;
		tree_point.x = curr_point(0);
		tree_point.y = curr_point(1);
		tree_point.z = curr_point(2);
		tree_points.push_back(tree_point);
		out_point_cloud.push_back(tree_point);
	}
	std::vector<double> tree_parameters;
	double final_error;

    //if success, update fit_para_
	if (Forestry_SLAM::FittingCylinderModel(tree_points, tree_parameters, final_error) == 1)
	{
		fit_para_.x = Eigen::Vector3d(tree_parameters[0], tree_parameters[1], tree_parameters[2]);
		fit_para_.n = Eigen::Vector3d(tree_parameters[3], tree_parameters[4], tree_parameters[5]);
		fit_para_.r = tree_parameters[6];
		final_fit_error_ = final_error;
		tree_parameters.push_back(final_fit_error_);
	}
	else
	{
		// fit_para_.r = -1;
		fit_tree_failed_ = true;
		tree_parameters.clear();
	}
	// update trans matrix
	derive_trans_matrix(FIT_CYLINDER);
	return tree_parameters;
}

/**
 * @brief Fits a cylinder model to all tree points in the current MapTree.
 *
 * Calculates the best-fit cylinder parameters for the aggregated points across scans. 
 * Updates the `fit_para_` member if successful.
 */

void MapTree::CalculateAllPointsFitCylinderParameters()
{
	// update parameter
	Eigen::Vector3d cur_center = Eigen::Vector3d::Zero();
	std::vector<PointType> tree_points;
	for (int nT = 0; nT < visibleInfo.size(); nT++)
	{
		int iscanId = visibleInfo[nT].first;
		int treeId = visibleInfo[nT].second;
		int num_points_cand_tree = pMap->mvpIScans[iscanId]->vTreePointMapping[treeId].size();
		for (int nP = 0; nP < num_points_cand_tree; nP++)
		{
			// points in mapping
			Eigen::Vector3d curr_point(pMap->mvpIScans[iscanId]->vTreePointMapping[treeId][nP].x,
				pMap->mvpIScans[iscanId]->vTreePointMapping[treeId][nP].y,
				pMap->mvpIScans[iscanId]->vTreePointMapping[treeId][nP].z);

			cur_center += curr_point;
			PointType tree_point;
			tree_point.x = curr_point(0);
			tree_point.y = curr_point(1);
			tree_point.z = curr_point(2);
			tree_points.push_back(tree_point);
		}
	}
	std::vector<double> tree_parameters;
	double final_error;
	if (Forestry_SLAM::FittingCylinderModel(tree_points, tree_parameters, final_error) == 1)
	{
		fit_para_.x = Eigen::Vector3d(tree_parameters[0], tree_parameters[1], tree_parameters[2]);
		fit_para_.n = Eigen::Vector3d(tree_parameters[3], tree_parameters[4], tree_parameters[5]);
		fit_para_.r = tree_parameters[6];
		final_fit_error_ = final_error;
	}
	else
	{
		fit_tree_failed_ = true;
		tree_parameters.clear();
	}
	// update trans matrix
	derive_trans_matrix(FIT_CYLINDER);
}


/*****************************************************************
 * Update parameters, candidate Iscan trees, #of point, and R_trans
 * ***************************************************************/
/**
 * @brief Updates the MapTree with a newly matched IScan tree.
 *
 * Updates tree parameters, center, number of points, and transformation matrix.
 * Optionally recomputes a global cylinder fit if enabled.
 */
void MapTree::update()
{
    // update parameter
    para = cand_para;
#ifdef FIT_CYLINDER_TREE
	para = fit_para_;
#endif
    visibleInfo.push_back(candTreeInfo);
    int iscanId = candTreeInfo.first;
    int treeId = candTreeInfo.second;

    int num_points_cand_tree = pMap->mvpIScans[iscanId]->vTreePointMapping[treeId].size();

    Eigen::Vector3d cur_center = Eigen::Vector3d::Zero();
    for (int nP = 0; nP < pMap->mvpIScans[iscanId]->vTreePointMapping[treeId].size(); nP++)
    {
        // points in mapping
        Eigen::Vector3d curr_point(pMap->mvpIScans[iscanId]->vTreePointMapping[treeId][nP].x,
                                   pMap->mvpIScans[iscanId]->vTreePointMapping[treeId][nP].y,
                                   pMap->mvpIScans[iscanId]->vTreePointMapping[treeId][nP].z);

        cur_center += curr_point;
    }
    center = (center * numPoint + cur_center)/double(numPoint + num_points_cand_tree);
    numPoint += num_points_cand_tree;

    // update trans matrix
#ifdef FIT_CYLINDER_TREE
	if (visibleInfo.size() > 1)
	{
		CalculateAllPointsFitCylinderParameters();
		derive_trans_matrix(FIT_CYLINDER);
		para = fit_para_;
	}
#endif
	derive_trans_matrix();
    if(visibleInfo.size() == 1)
    {
        start_iscan_index_ = iscanId;
    }
    end_iscan_index_ = iscanId;

}

/**
 * @brief Computes residual error between a candidate iscan tree and the current MapTree model.
 *
 * @param cand_tree_info Pair containing scan and tree IDs.
 * @param paraVer Which set of parameters to use for computing residuals (CANDIDATE, FIT_CYLINDER, etc.).
 * @return Root Mean Square Error (RMSE) of point-to-cylinder distances.
 */

double MapTree::compute_iscan_tree_residual(pair<int, int> cand_tree_info, eParaVer paraVer)
{
    double dis_rmse = 0.0;
    int iscanId = cand_tree_info.first;
    int treeId = cand_tree_info.second;
    for (int nP = 0; nP < pMap->mvpIScans[iscanId]->vTreePointMapping[treeId].size(); nP++)
    {
        // points in mapping
        Eigen::Vector3d curr_point(pMap->mvpIScans[iscanId]->vTreePointMapping[treeId][nP].x,
                                   pMap->mvpIScans[iscanId]->vTreePointMapping[treeId][nP].y,
                                   pMap->mvpIScans[iscanId]->vTreePointMapping[treeId][nP].z);
        double dis = compute_dis(curr_point, paraVer);
        dis_rmse = dis_rmse + dis * dis;
    }

    dis_rmse = sqrt(dis_rmse / double(pMap->mvpIScans[iscanId]->vTreePointMapping[treeId].size()));
    return dis_rmse;
}

/**
 * @brief Computes residual error for all points currently assigned to this map tree.
 *
 * @param paraVer Parameter set to use for residual computation.
 * @return Root Mean Square Error (RMSE) of distances from points to current tree cylinder model.
 */

double MapTree::compute_map_tree_residual(eParaVer paraVer)
{
    double dis_rmse = 0.0;
    if(visibleInfo.size() == 0)
        return dis_rmse; 

    for (int nT = 0; nT < visibleInfo.size(); nT++)
    {
        int iscanId = visibleInfo[nT].first;
        int treeId = visibleInfo[nT].second;
        // points in mapping
        for (int nP = 0; nP < pMap->mvpIScans[iscanId]->vTreePointMapping[treeId].size(); nP++)
        {
            Eigen::Vector3d curr_point(pMap->mvpIScans[iscanId]->vTreePointMapping[treeId][nP].x,
                                       pMap->mvpIScans[iscanId]->vTreePointMapping[treeId][nP].y,
                                       pMap->mvpIScans[iscanId]->vTreePointMapping[treeId][nP].z);
            double dis = compute_dis(curr_point, paraVer);
            dis_rmse = dis_rmse + dis * dis;
        }
    }
    dis_rmse = sqrt(dis_rmse / double(numPoint));
    return dis_rmse;
}
