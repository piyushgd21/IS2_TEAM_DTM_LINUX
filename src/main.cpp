#include "../header/featureExtraction.h"
#include "../header/Mapping.h"
#include "../header/CTrajectory.h"
#include "../header/utility.h"

int execute(std::string folder_path,std::string setting_file_name,std::string output_folder_path, std::string exe_name)
{
    //system("read -p 'Press Enter to continue...' var");
    SCAN_DURATION = 100.0;

    std::string dataset_folder = folder_path;
    dataset_folder = dataset_folder + "/";

    input_folder = dataset_folder + "sequence/";
	output_folder = output_folder_path;//dataset_folder + "result/";
    boost::filesystem::create_directories(output_folder);

    output_folder_odometry = output_folder + "odometry/";
    boost::filesystem::create_directories(output_folder_odometry);

    output_folder_loopclouse = output_folder + "loopclosure/";
    boost::filesystem::create_directories(output_folder_loopclouse);

    //-------------Load Parameter ----------------
	std::string paraPath = input_folder+setting_file_name;//input_folder + "setting.txt";
	std::string copy_to_file = output_folder + setting_file_name;
	CopyFile(wstring(paraPath.begin(), paraPath.end()).c_str(), wstring(copy_to_file.begin(), copy_to_file.end()).c_str(),false); 
    SettingPara sPara;
    loadSettingPara(paraPath, sPara);
    cout << sPara.nChannel << "\t" << sPara.treeAngleThreshold << "\t" << sPara.groundBufferTree << endl;

    //-------------Load Time tag------------------
    std::vector<double> vTime;     // in ms
    std::vector<double> vDuration; // in ms
	std::vector<int> vPointsNumber;
    std::string timePass = input_folder + "times.txt";
    ifstream fTime;
    fTime.open(timePass);

    if (!fTime)
    {
        cout << "failed to open" << timePass << endl;
        throw std::runtime_error("Wrong file name");
    }
    std::string str1, str2, str3, str4;

    while (fTime >> str1 >> str2 >> str3 >> str4)
    {
        vTime.push_back(stod(str1));
        vDuration.push_back(stod(str2));
		vPointsNumber.push_back(stoi(str4));
    }
    // vTime.push_back(vTime.back() + vTime.back() - vTime[vTime.size() - 2]);
    cout << vTime[0] << "\t" << vTime[vTime.size() - 2] << "\t" << vTime.back() << "\t" << vTime.size() << endl;
    cout << vTime[0] << "\t" << vTime[vTime.size() - 2] << "\t" << vTime.back() << "\t" << vTime.size() << endl;

    //-------------Load Trajectory------------------
    CTrajectory *cTrajectory = NULL;
    //if (argc == 3)
    {
        cTrajectory = new CTrajectory();
        std::string lidarParaPass = input_folder + "lidar_para.txt";
        cTrajectory->loadLidarPara(lidarParaPass);
		//(TODO)
        std::string trajPass = input_folder + "trajectory.txt";
        cTrajectory->loadTraj(trajPass);
        //export EOP 
        cTrajectory->computeEop();
        cout << "Number of trajectory event: " << cTrajectory->bopList.size() << endl;
    }
    //system("read -p 'Press Enter to continue...' var");

    std::string outPass =  output_folder_odometry + "trajectory_info.txt";
    std::string outPass2 = output_folder_odometry + "mapped_lidar.txt";
    std::string outPass3 = output_folder_odometry + "_Odometry_LOG.txt";
    std::string outPass4 = output_folder + "_Debug.txt";
    std::string outPass6 = output_folder_odometry + "trajectory_ref.txt";
    std::string outPass7 = output_folder_odometry + "trajectory_residual.txt";

    std::string outPassResPb = output_folder_odometry + "trajectory_res_pb.txt";
    std::string outPassResFb = output_folder_odometry + "trajectory_res_fb.txt";

    f_odometry_debug.open(output_folder_odometry + "_Debug_Odo.txt", std::ifstream::out);
    f_mapping_debug.open(output_folder+"_Debug_Mapping.txt", std::ifstream::out);
    fDebug.open(outPass4, std::ifstream::out);
	fDebug<< "exe: " << exe_name <<endl;

#ifdef EXPORT_TRAJECTORY
    fMapping.open(outPass, std::ifstream::out);
    fMapping << fixed << std::setprecision(4);
    fMapping << "ScanId\tX\tY\tZ\tOmega\tPhi\tKappa\tNumTree\tOdometryFlag" <<endl;
    if (cTrajectory)
    {
        fTrajectoryRef.open(outPass6, std::ifstream::out);
        fTrajectoryRef << fixed << std::setprecision(4);

        fTrajectoryRes.open(outPass7, std::ifstream::out);
        fTrajectoryRes << fixed << std::setprecision(4);

        fTrajectoryResPb.open(outPassResPb, std::ifstream::out);
        fTrajectoryResPb << fixed << std::setprecision(4);

        fTrajectoryResFb.open(outPassResFb, std::ifstream::out);
        fTrajectoryResFb << fixed << std::setprecision(4);
    }
#endif

#ifdef EXPORT_LOG
    fLog.open(outPass3, std::ifstream::out);
	fLog <<"exe: "<< exe_name <<endl;
    fLog << fixed << std::setprecision(4);


#endif

#ifdef EXPORT_RESULT
    int subCount = 1;
    if (cTrajectory)
    {
        std::string outPass5 = output_folder + "mapped_lidar_ref_" + to_string(subCount) + ".txt";
        fMappedRef.open(outPass5);
        fMappedRef << fixed << std::setprecision(4);
    }
#endif

    TicToc t;
    int trackingID = 0;
    int trackedNum = 0;
    vector<int> vecTrackedNum;
    LidarScan *prevScan = NULL;

    int scanCount = 0;

    Mapping *pMap = NULL;
#ifdef MAP
    pMap = new Mapping(sPara, cTrajectory);
    std::thread tMap(&Mapping::buildMap, pMap);
#endif 

    //************ Main loop for lidar odometry tracking ***********************
    for (size_t scanIndex = sPara.initScan; scanIndex <= sPara.endScan; scanIndex++) // index for scan  //natural, 300 - 4800 //ouster plantation: 400-28600
    {
        std::stringstream lidar_data_path;
        lidar_data_path << input_folder << "vlp/" << scanIndex << ".bin";

		std::stringstream tree_feature_file;
		tree_feature_file << input_folder << "tree_trunk/" << scanIndex << ".txt";

		std::stringstream ground_feature_file;
		ground_feature_file << input_folder << "ground/" << scanIndex << ".txt";

        // check if tracking to previous scan and if stop tracking due to large gap
        // two types of gap: 1. Gap in LiDAR scan. 2. Gap in the timestamp. end of current scan v.s. end of previous scan.
        if (prevScan)
        {
            // lose track for long (5 scan or 0.5 second), initialize a new tracking.
            if (scanIndex - prevScan->id > 5 || (vTime[scanIndex] + vDuration[scanIndex] - prevScan->mTimeEnd) > 3.0 * 1000.0) //@@modify
            {
                prevScan = NULL;
                trackingID++; // number of tracked
                vecTrackedNum.push_back(trackedNum);
                trackedNum = 0;
            }
        }

        //// Scans.emplace_back(lidar_data_path.str(), index);
        //// duration = duration > 150.0 ? SCAN_DURATION : duration;
        //// map->mpVecScan.push_back(tempScan);
        //// Scans.push_back(tempScan);

        // check the duration of a scan, if the duration is too short -> not a complete scan, ignore
        double duration = vDuration[scanIndex];
        if (duration < 0.8 * SCAN_DURATION)
            continue;

        // call odometry thread for this scan

        LidarScan *tempScan = new LidarScan(lidar_data_path.str(), scanIndex, prevScan, trackingID, vTime[scanIndex], duration, sPara, tree_feature_file.str(),
			ground_feature_file.str(), vPointsNumber[scanIndex],cTrajectory, pMap);
        scanCount++;

        // if current scan is valid, used as reference scan for the next scan
        if (tempScan->valid_odo_flag_)
        {
            delete prevScan;
            prevScan = tempScan;
            trackedNum++;
        }
        else
        {
            delete tempScan;
            tempScan = nullptr;
        }

#ifdef EXPORT_RESULT
        int subCount = 1;
        if (cTrajectory)
        {
            if (scanCount == 600)
            {
                scanCount = 0;
                subCount++;
                fMappedRef.close();
                std::string outPass5 = output_folder + "mapped_lidar_ref_" + to_string(subCount) + ".txt";
                fMappedRef.open(outPass5);
                fMappedRef << fixed << std::setprecision(4);
            }
        }
#endif      
    }

    vecTrackedNum.push_back(trackedNum);
    cout << t.toc() / 1000.0 << endl;
    cout << vecTrackedNum.size() << endl;

    for (int i = 0; i < vecTrackedNum.size(); i++)
    {
        if (vecTrackedNum[i] > 1)
            cout << i << " - " << vecTrackedNum[i] << endl;
    }
    cout << "Loading: " << TimeLoading / 1000.0 << endl;
    cout << "Segmentation: " << TimeSegment / 1000.0 << endl;
    cout << "Pb: " << TimePoint / 1000.0 << endl;
    cout << "Fb: " << TimeFeature / 1000.0 << endl;
    cout << "Ground extraction: " << Tground / 1000.0 << endl;
    cout << "Tree extraction: " << Ttree / 1000.0 << endl;
    cout << "Topt: " << Topt / 1000.0 << endl;
    cout << "Ttemp: " << Tcheck / 1000.0 << endl;

    ofstream fSummary(output_folder_odometry + "Odometry_Summary.txt");
    fSummary << "Utilized scan index: " << sPara.initScan << " -> " << sPara.endScan << endl;
    fSummary << "Odometry setting: " << endl;
    fSummary << "\tFeature-based: " << sPara.bFeatureBased << ", level information: " << sPara.bLevel << endl;
    fSummary << "\tPoint-based: " << sPara.bPointBased << ", distortion: " << sPara.bDistortion << endl;

    fSummary << endl;
    fSummary << "Processing time: " <<  t.toc() / 1000.0 << endl;
    fSummary << "\tLoading: " << TimeLoading / 1000.0 << endl;
    fSummary << "\tSegmentation: " << TimeSegment / 1000.0 << endl;

    if (sPara.bFeatureBased)
    {
        fSummary << "\tFb: " << TimeFeature / 1000.0 << endl;
        fSummary << "\t\tGround extraction: " << Tground / 1000.0 << endl;
        fSummary << "\t\tTree extraction: " << Ttree / 1000.0 << endl;
        fSummary << "\t\tTopt: " << Topt / 1000.0 << endl;
    }
    if (sPara.bPointBased)
    {
        fSummary << "\tPb: " << TimePoint / 1000.0 << endl;
        fSummary << "\t\tPoint extraction: " << Tp_extraction / 1000.0 << endl;
        fSummary << "\t\tMatching Optimization : " << Tp_opt / 1000.0 << endl;
    }

#ifdef MAP
    tMap.join();
#endif

    // outPass = output_folder + "trajectory_info2.txt";
    // cout << outPass << endl;
    // ofstream fNew;
    // fNew.open(outPass, std::ifstream::out);
    // fNew << fixed << std::setprecision(4);
    // cout << map->mVecScan.size() << endl;

    // for (int i = 0; i < map->mVecScan.size(); i++)
    // {
    //     fNew << map->mVecScan[i].scanID << "\t" << map->mVecScan[i].eop.XO << "\t" << map->mVecScan[i].eop.YO << "\t" << map->mVecScan[i].eop.ZO
    //          << "\t" << rad2deg(map->mVecScan[i].eop.omega) << "\t" << rad2deg(map->mVecScan[i].eop.phi) << "\t" << rad2deg(map->mVecScan[i].eop.kappa)
    //          << "\t" << map->mVecScan[i].partID << endl;
    // }
    // fNew.close();
    // fMapping.close();

    // outPass = output_folder + "trees.txt";

    // fMapping.open(outPass, std::ifstream::out);
    // fMapping << fixed << std::setprecision(4);
    // for (int i = 0; i < map->mVecTrees.size(); i++)
    // {
    //     fMapping << map->mVecTrees[i].id << ": ";
    //     for (int j = 0; j < map->mVecTrees[i].vecScanId.size(); j++)
    //         fMapping << map->mVecTrees[i].vecScanId[j] << " ";
    //     fMapping << endl;
    // }
    // fMapping.close();
    // fMappedPoints.close();
    // outPass = output_folder + "treesPoints.txt";

    // fMappedPoints.open(outPass, std::ifstream::out);
    // fMappedPoints << fixed << std::setprecision(4);
    // map->exportTree(true);
    // fMappedPoints.close();

#ifdef EXPORT_LOG
    fLog << endl;
    fLog << "----------------------------" << endl;
    fLog << "Processing Time " << t.toc() / 1000.0 << " s" << endl;
    fLog.close();
#endif
	f_odometry_debug.close();
	f_mapping_debug.close();
	fTrajectoryRef.close();
	fTrajectoryRes.close();
	fTrajectoryResPb.close();
	fTrajectoryResFb.close();
	fMappingTraj.close();
	fDebug.close();
	fMapping.close();


    //system("read -p 'Press Enter to continue...' var");
    delete prevScan;
    //system("read -p 'Press Enter to continue...' var");
	delete cTrajectory;
    delete pMap;
    //system("read -p 'Press Enter to continue...' var");
	return 1;

}

bool ReadBatchProcessFile(std::string file_name,std::vector<std::string>& batch_input_folders,std::vector<string>& batch_setting_file_names)
{
	std::ifstream file_stream(file_name);
	if (!file_stream)
	{
		std::cout << "Check file path please." << std::endl;
		return false;
	}
	std::string field_name, value;

	while (file_stream >> field_name >> value)
	{

		if (!field_name.compare("-folder"))
		{
			batch_input_folders.push_back(value);
		}
		if (!field_name.compare("-setting"))
		{
			batch_setting_file_names.push_back(value);
		}
	}
	return true;
}

std::string GetPureFileName(std::string *file_name)
{
	int slash_pos = (*file_name).find_last_of('/');
	int dot_pos = (*file_name).find_last_of('.');
	std::string suffix((*file_name).substr(dot_pos, (*file_name).length() - dot_pos));
	std::string pure_name((*file_name).substr(slash_pos + 1, (dot_pos - slash_pos - 1)));
	std::string dir_path((*file_name).substr(0, slash_pos + 1));

	return pure_name;
	//unsigned fileIndex = 1;
	//while (_access((*file_name).c_str(), 0) != -1)
	//{
	//	//printf( "File ACCESS.C exists " );
	//	(*file_name) = dir_path + pure_name + "_" + boost::to_string(fileIndex) + suffix;
	//	fileIndex++;
	//}

	//std::cout << "file_name : " << (*file_name) << std::endl;
}

int main(int argc, char **argv)
{
	std::cout << "this is first line for this program." << std::endl;
	std::string batch_process_file_name = argv[1];
	std::string exe_name= argv[0];
	std::vector<std::string> batch_input_folders;
	std::vector<string> batch_setting_file_names;
	ReadBatchProcessFile(batch_process_file_name, batch_input_folders, batch_setting_file_names);
	if (batch_input_folders.size() != batch_setting_file_names.size())
	{
		std::cout << "check input folders size and setting file size: "<< batch_input_folders.size()<<","<< batch_setting_file_names.size() << endl;
		throw std::runtime_error("Wrong file size");
	}

	for (int i = 0; i < batch_input_folders.size(); i++)
	{
		TimeFeature = 0.0;
		TimePoint = 0.0;
		TimeLoading = 0.0;
		TimeSegment = 0.0;
		Tground = 0.0;
		Ttree = 0.0;
		Topt = 0.0;
		Tcheck = 0.0;
		Tp_extraction = 0.0;
		Tp_opt = 0.0;
		std::string setting_pure_file_name = GetPureFileName(&batch_setting_file_names[i]);
		std::string output_folder_path = batch_input_folders[i] + "/" + setting_pure_file_name + "_result/";
		std::cout << "processing: " << output_folder_path << std::endl;
		int result=execute(batch_input_folders[i], batch_setting_file_names[i], output_folder_path, exe_name);

		std::cout << "finished: " << output_folder_path << " result:"<< result<<std::endl;
	}
	return 1; 

}