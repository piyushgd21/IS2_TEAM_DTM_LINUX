#include "../header/IntegratedScan.h"


/*****************************************************************
 * Transform the points in each individual scan to the start of integrated scan (local frame). 
 * Input: v_R_local/v_r_local, vTreePointRaw/vGroundPointRaw
 * Output: vTreePointLocal, vGroundPointLocal, pGroundPointLocalDs
 * downsample_distance: distance for downsampling
 *  option (default: 1):
 *      1.both ground points and tree points.
 *      2.tree points only
 *      3.ground points only
 * ***************************************************************/
void IntegratedScan::transform2Start(double downsample_distance, int option)
{
    if (option == 1 || option == 2)
    {
        vTreePointLocal.clear();
        for (int nTree = 0; nTree < vTreePointRaw.size(); nTree++)
        {
            vector<PointType> tempPoints;
            // for each tree point
            for (int nP = 0; nP < vTreePointRaw[nTree].size(); nP++)
            {
                int scanId = int(vTreePointRaw[nTree][nP].intensity / 100);
                Eigen::Vector3d curr_point(vTreePointRaw[nTree][nP].x, vTreePointRaw[nTree][nP].y, vTreePointRaw[nTree][nP].z);
                double t_ratio = double(vTreePointRaw[nTree][nP].intensity - int(vTreePointRaw[nTree][nP].intensity));

                Eigen::Quaterniond q_prev{v_R_local[scanId]};
                Eigen::Quaterniond q_cur_end{v_R_local[scanId + 1]};
                Eigen::Quaterniond q_cur = q_prev.slerp(t_ratio, q_cur_end);
                Eigen::Vector3d t_cur = (1.0 - t_ratio) * v_r_local[scanId] + t_ratio * v_r_local[scanId + 1];
                // Eigen::Vector3d rIm = q_cur * curr_point + t_cur;
                Eigen::Vector3d rIm = v_R_local[scanId] * curr_point +  v_r_local[scanId] ;

                PointType pm;
                pm.x = float(rIm(0));
                pm.y = float(rIm(1));
                pm.z = float(rIm(2));
                pm.intensity = vTreePointRaw[nTree][nP].intensity;
                tempPoints.push_back(pm);
            }
            vTreePointLocal.push_back(tempPoints);
        }
    }

    if (option == 1 || option == 3)
    {
        vGroundPointLocal.clear();
        pcl::PointCloud<PointType>::Ptr groundPoints(new pcl::PointCloud<PointType>);

        // for each ground point
        for (int nP = 0; nP < vGroundPointRaw.size(); nP++)
        {
            int scanId = int(vGroundPointRaw[nP].intensity / 100);
            Eigen::Vector3d curr_point(vGroundPointRaw[nP].x, vGroundPointRaw[nP].y, vGroundPointRaw[nP].z);
            double t_ratio = double(vGroundPointRaw[nP].intensity - int(vGroundPointRaw[nP].intensity));

            Eigen::Quaterniond q_prev{v_R_local[scanId]};
            Eigen::Quaterniond q_cur_end{v_R_local[scanId + 1]};
            Eigen::Quaterniond q_cur = q_prev.slerp(t_ratio, q_cur_end);
            Eigen::Vector3d t_cur = (1.0 - t_ratio) * v_r_local[scanId] + t_ratio * v_r_local[scanId + 1];
            // Eigen::Vector3d rIm = q_cur * curr_point + t_cur;
            Eigen::Vector3d rIm = v_R_local[scanId] * curr_point +  v_r_local[scanId];

            PointType pm;
            pm.x = float(rIm(0));
            pm.y = float(rIm(1));
            pm.z = float(rIm(2));
            pm.intensity = vGroundPointRaw[nP].intensity;
            vGroundPointLocal.push_back(pm);

            pm.intensity = float(index);
            groundPoints->push_back(pm);
        }

        pGroundPointLocalDs.reset(new pcl::PointCloud<PointType>);
        // pcl::VoxelGrid<PointType> downSizeFilter;
        // downSizeFilter.setInputCloud(groundPoints);
        // downSizeFilter.setLeafSize(downsample_distance, downsample_distance, downsample_distance);
        // downSizeFilter.filter(*pGroundPointLocalDs);

        //(TODO:pcl have integrater overflow problem)
        *pGroundPointLocalDs=DownSamplePointCloudBasedOnDistance(groundPoints,downsample_distance);
    }
}

/*****************************************************************
 * Compute the r_lu_m for each individual scan to based on r/R_lu_m_updated
 * and v_r/R_local:
 * Output:
 *  v_r_mapping, v_R_mapping
 * ***************************************************************/
void IntegratedScan::computeIndividualPoseMapping()
{
    v_r_mapping.clear();
    v_R_mapping.clear();
    for (int i = 0; i < v_r_local.size(); i++)
    {
        Eigen::Matrix3d R_map;
         R_map= R_local_m_updated * v_R_local[i];
        Eigen::Vector3d r_map;
        r_map = R_local_m_updated * v_r_local[i] + r_local_m_updated;

        //Eigen::Matrix3d R_map = R_local_m_updated * v_R_local[i];
        //Eigen::Vector3d r_map = R_local_m_updated * v_r_local[i] + r_local_m_updated;

        v_r_mapping.push_back(r_map);
        v_R_mapping.push_back(R_map);
    }
}

/*****************************************************************
 * Compute the mapping coordinates of raw ground points using v_r/R_mapping
 * vGroundPointRaw -> vGroundPointMapping
 * ***************************************************************/
void IntegratedScan::computeRawGroundPointToMapping()
{
    vector<PointType>().swap(vGroundPointMapping);

    // for each ground point
    for (int nP = 0; nP < vGroundPointRaw.size(); nP++)
    {
        int scanId = int(vGroundPointRaw[nP].intensity / 100);
        Eigen::Vector3d curr_point(vGroundPointRaw[nP].x, vGroundPointRaw[nP].y, vGroundPointRaw[nP].z);
        double t_ratio = double(vGroundPointRaw[nP].intensity - int(vGroundPointRaw[nP].intensity));

        Eigen::Quaterniond q_prev{v_R_mapping[scanId]};
        Eigen::Quaterniond q_cur_end{v_R_mapping[scanId + 1]};
        Eigen::Quaterniond q_cur = q_prev.slerp(t_ratio, q_cur_end);
        Eigen::Vector3d t_cur = (1.0 - t_ratio) * v_r_mapping[scanId] + t_ratio * v_r_mapping[scanId + 1];
        // Eigen::Vector3d rIm = q_cur * curr_point + t_cur;
        Eigen::Vector3d rIm = v_R_mapping[scanId] * curr_point +  v_r_mapping[scanId];

        PointType pm;
        pm.x = float(rIm(0));
        pm.y = float(rIm(1));
        pm.z = float(rIm(2));
        pm.intensity = vGroundPointRaw[nP].intensity;
        vGroundPointMapping.push_back(pm);
        
    }
}

/*****************************************************************
 * Compute the mapping coordinates of raw tree points using v_r/R_mapping
 * vTreePointRaw -> vTreePointMapping
 * ***************************************************************/
void IntegratedScan::computeRawTreePointToMapping()
{
    vector<vector<PointType>>().swap(vTreePointMapping);

    for (int nTree = 0; nTree < vTreePointRaw.size(); nTree++)
    {
        vector<PointType> tempPoints;
        // for each tree point
        for (int nP = 0; nP < vTreePointRaw[nTree].size(); nP++)
        {
            int scanId = int(vTreePointRaw[nTree][nP].intensity / 100);
            Eigen::Vector3d curr_point(vTreePointRaw[nTree][nP].x, vTreePointRaw[nTree][nP].y, vTreePointRaw[nTree][nP].z);
            double t_ratio = double(vTreePointRaw[nTree][nP].intensity - int(vTreePointRaw[nTree][nP].intensity));

            Eigen::Quaterniond q_prev{v_R_mapping[scanId]};
            Eigen::Quaterniond q_cur_end{v_R_mapping[scanId + 1]};
            Eigen::Quaterniond q_cur = q_prev.slerp(t_ratio, q_cur_end);
            Eigen::Vector3d t_cur = (1.0 - t_ratio) * v_r_mapping[scanId] + t_ratio * v_r_mapping[scanId + 1];
            // Eigen::Vector3d rIm = q_cur * curr_point + t_cur;
            Eigen::Vector3d rIm = v_R_mapping[scanId] * curr_point +  v_r_mapping[scanId];

            PointType pm;
            pm.x = float(rIm(0));
            pm.y = float(rIm(1));
            pm.z = float(rIm(2));
            pm.intensity = vTreePointRaw[nTree][nP].intensity;
            tempPoints.push_back(pm);
        }
        vTreePointMapping.push_back(tempPoints);
    }
}

/*****************************************************************
 * Compute the mapping coordinates of raw tree points using v_r/R_mapping
 * vTreePointRaw -> vTreePointMapping
 * ***************************************************************/
void IntegratedScan::computeBackupTreePointToMapping()
{
    vector<vector<vector<PointType>>>().swap(raw_tree_points_mapping);
    int count = 0;
    for (int nScan = 0; nScan < raw_tree_points.size(); nScan++)
    {
        vector<vector<PointType>> tree_points_per_scan;
        for (int nTree = 0; nTree < raw_tree_points[nScan].size(); nTree++)
        {
            vector<PointType> tempPoints;
            // for each tree point
            for (int nP = 0; nP < raw_tree_points[nScan][nTree].size(); nP++)
            {
                int scanId = nScan;
                Eigen::Vector3d curr_point(raw_tree_points[nScan][nTree][nP].x, raw_tree_points[nScan][nTree][nP].y, raw_tree_points[nScan][nTree][nP].z);
                double t_ratio = double(raw_tree_points[nScan][nTree][nP].intensity - int(raw_tree_points[nScan][nTree][nP].intensity));

                Eigen::Quaterniond q_prev{v_R_mapping[scanId]};
                Eigen::Quaterniond q_cur_end{v_R_mapping[scanId + 1]};
                Eigen::Quaterniond q_cur = q_prev.slerp(t_ratio, q_cur_end);
                Eigen::Vector3d t_cur = (1.0 - t_ratio) * v_r_mapping[scanId] + t_ratio * v_r_mapping[scanId + 1];
                // Eigen::Vector3d rIm = q_cur * curr_point + t_cur;
                Eigen::Vector3d rIm = v_R_mapping[scanId] * curr_point +  v_r_mapping[scanId];

                PointType pm;
                pm.x = float(rIm(0));
                pm.y = float(rIm(1));
                pm.z = float(rIm(2));
                pm.intensity = raw_tree_points[nScan][nTree][nP].intensity;
                tempPoints.push_back(pm);
                count++;
            }
            tree_points_per_scan.push_back(tempPoints);
        }
        raw_tree_points_mapping.push_back(tree_points_per_scan);
    }
    //cout << count <<endl;
}

/*****************************************************************
 * Compute the mapping coordinates of raw tree points using v_r/R_mapping
 * raw_non_ground_points -> raw_non_ground_points_mapping
 * ***************************************************************/
void IntegratedScan::computeNonGroundPointsToMapping()
{
    vector<vector<PointType>>().swap(raw_non_ground_points_mapping);
    int count = 0;
    for (int nScan = 0; nScan < raw_non_ground_points.size(); nScan++)
    {
        vector<PointType> points_per_scan;
        for (int nP = 0; nP < raw_non_ground_points[nScan].size(); nP++)
        {

            int scanId = nScan;
            Eigen::Vector3d curr_point(raw_non_ground_points[nScan][nP].x, raw_non_ground_points[nScan][nP].y, raw_non_ground_points[nScan][nP].z);
            double t_ratio = double(raw_non_ground_points[nScan][nP].intensity - int(raw_non_ground_points[nScan][nP].intensity));

            Eigen::Quaterniond q_prev{v_R_mapping[scanId]};
            Eigen::Quaterniond q_cur_end{v_R_mapping[scanId + 1]};
            Eigen::Quaterniond q_cur = q_prev.slerp(t_ratio, q_cur_end);
            Eigen::Vector3d t_cur = (1.0 - t_ratio) * v_r_mapping[scanId] + t_ratio * v_r_mapping[scanId + 1];
            // Eigen::Vector3d rIm = q_cur * curr_point + t_cur;
            Eigen::Vector3d rIm = v_R_mapping[scanId] * curr_point +  v_r_mapping[scanId];

            PointType pm;
            pm.x = float(rIm(0));
            pm.y = float(rIm(1));
            pm.z = float(rIm(2));
            pm.intensity = raw_non_ground_points[nScan][nP].intensity;
            points_per_scan.push_back(pm);
            count++;
        }
        raw_non_ground_points_mapping.push_back(points_per_scan);
    }
}


/*****************************************************************
 * Compute the mapping coordinates of local ground points using R_local_m_updated
 *  1. T_local_m_ini
 *  2. T_local_m_updated
 * pGroundPointLocalDs -> pGroundPointMappingDs
 * ***************************************************************/
void IntegratedScan::computeMapGroundPoints(const eUpdateFlag flag)
{
    pGroundPointMappingDs.reset(new pcl::PointCloud<PointType>);
    for (std::size_t nP = 0; nP < pGroundPointLocalDs->size(); nP++)
    {
        PointType p = pGroundPointLocalDs->points[nP];
        Eigen::Vector3f ril(p.x, p.y, p.z);
        Eigen::Vector3f rim;
        if (flag == INIT)
        {
            rim = R_local_m_ini.cast<float>() * ril + r_local_m_ini.cast<float>();
        }
        else if (flag == REFINED)
        {
            rim = R_local_m_updated.cast<float>() * ril + r_local_m_updated.cast<float>();
        }

        p.x = float(rim(0));
        p.y = float(rim(1));
        p.z = float(rim(2));
        pGroundPointMappingDs->push_back(p);
    }
}

/*****************************************************************
 * Compute the mapping coordinates of local ground points using R_local_m_updated
 *  1. T_local_m_ini
 *  2. T_local_m_updated
 * vTreePointLocal -> vTreePointMapping
 * ***************************************************************/
void IntegratedScan::computeMapTreePoints(const eUpdateFlag flag)
{
    vTreePointMapping.clear();
    for (int nT = 0; nT < vTreePointLocal.size(); nT++)
    {
        vector<PointType> tempPoints;
        // for each tree point
        for (int nP = 0; nP < vTreePointLocal[nT].size(); nP++)
        {
            PointType p = vTreePointLocal[nT][nP];
            Eigen::Vector3f rIl(p.x, p.y, p.z);
            Eigen::Vector3f rIm;
            if (flag == INIT)
            {
                rIm = R_local_m_ini.cast<float>() * rIl + r_local_m_ini.cast<float>();
            }
            else if (flag == REFINED)
            {
                rIm = R_local_m_updated.cast<float>() * rIl + r_local_m_updated.cast<float>();
            }

            PointType pm;
            pm.x = rIm(0);
            pm.y = rIm(1);
            pm.z = rIm(2);
            pm.intensity = p.intensity;
            tempPoints.push_back(pm);
        }
        vTreePointMapping.push_back(tempPoints);
    }
}

/*****************************************************************
 * Compute the tree parameters in mapping using
 *  1. T_local_m_ini
 *  2. T_local_m_updated
 * Input:
 *  vTreeParamLocal
 * Ouput:
 *  vTreeParamMapping
 * ***************************************************************/
void IntegratedScan::computeMapTreeParameter(const eUpdateFlag flag)
{
    vector<CylinderPara> vTreePara;
    CylinderPara tempPara;

    for (int nTree = 0; nTree < vTreeParamLocal.size(); nTree++)
    {
        if (flag == INIT)
        {
            // force nz = 1
            Eigen::Vector3d n_trans = R_local_m_ini * vTreeParamLocal[nTree].n;
            tempPara.n = n_trans / n_trans(2);

            tempPara.x = R_local_m_ini * vTreeParamLocal[nTree].x + r_local_m_ini;
        }
        else if (flag == REFINED)
        {
            Eigen::Vector3d n_trans = R_local_m_updated * vTreeParamLocal[nTree].n;
            tempPara.n = n_trans / n_trans(2);
            tempPara.x = R_local_m_updated * vTreeParamLocal[nTree].x + r_local_m_updated;
        }
        tempPara.r = vTreeParamLocal[nTree].r;
        vTreePara.push_back(tempPara);
    }
    vTreeParamMapping = vTreePara;
}

/*****************************************************************
 * Update the tree parameters for those who are not included in the scan integration
 * Input:
 *  vValidTree, vTreePointLocal
 * Output:
 *  vTreeParamLocal
 * **************************************************************/
void IntegratedScan::updateTreeParamLocal()
{
    for (int nTree = 0; nTree < vTreePointRaw.size(); nTree++)
    {
        // refined in optimization
        // if (vValidTree[nTree])
        //     continue;

        Eigen::Vector3f center = Eigen::Vector3f::Zero();
        vector<PointType> tempPoints;
        // for each tree point
        for (int nP = 0; nP < vTreePointLocal[nTree].size(); nP++)
        {
            Eigen::Vector3f r_I_local(vTreePointLocal[nTree][nP].x, vTreePointLocal[nTree][nP].y, vTreePointLocal[nTree][nP].z);
            center += r_I_local;
        }
        center = center / float(vTreePointLocal[nTree].size());
        CylinderPara treePara;
        treePara.x = center.cast<double>();
        treePara.n << 0.0, 0.0, 1.0;
        treePara.r = 0.0;
        vTreeParamLocal[nTree] = treePara;
    }
}

/*****************************************************************
 * Compute the tree parameters in mapping using
 *  1. T_local_m_ini
 *  2. T_local_m_updated
 * Input:
 *  vTreeParamLocal
 * Ouput:
 *  vTreeParamMapping
 * ***************************************************************/
void IntegratedScan::show_number_PointType()
{
    // tree points
    fDebug << endl;
    fDebug << "----------------------------------------" << endl;
    if(indScans.empty())
    {
        fDebug << "Problem----------------------------" <<endl;
        return;
    }
    fDebug << indScans[0].scanID << endl;
    fDebug << "Tree points " << endl;
    int count_all = 0;
    int count_temp = 0;

    for (int i = 0; i < vTreePointRaw.size(); i++)
        count_temp += vTreePointRaw[i].size();
    count_all += count_temp;
    fDebug << "\tRaw " << count_temp << endl;

    count_temp = 0;
    for (int i = 0; i < vTreePointLocal.size(); i++)
        count_temp += vTreePointLocal[i].size();
    count_all += count_temp;
    fDebug << "\tLocal " << count_temp << endl;

    count_temp = 0;
    for (int i = 0; i < vTreePointMapping.size(); i++)
        count_temp += vTreePointMapping[i].size();
    count_all += count_temp;
    fDebug << "\tMapping " << count_temp << endl;

    fDebug << "Ground points " << endl;

    count_temp = vGroundPointRaw.size();
    fDebug << "\tRaw " << count_temp << endl;
    count_all += count_temp;

    count_temp = vGroundPointLocal.size();
    fDebug << "\tLocal " << count_temp << endl;
    count_all += count_temp;

    count_temp = pGroundPointLocalDs->size();
    fDebug << "\tDS Local " << count_temp << endl;
    count_all += count_temp;

    count_temp = pGroundPointMappingDs->size();
    fDebug << "\tDS Map " << count_temp << endl;
    count_all += count_temp;

    fDebug << "All number of PointXYZI " << count_all << endl;

    fDebug << "Number of scans " << indScans.size() << endl;
    //count_all = 0;

    for (int i = 0; i < indScans.size(); i++)
    {
        count_temp = 0;
        count_temp += indScans[i].vGroundRawPoints.size();
        count_temp += indScans[i].vGroundPlanarRawPoints.size();
        for (int j = 0; j < indScans[i].vTreeRawPoints.size(); j++)
        {
            count_temp += indScans[i].vTreeRawPoints[j].size();
        }
        fDebug << "\tScan " << i << ": " << count_temp << endl;
        count_all += count_temp;
    }
    fDebug << "Including Points in ScanInfo " << count_all << endl;
}

/*****************************************************************
 * Form the new vTreePointRaw using:
 *  1. raw_tree_points
 *  2. matched_map_tree_id
 * Input:
 *  raw_tree_points, matched_map_tree_id
 * Ouput:
 *  vTreePointRaw
 * ***************************************************************/
void IntegratedScan::update_vTreePointRaw()
{
    vector<vector<PointType>>().swap( vTreePointRaw);
    vector<int> map_tree_ids;
    int map_tree_id = -1;
    int num_tree_points = 0;
    for (int nScan = 0; nScan < matched_map_tree_id.size(); nScan++)
    {
        for (int nTree = 0; nTree < matched_map_tree_id[nScan].size(); nTree++)
        {
            map_tree_id= matched_map_tree_id[nScan][nTree];
            if (map_tree_id == -1)
            {
                continue;
            }

            //find if the tree has been created, find the index correspond to  map_tree_ids/vTreePointRaw
            int index = -1;
            for(int i = 0; i <map_tree_ids.size(); i++)
            {
                if(map_tree_id == map_tree_ids[i])
                {
                    index = i;
                    break;
                }
            }
            //if not, create a new one
            if(index == -1)
            {
                index = map_tree_ids.size();
                map_tree_ids.push_back(map_tree_id);
                vTreePointRaw.push_back(vector<PointType>());
            }
            num_tree_points+= raw_tree_points[nScan][nTree].size();
            //add the points in raw_tree_points to vTreePointRaw
            for (int nP = 0; nP < raw_tree_points[nScan][nTree].size(); nP++)
            {
                PointType p = raw_tree_points[nScan][nTree][nP];
                p.intensity += float(100 * nScan);
                vTreePointRaw[index].push_back(p);
            }
        }
    }
    corresponding_map_tree_ids = map_tree_ids;
    fDebug<<"Iscan " << index << ": " <<  num_tree_points <<" points, " << map_tree_ids.size() << "trees" <<endl;
    computeRawTreePointToMapping();

}

