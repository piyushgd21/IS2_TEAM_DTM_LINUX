#include "../header/utility.h"
/** Output file streams for various logs and data exports. */
std::ofstream fMapping;
std::ofstream fMappingTraj;


std::ofstream fTrajectoryRef;
std::ofstream fTrajectoryRes;
std::ofstream fTrajectoryResPb;
std::ofstream fTrajectoryResFb;
std::ofstream f_odometry_debug;

std::ofstream fLog;
std::ofstream fMappedPoints;
std::ofstream f_mapping_debug;
std::ofstream fDebug;
std::ofstream fMappedRef;

/** Input and output folder paths. */
std::string output_folder;
std::string output_folder_odometry;
std::string output_folder_loopclouse;
std::string output_folder_iscan;
std::string output_folder_iscan_map;



std::string input_folder;


/** Timing metrics for profiling stages of SLAM. */
double TimeFeature = 0.0;
double TimePoint = 0.0;
double TimeLoading = 0.0;
double TimeSegment = 0.0;
double Tground = 0.0, Ttree = 0.0, Topt = 0.0, Tcheck = 0.0, Tp_extraction = 0.0, Tp_opt = 0.0;
/** Nominal LiDAR scan duration in milliseconds. */
double SCAN_DURATION = 100.0;

/** Mutex for synchronized transformations. */
std::mutex mTransMutex;

/**
 * @brief Computes the median of a float vector.
 * 
 * @param v Vector of float values.
 * @return Median value.
 */
float median(vector<float> &v)
{
    size_t n = v.size() / 2;
    nth_element(v.begin(), v.begin() + n, v.end());
    return v[n];
}

/* convert the plane parameters by setting the largest value of normal vector to 1
Input: plane parameters
Output: updated plane parameters, index of the component to 1
*/
/**
 * @brief Normalizes plane parameters so that the largest absolute component of the normal is 1.
 * 
 * @param params Plane parameters (a, b, c, d).
 * @return Index of the fixed component.
 */
int planeParamTrans(Eigen::Vector4d &params)
{
    double maxElement = max(max(abs(params(0)), abs(params(1))), abs(params(2)));

    double ratio;
    int fixIndex = -1;

    if (abs(params(0)) == maxElement)
    {
        fixIndex = 0;
        ratio = 1.0 / params(0);
    }
    else if (abs(params(1)) == maxElement)
    {
        fixIndex = 1;
        ratio = 1.0 / params(1);
    }
    else
    {
        fixIndex = 2;
        ratio = 1.0 / params(2);
    }

    params = params*ratio;
    return(fixIndex);
}

/**
 * @brief Fits a line to 3D points using PCA.
 * 
 * @param points Input 3D points.
 * @param normal Output direction vector of the line.
 * @param centerPoint Output centroid of the line.
 * @return True if line is valid.
 */
bool lineFitting(const std::vector<Eigen::Vector3f> points, Eigen::Vector3f &normal, Eigen::Vector3f &centerPoint)
{
    int numPoint = points.size();
    Eigen::Vector3f center(0, 0, 0);
    for (int j = 0; j < numPoint; j++)
    {
        center = center + points[j];
    }
    center = center / float(numPoint);

    Eigen::Matrix3f covMat = Eigen::Matrix3f::Zero();
    for (int j = 0; j < numPoint; j++)
    {
        Eigen::Matrix<float, 3, 1> tmpZeroMean = points[j] - center;
        covMat = covMat + tmpZeroMean * tmpZeroMean.transpose();
    }

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> saes(covMat);
    normal = saes.eigenvectors().col(2);
    normal = normal(2) < 0 ? -normal : normal;
    centerPoint = center;
    // if is indeed line feature
    // note Eigen library sort eigenvalues in increasing order
    if (saes.eigenvalues()[2] > 3 * saes.eigenvalues()[1])
    {
        return true;
    }
    else
        return false;
}

/*return the plane parameters, first 3 elements is the normal vector.
if pass the pca analysis, return yes*/
/**
 * @brief Fits a plane to 3D points using PCA and least squares.
 * 
 * @param points Input 3D points.
 * @param params Output plane coefficients (a, b, c, d).
 * @param centerPoint Output centroid of the plane.
 * @return True if the plane passes shape criteria.
 */
bool planeFitting(const std::vector<Eigen::Vector3f> points, Eigen::Vector4f &params, Eigen::Vector3f &centerPoint)
{
    int numPoint = points.size();
    // cout<< numPoint <<endl;
    const int n = numPoint;
    Eigen::Vector3f center(0, 0, 0);

    for (int j = 0; j < numPoint; j++)
    {
        center = center + points[j];
    }
    center = center / float(numPoint);
    centerPoint = center;

    Eigen::Matrix3f covMat = Eigen::Matrix3f::Zero();
    for (int j = 0; j < numPoint; j++)
    {
        Eigen::Matrix<float, 3, 1> tmpZeroMean = points[j] - center;
        covMat = covMat + tmpZeroMean * tmpZeroMean.transpose();
    }
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> saes(covMat);
    // cout << saes.eigenvalues().transpose() << endl;

    // large -> small
    float lam1 = saes.eigenvalues()(2);
    float lam2 = saes.eigenvalues()(1);
    float lam3 = saes.eigenvalues()(0);

    float a1 = (sqrt(lam1) - sqrt(lam2)) / sqrt(lam1);
    float a2 = (sqrt(lam2) - sqrt(lam3)) / sqrt(lam1);
    float a3 = sqrt(lam3) / sqrt(lam1);

    // system("read -p 'Press Enter to continue...' var");
    bool planeValid = true;
    if (!(a2 > 2.0 * a3 && a2 * 3.0 > a1))
        planeValid = false;
    // if (!(a2 > a1 && a2 > a3))
    //     planeValid = false;
    //  cout << a1 << "\t" << a2 << "\t" << a3 << "\t" << planeValid << endl;

    Eigen::MatrixXf matA0(n, 3);
    Eigen::MatrixXf matB0(n, 1);
    // Eigen::Matrix<float, numPoint, 3> matA0;
    // Eigen::Matrix<float, numPoint, 1> matB0 = -1 * Eigen::Matrix<float, numPoint, 1>::Ones();

    for (int j = 0; j < numPoint; j++)
    {
        matA0(j, 0) = points[j](0);
        matA0(j, 1) = points[j](1);
        matA0(j, 2) = points[j](2);
        matB0(j) = -1.0;
    }
    // find the norm of plane
    Eigen::Vector3f norm = matA0.colPivHouseholderQr().solve(matB0);
    float negative_OA_dot_norm = 1 / norm.norm();
    norm.normalize();

    // // Here n(pa, pb, pc) is unit norm of plane
    // for (int j = 0; j < numPoint; j++)
    // {
    //     // if OX * n > 0.2, then plane is not fit well
    //     // cout << fabs(norm.transpose() * points[j] + negative_OA_dot_norm) << "\t";
    //     if (fabs(norm.transpose() * points[j] + negative_OA_dot_norm) > 0.2)
    //     {
    //         planeValid = false;
    //         break;
    //     }
    // }

    if (norm(2) < 0)
        params << norm * (-1.0), -negative_OA_dot_norm;
    else
        params << norm, negative_OA_dot_norm;

    // if is indeed line feature
    // note Eigen library sort eigenvalues in increasing order
    if (planeValid)
        return true;
    else
        return false;
}

/**
 * @brief Fits a plane to 3D points and checks for outliers.
 * 
 * @param points Input 3D points.
 * @param params Output plane coefficients.
 * @param centerPoint Output centroid.
 * @param fail_type Reason for failure (1 = eigen check, 2 = too many outliers).
 * @return True if fitting succeeds.
 */
bool planeFitting_outlier_check(const std::vector<Eigen::Vector3f> points, Eigen::Vector4f &params, Eigen::Vector3f &centerPoint, int &fail_type)
{
    int numPoint = points.size();
    // cout<< numPoint <<endl;
    const int n = numPoint;
    Eigen::Vector3f center(0, 0, 0);

    for (int j = 0; j < numPoint; j++)
    {
        center = center + points[j];
    }
    center = center / float(numPoint);
    centerPoint = center;

    Eigen::Matrix3f covMat = Eigen::Matrix3f::Zero();
    for (int j = 0; j < numPoint; j++)
    {
        Eigen::Matrix<float, 3, 1> tmpZeroMean = points[j] - center;
        covMat = covMat + tmpZeroMean * tmpZeroMean.transpose();
    }
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> saes(covMat);
    // cout << saes.eigenvalues().transpose() << endl;

    // large -> small
    float lam1 = saes.eigenvalues()(2);
    float lam2 = saes.eigenvalues()(1);
    float lam3 = saes.eigenvalues()(0);

    float a1 = (sqrt(lam1) - sqrt(lam2)) / sqrt(lam1);
    float a2 = (sqrt(lam2) - sqrt(lam3)) / sqrt(lam1);
    float a3 = sqrt(lam3) / sqrt(lam1);

    // system("read -p 'Press Enter to continue...' var");
    bool planeValid = true;
    //if (!(a2 > 2.0 * a3 && a2 * 3.0 > a1))
    //if (!(a2 >  a3 && a2 * 3.0 > a1))
    if(a3 > a1 || a3 > a2)
    {
        planeValid = false;
        fail_type = 1;
        return false;
    }
    // if (!(a2 > a1 && a2 > a3))
    //     planeValid = false;
    //  cout << a1 << "\t" << a2 << "\t" << a3 << "\t" << planeValid << endl;

    Eigen::MatrixXf matA0(n, 3);
    Eigen::MatrixXf matB0(n, 1);
    // Eigen::Matrix<float, numPoint, 3> matA0;
    // Eigen::Matrix<float, numPoint, 1> matB0 = -1 * Eigen::Matrix<float, numPoint, 1>::Ones();

    for (int j = 0; j < numPoint; j++)
    {
        matA0(j, 0) = points[j](0);
        matA0(j, 1) = points[j](1);
        matA0(j, 2) = points[j](2);
        matB0(j) = -1.0;
    }
    // find the norm of plane
    Eigen::Vector3f norm = matA0.colPivHouseholderQr().solve(matB0);
    float negative_OA_dot_norm = 1 / norm.norm();
    norm.normalize();

    //check number of valid points
    int outlier_points = 0;
    for (int j = 0; j < numPoint; j++)
    {
        if (fabs(norm.transpose() * points[j] + negative_OA_dot_norm) > 0.1)
        {
            ++outlier_points;
            if( float(outlier_points)/float(numPoint) > 0.2)
            {
                planeValid = false;
                fail_type = 2;
                break;
            }
        }
    }

    // // Here n(pa, pb, pc) is unit norm of plane
    // for (int j = 0; j < numPoint; j++)
    // {
    //     // if OX * n > 0.2, then plane is not fit well
    //     // cout << fabs(norm.transpose() * points[j] + negative_OA_dot_norm) << "\t";
    //     if (fabs(norm.transpose() * points[j] + negative_OA_dot_norm) > 0.2)
    //     {
    //         planeValid = false;
    //         break;
    //     }
    // }

    if (norm(2) < 0)
        params << norm * (-1.0), -negative_OA_dot_norm;
    else
        params << norm, negative_OA_dot_norm;

    // if is indeed line feature
    // note Eigen library sort eigenvalues in increasing order
    if (planeValid)
        return true;
    else
        return false;
}

/**
 * @brief Computes 2D similarity transform (rotation + translation) between matched 2D points.
 * 
 * @param points1 Set 1 of 2D points.
 * @param points2 Set 2 of 2D points.
 * @param p Optional ID used for debug file naming.
 * @return Vector with [theta, tx, ty].
 */
Eigen::Vector3f compute2dSimilarity(const std::vector<Eigen::Vector3f> points1, const std::vector<Eigen::Vector3f> points2, int p)
{
    std::vector<pair<int, int>> pointPairs;
    int numP1 = points1.size(), numP2 = points2.size();
    std::vector<Eigen::Vector3f> m1;
    std::vector<Eigen::Vector3f> m2;

    for (int i = 0; i < numP1; i++)
    {
        float x1 = points1[i](0), y1 = points1[i](1);
        for (int j = 0; j < numP2; j++)
        {
            float x2 = points2[j](0), y2 = points2[j](1);
            float d_sqr = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);

            if (d_sqr < 1.0)
            {
                pointPairs.push_back(make_pair(i, j));
                m1.push_back(points1[i]);
                m2.push_back(points2[j]);
                //                cout << d_sqr << "\t" << i << "\t" << j << endl;
                break;
            }
        }
    }

    int numPair = pointPairs.size();
    const int n = 2 * numPair;

    Eigen::MatrixXf A(n, 4); // cos, sin, tx, ty. x' = x*cos + (-y)*sin + tx   y' = x*sin + y*cos+ty;
    Eigen::MatrixXf y(n, 1);
    for (int i = 0; i < numPair; i++)
    {
        A(i * 2, 0) = m1[i](0);
        A(i * 2, 1) = -m1[i](1);
        A(i * 2, 2) = 1;
        A(i * 2, 3) = 0;
        A(i * 2 + 1, 0) = m1[i](1);
        A(i * 2 + 1, 1) = m1[i](0);
        A(i * 2 + 1, 2) = 0;
        A(i * 2 + 1, 3) = 1;

        y(i * 2) = m2[i](0);
        y(i * 2 + 1) = m2[i](1);
    }
    Eigen::MatrixXf res1, res2;

    Eigen::Vector4f nonmial(1, 0, 0, 0);
    res1 = y - A * nonmial;

    Eigen::Vector4f result = A.colPivHouseholderQr().solve(y);
    res2 = y - A * result;

    cout << "Number of points: " << numP1 << "/" << numP2 << "\t" << numPair << endl;

    cout << res1.norm() << " -> " << res2.norm() << endl;
    // cout << result.transpose() << endl;
    float scale = sqrt(result(0) * result(0) + result(1) * result(1));
    float theta = atan2(result(1), result(0));
    float tx = result(2);
    float ty = result(3);
    cout << "Result: " << scale << "\t" << theta << "\t" << tx << "\t" << ty << endl;

#ifdef EXPORT_RESULT
    fLog << "M:"
         << "\t" << numP1 << "/" << numP2 << "\t" << numPair << endl;
    fLog << "\t" << scale << "\t" << rad2deg(theta) << "\t" << tx << "\t" << ty << endl;

    // if(numPair <5)
    // fDebug << id<<"-Number of tree pairs: " << numPair;
#endif

    if (0)
    {
        string out("/home/ztian/catkin_ws/src/A-LOAM/data/DPRG/Ouster_backpack/20210927_martell/sequence/2d_" + to_string(p) + ".txt");

        std::ofstream fSegFile(out, std::ifstream::out);
        fSegFile << fixed << std::setprecision(4);

        float tx_ = 0.0, ty_ = 0.0;
        for (int i = 0; i < numPair; i++)
        {
            float x1_ = m1[i](0) * cos(theta) + (-m1[i](1)) * sin(theta);
            float y1_ = m1[i](0) * sin(theta) + m1[i](1) * cos(theta);

            tx_ = tx_ + m2[i](0) - x1_;
            ty_ = ty_ + m2[i](1) - y1_;
        }
        tx_ = tx_ / float(numPair);
        ty_ = ty_ / float(numPair);
        cout << "Result 1: "
             << "\t" << tx_ << "\t" << ty_ << endl;

        for (int i = 0; i < numPair; i++)
        {
            float x1_ = m1[i](0) * cos(theta) + (-m1[i](1)) * sin(theta) + tx;
            float y1_ = m1[i](0) * sin(theta) + m1[i](1) * cos(theta) + ty;

            fSegFile << x1_ << "\t" << y1_ << "\t" << 0 << "\t" << i << "\t" << 1 << endl;
            // fSegFile << m1[i](0) << "\t" << m1[i](1) << "\t" << 0 << "\t" << i << "\t" << 1 << endl;

            fSegFile << m2[i](0) << "\t" << m2[i](1) << "\t" << 0 << "\t" << i << "\t" << 2 << endl;
        }
        fSegFile.close();
        // system("read -p 'Press Enter to continue...' var");
    }
    Eigen::Vector3f para;
    para << theta, tx, ty;
    return para;
}

/**
 * @brief Loads settings from a configuration file into a parameter structure.
 * 
 * @param settingFile Path to the settings file.
 * @param para Structure to store the parsed parameters.
 */
void loadSettingPara(string settingFile, SettingPara &para)
{
    // event file
    ifstream fSetting;
    fSetting.open(settingFile);

    if (!fSetting)
    {
        cout << "failed to open" << settingFile << endl;
        throw std::runtime_error("Wrong Lidar para file");
        exit(0);
    }
    string str1, str2;

    while (fSetting >> str1 >> str2)
    {

        if (!str1.compare("-flagFeature"))
        {
            para.bFeatureBased = stoi(str2);
        }
        if (!str1.compare("-flagLevelInfo"))
        {
            para.bLevel = stoi(str2);
        }
        if (!str1.compare("-flagPoint"))
        {
            para.bPointBased = stoi(str2);
        }
        if (!str1.compare("-flagDistortion"))
        {
            para.bDistortion = stoi(str2);
        }
        if (!str1.compare("-flagSimul"))
        {
            para.bFPSimultanesouly = stoi(str2);
        }
        if (!str1.compare("-initScanIndex"))
        {
            para.initScan = stoi(str2);
        }
        if (!str1.compare("-endScanIndex"))
        {
            para.endScan = stoi(str2);
        }

        if (!str1.compare("-channel"))
        {
            para.nChannel = stoi(str2);
        }

        if (!str1.compare("-minRange"))
        {
            para.minRange = stod(str2);
        }
        if (!str1.compare("-maxRange"))
        {
            para.maxRange = stod(str2);
        }
        if (!str1.compare("-treeAngleRange"))
        {
            para.treeAngleThreshold = stod(str2);
        }
        if (!str1.compare("-maxTreeInitHeight"))
        {
            para.groundBufferTree = stod(str2);
        }
        if (!str1.compare("-segmentLength"))
        {
            para.surfaceLength = stoi(str2);
        }
        if (!str1.compare("-segmentMinLength"))
        {
            para.minSegLength = stoi(str2);
        }
        if (!str1.compare("-maxSegDistance"))
        {
            para.maxTreeSegDistance = stod(str2);
        }
        if (!str1.compare("-minTreeDistance"))
        {
            para.minTreeDistance = stod(str2);
        }
        if (!str1.compare("-minNumTreeSeg"))
        {
            para.minNumSegTree = stoi(str2);
        }
        if (!str1.compare("-minTreePair"))
        {
            para.minTreePair = stoi(str2);
        }
        if(!str1.compare("-groundDownsampleDistanceScan"))
        {
            para.ground_downsample_distance_per_scan = stod(str2);
        }

        if (!str1.compare("-flagOdoIntegration"))
        {
            para.bOdometry = stoi(str2);
        }

        if (!str1.compare("-numInitIntegratedScan"))
        {
            para.numInitIntegratedScan = stoi(str2);
        }
        if (!str1.compare("-numIntegratedScan"))
        {
            para.numIntegratedScan = stoi(str2);
        }
        if(!str1.compare("-groundDownsampleDistanceScanIntegration"))
        {
            para.ground_downsample_distance_scan_integration = stod(str2);
        }

        if(!str1.compare("-numIscanLoopClosure"))
        {
            para.loop_closure_iscan_num = stoi(str2);
        }

        if(!str1.compare("-patch_size_lc"))
        {
            para.patch_size_lc = stod(str2);
        }


        if (!str1.compare("-globalGroundMap"))
        {
            para.global_map_flag = true;
            para.global_DTM = true;
            para.global_ground_map_pass = input_folder + str2;
        }

        //(TODO:dtm only)
        if (!str1.compare("-globalTreeMap"))
        {
            para.global_map_flag = true;
            para.global_tree_locations = true;
            para.global_tree_map_pass = input_folder + str2;
        }

        if(!str1.compare("-export_inter_result_odo"))
        {
            para.export_intermediate_result_odometry = stoi(str2);
        }

        if(!str1.compare("-export_inter_result_map"))
        {
            para.export_intermediate_result_mapping = stoi(str2);
        }

        if(!str1.compare("-odometryFromTrajectory"))
        {
            para.odo_from_trajectory_flag = stoi(str2);
        }

        if(!str1.compare("-input_intensity_flag"))
        {
            para.intensity_flag = stoi(str2);
        }

        if(!str1.compare("-iscan_tree_pts_per_scan"))
        {
            para.minNumTreePointPerScan = stod(str2);
        }
        if(!str1.compare("-iscan_pts_per_tree"))
        {
            para.minNumPtsPerIscanTree = stod(str2);
        }
        if(!str1.compare("-iscan_min_tree"))
        {
            para.minimum_tree_iscan = stoi(str2);
        }
        if(!str1.compare("-iscan_opt_flag"))
        {
            para.iscan_opt_flag = stoi(str2);
        }
        if(!str1.compare("-lc_iter"))
        {
            para.iscan_lc_iterative_flag = stoi(str2);
        }
        if(!str1.compare("-lc_end_iter"))
        {
            para.perform_iscan_lc_at_last_flag = stoi(str2);
        }
        if(!str1.compare("-raw_scan_lc_iterations"))
        {
            para.number_iterations_raw_scan_lc = stoi(str2);
        }
        if(!str1.compare("-reoptimize_tree_feature"))
        {
            para.reoptimize_fetch_tree_points = stoi(str2);
        }
        if(!str1.compare("-reoptimize_raw_point"))
        {
            para.reoptimize_fetch_raw_points = stoi(str2);
        }

        if(!str1.compare("-reoptimize_iterations"))
        {
            para.reoptimize_iterations = stoi(str2);
        }
        if(!str1.compare("-reoptimize_distance_threshold"))
        {
            para.distance_threshold_fetch_points = stod(str2);
        }
    }

    if(para.global_DTM && (!para.global_tree_locations))
    {
        para.global_DTM_only = true;
    }

    std::string outSetting = output_folder + "SettingInfo.txt";
    ofstream fSettingOut;

    fSettingOut.open(outSetting, std::ifstream::out);
    fSettingOut << fixed << std::setprecision(4);

    fSettingOut << "------------- General Setting ----------------" << endl;
    fSettingOut << "Start and End Scan Index: " << para.initScan << " -> " << para.endScan << endl;
    fSettingOut << endl;

    fSettingOut << "------------- Odometry Setting ----------------" << endl;
    fSettingOut << "Distance for ground point downsampling: " << para.ground_downsample_distance_per_scan << " m" << endl;
    fSettingOut << endl;

    fSettingOut << "------------- Mapping Setting ----------------" << endl;
    if (para.global_map_flag)
    {
        fSettingOut << "Reference point cloud is included:" << endl;
        if (!(para.global_ground_map_pass.empty()))
            fSettingOut << "Global ground map : " << para.global_ground_map_pass << endl;
        if (!(para.global_tree_map_pass.empty()))
            fSettingOut << "Global tree map : " << para.global_tree_map_pass << endl;
        fSettingOut << "Accuracy of 2d location of the global tree: " << para.global_tree_loc_std << endl;
    }

    if (para.odo_from_trajectory_flag)
    {
        fSettingOut << "Trajectory information is used in the mapping thread" << endl;
        fSettingOut << "\tPositional/rotational accuracy from trajectory: " << para.odoPositionStdTraj << " m, " << para.odoOrientationStdTraj << " deg" << endl;
    }

    fSettingOut << "Scan integration:" << endl;
    if (!para.global_map_flag)
    {
        fSettingOut << "\tNumber of scans for initialization integration : " << para.numInitIntegratedScan << endl;
    }
    if(!para.iscan_opt_flag)
    {
        fSettingOut << "\tPose are not refinement in Iscan" <<endl;
    }
    fSettingOut << "\tNumber of scans for integration : " << para.numIntegratedScan << endl;
    fSettingOut << "\tWhether including odometry result in scan itegration: " << para.bOdometry << endl;
    if (para.bOdometry)
    {
        fSettingOut << "\t\tPositional/rotational accuracy: " << para.odoPositionStd << " m, " << para.odoOrientationStd << " deg" << endl;
    }
    fSettingOut << "\tTree/ground points expected accuarcy in scan integration: " << para.iscanTreeStd << " / " << para.iscanGroundStd << " m" << endl;
    fSettingOut << "\tNumber of minimum points per scan for tree to be included in optimization: " << para.minNumTreePointPerScan <<endl;
    fSettingOut << "\tNumber of minimum points for a tree in Iscan: " << para.minNumPtsPerIscanTree << endl;
    fSettingOut << "\tDistance for ground point downsampling in scan integration: " << para.ground_downsample_distance_scan_integration << endl;
    fSettingOut << "\tMinimum number of trees per scan or scan to map matching: " << para.minimum_tree_iscan << endl;


    fSettingOut << "Iscan to map estimation" << endl;
    fSettingOut << "\tTree/ground points expected accuarcy in iscan to mapping: " << para.mapTreeStd << " / " << para.mapGroundStd << " m" << endl;
    fSettingOut << "\tEstimating tree radius: " << para.bMapTreeRadius << endl;
    fSettingOut << "\tEstimating tree direction: " << para.bMapTreeDirec << endl;

    fSettingOut << "Loop closure" << endl;
    fSettingOut << "\tNumber of Iscan for loop closure: " << para.loop_closure_iscan_num << endl;
    fSettingOut << "\tPatch size for loop closure: " << para.patch_size_lc << " m" << endl;
    if(para.iscan_lc_iterative_flag)
        fSettingOut << "\tConduct iterative iscan level LC until number of merge tree is few"<< endl;
    if(para.perform_iscan_lc_at_last_flag)
        fSettingOut << "\tConduct iscan level LC at last before raw scan level"<< endl;
    fSettingOut << "\tNumber of iterations for final raw scan level LC: " << para.number_iterations_raw_scan_lc << " m" << endl;

    if(para.reoptimize_fetch_tree_points || para.reoptimize_fetch_raw_points)
    {
        if(para.reoptimize_fetch_tree_points)
            fSettingOut << "Fetch tree points from tree features of individual scan and reoptimize" << endl;
        if(para.reoptimize_fetch_raw_points)
            fSettingOut << "Fetch tree points from non-ground points of individual scan and reoptimize" << endl;
        
        fSettingOut << "\tNumber of iterations for reoptimization: " << para.reoptimize_iterations << endl;
        fSettingOut << "\tDistance threshold for finding neighboring trees: " << para.distance_threshold_fetch_points << " m" << endl;      
    }
        

    fSettingOut.close();
}

//*********************************************************
// Function to find fixed parameter to ceres
// input:
//*********************************************************
/**
 * @brief Checks and applies fixed and weighted parameter constraints in a Ceres optimization problem.
 * 
 * @param problem Pointer to the Ceres problem.
 * @param X Parameters to optimize.
 * @param stdX Standard deviations of parameters.
 * @param blocksize Size of each parameter block.
 * @param nblocks Number of blocks.
 */
void Checkforconstantparams(ceres::Problem *problem, double *X, double *stdX, const int blocksize, const int nblocks)
{

    vector<int> constindices;
    int n_fixed_params = 0;
    int n_weighted_params = 0;

    ceres::SubsetParameterization *subsetparam;
    for (int i = 0; i < nblocks; i++)
    {
        double *block = X + i * blocksize; // first element
        double *stdblock = stdX + i * blocksize;

        int n_weighted_param_inblock = 0;
        for (int j = 0; j < blocksize; j++)
        {

            if (stdblock[j] <= SigmaFixed) // fixed
            {

                constindices.push_back(j);
                n_fixed_params++;
            }
            else if (stdblock[j] > SigmaFixed && stdblock[j] < SigmaFree) // weighted
            {
                n_weighted_params++;
                n_weighted_param_inblock++;
            }
        }

        // n_fixed_params += constindices.size();
        if (constindices.size() == blocksize) // all of the block is constant
        {
            problem->SetParameterBlockConstant(block);
        }
        else if (constindices.size() > 0 && constindices.size() < blocksize) // a subset of block is constant
        {
            subsetparam = new ceres::SubsetParameterization(blocksize, constindices);
            problem->SetParameterization(block, subsetparam);
        }
        // constindices.clear();
        vector<int>().swap(constindices);
    }
    // cout << "Number of unknowns: " << blocksize*nblocks << ", number of fixed: " << n_fixed_params << ", number of weighted: " << n_weighted_params << endl;
}

/**
 * @brief Computes rotation matrix from Euler angles.
 */
void Compute_Rotation(double om, double phi, double kap, Eigen::Matrix3d &R)
{
	// radius
	double c1, s1, c2, s2, c3, s3;

	c1 = cos(om);
	s1 = sin(om);
	c2 = cos(phi);
	s2 = sin(phi);
	c3 = cos(kap);
	s3 = sin(kap);

	R << c2 * c3, -c2 * s3, s2,
		c1 * s3 + s1 * s2 * c3, c1 * c3 - s1 * s2 * s3, -s1 * c2,
		s1 * s3 - c1 * s2 * c3, s1 * c3 + c1 * s2 * s3, c1 * c2;
}

Eigen::Matrix3d Compute_Rotation(double om, double phi, double kap)
{
	// radius
	double c1, s1, c2, s2, c3, s3;

	c1 = cos(om);
	s1 = sin(om);
	c2 = cos(phi);
	s2 = sin(phi);
	c3 = cos(kap);
	s3 = sin(kap);

	Eigen::Matrix3d R;
	R << c2 * c3, -c2 * s3, s2,
		c1 * s3 + s1 * s2 * c3, c1 * c3 - s1 * s2 * s3, -s1 * c2,
		s1 * s3 - c1 * s2 * c3, s1 * c3 + c1 * s2 * s3, c1 * c2;

	return R;
}

Eigen::Matrix3f Compute_Rotation(float om, float phi, float kap)
{
	// radius
	float c1, s1, c2, s2, c3, s3;

	c1 = cos(om);
	s1 = sin(om);
	c2 = cos(phi);
	s2 = sin(phi);
	c3 = cos(kap);
	s3 = sin(kap);

	Eigen::Matrix3f R;
	R << c2 * c3, -c2 * s3, s2,
		c1 * s3 + s1 * s2 * c3, c1 * c3 - s1 * s2 * s3, -s1 * c2,
		s1 * s3 - c1 * s2 * c3, s1 * c3 + c1 * s2 * s3, c1 * c2;

	return R;
}
/**
 * @brief Recovers Euler angles from a rotation matrix.
 */
void Find_Rotation(Eigen::Matrix3d R, double &ome, double &phi, double &kap)
{
	// radius
	double r11, r12, r13, r23, r33;
	r11 = R(0, 0);
	r12 = R(0, 1);
	r13 = R(0, 2);
	r23 = R(1, 2);
	r33 = R(2, 2);

	phi = asin(r13);
	kap = atan2(-r12 / cos(phi), r11 / cos(phi));
	ome = atan2(-r23 / cos(phi), r33 / cos(phi));
}

Eigen::Vector3d Find_Rotation(Eigen::Matrix3d R)
{
	// radius
	double r11, r12, r13, r23, r33;
	r11 = R(0, 0);
	r12 = R(0, 1);
	r13 = R(0, 2);
	r23 = R(1, 2);
	r33 = R(2, 2);

	double phi = asin(r13);
	double kap = atan2(-r12 / cos(phi), r11 / cos(phi));
	double ome = atan2(-r23 / cos(phi), r33 / cos(phi));

	Eigen::Vector3d result{ome, phi, kap};
	return result;
}

Eigen::Vector3f Find_Rotation(Eigen::Matrix3f R)
{
	// radius
	float r11, r12, r13, r23, r33;
	r11 = R(0, 0);
	r12 = R(0, 1);
	r13 = R(0, 2);
	r23 = R(1, 2);
	r33 = R(2, 2);

	float phi = asin(r13);
	float kap = atan2(-r12 / cos(phi), r11 / cos(phi));
	float ome = atan2(-r23 / cos(phi), r33 / cos(phi));

	Eigen::Vector3f result{ome, phi, kap};
	return result;
}

/**
 * @brief Downsamples a point cloud using octree-based radius filtering.
 * 
 * @param original_pc Input point cloud.
 * @param distance Minimum allowed distance between retained points.
 * @param downsampeld_pc Output downsampled cloud.
 */
void downsamplePointCloudDistance(pcl::PointCloud<PointType>::Ptr original_pc, double distance, pcl::PointCloud<PointType> &downsampeld_pc )
{

	//octree
    int num_original_pts = original_pc->size();
	float resolution = 0.05;
	pcl::octree::OctreePointCloudSearch<PointType> referenceCloudOctree(resolution);
	referenceCloudOctree.setInputCloud(original_pc);
	referenceCloudOctree.addPointsFromInputCloud();

	vector<int> flag(num_original_pts, 1);

	for (int nPts = 0; nPts < num_original_pts; nPts++) //for each point, find its neighboring
	{
		if (flag[nPts] == 0)
			continue;

		downsampeld_pc.push_back(original_pc->points[nPts]);

		PointType searchPt;
		vector<int> pointIdxNKNSearchNeighbors;
		vector<float> pointNKNSquaredDistanceNeighbors;

		searchPt.x = (float)((original_pc->points[nPts].x));
		searchPt.y = (float)((original_pc->points[nPts].y));
		searchPt.z = (float)((original_pc->points[nPts].z));

		int numNeighborPoints = 0;
		numNeighborPoints = referenceCloudOctree.radiusSearch(searchPt, distance , pointIdxNKNSearchNeighbors, pointNKNSquaredDistanceNeighbors);

		if (numNeighborPoints)
		{
			for (int nNeibor = 0; nNeibor < numNeighborPoints; nNeibor++)
			{
				flag[pointIdxNKNSearchNeighbors[nNeibor]] = 0; //remove the point
			}
		}
	}

}

/**
 * @brief Downsamples point cloud using brute-force radius rejection (more precise).
 * 
 * @param original_pc Input cloud.
 * @param distance Minimum allowed separation between points.
 * @return Downsampled point cloud.
 */
pcl::PointCloud<PointType> DownSamplePointCloudBasedOnDistance(pcl::PointCloud<PointType>::Ptr original_pc, double distance)
{
	pcl::PointCloud<PointType> downsampeld_pc;
	pcl::KdTreeFLANN<PointType> kdtree;
	kdtree.setInputCloud(original_pc);
	double radius = distance;
	int original_point_cloud_size = original_pc->size();
	std::vector<bool> seeds(original_point_cloud_size, true);
	std::cout << "before downsample: " << original_point_cloud_size << std::endl;
	for (int i = 0; i < original_point_cloud_size; i++)
	{
		if (seeds[i] == false)
		{
			continue;
		}
		PointType current_p = original_pc->points[i];
		std::vector<int> indices;
		std::vector<float> distances;
		kdtree.radiusSearch(current_p, radius, indices, distances);
		for (int idx = 0; idx < indices.size(); idx++)
		{
			seeds[indices[idx]] = false;
		}
		seeds[i] = true;
	}
	for (int i = 0; i < original_point_cloud_size; i++)
	{
		if (seeds[i] == false)
		{
			continue;
		}
		PointType current_p = original_pc->points[i];
		downsampeld_pc.push_back(current_p);
	}

	std::cout << "after downsample: " << downsampeld_pc.size() << std::endl;
	return downsampeld_pc;
}

// class TicToc
// {
//   public:
//     TicToc()
//     {
//         tic();
//     }

//     void tic()
//     {
//         start = std::chrono::system_clock::now();
//     }

//     double toc()
//     {
//         end = std::chrono::system_clock::now();
//         std::chrono::duration<double> elapsed_seconds = end - start;
//         return elapsed_seconds.count() * 1000;
//     }

//   private:
//     std::chrono::time_point<std::chrono::system_clock> start, end;
// };
