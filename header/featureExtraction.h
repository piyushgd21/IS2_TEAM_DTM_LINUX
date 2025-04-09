#ifndef _FEATURE_EXTRACTION_H_
#define _FEATURE_EXTRACTION_H_

#include "../header/utility.h"
#include "../header/CTrajectory.h"
#include "../header/Mapping.h"
//#include "../header/optimization.h"

// pcl::PointCloud<PointType>::Ptr laserCloudIn;
class LidarScan;

class LidarPoint
{
public:
    // compute in computeAttribute();
    PointType point;
    float range;
    int channel;
    int firing;
    bool bValid; //if the range in within the threshold
    double angHori;
    double angVert;

    // for continuity
    Eigen::Vector3f vecNext; // vector from current to next
    Eigen::Vector3f vecPrev; // vector from previous to current
    int neighborPrevious = -1; //index of previous point
    int neighborNext = -1;  //index of next point

    // compute from smooth
    float smooth = -1.0; // -2, without enough neighboring points, -1, bValid = false, so the smooth is not computed

    // compute from segment extraction
    int segID = -1; //segment id in a scan line
    int seg_init_type = -1; // segment type, 1 for long, 2 for short, currently not used

    // cluster type, cluster id
    int clusterID = -1;
    int clusterType = -1; // 1 ground, 2 tree

    LidarPoint(PointType point_, int channel_, int firing_) : point(point_), channel(channel_), firing(firing_)
    {
        range = sqrt(point.x * point.x + point.y * point.y + point.z * point.z);
        angHori = -atan2(point.y, point.x);
        angVert = atan2(point.z, sqrt(point.x * point.x + point.y * point.y));
    }

    // temp test
    float smooth2 = -1.0; // used in computing smoothness based on extracted ground segments
    int point_type_id = -1;   // planar/edge point: 5, ground planar points, 6, less flat ground planar points in feature based approach
};

struct SubSegment
{
    int s;   // start
    int e;   // end
    int len; // length
    int type = -1;
};

// Cluster: a group of segment for tree or ground
class ClusterSegment
{
public:
    //-------- from the extraction ------------
    int type = -1;                         // 1 for ground, 2 for tree
    std::vector<pair<int, int>> vSegments; // vec of the segID
    int len = 0;

    //-------- from plane fitting --------
    Eigen::Vector3f vCenterPlane; // this is based on the plane points of a segment

    //-------- from levelling --------
    Eigen::Vector3f vCenterPlaneTrans;

    //------- Tree cluster -----------
    vector<Eigen::Vector3f> segCenters; // save the leveled center of segments
    Eigen::Vector3f vCenterTrans;       // center of a tree cluster (leveled)
    Eigen::Vector3f vNorm;              // derived from fitLine(), norm in levelled, or from end points
    bool valid = true;                  // used in removing clost tree cluster
    // Eigen::Vector3f
    float uDistance; // distance between two segments

    //-------- Delivery for mapping -------------
    std::vector<PointType> mvRawPoints;        // points
    std::vector<PointType> mvRawPointsLeveled; // points

    int featureId = -1; // tracked tree id if map available

    ClusterSegment(int type_) : type(type_)
    {
        vCenterTrans << 0.0, 0.0, 0.0;
        vNorm << 0.0, 0.0, 1.0; // for tree segment, normal
        uDistance = 0.0;
    }

    // add segment to a cluster
    void addSegment(pair<int, int> segInfo, Eigen::Vector3f center)
    {
        len++;
        vSegments.push_back(segInfo);
        segCenters.push_back(center);
        vCenterTrans = vCenterTrans + (center - vCenterTrans) / float(len);
        if (len >= 3)
            vNorm = (segCenters.back() - segCenters[0]).normalized();
    }

    void addSegment(pair<int, int> segInfo, Eigen::Vector3f center, float distance)
    {

        len++;
        vSegments.push_back(segInfo);
        segCenters.push_back(center);
        vCenterTrans = vCenterTrans + (center - vCenterTrans) / float(len);

        // update mean different
        uDistance = uDistance + (distance - uDistance) / float(len - 1);

        if (len >= 3)
            vNorm = (segCenters.back() - segCenters[0]).normalized();
    }

    //update vCenterTrans, length, segCenters, vSegments, vNorm
    void mergeCluster(ClusterSegment cluster_2)
    {
        vCenterTrans = vCenterTrans * len + cluster_2.vCenterTrans * cluster_2.len;
        len = len + cluster_2.len;
        vCenterTrans = vCenterTrans / len;
        
        segCenters.insert(segCenters.end(), cluster_2.segCenters.begin(), cluster_2.segCenters.end());
        vSegments.insert(vSegments.end(), cluster_2.vSegments.begin(), cluster_2.vSegments.end());

        Eigen::Vector3f lowest_center = segCenters[0], highest_center = segCenters[0];
        for (int n_seg = 0; n_seg < segCenters.size(); ++ n_seg)
        {
            if(segCenters[n_seg](2) > highest_center(2))
            {
                highest_center =    segCenters[n_seg];
            }
            if(segCenters[n_seg](2) < lowest_center(2))
            {
                lowest_center =    segCenters[n_seg];
            }
        }
        vNorm = (highest_center- lowest_center).normalized();

    }

    float compute_height_range()
    {
        Eigen::Vector3f lowest_center = segCenters[0], highest_center = segCenters[0];
        for (int n_seg = 0; n_seg < segCenters.size(); ++ n_seg)
        {
            if(segCenters[n_seg](2) > highest_center(2))
            {
                highest_center =    segCenters[n_seg];
            }
            if(segCenters[n_seg](2) < lowest_center(2))
            {
                lowest_center =    segCenters[n_seg];
            }
        }
        return highest_center(2) - lowest_center(2);
    }

    // for line based cluster:
    bool fitLine()
    {
        Eigen::Vector3f c;
        return (lineFitting(segCenters, vNorm, c));
    }
};

class LidarSegment
{
public:
    //----------- threshold -------------
    int Length_Surface = 20;
    // int Length_ShortSegment = 5;
    // float Ang_Threshold = deg2rad(40.0); // mean angle of the line seg, for class 3, currently not used
    int Min_Length = 2;

    //----------- initialize ---------------
    int channel;           // channel id
    LidarScan *pLidarscan; // pointer to the scan

    //----------- filling in lidar scan --------
    std::vector<int> vecIndex; // coresponding to pLidarscan

    //----------- computed attribute -----------
    std::size_t length;                     // number of point in the segment
    std::vector<Eigen::Vector3f> vecPoints; // list of all point
    std::vector<float> vecAngle;            // angle between successive 3 points
    float uAngle;                           // mean of the angles
    Eigen::Vector3f centerPoint;            // center point of all points
    Eigen::Vector3f direction;              // difference between first and last point
    // segment type: 1. Large flat area. 2. Short surface. //optional 3. ShortSegment with undefined objects
    int initType = 0; // based on the classify(), 1 for long, 2 shor short
    // int type = -1;    // 1 for ground segments after ground extraction
    // firing ID (beginning and end)
    int maxFiring;
    int minFiring;

    // ---------- for ground cluster, plane fitting ------------------
    // for a segment, get the smooth point
    float Smooth_Threshold = 0.08;
    std::vector<Eigen::Vector3f> vecPlanePoints; // list of all plane point (based on smooth threshold)
    std::size_t planeLength;                     // number of points
    Eigen::Vector3f planeCenterPoint;            // center point of all points, currently not used

    //----------- from deriving cluster --------
    //int clusterID = -1; 
    int ground_cluster_id = -1;
    int outlier_removal_flag = 0; //0: insufficient number of segment, 1 insufficient number of points, 2 fail in plane fitting, 3, fail in normal direction check 4 fail in checking height info, 10 valid
    int tree_cluster_id = -1;
    int tree_outlier_removal_flag = -1;  //0: insufficient number of segment, 1. fail in line fitting, 3. fail in height range  10 valid
    int final_feature_id = -1; // after removing cluster outliers, final id (ground extraction/tree extraction)
    int final_feature_type = 0;       // cluster type: 1. ground, 2 tree

    //----------  from tree extraction ----------
    Eigen::Vector3f centerPointTrans; // center point after transformation

    //-------------------- function ---------------------------------
    // initialization
    LidarSegment(int channel_, LidarScan *pLidarscan_) : channel(channel_), pLidarscan(pLidarscan_)
    {
        // direVec = Eigen::Vector3f::Zero();
    }
    LidarSegment(int channel_, int min_length, int length_threshold, LidarScan *pLidarscan_) : Length_Surface(length_threshold), Min_Length(min_length), channel(channel_), pLidarscan(pLidarscan_)
    {

        // direVec = Eigen::Vector3f::Zero();
    }

    void computeSegAttribute();

    void computePlaneAttribute(float thSmooth);

    // based on the length and mean angle
    void classify()
    {
        if (length > Length_Surface)
            initType = 1;
        else if (length <= Length_Surface) // && uAngle < Ang_Threshold)
            initType = 2;
        // else if (length <= Length_Surface && uAngle >= Ang_Threshold)
        //     type = 3;
    }

    // valid based on the length
    bool valid()
    {
        if (length < Min_Length)
            return false;
        else
            return true;
    }

    // temp  (not using)
    //  int N_Smooth = 4;
    //  float *pCurv; // curvature of all point. type = 1 only
    //  std::vector<SubSegment> vecSubsegment;
    //  Eigen::Vector3f direVec;
    //  void computeNormDirection();
};

class LidarScan
{
private:

    bool flag_export_intermediate_result = false;

    // ------------ from initialize ------------
    SettingPara mPara;
    int mNScan;                             // channel
    int mNFiring;                           // firing per scan
    size_t mNPoints;                        // number of points
    float mScanDuration;                    // duration of a scan, ms
    float mTimeTracking;                    // duration to previous scan, = mTimeEnd - mpPrev->mTimeEnd. theoretically: mTimeTracking = mScanDuration
    std::pair<float, float> mRangeTreshold; // for computing valid points
    // for computing continuity
    int mTolerateGap;                            // gap: 1 -> mTolerateGap
    std::vector<pair<float, float>> mvContRange; // based on the tolerate angle
    // lidar point cloud. Order: firing(channel 1-32)
    pcl::PointCloud<PointType>::Ptr mpLaserCloud;
    std::ofstream fInfo; // output
    // expected transformation -> r_lut2_lut1, R_lut2_lut1
    Eigen::Quaterniond q_ini = Eigen::Quaterniond::Identity(); // for eop (without considering the rotation)
    Eigen::Vector3d t_ini = Eigen::Vector3d::Zero();
    Eigen::Quaterniond q_ini_level = Eigen::Quaterniond::Identity(); // for eop (without considering the rotation)
    Eigen::Vector3d t_ini_level = Eigen::Vector3d::Zero();

    // ------------ computing attribute ------------
    // smoothness
    int N_Smooth = 4; // half of it
    float Smooth_Threshold = 0.08;
    // segments
    int mMinSegLength = 2;
    int mSegLength = 20;
    std::vector<std::vector<LidarSegment>> mvSegments; // scanline ID, segment ID

    //--------------------- Ground extraction------------------
    int MinNumSegGround = 3;
    float mGroundBuffer = 0.5;             // threshold for removing outliers
    std::vector<ClusterSegment> ground_cluster_; // ground clusters..
    Eigen::Matrix3f R_lu_lup;              // levelling: mR_lu_lulevel * r_i_lu = r_i_lu with flat point cloud
    float mGroundH;                        // ground height of the leveling scan
    float mOme, mPhi;

    //--------------------- Tree extraction------------------
    int MinNumSegTree = 3;
    float MaxTreeSegDistance = 1.0;
    float MinTreeDistance = 2.0;       // for removing trees close to each other
    float TreeAngThreshold;            // deg2rad(20.0);    // mean angle of the line seg
    int ChannelContinous = 2;          // the most tolerate gap in channel of a tree cluster
    float mGroundBufferTreeInit = 2.5; // threshold for height of segment to initialize a tree cluster
    // std::vector<Eigen::Vector3f> mTreeLoc;     // XYZ of tree locations
    std::vector<ClusterSegment> mvClusterTree; // tree clusters


    //-------------------- Functions ---------------------
    void load_raw_data(const string pass);
    void init();
    void computeAttribute();
    bool conduct_feature_based_odo();
    bool conduct_point_based_odo();
    bool extractGroundCluster(bool bFinal);

    bool extractTreeCluster();
    bool extract_ground_segments(bool ground_plane_model_flag);

    void extractGround();
    bool matchTrees2d(bool bSearch);
    bool matchTrees3d(bool bSearch);
    void computeTransformation();
    void computeLevelCenters();

    void fetch_previous_ground_info();
    void update_previous_ground_info();

public:
    LidarScan(){};
	~LidarScan();/*{
	 mpSurfPointsFlatGround.reset();
	 pSurfPointsLessFlatGround_ds_.reset();
	 pSurfPointsLessFlatGround_.reset();
	 mpCornerPointsSharp.reset();
	 mpCornerPointsLessSharp.reset();
	 mpSurfPointsFlat.reset();
	 mpSurfPointsLessFlat.reset();
	};*/
    LidarScan(const std::string lidar_data_path, const int index, LidarScan *pPrev, int pid, double timetag, double timeduration,
              SettingPara setting, std::string single_scan_tree_file, std::string single_scan_ground_file,int points_number,
		CTrajectory *pTraj = NULL, Mapping *pMap = NULL);

    // flag for successful tracking: control flow
    //bool mbTracked = true;          // tracked from feature based or point based
    //bool mbTrackFeature = true;     // false: 1) zero ground cluster, 2) tree cluster < 3, 3) number of tree macthes < 3
    //bool mbExtractedFeature = true; // mbExtractedFeature false -> mbTrackFeature false. But mbExtractedFeature true, mbTrackFeature false, mbTrackPoint = true happens.
    //1. Feature-based
    bool valid_ground_info_flag_ = false;      // if ground plane information is valid: mGroundH and R_lu_lup. If invalid, R_lu_lup invalid
    bool valid_feature_flag_ = false;   //if ground info and tree matches are extracted
    bool valid_feature_based_odo_flag_ = false; //if feature based odometry is successful 

    //2. Point-based
    bool point_extraction_flag_ = false; //if the points are extracted from this scan.
    bool valid_point_flag_ = false; //if sufficient points are extracted from this scan.
    bool valid_point_based_odo_flag_ = false; //if point based odometry is successful 
    
    //3. Overall    
    bool valid_odo_flag_ = false; //if any of the odometry is successful

    //bool mbTrackPoint = true;
    //bool mbPointExtraction = false;

    bool mbInit = false; // first scan of a new tracking
    bool mbTraj = false; // for exporting reference
    bool mbMap = false;  // mapping or not

    // setting for odometry
    bool mbSimultaneously = false; // conduct both odometry regardless of failure
    bool mbPointBased = true;      //pb
    bool mLidarDistortion = false; //pb
    bool mbFeatureBased = true;
    bool mbLevel = false; // with leveling information or not, true for yes in cases when ground is not pronounced, false for estimating from the scan


    // ------------ from input ------------
    int mPartID = -2; // id of tracked part
    // double mTime;     // timetag of the scan, ms @@@
    double mTimeInit; // beginning of a scan, ms, theoretically = mpPrev->mTimeEnd
    double mTimeEnd;  // end of a scan, ms, = mTimeInit + mScanDuration
    int id;           // scan ID
    // pointer to other classes
    LidarScan *mpPrev;                // previous scan
    CTrajectory *mpTraj;              // reference trajectory
	// CTrajectory *mGNSSINS_Traj;              // GNSS_INS trajectory
    Mapping *mpMap;                   // map thread
    std::vector<LidarPoint> mvPoints; // lidar point cloud. Order: firing(channel 1-32)

    // ------------ from compute transformation ----------
    // t1: end of mpTraj. t2: end of this scan
    // final result for odometry
    Eigen::Vector3d r_lut2_lut1 = Eigen::Vector3d::Zero();
    Eigen::Matrix3d R_lut2_lut1 = Eigen::Matrix3d::Identity();
    // final accumulated trajectory from odometry, current EOP: R_lu_m & r_lu_m;
    Eigen::Matrix3d R_lu_m = Eigen::Matrix3d::Identity();
    Eigen::Vector3d r_lu_m = Eigen::Vector3d::Zero();
    Eigen::Isometry3d T_lu_m = Eigen::Isometry3d::Identity(); // Transformation 3*4
    EOP mEop;

    // valid only when mpTraj
    float ground_height_ref; //nominal ground height
    Eigen::Vector3d r_lu_b_ref;
    Eigen::Matrix3d R_lu_b_ref;
    Eigen::Vector3d r_lut2_lut1_ref = Eigen::Vector3d::Zero();
    Eigen::Matrix3d R_lut2_lut1_ref = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d R_lu_m_ref; // with same origin as the computed
    Eigen::Vector3d r_lu_m_ref;
    //(NEW 6.18)
    Eigen::Vector3d r_lu_m_t_ini_;
    Eigen::Matrix3d R_lu_m_t_ini_;
    std::vector<PointType> GetGroundFeatures(std::string single_scan_ground_features_file, Eigen::Matrix3f R_lu_lup, double& ground_height);
	std::vector<PointType> GetGroundFeaturesLoadOnce(std::string single_scan_ground_features_file, Eigen::Matrix3f R_lu_lup, double& ground_height);

    // ----------- functions -----------------


    // ----------- for extracting planar point in extracted plane segments
    bool extract_planar_points_from_ground_segment();
    pcl::PointCloud<PointType>::Ptr mpSurfPointsFlatGround;  //****saved point is in the leveled frame 
    pcl::PointCloud<PointType>::Ptr pSurfPointsLessFlatGround_ds_; //saved in the lu frame, downsampled for each scan line
    pcl::PointCloud<PointType>::Ptr pSurfPointsLessFlatGround_; //saved in the lu frame


    //------------ for tree matching in 3d
    std::vector<pair<int, int>> mVecTreePairs; //<prev, cur> order of current
    int MinTreePair;

    // ----------- for computing odo based on tree and ground planar points
    bool compute2dSimilarityLeveled();

    bool computeOdometryLeveled();
    Eigen::Vector3d r_lut2_lut1_gt; // t2: current
    Eigen::Matrix3d R_lut2_lut1_gt;

    // Eigen::Vector3f r_blut2_blut1;
    // Eigen::Matrix3f R_blut2_blut1;
    // Eigen::Vector3d r_blut2_blut1_ref;
    // Eigen::Matrix3d R_blut2_blut1_ref;

    // -----------test
    // float *mpCurv; // curvature array of all point. type = 1 only

    // 1. compute smooth
    int MaxSubArea = 6;
    int MinSubPoints = 80;
    int NumSubEdge = 2;
    int NumSubEdgeLess = 20;
    int NumSubPlanar = 4;
    pcl::PointCloud<PointType>::Ptr mpCornerPointsSharp;
    pcl::PointCloud<PointType>::Ptr mpCornerPointsLessSharp;
    pcl::PointCloud<PointType>::Ptr mpSurfPointsFlat;
    pcl::PointCloud<PointType>::Ptr mpSurfPointsLessFlat;
    bool extract_edge_planar_points();

    // 2. odo
    float DISTANCE_SQ_THRESHOLD = 25.0;
    int NEARBY_SCAN = 2;
    bool compute_point_based_odometry();
    Eigen::Vector3d r_lut2_lut1_pb; // t2: current
    Eigen::Matrix3d R_lut2_lut1_pb;
    Eigen::Quaterniond qTem = Eigen::Quaterniond::Identity();
    Eigen::Vector3d tTem = Eigen::Vector3d::Zero();
    // undistort lidar point to previous scan
    void TransformToStart(PointType const *const pi, PointType *const po)
    {
        double s;
        if (mLidarDistortion)
            s = (pi->intensity - int(pi->intensity));
        else
            s = 1.0;
        Eigen::Quaterniond q_point_last = Eigen::Quaterniond::Identity().slerp(s, qTem);
        Eigen::Vector3d t_point_last = s * tTem;
        Eigen::Vector3d point(pi->x, pi->y, pi->z);
        Eigen::Vector3d un_point = q_point_last * point + t_point_last;

        po->x = un_point.x();
        po->y = un_point.y();
        po->z = un_point.z();
        po->intensity = pi->intensity;
    }
    void TransformToEnd(PointType const *const pi, PointType *const po)
    {
        // undistort point first
        pcl::PointXYZI un_point_tmp;
        TransformToStart(pi, &un_point_tmp);

        Eigen::Vector3d un_point(un_point_tmp.x, un_point_tmp.y, un_point_tmp.z);
        Eigen::Vector3d point_end = qTem.inverse() * (un_point - tTem);

        po->x = point_end.x();
        po->y = point_end.y();
        po->z = point_end.z();

        // Remove distortion time info
        po->intensity = int(pi->intensity);
    }
    void TransformToLevel(PointType const *const pi, PointType *const po)
    {
        Eigen::Vector3f point(pi->x, pi->y, pi->z);
        Eigen::Vector3f un_point = R_lu_lup * point;

        po->x = un_point.x();
        po->y = un_point.y();
        po->z = un_point.z();
        po->intensity = pi->intensity;
    }

    /*Derive t_ini_level & q_ini_level using the R_lu_lup also from previous*/
    void PredictLevelInit()
    {
        Eigen::Matrix3d R_lu_lup_t1 = mpPrev->R_lu_lup.cast<double>();
        Eigen::Matrix3d R_lu_lup_t2 = R_lu_lup.cast<double>();
        t_ini_level = R_lu_lup_t1 * t_ini;
        q_ini_level = R_lu_lup_t1 * q_ini * R_lu_lup_t2.inverse();
    }

    // 3. add tree to map
    void AddTreetoMap();
    void AddScantoMap();
	void AddScanToMap(std::string single_scan_features_file,std::string single_scan_ground_file, int single_scan_id);
	int single_scan_points_number_;
    
    // -----------  Export -----------------------
    void export_scan_point(const std::string outPass);
    void export_scan_segment(const std::string outPass, bool flag_level = false);
    void export_final_planar_points(const std::string outPass);







    void exportGroundPoints(const std::string outPass, bool flagLevel);
    void exportGroundPoints2(const std::string outPass, bool flagLevel);

    void exportFeature(const std::string outPass);
    void exportScan();
    void exportPointMap();

    void exportRefScan();
    // void plotMatrix(cv::Mat *matrix, const std::string outPass);
    void exportPoint(const std::string outPass);
    void exportPointChannel(const std::string outPass);
    void exportPointPerChannel(const std::string folderName);
    void exportSegmentPerChannel(const std::string folderName);
    void exportSegmentInfoPerChannel(const std::string folderName);
    void exportPointPerChannelTransformed(const std::string outPass);
    void exportPointChannelTransformed(const std::string outPass);
    void exportSegments(const std::string outPass);
    void exportGroundSegments(const std::string outPass);
    void exportClusterInfo(const std::string outpass);
    void exportTreeClusterInfo(const std::string outPass);
    void exportSegmentInfo(const std::string outPass);
    void exportSegmentTransInfo(const std::string outPass);
    void exportFinal(const std::string outPass);
    void exportFinalMapping(const std::string outPass);
};

#endif
