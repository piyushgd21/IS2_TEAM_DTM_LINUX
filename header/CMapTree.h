#ifndef _MAPTREE_H_
#define _MAPTREE_H_

#include "../header/utility.h"
#include "../header/IntegratedScan.h"

class Mapping;

class MapTree
{
private:
    double Max_Nx_Ny = 1.0;
    double Max_Radius = 1.0;
    int Num_Bin = 20; //for computing the ratio
    double Min_Surface_Ratio = 0.6;
    int Min_Number_Pts_Per_Bin = 20;
    
    int Id = -1;
    Mapping *pMap; // pointer to map
    
public:
    // initilization: tree from global map
    MapTree(int id_, CylinderPara para_, Mapping *pMap_) : Id(id_),pMap(pMap_)
    {
        
        map_tree_flag = true;
        //(24/10/10) There should always be false when using DTM
        map_tree_flag = false;
        status = GLOBAL_INIT;
        //used to provide constraint in the X and Y 
        para_ref_ = para_;
        
        //shift x0, y0, z0 to avoid failure in 
        para = para_;
        para.x(0) = para.x(0) - 0.1;
        para.x(1) = para.x(1) - 0.1;
        para.x(2) = para.x(2) - 0.1;

        
        derive_trans_matrix();
    }

    // initialization: tree from SLAM
    MapTree(int id_, CylinderPara para_,pair<int, int> tree_info, Mapping *pMap_) : Id(id_),pMap(pMap_)
    {
        //lastValidPara = para;
        cand_para = para_;
        candTreeInfo = tree_info;
        status = INIT;
#ifdef FIT_CYLINDER_TREE
        //derive tree parameters for the new iscan tree

        //only if fitted & not poor state, add the tree
		CalculateIncomingPointsFitCylinderParameters();
		if (!fit_tree_failed_)
		{
			if (IsPoorStateTree())
			{
				poor_state_tree_ = true;
			}
			else
			{
				update();
			}
		}
#else
        update();
#endif

    }


    //check based on the estimated cylinder parameters
    bool poor_state_tree_ = false;
	bool fit_tree_failed_ = false;
	bool IsPoorStateTree()
	{
		return (fit_para_.r > 1 || final_fit_error_ > 0.1 || fit_para_.n(2) < 0.9);
	}


    // Input sensor
    enum eStatus
    {
        DEACTIVATE =-2,    // merged to other tree
        GLOBAL_INIT = -1, // initialize with map tree, without point
        INIT = 0,         // just initialized, not included in the iscan to map yet
        TBD = 1,          // included in the iscan to map, but not enough points yet
        FITTED = 2,       // enough points, but with large radius but small surface ratio
        ESTABLISHED = 3,  // well defined tree trunk with certain radius
        SOLID = 4        // well defined tree trunk with small radius
    };

    enum eParaVer{
        VALID=0, //reliable parameter
        CANDIDATE=1, //candidate parameters
		FIT_CYLINDER = 2
    };

    CylinderPara para;
    CylinderPara cand_para;
    CylinderPara para_ref_;
	CylinderPara fit_para_;
	double final_fit_error_ = 0;
	std::vector<double> CalculateIncomingPointsFitCylinderParameters();
	void CalculateAllPointsFitCylinderParameters();
   // CylinderPara lastValidPara; //backup of reliable tree para
    Eigen::Matrix3d R_trans; //rotation matrix to transform n to (0, 0, 1)
    Eigen::Matrix3d R_trans_cand; //rotation matrix to transform n to (0, 0, 1)


    int numPoint = 0; //number of tree points
    vector<pair<int,int>> visibleInfo; //(index of integrated scan, tree id)
    pair<int,int> candTreeInfo;

    Eigen::Vector3d center = Eigen::Vector3d::Zero(); //center in mapping for all points
    double surfaceRatio = -1.0; //whether the definition of cylinder is complete
    eStatus status = INIT;
    bool map_tree_flag = false;
    bool update_flag = false; //used for loop closure

    int start_iscan_index_ = -1;
    int end_iscan_index_ = -1;

    double rmse = -1.0; // rmse values of points in mapTree
    double candRmse = -1.0; //rmse values of points in Iscan tree
    vector<double> residual_distribution;

    Eigen::Vector3d ground_point;



    //in case add a Iscan tree, update tree points and tree parameter
    bool add_iscan_tree(CylinderPara updatedPara, pair<int,int> info);
    //in case no change in the tree points, purely the parameter
    bool UpdateTreeParam(CylinderPara updatedPara);
    //add another map tree to this map tree
    void AddMapTree(MapTree* p_cand_tree);

    
    void ResetPara();

    void update();
    
    bool is_valid_para(CylinderPara updatedPara); //check if a new parameter is valid


    //related to attributes
    void ComputeAttributes();
    void derive_trans_matrix(eParaVer paraVer = VALID);  // derive the trans matrix for computing surface ratio
    void compute_surface_ratio(eParaVer paraVer = VALID); //compute surface ratio
    void is_established(); // check the status of a tree based on the parameter and its surface ratio
    void ComputeCentroid(); // compute the center of current map tree
    void count_point();//count number of points

    //compute attributeds from other objects to current tree
    double compute_dis(Eigen::Vector3d p, eParaVer paraVer = VALID );
    double compute_angle(Eigen::Vector3d p, eParaVer paraVer = VALID);
    double compute_iscan_tree_residual(pair<int, int> tree_info, eParaVer paraVer = VALID);
    double compute_map_tree_residual( eParaVer paraVer = VALID);

    
    //relate to update the observation
    void clean_obs(){ vector<pair<int,int>>().swap(visibleInfo); }
    void update_status_new_obs();








    //bool update_para(CylinderPara updatedPara);

    //void reset_to_previous();



    //temp testing
    //int groupId = -1;
    //bool bVisit = false;


};
#endif