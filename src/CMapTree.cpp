#include "../header/CMapTree.h"

/*****************************************************************
 * Compute distance from point p to the cylinder
 * ***************************************************************/

/**
 * @brief Computes the shortest distance from a 3D point to the surface of a cylinder.
 *
 * Internally, the function computes the perpendicular distance from point `p` to the cylinder axis,
 * and then subtracts the cylinder radius to return the surface distance.
 *
 * @param[in] p The 3D point.
 * @param[in] paraVer The parameter version to use (VALID, CANDIDATE, FIT_CYLINDER).
 * @return The signed distance from the point to the cylinder surface.
 */

double MapTree::compute_dis(Eigen::Vector3d p, eParaVer paraVer)
{
    if (paraVer == VALID)
        return (computePoint2lineDistance(p, para.x, para.x + para.n) - para.r);
    else if (paraVer == CANDIDATE)
        return (computePoint2lineDistance(p, cand_para.x, cand_para.x + cand_para.n) - cand_para.r);
	else if (paraVer == FIT_CYLINDER)
		return (computePoint2lineDistance(p, fit_para_.x, fit_para_.x + fit_para_.n) - fit_para_.r);

    return 0.0;
}

/*****************************************************************
 * Compute horizontal angle from the point to the cylinder (apply R_trans)
 * ***************************************************************/

 /**
 * @brief Computes the horizontal angle from a 3D point to the cylinder axis after applying transformation.
 *
 * Projects the point onto the cylinder axis and computes the angle in the horizontal plane using atan2.
 * Rotation matrix used depends on the parameter version selected.
 *
 * @param[in] p The 3D point.
 * @param[in] paraVer The parameter version to use (VALID, CANDIDATE, FIT_CYLINDER).
 * @return The horizontal angle in radians.
 */

double MapTree::compute_angle(Eigen::Vector3d p, eParaVer paraVer)
{
    Eigen::Vector3d v_trans;
    if (paraVer == VALID)
    {
        // projection on the axis
        Eigen::Vector3d p_projection = computePoint2LineProjection(p, para.x, para.x + para.n);
        // vector from the projected point to p
        v_trans = R_trans * (p - p_projection);
    }
    else if (paraVer == CANDIDATE)
    {
        // projection on the axis
        Eigen::Vector3d p_projection = computePoint2LineProjection(p, cand_para.x, cand_para.x + cand_para.n);
        // vector from the projected point to p
        v_trans = R_trans_cand * (p - p_projection);
    }
	else if (paraVer == FIT_CYLINDER)
	{
		// projection on the axis
		Eigen::Vector3d p_projection = computePoint2LineProjection(p, fit_para_.x, fit_para_.x + fit_para_.n);
		// vector from the projected point to p
		v_trans = R_trans_cand * (p - p_projection);
	}
    return atan2(v_trans(0), v_trans(1));
}


/*****************************************************************
 * Check if a parameter is valid or not 
 * Input: updatedPara
 * ***************************************************************/

/**
 * @brief Validates whether a given cylinder parameter is within acceptable limits.
 *
 * Checks normal vector components and radius for validity.
 *
 * @param[in] updatedPara The cylinder parameter to validate.
 * @return true if valid, false otherwise.
 */

bool MapTree::is_valid_para(CylinderPara updatedPara)
{
    if (abs(updatedPara.n(0)) > Max_Nx_Ny || abs(updatedPara.n(1)) > Max_Nx_Ny || updatedPara.r > Max_Radius)
        return false;
    else
    return true;

}

/*****************************************************************
 * Based on the para, derive the rotation matrix to level the axis
 * para -> R_trans
 * ***************************************************************/
/**
 * @brief Derives the rotation matrix to align the cylinder axis with the vertical (Z) axis.
 *
 * Based on the cylinder's normal vector, computes roll and pitch angles to produce
 * a transformation that levels the cylinder upright.
 *
 * @param[in] paraVer The parameter version to derive rotation for (VALID, CANDIDATE, FIT_CYLINDER).
 */

void MapTree::derive_trans_matrix(eParaVer paraVer)
{
    if (paraVer == VALID)
    {
        Eigen::Vector3d n_ = para.n.normalized();
        double phi = asin(n_(0));
        double ome = atan2(-n_(1), n_(2));
        Eigen::Matrix3d R;
        R = Compute_Rotation(ome, phi, double(0.0));
        R_trans = R.transpose();
    }
    else if (paraVer == CANDIDATE)
    {
        Eigen::Vector3d n_ = cand_para.n.normalized();
        double phi = asin(n_(0));
        double ome = atan2(-n_(1), n_(2));
        Eigen::Matrix3d R;
        R = Compute_Rotation(ome, phi, double(0.0));
        R_trans_cand = R.transpose();
    }
	else if (paraVer == FIT_CYLINDER)
	{
		Eigen::Vector3d n_ = fit_para_.n.normalized();
		double phi = asin(n_(0));
		double ome = atan2(-n_(1), n_(2));
		Eigen::Matrix3d R;
		R = Compute_Rotation(ome, phi, double(0.0));
		R_trans_cand = R.transpose();
	}

}


/*****************************************************************
 * check the status of a tree fitted from Iscan to Map optimization
 * ***************************************************************/
/**
 * @brief Determines the status of a tree based on geometry and quality criteria.
 *
 * Uses heuristics based on number of points, radius, and surface ratio to
 * classify the tree as ESTABLISHED, FITTED, SOLID, or TBD.
 */

void MapTree::is_established()
{  
    // //if the status is already established, return
    // if(status == ESTABLISHED || status == SOLID )
    // {
    //     return;
    // }

    // if enough number of points 
#ifdef FIT_CYLINDER_TREE
    if(final_fit_error_>0.1)
    {
        status = TBD;
        return;
    }
#endif
    if (numPoint > 1000)
    {
        // if radius is not too small
        if (para.r > 0.04)
        {
            // compute surface ratio
            if (surfaceRatio >= 0.6)
            {
                status = ESTABLISHED;
            }
            else
            {
                status = FITTED;
                //CHECK
                //ResetPara();
            }
        }
        else
        {
            status = SOLID;
        }
    }
    else
        {
            status = TBD;
           // ResetPara();
        }
}

/**
 * @brief Resets cylinder parameters to default based on the centroid.
 */
void MapTree::ResetPara()
{
    para.x = center;
    para.n(0) = 0.0;para.n(1) = 0.0;para.n(2) = 1.0;
    para.r = 0.0;
}

//once a new parameter is estimated for a map tree, if it's valid, update_flag-> true
/**
 * @brief Updates tree parameters if the new estimate is valid.
 *
 * Calls transformation matrix update as well.
 *
 * @param[in] updatedPara Newly estimated cylinder parameters.
 * @return true if update was successful, false if parameters were invalid.
 */

bool MapTree::UpdateTreeParam(CylinderPara updatedPara)
{
    //if the parameter not valid, return false
    if(!is_valid_para(updatedPara))
    {
        if(numPoint > 1000)
        {
            status = FITTED;
        }
        else
        {
            status = TBD;
        }
        return false;
    }
    para = updatedPara;
    derive_trans_matrix();
    //update_flag = true;
    return true;

    // // should wait until tree points are updated.
    // derive_trans_matrix();
    // compute_surface_ratio();
    // is_established();

}

//after vTreePointMapping is updated, update the attributes
//if parameters are updated, compute surface ratio, centroid, type
//otherwise only the centroid
/**
 * @brief Computes geometric attributes of the tree based on assigned points.
 *
 * Recomputes the centroid and, if parameters are stable, also updates surface ratio and status.
 */
void MapTree::ComputeAttributes()
{
    if (status < 0)
        return;
        
    ComputeCentroid();

    if (status > 0)
    {
        compute_surface_ratio();
        is_established();
    }
    
}

/**
 * @brief Merges another tree candidate into the current one.
 *
 * Transfers visibility info, updates centroid and point count, 
 * and updates tracking indices.
 *
 * @param[in] p_cand_tree Pointer to the tree being merged.
 */

void MapTree::AddMapTree(MapTree* p_cand_tree)
{
    //include the information from previous to current test
    visibleInfo.insert(visibleInfo.end(), p_cand_tree->visibleInfo.begin(), p_cand_tree->visibleInfo.end());
    center = (center * double(numPoint) + p_cand_tree->center * double(p_cand_tree->numPoint)) /(double(numPoint + p_cand_tree->numPoint));
    numPoint += p_cand_tree->numPoint;

    p_cand_tree->status = DEACTIVATE;
    p_cand_tree->update_flag = false;

    start_iscan_index_ = min(start_iscan_index_, p_cand_tree->start_iscan_index_);
    end_iscan_index_ = max(end_iscan_index_, p_cand_tree->end_iscan_index_);
    update_flag = true;
}

/**
 * @brief Updates the tree's classification status when new observations are added.
 *
 * Simplified rule: if newly observed but under threshold, it’s TBD. Otherwise, FITTED.
 */

void MapTree::update_status_new_obs()
{
    if(numPoint == 0)
    {
        status = DEACTIVATE;
    }
    else if (numPoint <= 1000)
    {
        status = TBD;
    }
    else
    {
        //for the ones with fitted, or solid or well estabilished, dont change
        if(status == INIT || status == TBD)
        {
            status = FITTED;
        }
    }
}
