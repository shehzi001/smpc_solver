/** 
 * @file
 * @author Alexander Sherikov
 */

#include <cmath> // sqrt

#include "WMG.h"
#include "footstep.h"



/**
 * @brief Returns index of the next SS.
 *
 * @param[in] start_ind start search from this index.
 * @param[in] type search for a footstep of certain type,
 *                 by default (FS_TYPE_AUTO) both left and right
 *                 are searched.
 *
 * @return index of the next SS.
 */
int WMG::getNextSS(const int start_ind, const fs_type type)
{
    int index = start_ind + 1;
    for (; index < (int) FS.size(); index++) 
    {
        if (FS[index].type != FS_TYPE_DS)
        {
            if (type == FS_TYPE_AUTO)
            {
                break;
            }
            else
            {
                if (FS[index].type == type)
                {
                    break;
                }
            }
        }
    }
    return (index);
}



/**
 * @brief Returns index of the previous SS.
 *
 * @param[in] start_ind start search from this index.
 * @param[in] type search for a footstep of certain type,
 *                 by default (FS_TYPE_AUTO) both left and right
 *                 are searched.
 *
 * @return index of the previous SS.
 */
int WMG::getPrevSS(const int start_ind, const fs_type type)
{
    int index = start_ind - 1;
    for (; index >= 0; index--)
    {
        if (FS[index].type != FS_TYPE_DS)
        {
            if (type == FS_TYPE_AUTO)
            {
                break;
            }
            else
            {
                if (FS[index].type == type)
                {
                    break;
                }
            }
        }
    }
    return (index);
}



/**
 * @brief Determine position and orientation of feet in DS
 *
 * @param[in] support_number number of the support
 * @param[out] left_foot_pos 4x4 homogeneous matrix, which represents position and orientation
 * @param[out] right_foot_pos 4x4 homogeneous matrix, which represents position and orientation
 */
void WMG::getDSFeetPositions (
        const int support_number,
        double *left_foot_pos,
        double *right_foot_pos)
{
    int left_ind, right_ind;

    left_ind = getNextSS (support_number);
    if (FS[left_ind].type == FS_TYPE_SS_L)
    {
        right_ind = getPrevSS (support_number);
    }
    else
    {
        right_ind = left_ind;
        left_ind = getPrevSS (support_number);
    }

    Matrix4d::Map(left_foot_pos) = FS[left_ind].posture.matrix();
    Matrix4d::Map(right_foot_pos) = FS[right_ind].posture.matrix();
}



/**
 * @brief Determine position and orientation of feet
 *
 * @param[in] support_number number of the support
 * @param[in] theta a number between 0 and 1, a fraction of support time that have passed 
 * @param[out] left_foot_pos 4x4 homogeneous matrix, which represents position and orientation
 * @param[out] right_foot_pos 4x4 homogeneous matrix, which represents position and orientation
 */
void WMG::getSSFeetPositions (
        const int support_number,
        const double theta,
        double *left_foot_pos,
        double *right_foot_pos)
{
    double *swing_foot_pos, *ref_foot_pos;
    int next_swing_ind, prev_swing_ind;
    footstep& current_step = FS[support_number];


    if (current_step.type == FS_TYPE_SS_L)
    {
        ref_foot_pos = left_foot_pos;
        swing_foot_pos = right_foot_pos;

        prev_swing_ind = getPrevSS (support_number, FS_TYPE_SS_R);
        next_swing_ind = getNextSS (support_number, FS_TYPE_SS_R);
    }
    else
    {
        ref_foot_pos = right_foot_pos;
        swing_foot_pos = left_foot_pos;

        prev_swing_ind = getPrevSS (support_number, FS_TYPE_SS_L);
        next_swing_ind = getNextSS (support_number, FS_TYPE_SS_L);
    }

    Matrix4d::Map(ref_foot_pos) = current_step.posture.matrix();


    double dx = FS[next_swing_ind].x() - FS[prev_swing_ind].x();
    double dy = FS[next_swing_ind].y() - FS[prev_swing_ind].y();
    double l = sqrt(dx*dx + dy*dy);


    double x[3] = {0.0, l/2, l};
    double b_coef = - (x[2]*x[2] /*- x[0]*x[0]*/)/(x[2] /*- x[0]*/);
    double a = step_height / (x[1]*x[1] /*- x[0]*x[0]*/ + b_coef*(x[1] /*- x[0]*/));
    double b = a * b_coef;
    //double c = - a*x[0]*x[0] - b*x[0];


    double dl = /*(1-theta)*x[0] +*/ theta * l;

    Matrix4d::Map(swing_foot_pos) = (
            FS[prev_swing_ind].posture
          * Translation<double, 3>(theta * dx, theta * dy, a*dl*dl + b*dl)
          * AngleAxisd(FS[next_swing_ind].angle - FS[prev_swing_ind].angle, Vector3d::UnitZ())
            ).matrix();
}



