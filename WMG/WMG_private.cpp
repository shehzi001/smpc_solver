/** 
 * @file
 * @author Alexander Sherikov
 */


#include "WMG.h"
#include "footstep.h"




/**
 * @brief Returns index of the next SS.
 *
 * @param[in] type search for a footstep of certain type,
 *                 by default (FS_TYPE_AUTO) both left and right
 *                 are searched.
 *
 * @return index of the next SS.
 */
int WMG::getNextSS(const fs_type type)
{
    return (getNextSS (current_step_number, type));
}



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
 * @param[in] type search for a footstep of certain type,
 *                 by default (FS_TYPE_AUTO) both left and right
 *                 are searched.
 *
 * @return index of the previous SS.
 */
int WMG::getPrevSS(const fs_type type)
{
    return (getPrevSS (current_step_number, type));
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
 * @param[out] left_foot_pos 3x1 vector of coordinates [x y z] + angle (orientation in x,y plane)
 * @param[out] right_foot_pos 3x1 vector of coordinates [x y z] + angle (orientation in x,y plane)
 */
void WMG::getDSFeetPositions (
        double *left_foot_pos,
        double *right_foot_pos)
{
    int left_ind, right_ind;

    left_ind = getNextSS ();
    if (FS[left_ind].type == FS_TYPE_SS_L)
    {
        right_ind = getPrevSS ();
    }
    else
    {
        right_ind = left_ind;
        left_ind = getPrevSS ();
    }

    left_foot_pos[0] = FS[left_ind].x;
    left_foot_pos[1] = FS[left_ind].y;
    left_foot_pos[2] = 0.0;
    left_foot_pos[3] = FS[left_ind].angle;

    right_foot_pos[0] = FS[right_ind].x;
    right_foot_pos[1] = FS[right_ind].y;
    right_foot_pos[2] = 0.0;
    right_foot_pos[3] = FS[right_ind].angle;
}



/**
 * @brief Determine position and orientation of feet
 *
 * @param[in] loops_per_preview_iter number of control loops per preview iteration
 * @param[in] loops_in_current_preview number of control loops passed in the current 
 *                                     preview iteration
 * @param[out] left_foot_pos 3x1 vector of coordinates [x y z] + angle (orientation in x,y plane)
 * @param[out] right_foot_pos 3x1 vector of coordinates [x y z] + angle (orientation in x,y plane)
 */
void WMG::getSSFeetPositions (
        const int loops_per_preview_iter,
        const int loops_in_current_preview,
        double *left_foot_pos,
        double *right_foot_pos)
{
    double *swing_foot_pos, *ref_foot_pos;
    int next_swing_ind, prev_swing_ind;


    if (FS[current_step_number].type == FS_TYPE_SS_L)
    {
        ref_foot_pos = left_foot_pos;
        swing_foot_pos = right_foot_pos;

        prev_swing_ind = getPrevSS (FS_TYPE_SS_R);
        next_swing_ind = getNextSS (FS_TYPE_SS_R);
    }
    else
    {
        ref_foot_pos = right_foot_pos;
        swing_foot_pos = left_foot_pos;

        prev_swing_ind = getPrevSS (FS_TYPE_SS_L);
        next_swing_ind = getNextSS (FS_TYPE_SS_L);
    }

    ref_foot_pos[0] = FS[current_step_number].x;
    ref_foot_pos[1] = FS[current_step_number].y;
    ref_foot_pos[2] = 0.0;
    ref_foot_pos[3] = FS[current_step_number].angle;


    int num_iter_in_ss = FS[current_step_number].repeat_times;
    int num_iter_in_ss_passed = num_iter_in_ss - FS[current_step_number].repeat_counter;
    double theta =
        (double) (loops_per_preview_iter * num_iter_in_ss_passed + loops_in_current_preview) /
        (loops_per_preview_iter * num_iter_in_ss);


    double x[3] = {
        FS[prev_swing_ind].x,
        (FS[prev_swing_ind].x + FS[next_swing_ind].x)/2,
        FS[next_swing_ind].x};

    double b_coef = - (x[2]*x[2] - x[0]*x[0])/(x[2] - x[0]);
    double a = step_height / (x[1]*x[1] - x[0]*x[0] + b_coef*(x[1] - x[0]));
    double b = a * b_coef;
    double c = - a*x[0]*x[0] - b*x[0];


    swing_foot_pos[0] = (1-theta)*x[0] + theta*x[2]; // linear equation
    swing_foot_pos[1] = FS[next_swing_ind].y;
    swing_foot_pos[2] = a * swing_foot_pos[0] * swing_foot_pos[0] + b * swing_foot_pos[0] + c;
    swing_foot_pos[3] = FS[next_swing_ind].angle;
}



