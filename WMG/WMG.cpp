/** 
 * @file
 * @author Alexander Sherikov
 * @date 27.09.2011 18:59:03 MSD
 */


#include <stdio.h>
#include <math.h> // cos, sin


#include "WMG.h"
#include "footstep.h"


/** 
 * @brief Initializes a WMG object.
 *
 * @param[in] N_ Number of sampling times in a preview window
 * @param[in] T_ Sampling time [ms.]
 * @param[in] step_height_ step height (for interpolation of feet movements) [meter]
 */
WMG::WMG (
        const unsigned int N_,
        const unsigned int T_, 
        const double step_height_)
{
    step_height = step_height_;
    N = N_;
    sampling_period = T_;

    T_ms = new unsigned int[N];
    for (unsigned int i = 0; i < N; i++)
    {
        T_ms[i] = 0;
    }
    

    current_step_number = 0;
    last_time_decrement = 0;
    first_preview_step = current_step_number;


    /// NAO constraint with safety margin.
    def_ss_constraint[0] = 0.09;
    def_ss_constraint[1] = 0.025;
    def_ss_constraint[2] = 0.03;
    def_ss_constraint[3] = 0.025;

    def_auto_ds_constraint[0] = 0.07;
    def_auto_ds_constraint[1] = 0.025;
    def_auto_ds_constraint[2] = 0.025;
    def_auto_ds_constraint[3] = 0.025;

    addstep_constraint[0] = def_ss_constraint[0];
    addstep_constraint[1] = def_ss_constraint[1];
    addstep_constraint[2] = def_ss_constraint[2];
    addstep_constraint[3] = def_ss_constraint[3];

    def_time_ms = 400;
    ds_time_ms = 0;
}


/** \brief Default destructor. */
WMG::~WMG()
{
    if (T_ms != NULL)
    {
        delete T_ms;
    }
}



/**
 * @brief Set default parameters of footsteps
 *
 * @param[in] def_num_ default time in ms (for SS or DS depending on the type of added footstep)
 * @param[in] ds_num_ number of DS to be generated automatically.
 * @param[in] constraints Vector of the PoS constraints (assumed to be [4 x 1], can be NULL).
 *
 * @note length of one automatically generated DS is assumed to be equal to sampling_period.
 */
void WMG::setFootstepDefaults(
        const unsigned int def_num_, 
        const unsigned int ds_num_,
        const double *constraints) 
{
    if (constraints != NULL)
    {
        addstep_constraint[0] = constraints[0];
        addstep_constraint[1] = constraints[1];
        addstep_constraint[2] = constraints[2];
        addstep_constraint[3] = constraints[3];
    }

    def_time_ms = def_num_ * sampling_period;
    ds_time_ms  = sampling_period; 
    ds_num      = ds_num_;
}



/**
 * @brief Set default parameters of footsteps
 *
 * @param[in] def_time_ms_ default time in ms (for SS or DS depending on the type of added footstep)
 * @param[in] ds_time_ms_ length of one automatically generated DS [ms.]
 * @param[in] ds_num_ number of DS to be generated automatically.
 * @param[in] constraints Vector of the PoS constraints (assumed to be [4 x 1], can be NULL).
 */
void WMG::setFootstepDefaults(
        const unsigned int def_time_ms_, 
        const unsigned int ds_time_ms_, 
        const unsigned int ds_num_,
        const double *constraints) 
{
    if (constraints != NULL)
    {
        addstep_constraint[0] = constraints[0];
        addstep_constraint[1] = constraints[1];
        addstep_constraint[2] = constraints[2];
        addstep_constraint[3] = constraints[3];
    }

    def_time_ms = def_time_ms_;
    ds_time_ms  = ds_time_ms_; 
    ds_num      = ds_num_;
}



/**
 * @brief Adds a footstep to FS.
 *
 * @param[in] x_relative x_relative X position [meter] relative to the previous footstep.
 * @param[in] y_relative y_relative Y position [meter] relative to the previous footstep.
 * @param[in] angle_relative angle_relative Angle [rad.] relative to the previous footstep.
 * @param[in] type (optional) type of the footstep.
 *
 * @note Coordinates and angle are treated as absolute for the first step in the preview window.
 */
void WMG::addFootstep(
        const double x_relative, 
        const double y_relative, 
        const double angle_relative, 
        fs_type type)
{
    Transform<double, 3>    posture (Translation<double, 3>(x_relative, y_relative, 0.0));
    Vector3d    zref_offset ((addstep_constraint[0] - addstep_constraint[2])/2, 0.0, 0.0);


    if (FS.size() == 0)
    {
        // this is the first ("virtual") step.
        if (type == FS_TYPE_AUTO)
        {
            type = FS_TYPE_DS;
        }

        posture *= AngleAxisd(angle_relative, Vector3d::UnitZ());
        // offset is the absolute position here.
        Vector3d zref_abs = posture * zref_offset;

        FS.push_back(
                footstep(
                    angle_relative, 
                    posture,
                    zref_abs,
                    def_time_ms, 
                    type,
                    addstep_constraint));
    }
    else
    {
        // determine type of the step
        if (type == FS_TYPE_AUTO)
        {
            switch (FS[getPrevSS(FS.size())].type)
            {
                case FS_TYPE_SS_L:
                    type = FS_TYPE_SS_R;
                    break;
                case FS_TYPE_SS_R:
                    type = FS_TYPE_SS_L;
                    break;
                case FS_TYPE_DS:
                default:
                    type = FS_TYPE_SS_R;
                    break;
            }
        }

        // Position of the next step
        posture = (*FS.back().posture) * posture * AngleAxisd(angle_relative, Vector3d::UnitZ());

        double prev_a = FS.back().angle;
        double next_a = prev_a + angle_relative;
        Vector3d next_zref = posture * zref_offset;


        // Add double support constraints that lie between the
        // newly added step and the previous step
        double theta = (double) 1/(ds_num + 1);
        double angle_shift = angle_relative * theta;
        double x_shift = theta*x_relative;
        double y_shift = theta*y_relative;
        Vector3d *ds_zref = &FS.back().ZMPref;
        for (unsigned int i = 0; i < ds_num; i++)
        {
            Transform<double, 3> ds_posture = (*FS.back().posture) 
                       * Translation<double, 3>(x_shift, y_shift, 0.0)
                       * AngleAxisd(angle_shift, Vector3d::UnitZ());

            if (i == ds_num / 2)
            {
                ds_zref = &next_zref;
            }

            FS.push_back(
                    footstep(
                        FS.back().angle + angle_shift,
                        ds_posture,
                        *ds_zref,
                        ds_time_ms, 
                        FS_TYPE_DS,
                        def_auto_ds_constraint));
        }


        // add the new step
        FS.push_back(
                footstep(
                    next_a, 
                    posture, 
                    next_zref,
                    def_time_ms, 
                    type,
                    addstep_constraint));
    }    
}



/**
 * @brief Determine position and orientation of feet
 *
 * @param[in] shift_from_current_ms a positive shift in time (ms.) from the current time
 *  (allows to get positions for the future supports)
 * @param[out] left_foot_pos 4x4 homogeneous matrix, which represents position and orientation
 * @param[out] right_foot_pos 4x4 homogeneous matrix, which represents position and orientation
 *
 * @attention This function requires the walking pattern to be started and finished
 * by single support.
 *
 * @attention Cannot be called on the first or last SS  =>  must be called after 
 * FormPreviewWindow().
 */
void WMG::getFeetPositions (
        const unsigned int shift_from_current_ms,
        double *left_foot_pos,
        double *right_foot_pos)
{
    unsigned int support_number = first_preview_step;
    // formPreviewWindow() have already decremented the time
    unsigned int step_time_left = FS[support_number].time_left + last_time_decrement;
    unsigned int shift_ms = shift_from_current_ms;


    while (shift_ms > step_time_left)
    {
        shift_ms -= step_time_left;
        ++support_number;
        if (support_number >= FS.size())
        {
            return;
        }
        step_time_left = FS[support_number].time_left;
    }


    if (FS[support_number].type == FS_TYPE_DS)
    {
        getDSFeetPositions (support_number, left_foot_pos, right_foot_pos);
    }
    else
    {
        double theta = (double) 
            ((FS[support_number].time_period - step_time_left) + shift_ms) 
            / FS[support_number].time_period;

        getSSFeetPositionsBezier (
                support_number,
                theta,
                left_foot_pos, 
                right_foot_pos);
    }
}



/**
 * @brief Checks if the support foot switch is needed.
 *
 * @return true if the support foot must be switched. 
 */
bool WMG::isSupportSwitchNeeded ()
{
    // current_step_number is the number of step, which will
    // be the first in the preview window, when formPreviewWindow()
    // is called.
    if (FS[current_step_number].type == FS_TYPE_DS)
    {
        return (false);
    }
    else // single support
    {
        if (// if we are not in the initial support
            (current_step_number != 0) &&
            // this is the first iteration in SS
            (FS[current_step_number].time_period == FS[current_step_number].time_left) &&
            // the previous SS was different
            (FS[getPrevSS(first_preview_step)].type != FS[current_step_number].type))
        {
            return (true);
        }
    }

    return (false);
}



/**
 * @brief Changes position of the next SS.
 *
 * @param[in] posture a 4x4 homogeneous matrix representing new position and orientation
 * @param[in] zero_z_coordinate set z coordinate to 0.0
 *
 * @todo DS must be adjusted as well.
 */
void WMG::changeNextSSPosition (const double* posture, const bool zero_z_coordinate)
{
    FS[getNextSS(first_preview_step)].changePosture(posture, zero_z_coordinate);
}



/**
 * @brief Forms a preview window.
 *
 * @return WMG_OK or WMG_HALT (simulation must be stopped)
 */
WMGret WMG::formPreviewWindow(smpc_parameters & par)
{
    WMGret retval = WMG_OK;
    unsigned int win_step_num = current_step_number;
    unsigned int step_time_left = FS[win_step_num].time_left;


    for (unsigned int i = 0; i < N;)
    {
        if (step_time_left > 0)
        {
            par.angle[i] = FS[win_step_num].angle;

            par.fp_x[i] = FS[win_step_num].x();
            par.fp_y[i] = FS[win_step_num].y();


            // ZMP reference coordinates
            par.zref_x[i] = FS[win_step_num].ZMPref.x();
            par.zref_y[i] = FS[win_step_num].ZMPref.y();


            par.lb[i*2] = -FS[win_step_num].d[2];
            par.ub[i*2] = FS[win_step_num].d[0];

            par.lb[i*2 + 1] = -FS[win_step_num].d[3];
            par.ub[i*2 + 1] = FS[win_step_num].d[1];


            unsigned int step_len_ms;
            if (T_ms[i] == 0)
            {
                if (sampling_period > step_time_left) 
                {
                    step_len_ms = step_time_left;
                }
                else
                {
                    step_len_ms = sampling_period;
                }
            }
            else
            {
                if (T_ms[i] > step_time_left) 
                {
                    retval = WMG_HALT;
                    break;
                }
                step_len_ms = T_ms[i];
            }
            par.T[i] = (double) step_len_ms / 1000;
            step_time_left -= step_len_ms;

            if (i == 0)
            {
                last_time_decrement = step_len_ms;
            }
            i++;
        }
        else
        {
            win_step_num++;
            if (win_step_num == FS.size())
            {
                retval = WMG_HALT;
                break;
            }
            step_time_left = FS[win_step_num].time_left;
        }
    }


    if (retval == WMG_OK)
    {
        while (FS[current_step_number].time_left == 0)
        {
            current_step_number++;
        }

        first_preview_step = current_step_number;
        FS[current_step_number].time_left -= last_time_decrement;
        if (FS[current_step_number].time_left == 0)
        {
            current_step_number++;
        }
    }

    return (retval);
}



/**
 * @brief Outputs the footsteps in FS to a file, that can be executed in
 * Matlab/Octave to get a figure of the steps.
 *
 * @param[in] filename output file name.
 * @param[in] plot_ds enable/disable plotting of double supports
 */
void WMG::FS2file(const std::string filename, const bool plot_ds)
{
    
    FILE *file_op = fopen(filename.c_str(), "w");
    
    if(!file_op)
    {
        fprintf(stderr, "Cannot open file (for writing)\n");
        return;
    }
    
    fprintf(file_op,"%%\n%% Footsteps generated using the c++ version of the WMG\n%%\n\n");
    fprintf(file_op,"cla;\n");
    fprintf(file_op,"clear FS;\n\n");
    
    int i;
    for (i=0; i< (int) FS.size(); i++ )
    {
        if ((plot_ds) || (FS[i].type != FS_TYPE_DS))
        {
            fprintf(file_op, "FS(%i).a = %f;\nFS(%i).p = [%f;%f];\nFS(%i).d = [%f;%f;%f;%f];\n", 
                    i+1, FS[i].angle, 
                    i+1, FS[i].x(), FS[i].y(), 
                    i+1, FS[i].d[0], FS[i].d[1], FS[i].d[2], FS[i].d[3]);

            fprintf(file_op, "FS(%i).D = [%f %f;%f %f;%f %f;%f %f];\n", 
                    i+1, FS[i].D[0], FS[i].D[4],
                         FS[i].D[1], FS[i].D[5],
                         FS[i].D[2], FS[i].D[6],
                         FS[i].D[3], FS[i].D[7]); 

            fprintf(file_op, "FS(%i).v = [%f %f; %f %f; %f %f; %f %f; %f %f];\n", 
                    i+1, FS[i].vert(0,0), FS[i].vert(0,1), 
                         FS[i].vert(1,0), FS[i].vert(1,1), 
                         FS[i].vert(2,0), FS[i].vert(2,1), 
                         FS[i].vert(3,0), FS[i].vert(3,1), 
                         FS[i].vert(0,0), FS[i].vert(0,1));

            if (FS[i].type == FS_TYPE_DS)
            {
                fprintf(file_op, "FS(%i).type = 1;\n\n", i+1);
            }
            if ((FS[i].type == FS_TYPE_SS_L) || (FS[i].type == FS_TYPE_SS_R))
            {
                fprintf(file_op, "FS(%i).type = 2;\n\n", i+1);
            }
        }
    }

    fprintf(file_op,"hold on\n");    
    fprintf(file_op,"for i=1:length(FS)\n");
    fprintf(file_op,"    if FS(i).type == 1;\n");
    fprintf(file_op,"        plot (FS(i).p(1),FS(i).p(2),'gs','MarkerFaceColor','r','MarkerSize',2)\n");
    fprintf(file_op,"        plot (FS(i).v(:,1), FS(i).v(:,2), 'c');\n");
    fprintf(file_op,"    end\n");
    fprintf(file_op,"    if FS(i).type == 2;\n");
    fprintf(file_op,"        plot (FS(i).p(1),FS(i).p(2),'gs','MarkerFaceColor','g','MarkerSize',4)\n");
    fprintf(file_op,"        plot (FS(i).v(:,1), FS(i).v(:,2), 'r');\n");
    fprintf(file_op,"    end\n");
    fprintf(file_op,"end\n");
    fprintf(file_op,"grid on; %%axis equal\n");
    fclose(file_op);  
}


/**
 * @brief Return coordinates of footstep reference points and rotation 
 * angles of footsteps (only for SS).
 *
 * @param[out] x_coord x coordinates
 * @param[out] y_coord y coordinates
 * @param[out] angle_rot angles
 */
void WMG::getFootsteps(
        std::vector<double> & x_coord,
        std::vector<double> & y_coord,
        std::vector<double> & angle_rot)
{
    for (unsigned int i = 0; i < FS.size(); i++)
    {
        if ((FS[i].type == FS_TYPE_SS_L) || (FS[i].type == FS_TYPE_SS_R))
        {
            x_coord.push_back(FS[i].x());
            y_coord.push_back(FS[i].y());
            angle_rot.push_back(FS[i].angle);
        }
    }
}
