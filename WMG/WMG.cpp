/** 
 * @file
 * @author Alexander Sherikov
 * @date 27.09.2011 18:59:03 MSD
 */


#include <stdio.h>
#include <math.h> // cos, sin


#include "WMG.h"
#include "point2d.h"
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
        T_ms[i] = T_;
    }
    

    current_step_number = 0;
    last_time_decrement = 0;
    first_preview_step = current_step_number;


    /// NAO constraint with safety margin.
    def_ss_constraint[0] = 0.09;
    def_ss_constraint[1] = 0.025;
    def_ss_constraint[2] = 0.03;
    def_ss_constraint[3] = 0.025;

    def_ds_constraint[0] = 0.07;
    def_ds_constraint[1] = 0.025;
    def_ds_constraint[2] = 0.025;
    def_ds_constraint[3] = 0.025;

    addstep_constraint[0] = def_ss_constraint[0];
    addstep_constraint[1] = def_ss_constraint[1];
    addstep_constraint[2] = def_ss_constraint[2];
    addstep_constraint[3] = def_ss_constraint[3];

    def_repeat_times = 4;
    def_ds_num = 0;
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
 * @brief Adds a footstep to FS; sets the default constraints, the total number of 
 * iterations and the number of iterations in single support.
 *
 * @param[in] x_relative x_relative X position [meter] relative to the previous footstep.
 * @param[in] y_relative y_relative Y position [meter] relative to the previous footstep.
 * @param[in] angle_relative angle_relative Angle [rad.] relative to the previous footstep.
 * @param[in] n_this Number of (preview window) iterations in the added step.
 * @param[in] n Total number of (preview window) iterations, i.e., nSS + nDS.
 * @param[in] d Vector of the PoS constraints (assumed to be [4 x 1]).
 * @param[in] type (optional) type of the footstep.
 *
 * @note Coordinates and angle are treated as absolute for the first step in the preview window.
 */
void WMG::AddFootstep(
        const double x_relative, 
        const double y_relative, 
        const double angle_relative, 
        const int n_this, 
        const int n, 
        const double *d, 
        const fs_type type)
{
    addstep_constraint[0] = d[0];
    addstep_constraint[1] = d[1];
    addstep_constraint[2] = d[2];
    addstep_constraint[3] = d[3];
    def_repeat_times = n_this;
    def_ds_num = n - n_this;
    AddFootstep(x_relative, y_relative, angle_relative, type);
}



/**
 * @brief Adds a footstep to FS; sets the default total number of iterations and the 
 * number of iterations in single support.
 *
 * @param[in] x_relative x_relative X position [meter] relative to the previous footstep.
 * @param[in] y_relative y_relative Y position [meter] relative to the previous footstep.
 * @param[in] angle_relative angle_relative Angle [rad.] relative to the previous footstep.
 * @param[in] n_this Number of (preview window) iterations in the added step.
 * @param[in] n Total number of (preview window) iterations, i.e., nSS + nDS.
 * @param[in] type (optional) type of the footstep.
 *
 * @note Coordinates and angle are treated as absolute for the first step in the preview window.
 * @note Default vector of the PoS constraints is used.
 */
void WMG::AddFootstep(
        const double x_relative, 
        const double y_relative, 
        const double angle_relative, 
        const int n_this, 
        const int n, 
        const fs_type type)
{
    def_repeat_times = n_this;
    def_ds_num = n - n_this;
    AddFootstep(x_relative, y_relative, angle_relative, type);
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
 * @note Default vector of the PoS constraintsi, number of iterations in single support and
 * total number of iterations are used.
 */
void WMG::AddFootstep(
        const double x_relative, 
        const double y_relative, 
        const double angle_relative, 
        fs_type type)
{
    Point2D zref_offset ((addstep_constraint[0] - addstep_constraint[2])/2, 0);
    Point2D offset (x_relative, y_relative);

    if (FS.size() == 0)
    {
        // this is the first ("virtual") step.
        if (type == FS_TYPE_AUTO)
        {
            type = FS_TYPE_DS;
        }

        // offset is the absolute position here.
        Point2D zref_abs (offset, zref_offset, sin(angle_relative), cos(angle_relative));

        FS.push_back(
                FootStep(
                    angle_relative, 
                    offset,
                    zref_abs,
                    def_repeat_times * sampling_period, 
                    type,
                    addstep_constraint));
    }
    else
    {
        // determine type of the step
        if (type == FS_TYPE_AUTO)
        {
            switch (FS.back().type)
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
        Point2D prev_p (FS.back());
        Point2D next_p (prev_p, offset, FS.back().sa, FS.back().ca);

        double prev_a = FS.back().angle;
        double next_a = prev_a + angle_relative;
        Point2D next_zref (next_p, zref_offset, sin(next_a), cos(next_a));

        // Add double support constraints that lie between the
        // newly added step and the previous step
        for (unsigned int i = 0; i < def_ds_num; i++)
        {
            double theta = (double) (i+1)/(def_ds_num + 1);
            double ds_a = prev_a + angle_relative * theta;

            Point2D position ((1-theta)*prev_p.x + theta*next_p.x, (1-theta)*prev_p.y + theta*next_p.y);
            Point2D *zmp_ref;
            if (i < def_ds_num/2)
            {
                zmp_ref = &FS.back().ZMPref;
            }
            else
            {
                zmp_ref = &next_zref;
            }
            FS.push_back(
                    FootStep(
                        ds_a, 
                        position,
                        *zmp_ref,
                        sampling_period, 
                        FS_TYPE_DS,
                        def_ds_constraint));
        }


        // add the new step
        FS.push_back(
                FootStep(
                    next_a, 
                    next_p, 
                    next_zref,
                    def_repeat_times * sampling_period, 
                    type,
                    addstep_constraint));
    }    
}



/**
 * @brief Determine position and orientation of feet
 *
 * @param[in] shift_from_current a positive shift from the current support (allows to get
 *  positions for the future supports)
 * @param[in] loops_per_preview_iter number of control loops per preview iteration
 * @param[in] loops_in_current_preview number of control loops passed in the current 
 *                                     preview iteration
 * @param[out] left_foot_pos 3x1 vector of coordinates [x y z] + angle (orientation in x,y plane)
 * @param[out] right_foot_pos 3x1 vector of coordinates [x y z] + angle (orientation in x,y plane)
 *
 * @attention This function requires the walking pattern to be started and finished
 * by single support.
 *
 * @attention Cannot be called on the first or last SS  =>  must be called after 
 * FormPreviewWindow().
 *
 * @note If loops_per_preview_iter is set to 1, then the function returns a position 
 * at the end of preview window with number loops_in_current_preview.
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

        getSSFeetPositions (
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
 * @brief Corrects position of the next SS.
 *
 * @param[in] pos_error 3x1 vector of error between expected and real position.
 *
 * @todo Following DS and the angle must be adjusted as well.
 */
void WMG::correctNextSSPosition (const double *pos_error)
{
    FS[getNextSS(first_preview_step)].correct(pos_error[0], pos_error[1]);
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

            par.fp_x[i] = FS[win_step_num].x;
            par.fp_y[i] = FS[win_step_num].y;


            // ZMP reference coordinates
            par.zref_x[i] = FS[win_step_num].ZMPref.x;
            par.zref_y[i] = FS[win_step_num].ZMPref.y;


            par.lb[i*2] = -FS[win_step_num].d[2];
            par.ub[i*2] = FS[win_step_num].d[0];

            par.lb[i*2 + 1] = -FS[win_step_num].d[3];
            par.ub[i*2 + 1] = FS[win_step_num].d[1];

            if (T_ms[i] > step_time_left) 
            {
                retval = WMG_HALT;
                break;
            }
            step_time_left -= T_ms[i];
            par.T[i] = (double) T_ms[i] / 1000;
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
        last_time_decrement = T_ms[0];
        FS[current_step_number].time_left -= T_ms[0];
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
                    i+1, FS[i].x, FS[i].y, 
                    i+1, FS[i].d[0], FS[i].d[1], FS[i].d[2], FS[i].d[3]);

            fprintf(file_op, "FS(%i).D = [%f %f;%f %f;%f %f;%f %f];\n", 
                    i+1, FS[i].D[0], FS[i].D[4],
                         FS[i].D[1], FS[i].D[5],
                         FS[i].D[2], FS[i].D[6],
                         FS[i].D[3], FS[i].D[7]); 

            fprintf(file_op, "FS(%i).v = [%f %f; %f %f; %f %f; %f %f; %f %f];\n", 
                    i+1, FS[i].vert[0].x, FS[i].vert[0].y, 
                         FS[i].vert[1].x, FS[i].vert[1].y, 
                         FS[i].vert[2].x, FS[i].vert[2].y, 
                         FS[i].vert[3].x, FS[i].vert[3].y, 
                         FS[i].vert[0].x, FS[i].vert[0].y);

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
            x_coord.push_back(FS[i].x);
            y_coord.push_back(FS[i].y);
            angle_rot.push_back(FS[i].angle);
        }
    }
}
