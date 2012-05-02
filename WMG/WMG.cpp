/** 
 * @file
 * @author Alexander Sherikov
 * @date 27.09.2011 18:59:03 MSD
 */


#include <stdio.h>
#include <math.h> // cos, sin


#include "WMG.h"
#include "footstep.h"

WMG::WMG (
        const unsigned int N_,
        const unsigned int T_, 
        const double step_height_,
        const double bezier_weight_1_,
        const double bezier_weight_2_,
        const double bezier_inclination_1_,
        const double bezier_inclination_2_,
        const bool use_fsr_constraints)
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

    if (use_fsr_constraints)
    {
        def_constraints.init_FSR();
    }
    else
    {
        def_constraints.init();
    }

    setFootstepParameters (4, 0, 0);

    for (int i = 0; i < 4; ++i)
    {
        // auto_ds constraints are safer
        user_constraints[i] = user_constraints_auto_ds[i] = def_constraints.auto_ds[i];
    }

    bezier_weight_1 = bezier_weight_1_;
    bezier_weight_2 = bezier_weight_2_;
    bezier_inclination_1 = bezier_inclination_1_;
    bezier_inclination_2 = bezier_inclination_2_;
}



WMG::~WMG()
{
    if (T_ms != NULL)
    {
        delete T_ms;
    }
}



void WMG::setFootstepParameters(
        const unsigned int def_periods, 
        const unsigned int ds_periods, 
        const unsigned int ds_number,
        bool use_user_constraints_)
{
    setFootstepParametersMS(
            def_periods * sampling_period,
            ds_periods * sampling_period,
            ds_number,
            use_user_constraints_);
}


void WMG::setFootstepParametersMS(
        const unsigned int def_time_ms_, 
        const unsigned int ds_time_ms_, 
        const unsigned int ds_number,
        bool use_user_constraints_)
{
    def_time_ms = def_time_ms_;
    ds_time_ms  = ds_time_ms_; 
    ds_num      = ds_number;
    use_user_constraints = use_user_constraints_;
}


void WMG::addFootstep(
        const double x_relative, 
        const double y_relative, 
        const double angle_relative, 
        fs_type type)
{
    const double *constraints;
    const double *constraints_auto_ds;

    // determine type of the step
    if (type == FS_TYPE_AUTO)
    {
        if (FS.size() == 0)
        {
            type = FS_TYPE_DS;
        }
        else
        {
            switch (FS[getPrevSS(FS.size())].type)
            {
                case FS_TYPE_SS_R:
                    type = FS_TYPE_SS_L;
                    break;
                case FS_TYPE_SS_L:
                case FS_TYPE_DS:
                default:
                    type = FS_TYPE_SS_R;
                    break;
            }
        }
    }

    if (use_user_constraints)
    {
        constraints = user_constraints;
        constraints_auto_ds = user_constraints_auto_ds;
    }
    else
    {
        constraints_auto_ds = def_constraints.auto_ds;
        switch (type)
        {
            case FS_TYPE_SS_R:
                constraints = def_constraints.ss_right;
                break;
            case FS_TYPE_SS_L:
                constraints = def_constraints.ss_left;
                break;
            case FS_TYPE_DS:
            default:
                constraints = def_constraints.ds;
                break;
        }
    }


    Transform<double, 3>    posture (Translation<double, 3>(x_relative, y_relative, 0.0));
    Vector3d    zref_offset ((constraints[0] - constraints[2])/2, 0.0, 0.0);


    if (FS.size() == 0)
    {
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
                    constraints));
    }
    else
    {
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
                        constraints_auto_ds));
        }


        // add the new step
        FS.push_back(
                footstep(
                    next_a, 
                    posture, 
                    next_zref,
                    def_time_ms, 
                    type,
                    constraints));
    }    
}



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



void WMG::changeNextSSPosition (const double* posture, const bool zero_z_coordinate)
{
    FS[getNextSS(first_preview_step)].changePosture(posture, zero_z_coordinate);
}


void WMG::repositionFootsteps (const double diff_x, const double diff_y)
{
    int ind = 0;
    fs_type fixed_fs_type;
    

    ind = first_preview_step;
    if (FS[first_preview_step].type == FS_TYPE_DS)
    {
        ind = getNextSS(ind);
    }
    fixed_fs_type = FS[ind].type;

    for (; (ind < (int) FS.size()) && (FS[ind].type != fixed_fs_type); ++ind);

    Translation<double, 3>  diff(diff_x, diff_y, 0.0);
   
    for (; ind < (int) FS.size(); ++ind)
    {
        *FS[ind].posture = diff * (*FS[ind].posture);
        FS[ind].rotate_translate(FS[ind].ca, FS[ind].sa, FS[ind].x(), FS[ind].y());
    }
}


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
