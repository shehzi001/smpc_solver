/** 
 * @file
 * @author Alexander Sherikov
 * @date 27.09.2011 18:59:03 MSD
 */


#include <stdio.h>


#include "WMG.h"
#include "point2d.h"
#include "footstep.h"


/** \brief Default constructor. */
WMG::WMG()
{
    X = NULL;

    T = NULL;
    h = NULL;

    angle = NULL;
    zref_x = NULL;
    zref_y = NULL;
    lb = NULL;
    ub = NULL;

    current_reference_foot = FS_TYPE_AUTO;
}


/** \brief Default destructor. */
WMG::~WMG()
{
    if (X != NULL)
    {
        delete [] X;
        X = NULL;
    }

    if (T != NULL)
    {
        delete T;
    }
    if (h != NULL)
    {
        delete h;
    }

    if (angle != NULL)
    {
        delete angle;
    }
    if (zref_x != NULL)
    {
        delete zref_x;
    }
    if (zref_y != NULL)
    {
        delete zref_y;
    }
    if (lb != NULL)
    {
        delete lb;
    }
    if (ub != NULL)
    {
        delete ub;
    }
}



/** \brief Initializes a WMG object.

    \param[in] _N Number of sampling times in a preview window
 */
void WMG::init(const int _N)
{
    current_step_number = 0;
    gravity = 9.81; // hard-coded

    N = _N;


    X = new double[SMPC_NUM_VAR*N];

    for(int i=0; i<SMPC_NUM_STATE_VAR; i++)
    {
        X_tilde[i] = 0.0;
    }

    T = new double[N];
    h = new double[N];

    angle = new double[N];
    zref_x = new double[N];
    zref_y = new double[N];
    lb = new double[2*N];
    ub = new double[2*N];



    /// NAO constraint with safety margin.
    def_constraint[0] = 0.09;
    def_constraint[1] = 0.025;
    def_constraint[2] = 0.03;
    def_constraint[3] = 0.025;
    def_repeat_times = 4;

    def_ds_constraint[0] = 0.07;
    def_ds_constraint[1] = 0.025;
    def_ds_constraint[2] = 0.025;
    def_ds_constraint[3] = 0.025;
    def_ds_num = 0;
}


/** 
 * @brief Initializes a WMG object.
 *
 * @param[in] T_ Sampling time (for the moment it is assumed to be constant) [sec.]
 * @param[in] hCoM_ Height of the Center of Mass [meter]
 * @param[in] step_height_ step hight (for interpolation of feet movements)
 */
void WMG::init_param(
        const double T_, 
        const double hCoM_,
        const double step_height_)
{
    hCoM = hCoM_;      

    for (int i = 0; i < N; i++)
    {
        T[i] = T_;
        h[i] = hCoM/gravity;
    }

    step_height = step_height_;
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
    def_constraint[0] = d[0];
    def_constraint[1] = d[1];
    def_constraint[2] = d[2];
    def_constraint[3] = d[3];
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
    if (FS.size() == 0)
    {
        // this is the first ("virtual") step.
        if (type == FS_TYPE_AUTO)
        {
            type = FS_TYPE_DS;
        }

        FS.push_back(
                FootStep(
                    angle_relative, 
                    x_relative,
                    y_relative, 
                    def_repeat_times, 
                    type,
                    def_constraint));
    }
    else
    {
        // Angle and position of the previous step
        double prev_a = FS.back().angle;
        Point2D prev_p(FS.back());        

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
        Point2D next_p(prev_p);        
        next_p.x += FS.back().ca*x_relative - FS.back().sa*y_relative;
        next_p.y += FS.back().sa*x_relative + FS.back().ca*y_relative;


        // Add double support constraints that lie between the
        // newly added step and the previous step
        for (int i = 0; i < def_ds_num; i++)
        {
            double theta = (double) (i+1)/(def_ds_num + 1);
            double ds_a = prev_a + angle_relative * theta;

            FS.push_back(
                    FootStep(
                        ds_a, 
                        (1-theta)*prev_p.x + theta*next_p.x,
                        (1-theta)*prev_p.y + theta*next_p.y,
                        1, 
                        FS_TYPE_DS,
                        def_ds_constraint));
        }


        // add the new step
        FS.push_back(
                FootStep(
                    prev_a + angle_relative, 
                    next_p.x, 
                    next_p.y, 
                    def_repeat_times, 
                    type,
                    def_constraint));
    }    
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
    int index = start_ind;
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
    int index = start_ind;
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
 * @brief Returns current feet configuration.
 *
 * @param[in,out] prev_swing_pos the previous position of the swing foot.
 * @param[in,out] next_swing_pos the next position of the swing foot.
 * @param[out] repeat_times number of SS iterations in the current support.
 * @param[out] repeated_times the number of iterations, which were already
 *      spent in the current support.
 *
 * @return type of the current support feet: left or right.
 *
 * @attention This function requires the walking pattern to be started and finished
 * by single support.
 *
 * @attention Cannot be called on the first or last SS  =>  must be called after 
 * FormPreviewWindow().
 *
 * @note If we are in the double support then prev_swing_pos = next_swing_pos, and
 * repeat_times = repeated_times = 0.
 */
fs_type WMG::getSwingFootNextPrevPos (
        double *prev_swing_pos,
        double *next_swing_pos,
        int *repeat_times,
        int *repeated_times)
{
    fs_type cur_fs_type = FS[current_step_number].type;
    int next_swing_ind, prev_swing_ind;
    int next_ss_ind, prev_ss_ind;

    switch (cur_fs_type)
    {
        case FS_TYPE_SS_L:
        case FS_TYPE_SS_R:
            prev_swing_ind = getPrevSS (
                    current_step_number, 
                    (current_reference_foot == FS_TYPE_SS_R) ? FS_TYPE_SS_L : FS_TYPE_SS_R);
            next_swing_ind = getNextSS (
                    current_step_number, 
                    (current_reference_foot == FS_TYPE_SS_R) ? FS_TYPE_SS_L : FS_TYPE_SS_R);

            *repeat_times = FS[current_step_number].repeat_times;
            *repeated_times = *repeat_times - FS[current_step_number].repeat_counter;
            break;

        case FS_TYPE_DS:
            next_ss_ind = getNextSS (current_step_number);
            prev_ss_ind = getPrevSS (current_step_number);
            if (FS[next_ss_ind].type == current_reference_foot)
            {
                next_swing_ind = prev_swing_ind = prev_ss_ind;
            }
            else
            {
                next_swing_ind = prev_swing_ind = next_ss_ind;
            }

            *repeat_times = 0;
            *repeated_times = 0;
            break;

        case FS_TYPE_AUTO:
        default:
            // should never happen, since type is checked when steps are added.
            return (cur_fs_type);
    }


    prev_swing_pos[0] = FS[prev_swing_ind].x;
    prev_swing_pos[1] = FS[prev_swing_ind].y;
    prev_swing_pos[2] = FS[prev_swing_ind].angle;

    next_swing_pos[0] = FS[next_swing_ind].x;
    next_swing_pos[1] = FS[next_swing_ind].y;
    next_swing_pos[2] = FS[next_swing_ind].angle;

    return (cur_fs_type);
}



/**
 * @brief Determine position and orientation of the swing foot
 *
 * @param[in] swing_type which method to use.
 * @param[in] loops_per_preview_iter number of control loops per preview iteration
 * @param[out] swing_foot_position 3x1 vector of coordinates [x y z]
 * @param[out] angle orientation of the foot in x axis
 *
 * @attention This function requires the walking pattern to be started and finished
 * by single support.
 *
 * @attention Cannot be called on the first or last SS  =>  must be called after 
 * FormPreviewWindow().
 */
void WMG::getSwingFootPosition (
        const swing_foot_pos_type swing_type,
        const int loops_per_preview_iter,
        const int loops_in_current_preview,
        double *swing_foot_pos,
        double *angle_x)
{
    if (swing_type == WMG_SWING_PARABOLA)
    {
        int num_iter_in_ss;
        int num_iter_in_ss_passed;
        double prev_swing_pos[3];
        double next_swing_pos[3];

        if (getSwingFootNextPrevPos(
                prev_swing_pos,
                next_swing_pos,
                &num_iter_in_ss,
                &num_iter_in_ss_passed) != FS_TYPE_DS)
        {
            double theta =
                (double) (loops_per_preview_iter * num_iter_in_ss_passed + loops_in_current_preview) /
                (loops_per_preview_iter * num_iter_in_ss);

            double x[3] = {
                prev_swing_pos[0],
                (prev_swing_pos[0] + next_swing_pos[0])/2,
                next_swing_pos[0]};

            double b_coef = - (x[2]*x[2] - x[0]*x[0])/(x[2] - x[0]);
            double a = step_height / (x[1]*x[1] - x[0]*x[0] + b_coef*(x[1] - x[0]));
            double b = a * b_coef;
            double c = - a*x[0]*x[0] - b*x[0];

            swing_foot_pos[0] = (1-theta)*x[2] + theta*x[0]; // linear equation
            swing_foot_pos[1] = next_swing_pos[1];
            swing_foot_pos[2] = a * swing_foot_pos[0] * swing_foot_pos[0] + b * swing_foot_pos[0] + c;
        }
        else
        {
            swing_foot_pos[0] = next_swing_pos[0];
            swing_foot_pos[1] = next_swing_pos[1];
            swing_foot_pos[2] = 0.0;
        }

        *angle_x = next_swing_pos[2];
    }
}




/**
 * @brief Forms a preview window.
 *
 * @return WMG_OK or WMG_HALT (simulation must be stopped)
 */
WMGret WMG::FormPreviewWindow()
{
    WMGret retval = WMG_OK;
    int win_step_num = current_step_number;
    int step_repeat_times = FS[win_step_num].repeat_counter;


    // If reference foot is unknown, this function is executed
    // for the first time: we need to use the next SS as the
    // reference foot, since the first SS must be fake.
    if (current_reference_foot == FS_TYPE_AUTO)
    {
        current_reference_foot = FS[getNextSS (win_step_num)].type;
    }

    // Indicate switch of support foot
    if (    // if we are not in the initial support
            (win_step_num != 0) && 
            // we are in DS
            (FS[win_step_num].type == FS_TYPE_DS) &&
            // this is the first iteration in DS
            (FS[win_step_num].repeat_counter == FS[win_step_num].repeat_times))
    {
        retval = WMG_SWITCH_REFERENCE_FOOT;
        current_reference_foot = FS[getNextSS (win_step_num)].type;
    }


    for (int i = 0; i < N;)
    {
        if (step_repeat_times > 0)
        {
            angle[i] = FS[win_step_num].angle;
            zref_x[i] = FS[win_step_num].x;
            zref_y[i] = FS[win_step_num].y;

            lb[i*2] = -FS[win_step_num].d[2];
            ub[i*2] = FS[win_step_num].d[0];

            lb[i*2 + 1] = -FS[win_step_num].d[3];
            ub[i*2 + 1] = FS[win_step_num].d[1];

            step_repeat_times--;
            i++;
        }
        else
        {
            win_step_num++;
            if (win_step_num == (int) FS.size())
            {
                retval = WMG_HALT;
                break;
            }
            step_repeat_times = FS[win_step_num].repeat_counter;
        }
    }

    FS[current_step_number].repeat_counter--;
    if (FS[current_step_number].repeat_counter <= 0)
    {
        current_step_number++;
    }

    return (retval);
}



/**
 * @brief Outputs the footsteps in FS to a file, that can be executed in
 * Matlab/Octave to get a figure of the steps.
 *
 * @param[in] filename output file name.
 */
void WMG::FS2file(const std::string filename)
{
    
    FILE *file_op = fopen(filename.c_str(), "w");
    
    if(!file_op)
    {
        fprintf(stderr, "Cannot open file (for writing)\n");
        return;
    }
    
    fprintf(file_op,"%%\n%% Footsteps generated using the c++ version of the WMG\n%%\n\n");
    
    int i;
    for (i=0; i< (int) FS.size(); i++ )
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
    fprintf(file_op,"grid on; axis equal\n");
    fclose(file_op);  
}
