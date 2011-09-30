/** 
 * @file
 * @author Alexander Sherikov
 * @date 27.09.2011 18:59:03 MSD
 */


#include <iostream>


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
    \param[in] _T Sampling time (for the moment it is assumed to be constant) [sec.]
    \param[in] _hCoM Height of the Center of Mass [meter]
 */
void WMG::init(const int _N, const double _T, const double _hCoM)
{
    current_step_number = 0;
    gravity = 9.81; // hard-coded

    N = _N;
    hCoM = _hCoM;      



    X = new double[NUM_VAR*N];

    for(int i=0; i<NUM_STATE_VAR; i++)
    {
        X_tilde[i] = 0.0;
    }

    T = new double[N];
    h = new double[N];
    for (int i = 0; i < N; i++)
    {
        T[i] = _T;
        h[i] = hCoM/gravity;
    }

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
 * @brief Adds a footstep to FS; sets the default constraints, the total number of 
 * iterations and the number of iterations in single support.
 *
 * @param[in] x_relative x_relative X position [meter] relative to the previous footstep.
 * @param[in] y_relative y_relative Y position [meter] relative to the previous footstep.
 * @param[in] angle_relative angle_relative Angle [rad.] relative to the previous footstep.
 * @param[in] nSS Number of (preview window) iterations in Single Support.
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
        const int nSS, 
        const int n, 
        const double *d, 
        const fs_type type)
{
    def_constraint[0] = d[0];
    def_constraint[1] = d[1];
    def_constraint[2] = d[2];
    def_constraint[3] = d[3];
    def_repeat_times = nSS;
    def_ds_num = n - nSS;
    AddFootstep(x_relative, y_relative, angle_relative, type);
}



/**
 * @brief Adds a footstep to FS; sets the default total number of iterations and the 
 * number of iterations in single support.
 *
 * @param[in] x_relative x_relative X position [meter] relative to the previous footstep.
 * @param[in] y_relative y_relative Y position [meter] relative to the previous footstep.
 * @param[in] angle_relative angle_relative Angle [rad.] relative to the previous footstep.
 * @param[in] nSS Number of (preview window) iterations in Single Support.
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
        const int nSS, 
        const int n, 
        const fs_type type)
{
    def_repeat_times = nSS;
    def_ds_num = n - nSS;
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
        const fs_type type)
{
    if (FS.size() == 0)
    {
        // this is the first ("virtual") step.
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
 * @brief Forms a preview window.
 *
 * @return WMG_OK or WMG_HALT (simulation must be stopped)
 */
WMGret WMG::FormPreviewWindow()
{
    WMGret retval = WMG_OK;
    int win_step_num = current_step_number;
    int step_repeat_times = FS[win_step_num].repeat_times;


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
            /// @todo "-1" is added to be compatible with the old version,
            /// this must be a bug.
            if (win_step_num == (int) FS.size() - 1)
            {
                retval = WMG_HALT;
                printf(" \n\n====================================\n ");        
                printf(" WARNING: not enough steps in FS.");
                printf(" \n====================================\n\n ");        
                break;
            }
            step_repeat_times = FS[win_step_num].repeat_times;
        }
    }

    FS[current_step_number].repeat_times--;
    if (FS[current_step_number].repeat_times == 0)
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
        std::cerr << "Cannot open file (for writing) " << std::endl;
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
        if (FS[i].type == FS_TYPE_SS)
        {
            fprintf(file_op, "FS(%i).type = 2;\n\n", i+1);
        }
    }

    fprintf(file_op,"hold on\n");    
    fprintf(file_op,"for i=1:length(FS)\n");
    fprintf(file_op,"    if FS(i).type == 1;");
    fprintf(file_op,"        plot (FS(i).p(1),FS(i).p(2),'gs','MarkerFaceColor','b','MarkerSize',2)\n");
    fprintf(file_op,"        plot (FS(i).v(:,1), FS(i).v(:,2), 'c');\n");
    fprintf(file_op,"    end\n");
    fprintf(file_op,"    if FS(i).type == 2;");
    fprintf(file_op,"        plot (FS(i).p(1),FS(i).p(2),'gs','MarkerFaceColor','g','MarkerSize',4)\n");
    fprintf(file_op,"        plot (FS(i).v(:,1), FS(i).v(:,2), 'r');\n");
    fprintf(file_op,"    end\n");
    fprintf(file_op,"end\n");
    fprintf(file_op,"grid on; axis equal\n");

    fprintf(file_op,"\n%% Note: the -FS(i).D is because Constraints2Vert expects constraints of the form D*z+d>=0.\n");
    
    fclose(file_op);  
}
