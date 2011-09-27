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
    ZMP_ref = NULL;
    ind = NULL;

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
    if (ZMP_ref  != NULL)
    {
        delete [] ZMP_ref;
        ZMP_ref = NULL;
    }     

    if (ind  != NULL)
    {
        delete [] ind;
        ind = NULL;
    }

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
    /* ----------------------------- */
}



/** \brief Initializes a WMG object.

    \param[in] _N Number of sampling times in a preview window
    \param[in] _T Sampling time (for the moment it is assumed to be constant) [sec.]
    \param[in] _hCoM Height of the Center of Mass [meter]
 */
void WMG::init(int _N, double _T, double _hCoM)
{
    current_step_number = 0;
    counter = 0;
    halt = 0;
    gravity = 9.81; // hard-coded

    N = _N;
    hCoM = _hCoM;      



    ZMP_ref = new double[NUM_CONTROL_VAR*N];
    ind = new int[N];
    X = new double[NUM_VAR*N];


    for (int i=0; i<N; i++)
    {
        ind[i] = 0;
    }
    
    for(int i=0; i<NUM_STATE_VAR; i++)
    {
        X_tilde[i] = 0.0;
    }
    init_angle = 0;

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
}



/** \brief Adds a footstep to FS.

    \param[in] x_relative X position [meter] relative to the previous footstep (only for the first
    "virtual" step x_relative is with respect to the world frame.)

    \param[in] y_relative Y position [meter] relative to the previous footstep (only for the first
    "virtual" step y_relative is with respect to the world frame.)

    \param[in] angle_relative Angle [rad.] relative to the previous footstep (only for the first
    "virtual" step angle_relative is with respect to the world frame.)

    \param[in] _nSS Number of (preview window) iterations in Single Support.
    \param[in] _n Total number of (preview window) iterations, i.e., nSS + nDS.
    \param[in] _RL Right or left foot (1 for left, -1 for right).
    \param[in] _d Vector of the PoS constraints (assumed to be [4 x 1]).

    \note For now, it is assumed that there can not be a right support after right support
    (likewise for left support). In other words FS[i+1].RL = -FS[i].RL. Hence we can only set
    FootStep.RL for the first "virtual" step.
*/
void WMG::AddFootstep(double x_relative, double y_relative, double angle_relative, int _nSS, int _n, int _RL, double *_d)
{

    // if this is the first ("virtual") step. In this case only x_relative, y_relative, angle_relative are in the world frame
    if (FS.size() == 0)
    {
        FS.push_back(FootStep(angle_relative, Point2D(x_relative,y_relative), _nSS, _n, _RL, _d));
    }
    else
    {
        // Calculate (and set) the absolute position and orientation
        double a = FS.back().angle + angle_relative;

        Point2D p(FS.back().p.x, FS.back().p.y);        
        p.x += FS.back().ca*x_relative - FS.back().sa*y_relative;
        p.y += FS.back().sa*x_relative + FS.back().ca*y_relative;

        // consistency check 
        if (FS.back().RL == _RL)
        {
            printf("------------------------------------------------------------\n");
            printf(" WARNING: Two consecutive right leg (or left leg) supports. \n");
            printf("------------------------------------------------------------\n");
            return;
        }

        FS.push_back(FootStep(a, p, _nSS, _n, _RL, _d));
    }    
}



/** \brief Adds a footstep to FS.

    \param[in] x_relative X position [meter] relative to the previous footstep (only for the first
    "virtual" step x_relative is with respect to the world frame.)

    \param[in] y_relative Y position [meter] relative to the previous footstep (only for the first
    "virtual" step y_relative is with respect to the world frame.)

    \param[in] angle_relative Angle [rad.] relative to the previous footstep (only for the first
    "virtual" step angle_relative is with respect to the world frame.)

    \param[in] _nSS Number of (preview window) iterations in Single Support.
    \param[in] _n Total number of (preview window) iterations, i.e., nSS + nDS.

    \attention This function assumes that there is at least one footstep in FS.
*/
void WMG::AddFootstep(double x_relative, double y_relative, double angle_relative, int _nSS, int _n)
{
    // Calculate (and set) the absolute position and orientation
    double a = FS.back().angle + angle_relative;
    
    Point2D p(FS.back().p.x, FS.back().p.y);        
    p.x += FS.back().ca*x_relative - FS.back().sa*y_relative;
    p.y += FS.back().sa*x_relative + FS.back().ca*y_relative;
            
    FS.push_back(FootStep(a, p, _nSS, _n, -FS.back().RL, FS.back().ctr.d0));
}



/** \brief Adds a footstep to FS.

    \param[in] x_relative X position [meter] relative to the previous footstep (only for the first
    "virtual" step x_relative is with respect to the world frame.)

    \param[in] y_relative Y position [meter] relative to the previous footstep (only for the first
    "virtual" step y_relative is with respect to the world frame.)

    \param[in] angle_relative Angle [rad.] relative to the previous footstep (only for the first
    "virtual" step angle_relative is with respect to the world frame.)

    \attention This function assumes that there is at least one footstep in FS.
*/
void WMG::AddFootstep(double x_relative, double y_relative, double angle_relative)
{
    // Calculate (and set) the absolute position and orientation
    double a = FS.back().angle + angle_relative;
    
    Point2D p(FS.back().p.x, FS.back().p.y);        
    p.x += FS.back().ca*x_relative - FS.back().sa*y_relative;
    p.y += FS.back().sa*x_relative + FS.back().ca*y_relative;

    FS.push_back(FootStep(a, p, FS.back().nSS, FS.back().n, -FS.back().RL, FS.back().ctr.d0));    
}



/** \brief Forms the preview window. Using only single support!

    \note If the i-th step is old, then it will have FS[i].n = 0
 */
void WMG::FormPreviewWindow()
{
    int j, i=0, k = current_step_number;

    while (i < N && k < (int) FS.size()-1) 
    {        
        j = 1;
        while (j <= FS[k].n && i < N)
        {
            ind[i] = k;
            i++;
            j++;
        }  
        k++;   
    } 

    // if the reason for exiting the first while loop above is "k == FS.size()-1"
    if (i < N) // note that i is incremented once at the end
    {    
        halt = 1;
        printf(" \n\n====================================\n ");        
        printf(" WARNING: not enough steps in FS.");
        printf(" \n====================================\n\n ");        
    }


    for (i = 0; i < N; i++)
    {
        angle[i] = FS[ind[i]].angle;
        zref_x[i] = FS[ind[i]].p.x;
        zref_y[i] = FS[ind[i]].p.y;

        lb[i*2] = -FS[ind[i]].ctr.d[2];
        ub[i*2] = FS[ind[i]].ctr.d[0];

        lb[i*2 + 1] = -FS[ind[i]].ctr.d[3];
        ub[i*2 + 1] = FS[ind[i]].ctr.d[1];
    }

    this->slide();
}



/** \brief Increase #counter. Decrease #FS[#current_step_number].n or increase #current_step_number.  
        
    \note This function should be called each time before we update the preview window using
    #FormPreviewWindow (except for the first time).
*/
void WMG::slide()
{
    counter++;
    
    if (FS[current_step_number].n == 1)
    {
        FS[current_step_number].n--;
        current_step_number++;
    }
    else
    {
        FS[current_step_number].n--;
    }
}



/** \brief Print #ind in the terminal. */
void WMG::print_ind()
{
    printf("\n (counter = %i) ind = [",counter);
    for(int i=0; i<N; i++)
        printf(" %i ", ind[i]);
    printf("] (current_step_number = %i) (%i) \n", current_step_number, FS[current_step_number].n);    
}



/** \brief Outputs the footsteps in FS to a file (forms a Matlab structure FS). */
void WMG::FS2file()
{
    
    FILE *file_op = fopen("output/foot_steps_cpp.m", "w");
    
    if(!file_op)
    {
        std::cerr << "Cannot open file (for writing) " << std::endl;
        return;
    }
    
    fprintf(file_op,"%%\n%% Footsteps generated using the c++ version of the WMG\n%%\n\n");
    
    //fprintf(file_op,"clear;clc;\n\n");
    
    int i;
    for (i=0; i< (int) FS.size(); i++ )
    {
        fprintf(file_op, "FS(%i).a = %f;\nFS(%i).p = [%f;%f];\nFS(%i).Offset = [%f;%f];\nFS(%i).d = [%f;%f;%f;%f];\n", 
                i+1, FS[i].angle, 
                i+1, FS[i].p.x, FS[i].p.y, 
                i+1, FS[i].Offset.x, FS[i].Offset.y,
                i+1, FS[i].ctr.d[0], FS[i].ctr.d[1], FS[i].ctr.d[2], FS[i].ctr.d[3]);

        fprintf(file_op, "FS(%i).D = [%f %f;%f %f;%f %f;%f %f];\n", 
                i+1, FS[i].ctr.D[0], FS[i].ctr.D[4],
                     FS[i].ctr.D[1], FS[i].ctr.D[5],
                     FS[i].ctr.D[2], FS[i].ctr.D[6],
                     FS[i].ctr.D[3], FS[i].ctr.D[7]); 

        fprintf(file_op, "FS(%i).v = [%f %f;%f %f;%f %f;%f %f];\n\n", 
                i+1, FS[i].ctr.vert[0].x, FS[i].ctr.vert[0].y, 
                     FS[i].ctr.vert[1].x, FS[i].ctr.vert[1].y, 
                     FS[i].ctr.vert[2].x, FS[i].ctr.vert[2].y, 
                     FS[i].ctr.vert[3].x, FS[i].ctr.vert[3].y); 

    }

    fprintf(file_op,"hold on\n");    
    fprintf(file_op,"for i=1:length(FS)\n");
    fprintf(file_op,"    plot(FS(i).p(1),FS(i).p(2),'gs','MarkerFaceColor','g','MarkerSize',4)\n");
    fprintf(file_op,"%%    plot_ConvHull(FS(i).v, 'g')\n");
    fprintf(file_op,"%%    plot(FS(i).v(:,1),FS(i).v(:,2),'bs')\n"); 
    fprintf(file_op,"    plot_ConvHull(Constraints2Vert(-FS(i).D,FS(i).d), 'r--')\n");
    fprintf(file_op,"end\n");
    fprintf(file_op,"grid on; axis equal\n");

    fprintf(file_op,"\n%% Note: the -FS(i).D is because Constraints2Vert expects constraints of the form D*z+d>=0.\n");
    
    fclose(file_op);  
}
