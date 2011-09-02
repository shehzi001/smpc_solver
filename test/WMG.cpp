// --------------------
// Definition of: 
// --------------------
// class Point2D
// class RectangularConstraint_ZMP
// class FootStep 
// class WMG
// --------------------
//
// --------------------
// Notes: 
// --------------------
// 1. It is assumed that the arrays are stored column-wise (Fortran format)
//
// 2. The ZMP constraints here are D*z <= d, while in the Matlab code they are D*z + d >= 0
//
// 3. Currently no double support is considered (when forming the preview window).
//
// [2011/07/02] DD
// This file was originally written by Dimitar Dimitrov.

// --------------------
// TODO:
// --------------------
// 1. Fill in ZMP_ref (and use Offset)
//


#include "WMG.h" // for cos and sin


// ============================================================================================
// Utility functions 
// ============================================================================================

/** \brief Print a [m by n] matrix in the terminal

    \param[in] m Number of rows
    \param[in] n Number of columns
    \param[in] A Matrix of dimension [m by n]
    \param[in] description A short description of the array (it is printed above the array)
    
    \attention It is assumed that A is stored column-wise (Fortran format)
 */
void Matrix_print(int m, int n, double * A, const char * description)
{
 
    int  i, j;
    
    printf(" %s", description);
    for (i=0; i<m; i++)
    {
        printf("\n");
        for (j=0; j<n; j++)
        {
            printf("% f ", A[ j*m + i ]);
        }
    } 
    printf("\n");
    
}

/** \brief Outputs a double array in ASCII file

    \param[in] A A double array
    \param[in] row Number of rows (the array has)
    \param[in] col Number of columns (the array has)
    \param[in] output_file File name (where to output the array)
    \param[in] mode If mode = "w" - write. If mode = "a" - append 
    
    \attention It is assumed that the array is stored column-wise (Fortran format)
 */
void write_file(double *A, int row, int col, const char *output_file, const char * mode) 
{

    FILE *file_op = fopen(output_file, mode);

    if(!file_op)
    {
        std::cerr << "Cannot open file (for writing) " << output_file << std::endl;
        return;
    }
    
    int i, j;
    for (i=0 ; i<row ; i++ )
    {
        for ( j=0 ; j<col ; j++ )
        {
            fprintf(file_op, "%4.15f ", A[ j*row + i ]);
        }
            
        fprintf(file_op, "\n");
    } 

    fclose(file_op);
}

// ============================================================================================
// End utility functions 
// ============================================================================================

Point2D::Point2D(){}

Point2D::Point2D(double _x, double _y)  
{
    set(_x, _y);
}

void Point2D::print()
{
    printf(" p = [% f, % f] \n\n", x, y);
}

void Point2D::set(double _x, double _y)
{
    x = _x;
    y = _y;
}

// ============================================================================================


/** \brief Default constructor.
    \brief NAO constraint with safety margin.
 */
RectangularConstraint_ZMP::RectangularConstraint_ZMP()
{
    D[0] =  1.0; D[4] =  0.0; d0[0] = 0.09;  d[0] = 0.09; 
    D[1] =  0.0; D[5] =  1.0; d0[1] = 0.025; d[1] = 0.025;
    D[2] = -1.0; D[6] =  0.0; d0[2] = 0.03;  d[2] = 0.03;
    D[3] =  0.0; D[7] = -1.0; d0[3] = 0.025; d[3] = 0.025;

    Constraints2Vert();
}



/** \brief Set the foot size.
    \param[in] _d Vector of the PoS constraints (assumed to be [4 x 1]).
 */
void RectangularConstraint_ZMP::set_size(double * _d)
{
    d0[0] = _d[0]; 
    d0[1] = _d[1]; 
    d0[2] = _d[2]; 
    d0[3] = _d[3];
    
    Constraints2Vert();
}



/** \brief Prints the constraints. */
void RectangularConstraint_ZMP::print()
{
    printf("           D                 d \n");
    printf(" |% f % f|  % f\n", D[0], D[4], d[0]);
    printf(" |% f % f|  % f\n", D[1], D[5], d[1]);
    printf(" |% f % f|  % f\n", D[2], D[6], d[2]);
    printf(" |% f % f|  % f\n\n", D[3], D[7], d[3]);
    printf("          vert \n");
    for (int i=0; i<4; i++)
    {
        printf(" % f % f \n", vert[i].x, vert[i].y);
    }
}



/** \brief translates from [0;0] to "p" and rotates from 0 to "angle". The 0 initial angle implies
    that D = [eye(2); -eye(2)].

    \param[in] ca cos(angle)
    \param[in] sa sin(angle)
    \param[in] p 2D point 

    \note This is used when the constraints are initialized (only then the orientation can be set).
 */
void RectangularConstraint_ZMP::rotate_translate(double ca, double sa, Point2D p)
{
    // D = D*R'
    D[0] =  ca; D[4] =  sa;
    D[1] = -sa; D[5] =  ca;
    D[2] = -ca; D[6] = -sa;
    D[3] =  sa; D[7] = -ca;
    
    // d = d0 - D*p
    d[0] = d0[0] + D[0]*p.x + D[4]*p.y;
    d[1] = d0[1] + D[1]*p.x + D[5]*p.y;
    d[2] = d0[2] + D[2]*p.x + D[6]*p.y;
    d[3] = d0[3] + D[3]*p.x + D[7]*p.y;

    Constraints2Vert();
}



/** \brief Computes the vertices of a polygon from the constraints (D*x <= d)
        
    \note It is assumed that D is a [4 x 2] matrix stored column-wise.
    
    The numbering of the vertices is (see the definition of #d0)
    \verbatim
    3---------2
    |         |
    |         |
    4---------1
    \endverbatim
    
    \note I am too lazy to check whether det = 0 :).
*/
void RectangularConstraint_ZMP::Constraints2Vert()
{
    double det;
    vert.clear();
    
    // Indexes of ctr.D
    // 0 4
    // 1 5 
    // 2 6
    // 3 7
    
    // |0 4|   0     | 7 -4|
    // |3 7|   3     |-3  0|
    det = D[0]*D[7] - D[3]*D[4];
    vert.push_back(Point2D( D[7]/det*d[0] - D[4]/det*d[3],-D[3]/det*d[0] + D[0]/det*d[3])); 
    
    // |0 4|   0     | 5 -4|
    // |1 5|   1     |-1  0|
    det = D[0]*D[5] - D[4]*D[1]; 
    vert.push_back(Point2D( D[5]/det*d[0] - D[4]/det*d[1],-D[1]/det*d[0] + D[0]/det*d[1])); 
    
    // |1 5|   1     | 6 -5|
    // |2 6|   2     |-2  1|
    det = D[1]*D[6] - D[5]*D[2]; 
    vert.push_back(Point2D( D[6]/det*d[1] - D[5]/det*d[2],-D[2]/det*d[1] + D[1]/det*d[2])); 
    
    // |2 6|   2     | 7 -6|
    // |3 7|   3     |-3  2|
    det = D[2]*D[7] - D[3]*D[6]; 
    vert.push_back(Point2D( D[7]/det*d[2] - D[6]/det*d[3],-D[3]/det*d[2] + D[2]/det*d[3])); 
}

// ============================================================================================


/** \brief Default constructor. */
FootStep::FootStep()
{
    p.set(0.0, 0.0);
    angle = 0.0; ca = 0.0; sa = 0.0;
    Offset.set(0.0, 0.0);
    nSS = 0;
    n = 0;
    RL = 0;
}



/** \brief Defines a footstep at a given position with a given orientation. 
        
    \param[in] _p 2D position (of the "point of interest")
    \param[in] _angle orientation
    \param[in] _nSS Number of (preview window) iterations in Single Support.
    \param[in] _n Total number of (preview window) iterations, i.e., nSS + nDS.
    \param[in] _RL Right or left foot (1 for left, -1 for right).
 */
FootStep::FootStep(double _angle, Point2D _p, int _nSS, int _n, int _RL)
{
    set(_angle, _p, _nSS, _n, _RL);
    Offset.set(0.0, 0.0);
}



/** \brief Defines a footstep at a given position with a given orientation. 
        
    \param[in] _p 2D position (of the "point of interest")
    \param[in] _angle orientation
    \param[in] _nSS Number of (preview window) iterations in Single Support.
    \param[in] _n Total number of (preview window) iterations, i.e., nSS + nDS.
    \param[in] _RL Right or left foot (1 for left, -1 for right).
    \param[in] _d Vector of the PoS constraints (assumed to be [4 x 1])
 */
FootStep::FootStep(double _angle, Point2D _p, int _nSS, int _n, int _RL, double *_d)
{
    ctr.set_size(_d);
    set(_angle, _p, _nSS, _n, _RL);
    Offset.set(0.0, 0.0);
}



/** \brief Sets the position and orientation of a footstep.
        
    \param[in] _p 2D Position (of the point of interest)
    \param[in] _angle Orientation
    \param[in] _nSS Number of (preview window) iterations in Single Support.
    \param[in] _n Total number of (preview window) iterations, i.e., nSS + nDS.
    \param[in] _RL Right or left foot (1 for left, -1 for right).
 */
void FootStep::set(double _angle, Point2D _p, int _nSS, int _n, int _RL)  
{
    p.set(_p.x, _p.y);
    angle = _angle; ca = cos(angle); sa = sin(angle);
    ctr.rotate_translate(ca, sa, p);
    nSS = _nSS;
    n = _n;
    RL = _RL;
}



/** \brief Prints information about a footstep. */
void FootStep::print()
{
    printf("----------------------------------\n");
    printf(" nSS = %d, n = %d, RL = %d  \n", nSS, n, RL);
    printf(" angle = % f, \n", angle);
    printf(" Offset = [% f, %f] \n", Offset.x, Offset.y);
    p.print();
    ctr.print();
    printf("----------------------------------\n");
}

// ============================================================================================


/** \brief Default constructor (set some pointers to NULL). */
WMG::WMG()
{
    ZMP_ref = NULL;
    ind = NULL;
    FP_init = NULL;
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

    if (FP_init != NULL)
    {
        delete [] FP_init;
        FP_init = NULL;
    }

    if (T != NULL)
    {
        delete T;
    }
    if (h != NULL)
    {
        delete h;
    }
    if (dh != NULL)
    {
        delete dh;
    }
    /* ----------------------------- */
}



/** \brief Initializes a WMG object.

    \param[in] _N Number of sampling times in a preview window
    \param[in] _T Sampling time (for the moment it is assumed to be constant) [sec.]
    \param[in] _hCoM Height of the Center of Mass [meter]
 */
WMG::WMG(int _N, double _T, double _hCoM)
{
    current_step_number = 0;
    counter = 0;
    halt = 0;
    gravity = 9.81; // hard-coded

    N = _N;
    hCoM = _hCoM;      



    ZMP_ref = new double[NUM_CONTROL_VAR*N];
    ind = new int[N];
    FP_init = new double[NUM_VAR*N];


    for (int i=0; i<N; i++)
    {
        ind[i] = 0;
    }
    
    for(int i=0; i<NUM_STATE_VAR; i++)
    {
        X[i] = 0.0;
    }

    T = new double[N];
    h = new double[N];
    dh = new double[N-1]();
    for (int i = 0; i < N; i++)
    {
        T[i] = _T;
        h[i] = hCoM/gravity;
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
    
    if (counter == 0)
        first_angle_old = FS[current_step_number].angle;

    double ca = cos(first_angle_old);
    double sa = sin(first_angle_old);
 
    // update the tilde state (used when generating an initial feasible point in WMG.form_FP_init,
    // and when updating the variables wmg->CoM and wmg->ZMP in form_problem_eigen.solve_QP)
    X_tilde[0] = ca*X[0] - sa*X[3];
    X_tilde[1] = X[1];
    X_tilde[2] = X[2];
    X_tilde[3] = sa*X[0] + ca*X[3];
    X_tilde[4] = X[4];
    X_tilde[5] = X[5];
    
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
}



/** \brief Increase #counter. Decrease #FS[#current_step_number].n or increase #current_step_number.  
        
    \note This function should be called each time before we update the preview window using
    #FormPreviewWindow (except for the first time).
*/
void WMG::slide()
{
    counter++;
    
    first_angle_old = FS[current_step_number].angle;

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



/** \brief Generates an initial feasible point #FP_init. First we perform a change of variable to "tilde states",
    generate a feasible point, and then we go back to "bar states".

    \param[in] output_flag If output_flag = 1 outputs all initial guesses in a file output/FP.dat
*/
void WMG::form_FP_init(int output_flag)
{
    /** \brief State transition matrix. */
    double A[3*3];

    /** \brief Control matrix. */
    double B[3];

    /** \brief inv(Cp*B). This is a [2 x 2] diagonal matrix (which is invertible if #T^3/6-#h*#T is
     * not equal to zero). The two elements one the main diagonal are equal, and only one of them 
     * is stored, which is equal to
        \verbatim
        1/(T^3/6 - h*T)
        \endverbatim      
     */
    double iCpB;

    /** \brief inv(Cp*B)*Cp*A. This is a [2 x 6] matrix with the following structure
        \verbatim
        iCpB_CpA = [a b c 0 0 0;
                      0 0 0 a b c];

        a = iCpB
        b = iCpB*T
        c = iCpB*T^2/2      
        \endverbatim

     * Only a,b and c are stored.
     */      
    double iCpB_CpA[3];


    int k, k1, k2, k3=0;
    double tmp;

    k  = NUM_STATE_VAR*N; 
    k1 = k + 1;

    //------------------------------------
    /*A[0] = 1.0;*/  A[3]  = T[0];   A[6] = T[0]*T[0]/2;
    /*A[1] = 0.0;    A[4]  = 1.0;*/  A[7] = T[0];    
    /*A[2] = 0.0;    A[5]  = 0.0;    A[8] = 1.0;*/

    B[0] = T[0]*T[0]*T[0]/6 - h[0]*T[0]; 
    B[1] = T[0]*T[0]/2;         
    B[2] = T[0];             

    iCpB = 1/(T[0]*T[0]*T[0]/6 - h[0]*T[0]);

    iCpB_CpA[0] = iCpB;
    iCpB_CpA[1] = iCpB*T[0];
    iCpB_CpA[2] = iCpB*T[0]*T[0]/2;
    //------------------------------------
    

    //Vu(1:2)
    FP_init[k]  = -iCpB_CpA[0]*X_tilde[0] - iCpB_CpA[1]*X_tilde[1] - iCpB_CpA[2]*X_tilde[2] + iCpB*FS[ind[0]].p.x;
    FP_init[k1] = -iCpB_CpA[0]*X_tilde[3] - iCpB_CpA[1]*X_tilde[4] - iCpB_CpA[2]*X_tilde[5] + iCpB*FS[ind[0]].p.y;
 
    //Vc(1:6)
    FP_init[0] = X_tilde[0] + A[3]*X_tilde[1] + A[6]*X_tilde[2] + B[0]*FP_init[k];
    FP_init[1] =                   X_tilde[1] + A[7]*X_tilde[2] + B[1]*FP_init[k];
    FP_init[2] =                                     X_tilde[2] + B[2]*FP_init[k];
    FP_init[3] = X_tilde[3] + A[3]*X_tilde[4] + A[6]*X_tilde[5] + B[0]*FP_init[k1];
    FP_init[4] =                   X_tilde[4] + A[7]*X_tilde[5] + B[1]*FP_init[k1];
    FP_init[5] =                                     X_tilde[5] + B[2]*FP_init[k1];
    
    for (int i=1; i<N; i++)
    {
        k  = NUM_STATE_VAR*N + 2*i; 
        k1 = k + 1;
        k2 = NUM_STATE_VAR*i - NUM_STATE_VAR;
        k3 = NUM_STATE_VAR*i; //6*(i+1) - 6

        //------------------------------------
        /*A[0] = 1.0;*/  A[3]  = T[i];   A[6] = T[i]*T[i]/2;
        /*A[1] = 0.0;    A[4]  = 1.0;*/  A[7] = T[i];    
        /*A[2] = 0.0;    A[5]  = 0.0;    A[8] = 1.0;*/

        B[0] = T[i]*T[i]*T[i]/6 - h[i]*T[i]; 
        B[1] = T[i]*T[i]/2;         
        B[2] = T[i];             

        iCpB = 1/(T[i]*T[i]*T[i]/6 - h[i]*T[i]);

        iCpB_CpA[0] = iCpB;
        iCpB_CpA[1] = iCpB*T[i];
        iCpB_CpA[2] = iCpB*T[i]*T[i]/2;
        //------------------------------------


        FP_init[k]  = -iCpB_CpA[0]*FP_init[k2]   - iCpB_CpA[1]*FP_init[k2+1] - iCpB_CpA[2]*FP_init[k2+2] + iCpB*FS[ind[i]].p.x;
        FP_init[k1] = -iCpB_CpA[0]*FP_init[k2+3] - iCpB_CpA[1]*FP_init[k2+4] - iCpB_CpA[2]*FP_init[k2+5] + iCpB*FS[ind[i]].p.y;

        FP_init[k3]   = FP_init[k2]   + A[3]*FP_init[k2+1] + A[6]*FP_init[k2+2] + B[0]*FP_init[k];
        FP_init[k3+1] =                      FP_init[k2+1] + A[7]*FP_init[k2+2] + B[1]*FP_init[k];
        FP_init[k3+2] =                                           FP_init[k2+2] + B[2]*FP_init[k];
        FP_init[k3+3] = FP_init[k2+3] + A[3]*FP_init[k2+4] + A[6]*FP_init[k2+5] + B[0]*FP_init[k1];
        FP_init[k3+4] =                      FP_init[k2+4] + A[7]*FP_init[k2+5] + B[1]*FP_init[k1];
        FP_init[k3+5] =                                           FP_init[k2+5] + B[2]*FP_init[k1];

        // go back to bar states
        tmp           =  FS[ind[i-1]].ca*FP_init[k2] + FS[ind[i-1]].sa*FP_init[k2+3];
        FP_init[k2+3] = -FS[ind[i-1]].sa*FP_init[k2] + FS[ind[i-1]].ca*FP_init[k2+3];
        FP_init[k2]   = tmp;        
    }
    // go back to bar states (terminal state)
    tmp           =  FS[ind[N-1]].ca*FP_init[k3] + FS[ind[N-1]].sa*FP_init[k3+3];
    FP_init[k3+3] = -FS[ind[N-1]].sa*FP_init[k3] + FS[ind[N-1]].ca*FP_init[k3+3];   
    FP_init[k3]   = tmp; 


    /* Another version for going back to bar states
    for (int i=0; i<N; i++)
    {
        tmp                =  FS[ind[i]].ca*FP_init[6*(i+1)-6] + FS[ind[i]].sa*FP_init[6*(i+1)-3];
        FP_init[6*(i+1)-3] = -FS[ind[i]].sa*FP_init[6*(i+1)-6] + FS[ind[i]].ca*FP_init[6*(i+1)-3];        
        FP_init[6*(i+1)-6] =  tmp;        
    }    
    */

    if (output_flag)
    {
        if (counter == 0)
        {
            write_file(FP_init, NUM_VAR*N, 1, "output/FP.dat", "w"); 
        }
        else
        {
            write_file(FP_init, NUM_VAR*N, 1, "output/FP.dat", "a"); 
        }
    }
}



/** \brief Outputs the position of the CoM (on the floor) and ZMP (for plotting purposes). */
void WMG::output_CoM_ZMP()
{
    
    double tmp[2];
    
    if (counter == 0)
    {   
        tmp[0] = CoM.x; tmp[1] = CoM.y;
        write_file(tmp, 1, 2, "output/CoM.dat", "w"); 

        tmp[0] = ZMP.x; tmp[1] = ZMP.y;
        write_file(tmp, 1, 2, "output/ZMP.dat", "w"); 
    }
    else
    {
        tmp[0] = CoM.x; tmp[1] = CoM.y;
        write_file(tmp, 1, 2, "output/CoM.dat", "a"); 

        tmp[0] = ZMP.x; tmp[1] = ZMP.y;
        write_file(tmp, 1, 2, "output/ZMP.dat", "a");   
    }
}
