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


#ifndef WMG_H
#define WMG_H

#include <vector>
#include <iostream> // cerr, endl
#include <math.h> // for cos and sin


#define NUM_STATE_VAR 6
#define NUM_CONTROL_VAR 2
#define NUM_VAR 8


/** \brief Defines a 2D point. */
class Point2D
{
    public:
        Point2D();
        Point2D(double _x, double _y);
        void print();
        void set(double _x, double _y);



        /** \brief x position [meter]*/
        double x;

        /** \brief y position [meter] */
        double y;
};



/** \brief Defines rectangular constraints (of the form D*z <= d) for
    the ZMP.

    \note D is a [4 x 2] matrix stored column-wise (Fortran style). 
    D is always initialized as (later it could be rotated using D*R'):
    \verbatim
    D[0] =  1.0; D[4] =  0.0;
    D[1] =  0.0; D[5] =  1.0;
    D[2] = -1.0; D[6] =  0.0;
    D[3] =  0.0; D[7] = -1.0;
    \endverbatim
*/
class RectangularConstraint_ZMP
{
    public:
        RectangularConstraint_ZMP();
        void set_size(double * _d);
        void print();
        void rotate_translate(double ca, double sa, Point2D p);
        void Constraints2Vert();



        /** \brief Matrix of the constraints D*z <= d (where z is a 2D point). This variable is formed
            internaly.
            
            \note D is a [4 x 2] matrix stored column-wise (Fortran style). D takes the form 
            \verbatim
            D = [ 1  0;
                    0  1;
                   -1  0;
                    0 -1] * [cos(a) -sin(a)
                             sin(a)  cos(a)].
            \endverbatim
         */
        double D[4*2];

        /** \brief Size of the support polygon for a single support (USER INPUT)
                
            \verbatim
             ------------------------------------------------
             |                    |                         |
             |                    |d0(2)                    |
             |                    |                         |
             |      d0(3)         |          d0(1)          |
             |------------------- p ------------------------|
             |                    |                         |
             |                    |                         |
             |                    |d0(4)                    |
             |                    |                         |
             ------------------------------------------------
            \endverbatim
         */
        double d0[4];

        /** \brief Vector of the constraints D*z <= d (where z is a 2D point). This variable is formed internaly. 
                
            \note d is a [4 x 1] vector. d takes the form #d = #d0 - #D*p. p is a 2D reference point
            (expressed with respect to the world frame) for the rectangular constraint (see the sketch for
            #d0).
                
         */
        double d[4];
         
        std::vector<Point2D> vert;
};



/** \brief Defines a footstep. */
class FootStep
{
    public:
        FootStep();
        FootStep(double _angle, Point2D _p, int _nSS, int _n, int _RL);
        FootStep(double _angle, Point2D _p, int _nSS, int _n, int _RL, double *_d);
        void set(double _angle, Point2D _p, int _nSS, int _n, int _RL);
        void print();



        /** \brief Position (in the world frame) of a footstep [meter]. */
        Point2D p;

        /** \brief Angle (relative to the world frame) of a footstep [rad.]. */
        double angle;

        /** \brief cos(angle). */
        double ca;

        /** \brief sin(angle). */
        double sa;

        /** \brief Constraints defining the Polyogn of Support (PoS). */
        RectangularConstraint_ZMP ctr;

        /** \brief Offset from the "point of interest of a constraint" [meter] defined in the local frame. 
                
            \note This is used in order to define ZMP_ref
         */
        Point2D Offset;

        /** \brief If RL = 1 the left foot is in support. If RL = -1 the right foot is in support. */
        int RL;

        /** \brief Number of (preview window) iterations in Single Support. */
        int nSS;

        /** \brief Total number of (preview window) iterations, i.e., nSS + nDS. */
        int n;
};


/** \brief Defines the parameters of the Walking Pattern Generator. */
class WMG
{
    public:
        WMG();
        ~WMG();
        void init (int _N, double _T, double _hCoM);
        void AddFootstep(double x_relative, double y_relative, double angle_relative, int _nSS, int _n, int _RL, double *_d);
        void AddFootstep(double x_relative, double y_relative, double angle_relative, int _nSS, int _n);
        void AddFootstep(double x_relative, double y_relative, double angle_relative);
        void FormPreviewWindow();
        void slide();
        void print_ind();
        void FS2file();
        void output_CoM_ZMP();



        /** \brief A vector of footsteps. */
        std::vector<FootStep> FS; 


        /** \brief Number of iterations in a preview window. */
        int N;

        /** \brief Preview sampling time  */
        double *T;


        /** \brief Height of the CoM. */
        double hCoM;

        /** \brief Norm of the acceleration due to gravity. For the moment gravity is set equal to 9.81. */
        double gravity;
        
        /** \brief h = #hCoM/#gravity. */
        double *h;
        
        /** \brief Contains the reference profile for the ZMP. */
        double *ZMP_ref;

        /** \brief Initial feasible point with respect to the equality and inequality constraints. */
        double *X;

        /** \brief Initial state. */
        double X_tilde[NUM_STATE_VAR];

        /** \brief Initial angle */
        double init_angle;

        /** \brief Indexes of footsteps appearing in the current preview window
         */
        int *ind; 

        /** \brief Position of the CoM (on the floor). */
        Point2D CoM;

        /** \brief Position of the ZMP. */
        Point2D ZMP;


        /** \brief This is the step in FS that is at the start of the current preview window. I don't need
        to pop_front steps from FS and that is why I use stl vector and not stl deque). */
        int current_step_number;

        /** \brief Iteration counter. Starts from 0 and increases (by one) with each preview window. */
        int counter;

        /** \brief If halt = 1 then stop the execution (else keep going) */
        int halt;


        double *angle;
        double *zref_x;
        double *zref_y;
        double *lb;
        double *ub;
};


#endif
