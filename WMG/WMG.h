/**
 * @file
 * @author Alexander Sherikov
 * @date 27.09.2011 18:57:37 MSD
 */


#ifndef WMG_H
#define WMG_H


/****************************************
 * INCLUDES 
 ****************************************/

#include <vector>

#include "common_const.h"


/****************************************
 * TYPEDEFS 
 ****************************************/
class FootStep;


/// @addtogroup gWMG_API
/// @{

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

    private:
        void slide();
};

///@}

#endif /*WMG_H*/
