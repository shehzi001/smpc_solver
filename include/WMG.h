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

#include <string>
#include <vector>

#include "smpc_solver.h"

/****************************************
 * DEFINES
 ****************************************/



/****************************************
 * TYPEDEFS 
 ****************************************/
class FootStep;


/// @addtogroup gWMG_API
/// @{

enum WMGret
{
    WMG_OK,
    WMG_HALT
};

enum fs_type
{
    FS_TYPE_AUTO,   // let the library decide
    FS_TYPE_SS_L,   // left single support
    FS_TYPE_SS_R,   // right single support
    FS_TYPE_DS      // double support
};

enum swing_foot_pos_type
{
    WMG_SWING_2D_PARABOLA
};


/** \brief Defines the parameters of the Walking Pattern Generator. */
class WMG
{
    public:
// methods
        WMG();
        ~WMG();
        void init (const int);
        void init_param (
                const double, 
                const double,
                const double step_height_ = 0.0135);
        void AddFootstep(
                const double, 
                const double, 
                const double, 
                const int, 
                const int, 
                const double *, 
                const fs_type type = FS_TYPE_AUTO);
        void AddFootstep(
                const double, 
                const double, 
                const double, 
                const int, 
                const int, 
                const fs_type type = FS_TYPE_AUTO);
        void AddFootstep(
                const double, 
                const double, 
                const double, 
                fs_type type = FS_TYPE_AUTO);
        WMGret FormPreviewWindow(bool *switch_foot = NULL);
        void FS2file(const std::string, const bool plot_ds = true);

        void getSwingFootPosition (
                const swing_foot_pos_type,
                const int,
                const int,
                double *,
                double *);


        void initABMatrices (const double);
        void calculateNextState (const double *, double *);
        void initState (const double, const double, double *);

        void getFootsteps(
                std::vector<double> &,
                std::vector<double> &,
                std::vector<double> &);



// variables
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
        
        /// A chunk of memory allocated for solution.
        double *X;

        /** \brief Initial state. */
        double X_tilde[SMPC_NUM_STATE_VAR];

        /// Array of N absolute angles corresponding to supports in the preview window.
        double *angle;

        ///@{
        /// Coordinates of N points satisfying constraints, 
        double *fp_x;
        double *fp_y;
        ///@}
        
        ///@{
        /// 2*N bounds for coordinates of ZMP position.
        double *lb;
        double *ub;
        ///@}

        ///@{
        /// N reference coordinates of ZMP.
        double *zref_x;
        double *zref_y;
        ///@}


        /** This is the step in FS that is at the start of the current preview window. */
        int current_step_number;

        /// Reference foot (R/L), even in DS we need to choose a reference foot.
        fs_type current_reference_foot;

        /// The maximum height, that can be reached by a swing foot.
        double step_height;

    private:
        fs_type getSwingFootNextPrevPos (double *, double *, int *, int *);
        int getNextSS (const int, const fs_type type = FS_TYPE_AUTO);
        int getPrevSS (const int, const fs_type type = FS_TYPE_AUTO);

        double def_constraint[4];
        int def_repeat_times;

        double def_ds_constraint[4];
        int def_ds_num;

        
        ///@{
        /// State and control matrices, that can be used to estimate a
        /// @ref pX_tilde "state" based on the current state and controls.
        double *A;
        double *B;
        ///@}
};

///@}

#endif /*WMG_H*/
