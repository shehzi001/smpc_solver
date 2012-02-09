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
        bool isSupportSwitchNeeded ();
        WMGret formPreviewWindow ();
        void FS2file(const std::string, const bool plot_ds = true);

        void getFeetPositions (
                const int, 
                const int, 
                const int, 
                double *, 
                double *);


        void initABMatrices (const double);
        void calculateNextState (smpc::control&, smpc::state_orig&);

        void getFootsteps(
                std::vector<double> &,
                std::vector<double> &,
                std::vector<double> &);

        void correctNextSSPosition (const double *);


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

        /** Initial state. */
        smpc::state_tilde X_tilde;
        
        /** Initial state. */
        smpc::state_orig init_state;

        /// A storage for the controls, that must be applied to reach the next state.
        smpc::control next_control;

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

        /// The first step in the current preview window.
        int first_preview_step;

        /// The maximum height, that can be reached by a swing foot.
        double step_height;

        double def_ss_constraint[4];
        double def_ds_constraint[4];

    private:
        void getDSFeetPositions (const int, double *, double *);
        void getSSFeetPositions (const int, const double, double *, double *);
        int getNextSS (const int, const fs_type type = FS_TYPE_AUTO);
        int getPrevSS (const int, const fs_type type = FS_TYPE_AUTO);

        double addstep_constraint[4];
        int def_repeat_times;
        int def_ds_num;

        
        ///@{
        /// State and control matrices, that can be used to determine the next
        /// state based on the current state and the controls.
        double *A;
        double *B;
        ///@}
};

///@}

#endif /*WMG_H*/
