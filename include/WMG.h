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
class footstep;

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



/**
 * @brief A container for parameters of the SMPC solver.
 */
class smpc_parameters
{
    public:
        smpc_parameters (
                const unsigned int, 
                const double, 
                const double gravity_ = 9.81);
        ~smpc_parameters();



// variables
        double hCoM;    /// Height of the CoM.

        double gravity; /// Norm of the acceleration due to gravity.

        double *T;      /// Preview sampling time

        double *h;      /// h = #hCoM/#gravity for each preview step
        double h0;      /// initial h

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


        /** Initial state. */
        smpc::state_orig init_state;

        /// A chunk of memory allocated for solution.
        double *X;
};



/**
 * @brief Inverted pendulum model
 */
class IPM
{
    public:
        IPM (const double);
        ~IPM ();
        void calculateNextState (smpc::control&, smpc::state_orig&);


        ///@{
        /// State and control matrices, that can be used to determine the next
        /// state based on the current state and the controls.
        double *A;
        double *B;
        ///@}

        smpc::control control_vector;
        smpc::state_orig state_vector;
};



/** \brief Defines the parameters of the Walking Pattern Generator. */
class WMG
{
    public:
// methods
        WMG (   const unsigned int,
                const unsigned int, 
                const double step_height_ = 0.0135,
                double bezier_weight_1_ = 1.0,
                double bezier_weight_2_ = 2.0,
                double bezier_inclination_1_ = 0.01,
                double bezier_inclination_2_ = 0.008);
        ~WMG();


        void setFootstepDefaults(
                const unsigned int, 
                const unsigned int,
                const double * constraints = NULL);
        void setFootstepDefaults(
                const unsigned int, 
                const unsigned int, 
                const unsigned int,
                const double * constraints = NULL);

        void addFootstep(
                const double, 
                const double, 
                const double, 
                fs_type type = FS_TYPE_AUTO);


        bool isSupportSwitchNeeded ();

        WMGret formPreviewWindow (smpc_parameters &);

        void FS2file(const std::string, const bool plot_ds = true);

        void getFeetPositions (const unsigned int, double *, double *);

        void getFootsteps(
                std::vector<double> &,
                std::vector<double> &,
                std::vector<double> &);

        void changeNextSSPosition (const double *, const bool);


// variables
        /** \brief A vector of footsteps. */
        std::vector<footstep> FS; 


        /** \brief Number of iterations in a preview window. */
        unsigned int N;

        unsigned int *T_ms;
        unsigned int sampling_period;

        
        /** This is the step in FS that is at the start of the current preview window. */
        int current_step_number;

        /// The first step in the current preview window.
        int first_preview_step;

        /// The maximum height, that can be reached by a swing foot.
        double step_height;

        double def_ss_constraint[4];
        double def_auto_ds_constraint[4];

        /// @{
        /**
         * The foot trajectories, that are build using Bezier curve have four control
         * points (0,1,3,4). The first and the last points are defined by positions
         * of adjacent steps of a foot on the floor. The 1 and 2 points are computed
         * depending on the given parameters.
         *
         * The step height is fixed in the middle (theta = 0.5) of the step, note,
         * that this is not the geometrical middle.
         *
         * The control points have weights, the weights of the first and the last points
         * are always equal to 1.0. The weights of other two points can be given by the
         * user.
         *
         * The 1 and 2 points also have inclination, i.e. the shift along y axis (in 
         * the frame fixed in the step reference point). The values of inclinations are 
         * given in meters.
         */
        double bezier_weight_1;
        double bezier_weight_2;
        double bezier_inclination_1;
        double bezier_inclination_2;
        /// @}

    private:
        void getDSFeetPositions (const int, double *, double *);
        void getSSFeetPositions (const int, const double, double *, double *);
        void getSSFeetPositionsBezier (const int, const double, double *, double *);
        int getNextSS (const int, const fs_type type = FS_TYPE_AUTO);
        int getPrevSS (const int, const fs_type type = FS_TYPE_AUTO);


        double addstep_constraint[4];

        unsigned int def_time_ms;
        unsigned int ds_time_ms;
        unsigned int ds_num;

        unsigned int last_time_decrement;
};

///@}

#endif /*WMG_H*/
