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
//@{

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
        /**
         * @brief Allocate memory and initialize some of the parameters.
         *
         * @param[in] N preview window length
         * @param[in] hCoM_ Height of the Center of Mass [meter]
         * @param[in] gravity_ gravity [m/s^2]
         */
        smpc_parameters (
                const unsigned int, 
                const double, 
                const double gravity_ = 9.81);

        /**
         * @brief Default destructor
         */
        ~smpc_parameters();



// variables
        double hCoM;    /// Height of the CoM.

        double gravity; /// Norm of the acceleration due to gravity.

        double *T;      /// Preview sampling time

        double *h;      /// h = #hCoM/#gravity for each preview step
        double h0;      /// initial h

        /// Array of N absolute angles corresponding to supports in the preview window.
        double *angle;

        //@{
        /// Coordinates of N points satisfying constraints, 
        double *fp_x;
        double *fp_y;
        //@}
        
        //@{
        /// 2*N bounds for coordinates of ZMP position.
        double *lb;
        double *ub;
        //@}

        //@{
        /// N reference coordinates of ZMP.
        double *zref_x;
        double *zref_y;
        //@}


        /** Initial state. */
        smpc::state_orig init_state;

        /// A chunk of memory allocated for solution.
        double *X;
};



/**
 * @brief Default footstep constraints.
 */
class defConstraints
{
    public:
        /// Support constraints with safety margin.
        double ss[4];

        /// Both feet standing together.
        double ds[4];

        /// Automatically added DS.
        double auto_ds[4];

        /// Distance between reference points of the feet.
        double support_distance_y;


        void init()
        {
            // 2*HipOffsetY = 2*0.05
            support_distance_y = 0.1;

            // measured using ruler =)
            ss[0] = 0.09;
            ss[1] = 0.025;
            ss[2] = 0.03;
            ss[3] = 0.025;

            ds[0] = 0.09;
            ds[1] = 0.025 + support_distance_y/2;
            ds[2] = 0.03;
            ds[3] = 0.025 + support_distance_y/2;

            // SS constraints with increased safety margin
            auto_ds[0] = 0.07;
            auto_ds[1] = 0.025;
            auto_ds[2] = 0.025;
            auto_ds[3] = 0.025;
        };


        void init_FSR()
        {
            /*
               Positions of forse sensitve resistors (with respect to ankle frame
               of respective foot):

              Left foot                                              
                (0.07025 , 0.0299)  *---------* (0.07025 , -0.0231)   
                                    |         |                       
                                    |    ^    /
                                    |    |x  /
                                    | <--+   |
                                    |  y     |
                (-0.03025 , 0.0299) *--------* (-0.02965 , -0.0191)
                
                Rectangle:      
                ( 0.07  , 0.029)  ( 0.07  , -0.019)
                (-0.029 , 0.029)  (-0.029 , -0.019)


              Right foot                                                    
                (0.07025 , 0.0231)  *---------* (0.07025 , -0.0299)
                                    |         |                             
                                    \    ^    |
                                     \   |x   |
                                     |<--+    |
                                     | y      |
                (-0.03025 , 0.0191)  *--------* (-0.02965 , -0.0299)

                Rectangle:
                ( 0.07  , 0.019)  ( 0.07  , -0.029)
                (-0.029 , 0.019)  (-0.029 , -0.029)

              Both can be represented by the same rectangle:

                *-------------*
                |      ^ 0.07 |    Note that the distance between
                |      |      |    reference points must be changed
                |      |      |
                |      |      |
                |0.024 | 0.024|
                |<-----*----->|
                |      |      |
                |      V 0.029|
                *-------------*
             */


            // 2*HipOffsetY + 2*(0.024-0.019) = 2*(0.05+0.005)
            support_distance_y = 0.11;


            ss[0] = 0.07;
            ss[1] = 0.024;
            ss[2] = 0.029;
            ss[3] = 0.024;

            ds[0] = 0.07;
            ds[1] = 0.024 + support_distance_y/2;
            ds[2] = 0.029;
            ds[3] = 0.024 + support_distance_y/2;

            // SS constraints with increased safety margin
            auto_ds[0] = 0.065;
            auto_ds[1] = 0.020;
            auto_ds[2] = 0.025;
            auto_ds[3] = 0.020;
        };
};



/** 
 * @brief Defines the parameters of the Walking Pattern Generator. 
 */
class WMG
{
    public:

        // --------------------

        /**
         * @brief Initializes a WMG object.
         *
         * @param[in] N_ Number of sampling times in a preview window
         * @param[in] T_ Sampling time [ms.]
         * @param[in] step_height_ step height (for interpolation of feet movements) [meter]
         * @param[in] bezier_weight_1_ see #bezier_weight_1
         * @param[in] bezier_weight_2_ see #bezier_weight_2
         * @param[in] bezier_inclination_1_ see #bezier_inclination_1
         * @param[in] bezier_inclination_2_ see #bezier_inclination_2
         * @param[in] use_fsr_constraints Use constraints determined by positions of FSR sensors.
         *
         * @note The 'CenterOfPressure' value of the NAO robot can be computed correctly only 
         * when CoP lies within the polygon determined by force sensitive resistors (FSR) located 
         * on each foot. From the documentation: "If the real center of pressure go outside the 
         * 4 FSR parallelogram, the value of weight and centre of pressure could become bad, due 
         * to internal mechanical constraint". The rectanle that fits within the FSR polygon is
         * smaller than our default constraint rectangle, which was selected based on the 
         * measurements of the feet. See also #defConstraints.
         */
        WMG (   const unsigned int N_,
                const unsigned int T_, 
                const double step_height_ = 0.0135,
                double bezier_weight_1_ = 1.0,
                double bezier_weight_2_ = 2.0,
                double bezier_inclination_1_ = 0.01,
                double bezier_inclination_2_ = 0.008,
                bool use_fsr_constraints = false);

        /** @brief Destructor. */
        ~WMG();


        /** 
         * @brief Set default parameters of footsteps, a wrapper around not so safe 
         *  #setFootstepParametersMS function.
         *
         * @param[in] def_periods default number of sampling periods in a support
         *              (SS or DS depending on the type of added footstep)
         * @param[in] ds_periods default number of sampling periods in an 
         *              automatically generated DS.
         * @param[in] ds_number number of DS to be generated automatically.
         * @param[in] use_user_constraints_ use #user_constraints and #user_constraints_auto_ds 
         *              instead  of the default constraints.
         */
        void setFootstepParameters (
                const unsigned int def_preiods,
                const unsigned int ds_periods, 
                const unsigned int ds_number,
                bool use_user_constraints_ = false);


        /** 
         * @brief Set default parameters of footsteps.
         *
         * @param[in] def_time_ms_ default time spent in a support (SS or DS depending on the 
         *              type of added footstep)
         * @param[in] ds_time_ms_ default time spent in an automatically generated DS.
         * @param[in] ds_number number of DS to be generated automatically.
         * @param[in] use_user_constraints_ use #user_constraints and #user_constraints_auto_ds 
         *              instead  of the default constraints.
         */
        void setFootstepParametersMS (
                const unsigned int def_time_ms_,
                const unsigned int ds_time_ms_, 
                const unsigned int ds_number,
                bool use_user_constraints_ = false);


        /**
         * @brief Adds a footstep to FS.
         *
         * @param[in] x_relative x_relative X position [meter] relative to the previous footstep.
         * @param[in] y_relative y_relative Y position [meter] relative to the previous footstep.
         * @param[in] angle_relative angle_relative Angle [rad.] relative to the previous footstep.
         * @param[in] type (optional) type of the footstep.
         *
         * @note Coordinates and angle are treated as absolute for the first step in the preview window.
         */
        void addFootstep(
                const double, 
                const double, 
                const double, 
                fs_type type = FS_TYPE_AUTO);


        /**
         * @brief Forms a preview window.
         *
         * @return WMG_OK or WMG_HALT (simulation must be stopped)
         */
        WMGret formPreviewWindow (smpc_parameters &);


        // --------------------


        /**
         * @brief Determine position and orientation of feet
         *
         * @param[in] shift_from_current_ms a positive shift in time (ms.) from the current time
         *  (allows to get positions for the future supports)
         * @param[out] left_foot_pos 4x4 homogeneous matrix, which represents position and orientation
         * @param[out] right_foot_pos 4x4 homogeneous matrix, which represents position and orientation
         *
         * @attention This function requires the walking pattern to be started and finished
         * by single support.
         *
         * @attention Cannot be called on the first or last SS  =>  must be called after 
         * FormPreviewWindow().
         */
        void getFeetPositions (
                const unsigned int shift_from_current_ms, 
                double * left_foot_pos, 
                double * right_foot_pos);


        /**
         * @brief Checks if the support foot switch is needed.
         *
         * @return true if the support foot must be switched. 
         */
        bool isSupportSwitchNeeded ();


        /**
         * @brief Changes position of the next SS.
         *
         * @param[in] posture a 4x4 homogeneous matrix representing new position and orientation
         * @param[in] zero_z_coordinate set z coordinate to 0.0
         *
         * @todo DS must be adjusted as well.
         */
        void changeNextSSPosition (const double *posture, const bool zero_z_coordinate);


        // --------------------

        /**
         * @brief Outputs the footsteps in FS to a Matlab/Octave script for plotting.
         *
         * @param[in] filename output file name.
         * @param[in] plot_ds enable/disable plotting of double supports
         */
        void FS2file(const std::string filename, const bool plot_ds = true);


        /**
         * @brief Return coordinates of footstep reference points and rotation 
         * angles of footsteps (only for SS).
         *
         * @param[out] x_coord x coordinates
         * @param[out] y_coord y coordinates
         * @param[out] angle_rot angles
         */
        void getFootsteps(
                std::vector<double> &x_coord,
                std::vector<double> &y_coord,
                std::vector<double> &angle_rot);

        // --------------------


// variables
        /// A vector of footsteps. 
        std::vector<footstep> FS; 


        /// Number of iterations in a preview window.
        unsigned int N;

        unsigned int *T_ms;
        unsigned int sampling_period;

        
        /// This is the step in FS that is at the start of the current preview window.
        int current_step_number;

        /// The first step in the current preview window.
        int first_preview_step;

        /// The maximum height, that can be reached by a swing foot.
        double step_height;


        /// Default SS and DS constraints.
        defConstraints def_constraints;

        ///@{
        /// Constraints given by the user, initialized to the default values on initialization.
        double user_constraints[4];
        double user_constraints_auto_ds[4];
        bool use_user_constraints;
        ///@}



        //@{
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
        //@}

    private:
        void getDSFeetPositions (const int, double *, double *);
        void getSSFeetPositions (const int, const double, double *, double *);
        void getSSFeetPositionsBezier (const int, const double, double *, double *);
        int getNextSS (const int, const fs_type type = FS_TYPE_AUTO);
        int getPrevSS (const int, const fs_type type = FS_TYPE_AUTO);

        unsigned int def_time_ms;
        unsigned int ds_time_ms;
        unsigned int ds_num;

        unsigned int last_time_decrement;
};
//@}

#endif /*WMG_H*/
