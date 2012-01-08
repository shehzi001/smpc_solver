/**
 * @file
 * @brief Interface of the library.
 *
 * @author Alexander Sherikov
 * @date 02.09.2011 00:22:48 MSD
 */


#ifndef SMPC_SOLVER_H
#define SMPC_SOLVER_H

class qp_solver;


/// @addtogroup gAPI 
/// @{

/// Number of state variables.
#define SMPC_NUM_STATE_VAR 6
/// Number of control variables.
#define SMPC_NUM_CONTROL_VAR 2
/// Total number of variables.
#define SMPC_NUM_VAR 8

namespace smpc
{
    class solver;



    /**
     * @brief An abstract class representing state.
     */
    class state
    {
        public:
            /**
             *  \verbatim
                state[0] - x CoM or ZMP position [meter]
                state[1] - x CoM velocity [meter/s]
                state[2] - x CoM acceleration [meter/s^2]
                state[3] - y CoM or ZMP position [meter]
                state[4] - y CoM velocity [meter/s]
                state[5] - y CoM acceleration [meter/s^2]
                \endverbatim
             */
            double state_vector[SMPC_NUM_STATE_VAR];


            /**
             * @brief Default constructor (everything set to 0).
             */
            state();

            // -------------------------------


            /**
             * @brief X coordinate.
             */
            double &x  () {return state_vector[0];};
            /**
             * @brief Velocity along X axis.
             */
            double &vx () {return state_vector[1];};
            /**
             * @brief Acceleration along X axis.
             */
            double &ax () {return state_vector[2];};

            /**
             * @brief Y coordinate.
             */
            double &y  () {return state_vector[3];};
            /**
             * @brief Velocity along Y axis.
             */
            double &vy () {return state_vector[4];};
            /**
             * @brief Acceleration along Y axis.
             */
            double &ay () {return state_vector[5];};


            // -------------------------------
            
           
            /**
             * @brief Returns the next state.
             *  
             * @param[in] smpc_solver initialized solver
             */
            virtual void get_next_state (const solver &smpc_solver) = 0;


            /**
             * @brief The same as #get_next_state, but takes an additional 
             * parameter - index of the desired state in the state vector.
             *  
             * @param[in] smpc_solver initialized solver
             * @param[in] ind index of a state [0 : N-1].
             */
            virtual void get_state (const solver &smpc_solver, const int ind) = 0;


            // -------------------------------


            /**
             * @brief Set coordinates, velocities and accelerations 
             * are assumed to be zero.
             *
             * @param[in] x_ x CoM/ZMP position [meter]
             * @param[in] y_ y CoM/ZMP position [meter]
             */
            void set (double x_, double y_);


            /**
             * @brief Set state to given values.
             *
             * @param[in] x_  x CoM/ZMP position [meter]
             * @param[in] vx_ x CoM velocity [meter/s]
             * @param[in] ax_ x CoM acceleration [meter/s^2]
             * @param[in] y_  y CoM/ZMP position [meter]
             * @param[in] vy_ y CoM velocity [meter/s]
             * @param[in] ay_ y CoM acceleration [meter/s^2]
             */
            void set (double x_, double vx_, double ax_,
                     double y_, double vy_, double ay_);
    };


    /**
     * @brief A container for a state in the original form:
     *  \verbatim
        x CoM position [meter]
        x CoM velocity [meter/s]
        x CoM acceleration [meter/s^2]
        y CoM position [meter]
        y CoM velocity [meter/s]
        y CoM acceleration [meter/s^2]
        \endverbatim
     */
    class state_orig : public state
    {
        public:
            state_orig () : state () {};

            /// @{
            /// Refer to the base class state for description.
            void get_next_state (const solver &smpc_solver);
            void get_state (const solver &smpc_solver, const int ind);
            /// @}
    };

    
    /**
     * @brief A container for a state in @ref pX_tilde "tilde" form (after
     * variable substitution):
     *  \verbatim
        x ZMP position [meter]
        x CoM velocity [meter/s]
        x CoM acceleration [meter/s^2]
        y ZMP position [meter]
        y CoM velocity [meter/s]
        y CoM acceleration [meter/s^2]
        \endverbatim
     */
    class state_tilde : public state
    {
        public:
            state_tilde () : state () {};

            /// @{
            /// Refer to the base class state for description.
            void get_next_state (const solver &smpc_solver);
            void get_state (const solver &smpc_solver, const int ind);
            /// @}
    };


    /**
     * @brief A container for a control vector.
     */
    class control
    {
        public:
            /**
             *  \verbatim
                control[0] - jerk along x axis
                control[1] - jerk along y axis
                \endverbatim
             */
            double control_vector[SMPC_NUM_CONTROL_VAR];


            control();

            // -------------------------------


            /**
             * @brief Jerk along Y axis.
             */
            double &jx () {return control_vector[0];};
            /**
             * @brief Jerk along Y axis.
             */
            double &jy () {return control_vector[1];};


            // -------------------------------


            /**
             * @brief Returns the controls that must be applied to reach the 
             *  next state.
             *
             * @param[in] smpc_solver initialized solver
             */
            void get_first_controls (const solver &smpc_solver);


            /**
             * @brief The same as #get_first_controls, but takes an additional 
             * parameter - the index of the control inputs. Control[ind] is 
             * applied to State[ind-1] to reach State[ind].
             *
             * @param[in] smpc_solver initialized solver
             * @param[in] ind index of control inputs [0 : N-1].
             */
            void get_controls (const solver &smpc_solver, const int ind);
    };




    /**
     * @brief API of the sparse MPC solver.
     */
    class solver
    {
        public:


            // -------------------------------


            /** @brief Constructor: initialize an active set method solver.
             *
                @param[in] N Number of sampling times in a preview window
                @param[in] Alpha Velocity gain
                @param[in] Beta Position gain
                @param[in] Gamma Jerk gain
                @param[in] regularization regularization
                @param[in] tol tolerance
            */
            solver (
                    const int N, 
                    const double Alpha = 150.0, 
                    const double Beta = 2000.0, 
                    const double Gamma = 1.0,
                    const double regularization = 0.01,
                    const double tol = 1e-7);


            /** @brief Constructor: initialize an interior-point method solver.
             *
             * @param[in] N Number of sampling times in a preview window
             * @param[in] max_iter maximum number of internal loop iterations
             * @param[in] Alpha Velocity gain
             * @param[in] Beta Position gain
             * @param[in] Gamma Jerk gain
             * @param[in] regularization regularization
             * @param[in] tol tolerance (internal loop)
             * @param[in] tol_out tolerance of the outer loop, which resolves
             *                    the problem with new t (kappa) parameter.
             * @param[in] t logarithmic barrier parameter
             * @param[in] mu multiplier of t, >1.
             * @param[in] bs_alpha backtracking search parameter 0 < alpha < 0.5
             * @param[in] bs_beta  backtracking search parameter 0 < beta < 1
             *
             * @todo Which defaults are good?
             */
            solver (
                    const int N, 
                    const int max_iter, // no default to avoid ambiguity
                    const double Alpha = 150.0, 
                    const double Beta = 2000.0, 
                    const double Gamma = 1.0,
                    const double regularization = 0.01,
                    const double tol = 1e-3,
                    const double tol_out = 1e-2,
                    const double t = 100,
                    const double mu = 15,
                    const double bs_alpha = 0.01,
                    const double bs_beta = 0.5);


            ~solver();


            // -------------------------------


            /** @brief Initializes quadratic problem.

                @param[in] T sampling time for each time step [sec.]
                @param[in] h height of the center of mass divided by gravity for each time step
                @param[in] h_initial initial value of height of the center of mass divided by gravity
                @param[in] angle rotation angle for each state relative to the world frame
                @param[in] zref_x reference values of x coordinate of ZMP
                @param[in] zref_y reference values of y coordinate of ZMP
                @param[in] lb array of lower bounds for coordinates of ZMP
                @param[in] ub array of upper bounds for coordinates of ZMP
            */
            void set_parameters (
                    const double* T,
                    const double* h,
                    const double h_initial,
                    const double* angle,
                    const double* zref_x,
                    const double* zref_y,
                    const double* lb,
                    const double* ub);


            /** @brief Generates an initial feasible point. 

                @param[in] x_coord x coordinates of points satisfying constraints
                @param[in] y_coord y coordinates of points satisfying constraints
                @param[in] init_state initial state
                @param[in,out] X solution of optimization problem
             */
            void form_init_fp (
                    const double *x_coord,
                    const double *y_coord,
                    const state_orig &init_state,
                    double* X);


            /**
             * @brief Solve QP problem.
             *
             * @return A negative number on error. Number of activated constraints 
             *         for active set method and 0 for interior-point method on success.
             */
            int solve ();
       

            // -------------------------------

            /**
             * @brief Internal representation.
             */
            qp_solver *qp_sol;
    };
}
/// @}

#endif /*SMPC_SOLVER_H*/
