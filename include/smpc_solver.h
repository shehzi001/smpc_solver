/**
 * @file
 * @brief Interface of the library.
 *
 * @author Alexander Sherikov
 * @date 02.09.2011 00:22:48 MSD
 */


#ifndef SMPC_SOLVER_H
#define SMPC_SOLVER_H

class qp_as;


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
    class solver_as;



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
            virtual void get_next_state (const solver_as &smpc_solver) = 0;


            /**
             * @brief The same as #get_next_state, but takes an additional 
             * parameter - index of the desired state in the state vector.
             *  
             * @param[in] smpc_solver initialized solver
             * @param[in] ind index of a state [0 : N-1].
             */
            virtual void get_state (const solver_as &smpc_solver, const int ind) = 0;


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
            void get_next_state (const solver_as &smpc_solver);
            void get_state (const solver_as &smpc_solver, const int ind);
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
            void get_next_state (const solver_as &smpc_solver);
            void get_state (const solver_as &smpc_solver, const int ind);
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
            void get_first_controls (const solver_as &smpc_solver);


            /**
             * @brief The same as #get_first_controls, but takes an additional 
             * parameter - the index of the control inputs. Control[ind] is 
             * applied to State[ind-1] to reach State[ind].
             *
             * @param[in] smpc_solver initialized solver
             * @param[in] ind index of control inputs [0 : N-1].
             */
            void get_controls (const solver_as &smpc_solver, const int ind);
    };




    /**
     * @brief API of the sparse MPC solver.
     */
    class solver_as
    {
        public:

            /** @brief Constructor: initialize an active set method solver.
             *
                @param[in] N Number of sampling times in a preview window
                @param[in] Alpha Velocity gain
                @param[in] Beta Position gain
                @param[in] Gamma Jerk gain
                @param[in] regularization regularization
                @param[in] tol tolerance
            */
            solver_as (
                    const int N, 
                    const double Alpha = 150.0, 
                    const double Beta = 2000.0, 
                    const double Gamma = 1.0,
                    const double regularization = 0.01,
                    const double tol = 1e-7);


            ~solver_as();


            // -------------------------------

            
            /**
             * @brief Enable floating point exceptions: die rather than process 
             * incorrect data. SIGFPE is sent to the program, when an error occurs.
             * @note If feenableexcept() is not present on the system, the function 
             * does nothing.
             */
            void enable_fexceptions();


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


            /**
             * @brief Changes parameters which control the logic of the solver.
             *
             * @param[in] max_added_constraints_num limit the number of added constraints 
             *                                      (NOT the size of active set).
             * @param[in] constraint_removal_enabled enable/disable removal of activated constraints.
             *
             * @note These parameters affect the time required for solution. If the number
             * of added constraints is less than (length of preview window)*2 or constraint 
             * removal is disabled, the solution is approximate. How good is this approximation 
             * depends on the problem.
             */
            void set_limits (
                    const unsigned int max_added_constraints_num,
                    const bool constraint_removal_enabled);


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
             */
            void solve ();
       

            // -------------------------------

            /**
             * @brief Number of added constraints (the constraints, that were
             * removed are also counted).
             *
             * @note Updated by #solve function.
             */
            unsigned int added_constraints_num;

            /**
             * @brief Number of removed constraints.
             *
             * @note Updated by #solve function.
             */
            unsigned int removed_constraints_num;

            /**
             * @brief The final size of the active set.
             *
             * @note Updated by #solve function.
             */
            unsigned int active_set_size;
            

            /**
             * @brief Internal representation.
             */
            qp_as *qp_sol;
    };
}
/// @}

#endif /*SMPC_SOLVER_H*/
