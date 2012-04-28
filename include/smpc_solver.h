/**
 * @file
 * @brief Interface of the library.
 *
 * @author Alexander Sherikov
 * @date 02.09.2011 00:22:48 MSD
 */


#ifndef SMPC_SOLVER_H
#define SMPC_SOLVER_H

#include <vector>


class qp_as;
class qp_ip;


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
    // -------------------------------

    
    /**
     * @brief Enable floating point exceptions: die rather than process 
     * incorrect data. SIGFPE is sent to the program, when an error occurs.
     * @note If feenableexcept() is not present on the system, the function 
     * does nothing.
     */
    void enable_fexceptions();


    // -------------------------------


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
                @param[in] gain_position Position gain (Alpha)
                @param[in] gain_velocity Velocity gain (Beta)
                @param[in] gain_acceleration Acceleration gain (Gamma)
                @param[in] gain_jerk Jerk gain (Eta)
                @param[in] tol tolerance
                @param[in] obj_computation_on compute and keep values of the objective function
                @param[in] max_added_constraints_num limit the number of added constraints 
                        (NOT the size of active set), no limit if set to 0.
                @param[in] constraint_removal_on enable/disable removal of activated constraints.

              @note #max_added_constraints_num and #constraint_removal_on affect the time required 
              for solution. If the number of added constraints is less than (length of preview window)*2 
              or constraint removal is disabled, the solution is approximate. How good is this 
              approximation depends on the problem.
             */
            solver_as (
                    const int N, 
                    const double gain_position = 2000.0, 
                    const double gain_velocity = 150.0, 
                    const double gain_acceleration = 0.02,
                    const double gain_jerk = 1.0,
                    const double tol = 1e-7,
                    const bool obj_computation_on = false,
                    const unsigned int max_added_constraints_num = 0,
                    const bool constraint_removal_on = true);


            ~solver_as();


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



            ///@{
            /** @brief Generates an initial feasible point. 

                @param[in] x_coord x coordinates of points satisfying constraints
                @param[in] y_coord y coordinates of points satisfying constraints
                @param[in] init_state initial state (#state_orig or #state_tilde)
                @param[in,out] X solution of optimization problem
             */
            void form_init_fp (
                    const double *x_coord,
                    const double *y_coord,
                    const state_orig &init_state,
                    double* X);

            void form_init_fp (
                    const double *x_coord,
                    const double *y_coord,
                    const state_tilde &init_state,
                    double* X);
            ///@}


            /**
             * @brief Solve QP problem.
             */
            void solve ();
       

            // -------------------------------


            /// @{
            /**
             * @brief Returns the next state.
             *  
             * @param[out] s an output state.
             */
            void get_next_state (state_orig &s);
            void get_next_state (state_tilde &s);
            /// @}
            
            /// @{
            /**
             * @brief Returns a state with given index.
             *  
             * @param[out] s an output state.
             * @param[in] ind index of a state [0 : N-1].
             */
            void get_state (state_orig &s, const int ind);
            void get_state (state_tilde &s, const int ind);
            /// @}


            // -------------------------------


            /**
             * @brief Returns the controls that must be applied to reach the 
             *  next state.
             *
             * @param[out] c an output control vector
             */
            void get_first_controls (control &c);


            /**
             * @brief The same as #get_first_controls, but takes an additional 
             * parameter - the index of the control inputs. Control[ind] is 
             * applied to State[ind-1] to reach State[ind].
             *
             * @param[out] c an output control vector
             * @param[in] ind index of control inputs [0 : N-1].
             */
            void get_controls (control &c, const int ind);


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
             * @brief Contains values of objective function after each iteration,
             * the initial value is also included.
             *
             * @note Updated by #solve function (only if the respective flag is
             * set on initialization).
             */
            std::vector<double> objective_log;
           

            // -------------------------------


            /**
             * @brief Internal representation.
             */
            qp_as *qp_sol;
    };


    /**
     * @brief API of the sparse MPC solver.
     */
    class solver_ip
    {
        public:


            // -------------------------------

            /** @brief Constructor: initialize an interior-point method solver.
             *
             * @param[in] N Number of sampling times in a preview window
             * @param[in] max_iter maximum number of internal loop iterations
             * @param[in] gain_position Position gain (Alpha)
             * @param[in] gain_velocity Velocity gain (Beta)
             * @param[in] gain_acceleration Acceleration gain (Gamma)
             * @param[in] gain_jerk Jerk gain (Eta)
             * @param[in] tol tolerance (internal loop)
             * @param[in] tol_out tolerance of the outer loop, which resolves
             *          the problem with new t (kappa) parameter.
             * @param[in] t logarithmic barrier parameter
             * @param[in] mu multiplier of t, >1.
             * @param[in] bs_alpha backtracking search parameter 0 < alpha < 0.5
             * @param[in] bs_beta  backtracking search parameter 0 < beta < 1
             * @param[in] obj_computation_on enable computation of the objective function 
             *          (the results are kept in #objective_log)
             * @param[in] backtracking_search_on enable backtracking search, note that even
             *          if it is disabled, the 'bs_beta' parameter is still used.
             */
            solver_ip (
                    const int N, 
                    const int unsigned max_iter,
                    const double gain_position = 2000.0, 
                    const double gain_velocity = 150.0, 
                    const double gain_acceleration = 0.01,
                    const double gain_jerk = 1.0,
                    const double tol = 1e-3,
                    const double tol_out = 1e-2,
                    const double t = 100,
                    const double mu = 15,
                    const double bs_alpha = 0.01,
                    const double bs_beta = 0.5,
                    const bool obj_computation_on = false,
                    const bool backtracking_search_on = true);


            ~solver_ip();


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


            ///@{
            /** @brief Generates an initial feasible point. 

                @param[in] x_coord x coordinates of points satisfying constraints
                @param[in] y_coord y coordinates of points satisfying constraints
                @param[in] init_state initial state (#state_orig or #state_tilde)
                @param[in,out] X solution of optimization problem
             */
            void form_init_fp (
                    const double *x_coord,
                    const double *y_coord,
                    const state_orig &init_state,
                    double* X);

            void form_init_fp (
                    const double *x_coord,
                    const double *y_coord,
                    const state_tilde &init_state,
                    double* X);
            ///@}


            /**
             * @brief Solve QP problem.
             */
            void solve ();


            // -------------------------------


            /// @{
            /**
             * @brief Returns the next state.
             *  
             * @param[out] s an output state.
             */
            void get_next_state (state_orig &s);
            void get_next_state (state_tilde &s);
            /// @}


            /// @{
            /**
             * @brief Returns a state with given index.
             *  
             * @param[out] s an output state.
             * @param[in] ind index of a state [0 : N-1].
             */
            void get_state (state_orig &s, const int ind);
            void get_state (state_tilde &s, const int ind);
            /// @}


            // -------------------------------

            
            /**
             * @brief Returns the controls that must be applied to reach the 
             *  next state.
             *
             * @param[out] c an output control vector
             */
            void get_first_controls (control &c);


            /**
             * @brief The same as #get_first_controls, but takes an additional 
             * parameter - the index of the control inputs. Control[ind] is 
             * applied to State[ind-1] to reach State[ind].
             *
             * @param[out] c an output control vector
             * @param[in] ind index of control inputs [0 : N-1].
             */
            void get_controls (control &c, const int ind);


            // -------------------------------


            /**
             * @brief The number of iterations of the external loop.
             *
             * @note Updated by #solve function.
             */
            unsigned int ext_loop_iterations;

            /**
             * @brief The total number of iterations of the internal loop.
             *
             * @note Updated by #solve function.
             */
            unsigned int int_loop_iterations;

            /**
             * @brief The total number of iterations of backtracking search.
             *
             * @note Updated by #solve function.
             */
            unsigned int bt_search_iterations;

            /**
             * @brief Contains values of objective function after each iteration,
             * the initial value is also included.
             *
             * @note Updated by #solve function (only if the respective flag is
             * set on initialization).
             */
            std::vector<double> objective_log;


            // -------------------------------


            /**
             * @brief Internal representation.
             */
            qp_ip *qp_sol;
    };
}
/// @}

#endif /*SMPC_SOLVER_H*/
