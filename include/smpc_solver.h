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



/**
 * @brief API of the sparse MPC solver.
 */
class smpc_solver
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
        smpc_solver (
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
        smpc_solver (
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


        ~smpc_solver();


        // -------------------------------


        /** @brief Initializes quadratic problem.

            @param[in] T sampling time for each time step [sec.]
            @param[in] h height of the center of mass divided by gravity for each time step
            @param[in] angle rotation angle for each state relative to the world frame
            @param[in] zref_x reference values of x coordinate of ZMP
            @param[in] zref_y reference values of y coordinate of ZMP
            @param[in] lb array of lower bounds for coordinates of ZMP
            @param[in] ub array of upper bounds for coordinates of ZMP
        */
        void set_parameters (
                const double* T,
                const double* h,
                const double* angle,
                const double* zref_x,
                const double* zref_y,
                const double* lb,
                const double* ub);


        /** @brief Generates an initial feasible point. 

            First we perform a change of variable to @ref pX_tilde "X_tilde"
            generate a feasible point, and then we go back to @ref pX_bar "X_bar".
         
            @param[in] x_coord x coordinates of points satisfying constraints
            @param[in] y_coord y coordinates of points satisfying constraints
            @param[in] X_tilde current state (@ref pX_tilde "X_tilde")
            @param[in,out] X solution of optimization problem
         */
        void form_init_fp (
                const double *x_coord,
                const double *y_coord,
                const double *X_tilde,
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
         * @brief Returns the next state as @ref pX_tilde "X_tilde".
         *  
         * @param[in,out] state the state (#SMPC_NUM_STATE_VAR elements).
         *
         *  \verbatim
            state[0] - x ZMP position [meter]
            state[1] - x CoM velocity [meter/s]
            state[2] - x CoM acceleration [meter/s^2]
            state[3] - y ZMP position [meter]
            state[4] - y CoM velocity [meter/s]
            state[5] - y CoM acceleration [meter/s^2]
            \endverbatim
         */
        void get_next_state_tilde (double *state);


        /**
         * @brief Converts state from original variables to @ref pX_tilde "X_tilde".
         *
         * @param[in] h height of CoM divided by gravity.
         * @param[in,out] state see below:
         *
         * Input:
         *  \verbatim
            state[0] - x CoM position [meter]
            state[1] - x CoM velocity [meter/s]
            state[2] - x CoM acceleration [meter/s^2]
            state[3] - y CoM position [meter]
            state[4] - y CoM velocity [meter/s]
            state[5] - y CoM acceleration [meter/s^2]
            \endverbatim
         *   
         * Output:
         *  \verbatim
            state[0] - x ZMP position [meter]
            state[1] - x CoM velocity [meter/s]
            state[2] - x CoM acceleration [meter/s^2]
            state[3] - y ZMP position [meter]
            state[4] - y CoM velocity [meter/s]
            state[5] - y CoM acceleration [meter/s^2]
            \endverbatim
         */
        void convert_to_tilde (double h, double *state);


        /**
         * @brief Returns the next state represented by original variables.
         *  
         * @param[in,out] state the state (#SMPC_NUM_STATE_VAR elements).
         *
         *  \verbatim
            state[0] - x CoM position [meter]
            state[1] - x CoM velocity [meter/s]
            state[2] - x CoM acceleration [meter/s^2]
            state[3] - y CoM position [meter]
            state[4] - y CoM velocity [meter/s]
            state[5] - y CoM acceleration [meter/s^2]
            \endverbatim
         */
        void get_next_state (double *state);


        /**
         * @brief The same as #get_next_state, but takes additional parameter - 
         * index of the desired state in the state vector.
         *  
         * @param[in] ind index of a state [0, N-1].
         * @param[in,out] state the state (#SMPC_NUM_STATE_VAR elements).
         */
        void get_state (const int ind, double *state);

        /**
         * @brief Returns the controls that must be applied to reach the 
         *  next state.
         *
         * @param[in,out] controls the controls (#SMPC_NUM_CONTROL_VAR).
         * @verbatim
            control[0] - jerk along x axis
            control[1] - jerk along y axis
           @endverbatim
         */
        void get_first_controls (double *controls);

        // -------------------------------


    private:
        qp_solver *qp_sol;
};
/// @}

#endif /*SMPC_SOLVER_H*/
