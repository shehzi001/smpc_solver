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
/**
 * @brief API of the sparse MPC solver.
 */
class smpc_solver
{
    public:


        // -------------------------------


        /** @brief Constructor: initialization of the constant parameters

            @param[in] N Number of sampling times in a preview window
            @param[in] Alpha Velocity gain
            @param[in] Beta Position gain
            @param[in] Gamma Jerk gain
            @param[in] regularization regularization
            @param[in] tol tolerance
        */
        smpc_solver(
                int N, 
                double Alpha = 150.0, 
                double Beta = 2000.0, 
                double Gamma = 1.0,
                double regularization = 0.01,
                double tol = 1e-7);

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
            @param[in] X_tilde current state (@ref pX_tilde "X_tilde")
            @param[in,out] X solution of optimization problem
        */
        void init(
                double* T,
                double* h,
                double* angle,
                double* zref_x,
                double* zref_y,
                double* lb,
                double* ub,
                double* X_tilde,
                double* X);


        /**
         * @brief Solve QP problem.
         *
         * @return number of activated constraints
         */
        int solve ();
   

        // -------------------------------

        /**
         * @brief Returns the next state as @ref pX_tilde "X_tilde".
         *  
         * @param[in,out] state the state (#NUM_STATE_VAR elements).
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
         * @brief Returns the next state represented by original variables.
         *  
         * @param[in,out] state the state (#NUM_STATE_VAR elements).
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

        // -------------------------------


    private:
        qp_as *qpas_solver;
};
/// @}

#endif /*SMPC_SOLVER_H*/
