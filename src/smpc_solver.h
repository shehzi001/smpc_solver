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

            @param[in] T Sampling time (for the moment it is assumed to be constant) [sec.]
            @param[in] h Height of the Center of Mass divided by gravity
            @param[in] angle Rotation angle for each state in the preview window
            @param[in] zref_x reference values of z_x
            @param[in] zref_y reference values of z_y
            @param[in] lb array of lower bounds for z
            @param[in] ub array of upper bounds for z
            @param[in,out] X initial guess / solution of optimization problem
        */
        void init(
                double* T,
                double* h,
                double* angle,
                double* zref_x,
                double* zref_y,
                double* lb,
                double* ub,
                double* X);


        /**
         * @brief Solve QP problem.
         *
         * @return number of activated constraints
         */
        int solve ();
   

        // -------------------------------


        /**
         * @brief Determines coordinates of ZMP and CoM.
         *
         * @param[out] ZMP_x x coordinate of ZMP.
         * @param[out] ZMP_y y coordinate of ZMP.
         * @param[out] CoM_x x coordinate of CoM.
         * @param[out] CoM_y y coordinate of CoM.
         */
        void get_ZMP_CoM (double *ZMP_x, double *ZMP_y, double *CoM_x, double *CoM_y);


        // -------------------------------


    private:
        qp_as *qpas_solver;
};

#endif /*SMPC_SOLVER_H*/
