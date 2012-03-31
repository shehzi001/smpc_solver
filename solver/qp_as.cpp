/** 
 * @file
 * @author Alexander Sherikov
 * @date 19.07.2011 22:30:13 MSD
 */


/****************************************
 * INCLUDES 
 ****************************************/
#include "qp_as.h"
#include "state_handling.h"

#include <cmath> //cos,sin



/****************************************
 * FUNCTIONS
 ****************************************/

/** @brief Constructor: initialization of the constant parameters

    @param[in] N_ Number of sampling times in a preview window
    @param[in] Alpha Velocity gain
    @param[in] Beta Position gain
    @param[in] Gamma Jerk gain
    @param[in] regularization regularization
    @param[in] tol_ tolerance
*/
qp_as::qp_as(
        const int N_, 
        const double Alpha, 
        const double Beta, 
        const double Gamma, 
        const double regularization, 
        const double tol_) : 
    qp_solver (N_, Alpha, Beta, Gamma, regularization, tol_),
    chol (N_)
{
    iHg = new double[2*N];

    // there are no inequality constraints in the initial working set (no hot-starting for the moment)
    W = new int[2*N];
    W_sign = new int[2*N];

    constraints.resize(2*N);
}


/** Destructor */
qp_as::~qp_as()
{
    if (W  != NULL)
        delete W;

    if (iHg != NULL)
        delete iHg;

    if (W_sign != NULL)
        delete W_sign;
}


/** @brief Initializes quadratic problem.

    @param[in] T_ Sampling time (for the moment it is assumed to be constant) [sec.]
    @param[in] h_ Height of the Center of Mass divided by gravity
    @param[in] h_initial_ current h
    @param[in] angle Rotation angle for each state in the preview window
    @param[in] zref_x reference values of z_x
    @param[in] zref_y reference values of z_y
    @param[in] lb array of lower constraints for z_x and z_y
    @param[in] ub array of upper constraints for z_x and z_y
*/
void qp_as::set_parameters(
        const double* T_, 
        const double* h_, 
        const double h_initial_,
        const double* angle,
        const double* zref_x,
        const double* zref_y,
        const double* lb,
        const double* ub)
{
    nW = 0;

    h_initial = h_initial_;
    set_state_parameters (T_, h_, h_initial_);
    form_iHg (zref_x, zref_y);
    form_constraints(lb, ub, angle);
}



/**
 * @brief Forms inv(2*H) * g.
 *
 * @param[in] zref_x reference values of z_x
 * @param[in] zref_y reference values of z_y
 */
void qp_as::form_iHg(const double *zref_x, const double *zref_y)
{
    for (int i = 0; i < N; i++)
    {
        // inv (2*H) * Cp' * zref
        iHg[i*2] =     -i2Q[0] * zref_x[i] * gain_beta;
        iHg[i*2 + 1] = -i2Q[0] * zref_y[i] * gain_beta; 
    }
}



/**
 * @brief Forms the upper and lower constraints.
 *
 * @param[in] lb array of lower constraints for z
 * @param[in] ub array of upper constraints for z
 * @param[in] angle Rotation angle for each state in the preview window
 */
void qp_as::form_constraints(const double *lb, const double *ub, const double *angle)
{
    for (int i=0; i < N; i++)
    {
        double cosR = cos(angle[i]);
        double sinR = sin(angle[i]);

        constraints[i*2].set(
                SMPC_NUM_STATE_VAR*i, 
                cosR, sinR, 
                lb[i*2], ub[i*2], 
                false);
        constraints[i*2+1].set(
                SMPC_NUM_STATE_VAR*i, 
                -sinR, cosR, 
                lb[i*2+1], ub[i*2+1], 
                false);
    }
}



/**
 * @brief Checks for blocking constraints.
 *
 * @return sequential number of constraint to be added, -1 if no constraints.
 */
int qp_as::check_blocking_constraints()
{
    alpha = 1;

    /* Index to include in the working set, -1 if no constraint have to be included. */
    int activated_var_num = -1;


    for (int i = 0; i < 2*N; i++)
    {
        // Check only inactive constraints for violation. 
        // The constraints in the working set will not be violated regardless of 
        // the depth of descent
        if (!constraints[i].isActive)
        {
            int ind = constraints[i].ind;

            double coef_x = constraints[i].coef_x;
            double coef_y = constraints[i].coef_y;
            double constr = X[ind]*coef_x + X[ind+3]*coef_y;
            double d_constr = dX[ind]*coef_x + dX[ind+3]*coef_y;

            if ( d_constr < -tol )
            {
                double t = (constraints[i].lb - constr)/d_constr;
                if (t < alpha)
                {
                    alpha = t;
                    activated_var_num = i;
                    W_sign[nW] = -1;
                }
            }
            else if ( d_constr > tol )
            {
                double t = (constraints[i].ub - constr)/d_constr;
                if (t < alpha)
                {
                    alpha = t;
                    activated_var_num = i;
                    W_sign[nW] = 1;
                }
            }
        }
    }
    if (activated_var_num != -1)
    {
        W[nW] = activated_var_num;    
        nW++;
        constraints[activated_var_num].isActive = true;
    }

    return (activated_var_num);
}  


/**
 * @brief Selects a constraint for removal from active set.
 *
 * @param[in] lambda vector of Lagrange multipliers corresponding
 *  to inequality constraints.
 *
 * @return index of constraint in the active set, -1 if no constraint
 *  can be removed.
 *
 * @attention If a constraint for removal is selected, then it is 
 * removed from the active set (#W) and the number of constraints
 * in active set (#nW) is decremented.
 */
int qp_as::choose_excl_constr (const double *lambda)
{
    double min_lambda = -tol;
    int ind_exclude = -1;

    // find the constraint with the smallest lambda
    for (int i = 0; i < nW; i++)
    {
        if (lambda[i] * W_sign[i] < min_lambda)
        {
            min_lambda = lambda[i] * W_sign[i];
            ind_exclude = i;
        }
    }

    if (ind_exclude != -1)
    {
        constraints[W[ind_exclude]].isActive = false;
        for (int i = ind_exclude; i < nW-1; i++)
        {
            W[i] = W[i + 1];
            W_sign[i] = W_sign[i + 1];
        }
        nW--;
    }

    return (ind_exclude);
}


/**
 * @brief Solve QP problem.
 *
 * @return number of activated constraints
 */
int qp_as::solve ()
{
    // obtain dX
    chol.solve(*this, iHg, X, dX);

    for (;;)
    {
        int activated_var_num = check_blocking_constraints();

        // Move in the feasible descent direction
        for (int i = 0; i < N*SMPC_NUM_VAR ; i++)
        {
            X[i] += alpha * dX[i];
        }

        // no new inequality constraints
        if (activated_var_num == -1)
        {
            int ind_exclude = choose_excl_constr (chol.get_lambda(*this));
            if (ind_exclude != -1)
            {
                chol.down_resolve (*this, iHg, constraints, nW, W, ind_exclude, X, dX);
            }
            else
            {
                break;
            }
        }
        else
        {
            // add row to the L matrix and find new dX
            chol.up_resolve (*this, iHg, constraints, nW, W, X, dX);
        }
    }

    return (nW);
}
