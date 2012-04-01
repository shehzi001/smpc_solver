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

    constraints.resize(2*N);
}


/** Destructor */
qp_as::~qp_as()
{
    if (iHg != NULL)
        delete iHg;
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
    h_initial = h_initial_;
    set_state_parameters (T_, h_, h_initial_);


    active_set.clear();
    int cind = 0;

    // form inv(2*H) *g and initialize constraints
    // inv(2*H) * g  =  inv (2*(beta/2)) * beta * Cp' * zref  =  zref
    for (int i = 0; i < N; ++i)
    {
        double cosR = cos(angle[i]);
        double sinR = sin(angle[i]);

        constraints[cind].set(cind, cosR, sinR, lb[cind], ub[cind], false);
        iHg[cind] = -zref_x[i];
        ++cind;

        constraints[cind].set(cind, -sinR, cosR, lb[cind], ub[cind], false);
        iHg[cind] = -zref_y[i]; 
        ++cind;
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
    int sign = 0;


    for (int i = 0; i < 2*N; i++)
    {
        // Check only inactive constraints for violation. 
        // The constraints in the working set will not be violated regardless of 
        // the depth of descent
        if (!constraints[i].isActive)
        {
            constraint c = constraints[i];

            double constr = X[c.ind]*c.coef_x + X[c.ind+3]*c.coef_y;
            double d_constr = dX[c.ind]*c.coef_x + dX[c.ind+3]*c.coef_y;

            if ( d_constr < -tol )
            {
                double t = (c.lb - constr)/d_constr;
                if (t < alpha)
                {
                    alpha = t;
                    activated_var_num = i;
                    sign = -1;
                }
            }
            else if ( d_constr > tol )
            {
                double t = (c.ub - constr)/d_constr;
                if (t < alpha)
                {
                    alpha = t;
                    activated_var_num = i;
                    sign = 1;
                }
            }
        }
    }
    if (activated_var_num != -1)
    {
        constraints[activated_var_num].sign = sign;
        constraints[activated_var_num].isActive = true;
        active_set.push_back(constraints[activated_var_num]);
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
    for (unsigned int i = 0; i < active_set.size(); i++)
    {
        if (lambda[i] * active_set[i].sign < min_lambda)
        {
            min_lambda = lambda[i] * active_set[i].sign;
            ind_exclude = i;
        }
    }

    if (ind_exclude != -1)
    {
        active_set.erase(active_set.begin()+ind_exclude);
        constraints[ind_exclude].isActive = false;
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
                chol.down_resolve (*this, iHg, active_set, ind_exclude, X, dX);
            }
            else
            {
                break;
            }
        }
        else
        {
            // add row to the L matrix and find new dX
            chol.up_resolve (*this, iHg, active_set, X, dX);
        }
    }

    return (active_set.size());
}
