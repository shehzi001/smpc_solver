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
    problem_parameters (N_, Alpha, Beta, Gamma, regularization),
    chol (N_)
{
    tol = tol_;

    dX = new double[SMPC_NUM_VAR*N]();

    // there are no inequality constraints in the initial working set (no hot-starting for the moment)

    constraints.resize(2*N);
}


/** Destructor */
qp_as::~qp_as()
{
    if (dX != NULL)
        delete dX;
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
        const double* zref_x_,
        const double* zref_y_,
        const double* lb,
        const double* ub)
{
    h_initial = h_initial_;
    set_state_parameters (T_, h_, h_initial_);

    zref_x = (double *) zref_x_;
    zref_y = (double *) zref_y_;

    active_set.clear();

    // form inv(2*H) *g and initialize constraints
    // inv(2*H) * g  =  inv (2*(beta/2)) * beta * Cp' * zref  =  zref
    for (int i = 0, cind = 0; i < N; ++i)
    {
        double cosR = cos(angle[i]);
        double sinR = sin(angle[i]);
        double RTzref_x = (cosR*zref_x[i] + sinR*zref_y[i]);
        double RTzref_y = (-sinR*zref_x[i] + cosR*zref_y[i]);

        constraints[cind].set(
                cind, cosR, sinR, 
                lb[cind] - RTzref_x, 
                ub[cind] - RTzref_x, 
                false);
        ++cind;

        constraints[cind].set(
                cind, -sinR, cosR, 
                lb[cind] - RTzref_y, 
                ub[cind] - RTzref_y, 
                false);
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


    for (unsigned int i = 0; i < constraints.size(); ++i)
    {
        // Check only inactive constraints for violation. 
        // The constraints in the working set will not be violated regardless of 
        // the depth of descent
        if (!constraints[i].isActive)
        {
            const constraint c = constraints[i];

            const double constr = X[c.ind]*c.coef_x + X[c.ind+3]*c.coef_y;
            const double d_constr = dX[c.ind]*c.coef_x + dX[c.ind+3]*c.coef_y;

            if ( d_constr < -tol )
            {
                const double t = (c.lb - constr)/d_constr;
                if (t < alpha)
                {
                    alpha = t;
                    activated_var_num = i;
                    sign = -1;
                }
            }
            else if ( d_constr > tol )
            {
                const double t = (c.ub - constr)/d_constr;
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
 * removed from the active set (#active_set). 
 */
int qp_as::choose_excl_constr (const double *lambda)
{
    double min_lambda = -tol;
    int ind_exclude = -1;

    // find the constraint with the smallest lambda
    for (unsigned int i = 0; i < active_set.size(); ++i)
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
    for (int i = 0; i < N; ++i)
    {
        const int ind = i*SMPC_NUM_STATE_VAR;
        X[ind]   -= zref_x[i];
        X[ind+3] -= zref_y[i];
    }

    // obtain dX
    chol.solve(*this, X, dX);

    for (;;)
    {
        int activated_var_num = check_blocking_constraints();

        // Move in the feasible descent direction
        for (int i = 0; i < N*SMPC_NUM_VAR; i += SMPC_NUM_VAR)
        {
            X[i]   += alpha * dX[i];
            X[i+1] += alpha * dX[i+1];
            X[i+2] += alpha * dX[i+2];
            X[i+3] += alpha * dX[i+3];
            X[i+4] += alpha * dX[i+4];
            X[i+5] += alpha * dX[i+5];
            X[i+6] += alpha * dX[i+6];
            X[i+7] += alpha * dX[i+7];
        }

        // no new inequality constraints
        if (activated_var_num == -1)
        {
            int ind_exclude = choose_excl_constr (chol.get_lambda(*this));
            if (ind_exclude != -1)
            {
                chol.down_resolve (*this, active_set, ind_exclude, X, dX);
            }
            else
            {
                break;
            }
        }
        else
        {
            // add row to the L matrix and find new dX
            chol.up_resolve (*this, active_set, X, dX);
        }
    }

    for (int i = 0; i < N; ++i)
    {
        const int ind = i*SMPC_NUM_STATE_VAR;
        X[ind]   += zref_x[i];
        X[ind+3] += zref_y[i];
    }

    return (active_set.size());
}


/**
 * @brief Generates an initial feasible point. 
 * First we perform a change of variable to @ref pX_tilde "X_tilde"
 * and then generate a feasible point.
 *
 * @param[in] x_coord x coordinates of points satisfying constraints
 * @param[in] y_coord y coordinates of points satisfying constraints
 * @param[in] init_state current state
 * @param[in,out] X_ initial guess / solution of optimization problem
 */
void qp_as::form_init_fp (
        const double *x_coord, 
        const double *y_coord, 
        const double *init_state,
        double* X_)
{
    X = X_;

    double X_tilde[6] = {
        init_state[0], init_state[1], init_state[2],
        init_state[3], init_state[4], init_state[5]};
    double *control = &X[SMPC_NUM_STATE_VAR*N];
    double *cur_state = X;
    state_handling::orig_to_tilde (h_initial, X_tilde);
    const double *prev_state = X_tilde;

    
    for (int i=0; i<N; i++)
    {
        //------------------------------------
        /* inv(Cp*B). This is a [2 x 2] diagonal matrix (which is invertible if T^3/6-h*T is
         * not equal to zero). The two elements on the main diagonal are equal, and only one of them 
         * is stored, which is equal to
            1/(T^3/6 - h*T)
         */
        double iCpB = 1/(spar[i].B[0]);

        /* inv(Cp*B)*Cp*A. This is a [2 x 6] matrix with the following structure
            iCpB_CpA = [a b c 0 0 0;
                        0 0 0 a b c];

            a = iCpB
            b = iCpB*T
            c = iCpB*T^2/2
         * Only a,b and c are stored.
         */
        double iCpB_CpA[3] = {iCpB, iCpB*spar[i].A3, iCpB*spar[i].A6};
        //------------------------------------


        control[0] = -iCpB_CpA[0]*prev_state[0] - iCpB_CpA[1]*prev_state[1] - iCpB_CpA[2]*prev_state[2] + iCpB*x_coord[i];
        control[1] = -iCpB_CpA[0]*prev_state[3] - iCpB_CpA[1]*prev_state[4] - iCpB_CpA[2]*prev_state[5] + iCpB*y_coord[i];

        cur_state[0] = prev_state[0] + spar[i].A3*prev_state[1] + spar[i].A6*prev_state[2] + spar[i].B[0]*control[0];
        cur_state[1] =                            prev_state[1] + spar[i].A3*prev_state[2] + spar[i].B[1]*control[0];
        cur_state[2] =                                                       prev_state[2] + spar[i].B[2]*control[0];
        cur_state[3] = prev_state[3] + spar[i].A3*prev_state[4] + spar[i].A6*prev_state[5] + spar[i].B[0]*control[1];
        cur_state[4] =                            prev_state[4] + spar[i].A3*prev_state[5] + spar[i].B[1]*control[1];
        cur_state[5] =                                                       prev_state[5] + spar[i].B[2]*control[1];


        prev_state = &X[SMPC_NUM_STATE_VAR*i];
        cur_state = &X[SMPC_NUM_STATE_VAR*(i+1)];
        control = &control[SMPC_NUM_CONTROL_VAR];
    }
}
