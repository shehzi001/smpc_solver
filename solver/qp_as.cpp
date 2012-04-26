/** 
 * @file
 * @author Alexander Sherikov
 * @date 19.07.2011 22:30:13 MSD
 */


/****************************************
 * INCLUDES 
 ****************************************/
#include "qp.h"
#include "qp_as.h"
#include "state_handling.h"

#include <cmath> //cos,sin


using namespace AS;

/****************************************
 * FUNCTIONS
 ****************************************/

/** @brief Constructor: initialization of the constant parameters

    @param[in] N_ Number of sampling times in a preview window
    @param[in] gain_position Position gain
    @param[in] gain_velocity Velocity gain
    @param[in] gain_acceleration Acceleration gain
    @param[in] gain_jerk Jerk gain
    @param[in] tol_ tolerance
    @param[in] obj_computation_enabled_ enable computation of the objective function
*/
qp_as::qp_as(
        const int N_, 
        const double gain_position, 
        const double gain_velocity, 
        const double gain_acceleration, 
        const double gain_jerk, 
        const double tol_,
        const bool obj_computation_enabled_) : 
    problem_parameters (N_, gain_position, gain_velocity, gain_acceleration, gain_jerk),
    chol (N_),
    tol (tol_),
    obj_computation_enabled (obj_computation_enabled_)
{
    dX = new double[SMPC_NUM_VAR*N]();

    constraints.resize(2*N);

    max_added_constraints_num = N*2;
    constraint_removal_enabled = true;
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
    set_state_parameters (T_, h_, h_initial_);

    zref_x = zref_x_;
    zref_y = zref_y_;

    active_set.clear();

    added_constraints_num = 0;
    removed_constraints_num = 0;


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
 * @brief Generates an initial feasible point. 
 *
 * @param[in] x_coord x coordinates of points satisfying constraints
 * @param[in] y_coord y coordinates of points satisfying constraints
 * @param[in] init_state current state
 * @param[in] tilde_state if true the state is interpreted as @ref pX_tilde "X_tilde".
 * @param[in,out] X_ initial guess / solution of optimization problem
 */
void qp_as::form_init_fp (
        const double *x_coord, 
        const double *y_coord, 
        const double *init_state,
        const bool tilde_state,
        double* X_)
{
    X = X_;
    form_init_fp_tilde<problem_parameters>(*this, x_coord, y_coord, init_state, tilde_state, X);
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
 * @param[in,out] obj_log a vector of objective function values
 *
 * @return number of activated constraints
 */
void qp_as::solve (vector<double> &obj_log)
{
    for (int i = 0; i < N; ++i)
    {
        const int ind = i*SMPC_NUM_STATE_VAR;
        X[ind]   -= zref_x[i];
        X[ind+3] -= zref_y[i];
    }
    if (obj_computation_enabled)
    {
        obj_log.clear();
        obj_log.push_back(compute_obj());
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

        if (obj_computation_enabled)
        {
            obj_log.push_back(compute_obj());
        }

        if (added_constraints_num == max_added_constraints_num)
        {
            break;
        }

        if (activated_var_num != -1)
        {
            // add row to the L matrix and find new dX
            chol.up_resolve (*this, active_set, X, dX);
            ++added_constraints_num;
        }
        else if (constraint_removal_enabled)
        {
            // no new inequality constraints
            int ind_exclude = choose_excl_constr (chol.get_lambda(*this));
            if (ind_exclude == -1)
            {
                break;
            }

            chol.down_resolve (*this, active_set, ind_exclude, X, dX);
            ++removed_constraints_num;
        }
        else
        {
            break;
        }
    }

    for (int i = 0; i < N; ++i)
    {
        const int ind = i*SMPC_NUM_STATE_VAR;
        X[ind]   += zref_x[i];
        X[ind+3] += zref_y[i];
    }

    active_set_size = active_set.size();
}



/**
 * @brief Compute value of the objective function.
 *
 * @return value of the objective function.
 */
double qp_as::compute_obj()
{
    int i,j;
    double obj_pos = 0;
    double obj_vel = 0;
    double obj_acc = 0;
    double obj_jerk = 0;

    // phi_X = X'*H*X + g'*X
    for(i = 0, j = 0;
        i < N*SMPC_NUM_STATE_VAR;
        i += SMPC_NUM_STATE_VAR, j += 2)
    {
        const double X_copy[6] = {X[i], X[i+1], X[i+2], X[i+3], X[i+4], X[i+5]};

        // X'*H*X
        obj_pos += X_copy[0]*X_copy[0] + X_copy[3]*X_copy[3];
        obj_vel += X_copy[1]*X_copy[1] + X_copy[4]*X_copy[4];
        obj_acc += X_copy[2]*X_copy[2] + X_copy[5]*X_copy[5];
    }
    for (; i < N*SMPC_NUM_VAR; i += SMPC_NUM_CONTROL_VAR)
    {
        // X'*H*X
        obj_jerk += X[i] * X[i] + X[i+1] * X[i+1];
    }

    return (0.5*obj_pos/i2Q[0] + 0.5*obj_vel/i2Q[1] + 0.5*obj_acc/i2Q[2] + 0.5*obj_jerk/i2P);
}

