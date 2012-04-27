/** 
 * @file
 * @author Alexander Sherikov
 * @date 19.07.2011 22:30:13 MSD
 */


/****************************************
 * INCLUDES 
 ****************************************/
#include "qp_ip.h"
#include "state_handling.h"
#include "qp.h"


#include <cmath> // log

/****************************************
 * FUNCTIONS
 ****************************************/

using namespace IP;

//==============================================
// qp_ip

/** @brief Constructor: initialization of the constant parameters

    @param[in] N_ Number of sampling times in a preview window
    @param[in] gain_position (Alpha) Position gain
    @param[in] gain_velocity (Beta) Velocity gain
    @param[in] gain_acceleration (Gamma) Acceleration gain
    @param[in] gain_jerk (Eta) Jerk gain
    @param[in] tol_ tolerance
    @param[in] obj_computation_enabled_ enable computation of the objective function
*/
qp_ip::qp_ip(
        const int N_, 
        const double gain_position_,
        const double gain_velocity_,
        const double gain_acceleration_,
        const double gain_jerk_,
        const double tol_,
        const bool obj_computation_enabled_) : 
    problem_parameters (N_, gain_position_, gain_velocity_, gain_acceleration_, gain_jerk_),
    gain_position (gain_position_),
    tol (tol_),
    obj_computation_enabled (obj_computation_enabled_),
    chol (N_)
{
    dX = new double[SMPC_NUM_VAR*N]();
    g = new double[2*N];
    i2hess = new double[2*N];
    i2hess_grad = new double[N*SMPC_NUM_VAR];
    grad = new double[2*N];

    Q[0] = gain_position_/2;
    Q[1] = gain_velocity_/2;
    Q[2] = gain_acceleration_/2;
    P = gain_jerk_/2;
}


/** Destructor */
qp_ip::~qp_ip()
{
    if (g != NULL)
        delete g;
    if (i2hess != NULL)
        delete i2hess;
    if (i2hess_grad != NULL)
        delete i2hess_grad;
    if (grad != NULL)
        delete grad;
    if (dX  != NULL)
        delete dX;
}


/** @brief Initializes quadratic problem.

    @param[in] T Sampling time (for the moment it is assumed to be constant) [sec.]
    @param[in] h Height of the Center of Mass divided by gravity
    @param[in] h_initial_ current h
    @param[in] angle Rotation angle for each state in the preview window
    @param[in] zref_x reference values of z_x
    @param[in] zref_y reference values of z_y
    @param[in] lb_ array of lower bounds for z_x and z_y
    @param[in] ub_ array of upper bounds for z_x and z_y
*/
void qp_ip::set_parameters(
        const double* T, 
        const double* h, 
        const double h_initial_, 
        const double* angle,
        const double* zref_x_,
        const double* zref_y_,
        const double* lb_,
        const double* ub_)
{
    set_state_parameters (T, h, h_initial_, angle);

    lb = lb_;
    ub = ub_;

    zref_x = zref_x_;
    zref_y = zref_y_;

    form_g (zref_x, zref_y);
}



/**
 * @brief Forms vector @ref pg "g".
 *
 * @param[in] zref_x_ x coordinates of reference ZMP positions
 * @param[in] zref_y_ y coordinates of reference ZMP positions
 */
void qp_ip::form_g (const double *zref_x_, const double *zref_y_)
{
    double p0, p1;
    double cosA, sinA;

    for (int i = 0; i < N; i++)
    {
        cosA = spar[i].cos;
        sinA = spar[i].sin;

        // zref
        p0 = zref_x_[i];
        p1 = zref_y_[i];

        // inv (2*H) * R' * Cp' * zref
        g[i*2] = -(cosA*p0 + sinA*p1)*gain_position;
        g[i*2 + 1] = -(-sinA*p0 + cosA*p1)*gain_position; 
    }
}



/**
 * @brief Compute gradient of phi (partially), varying elements of 
 * i2hess, logarithmic barrier part of phi, i2hess_grad = -i2hess*grad.
 *
 * @param[in] kappa 1/t, a logarithmic barrier multiplicator.
 *
 * @return logarithmic barrier part of phi.
 */
double qp_ip::form_grad_i2hess_logbar (const double kappa)
{
    double phi_X_logbar = 1.0;

    // grad = H*X + g + kappa * b;
    // initialize inverted hessian
    // initialize logarithmic barrier in the function
    for (int i = 0; i < 2*N; i++)
    {
        const int j = 3*i;
        double lb_diff = -lb[i] + X[j];
        double ub_diff =  ub[i] - X[j];

        // logarithmic barrier
        phi_X_logbar *= lb_diff * ub_diff;

        lb_diff = 1/lb_diff;
        ub_diff = 1/ub_diff;

        // grad = H*X + g + kappa * (ub_diff - lb_diff)
        const double grad_el = X[j]*gain_position + g[i] + kappa * (ub_diff - lb_diff);
        grad[i] = grad_el;

        // only elements 1:3:N*SMPC_NUM_STATE_VAR on the diagonal of hessian 
        // can change
        // hess = 2H + kappa * (ub_diff^2 + lb_diff^2)
        const double i2hess_el = 1/(gain_position + kappa * (ub_diff*ub_diff + lb_diff*lb_diff));
        i2hess[i] = i2hess_el;

        i2hess_grad[j] = -grad_el * i2hess_el;
        i2hess_grad[j+1] = - X[j+1];  //grad[j+1] * i2Q[1]; 
        i2hess_grad[j+2] = - X[j+2];  //grad[j+2] * i2Q[2]; 
    }

    for (int i = N*SMPC_NUM_STATE_VAR; i < N*SMPC_NUM_VAR; i+= SMPC_NUM_CONTROL_VAR)
    {
        i2hess_grad[i]   = - X[i];    //grad[i]   * i2P;
        i2hess_grad[i+1] = - X[i+1];  //grad[i+1] * i2P;
    }

    return (-kappa * log(phi_X_logbar));
}



/**
 * @brief Compute phi_X for initial point, phi_X must already store
 *      logarithmic barrier term.
 */
double qp_ip::form_phi_X ()
{
    int i,j;
    double phi_X_pos = 0;
    double phi_X_vel = 0;
    double phi_X_acc = 0;
    double phi_X_jerk = 0;
    double phi_X_gX = 0;

    // phi_X = X'*H*X + g'*X
    for(i = 0, j = 0; 
        i < N*SMPC_NUM_STATE_VAR; 
        i += SMPC_NUM_STATE_VAR, j += 2)
    {
        const double X_copy[6] = {X[i], X[i+1], X[i+2], X[i+3], X[i+4], X[i+5]};

        // X'*H*X
        phi_X_pos += X_copy[0]*X_copy[0] + X_copy[3]*X_copy[3];
        phi_X_vel += X_copy[1]*X_copy[1] + X_copy[4]*X_copy[4];
        phi_X_acc += X_copy[2]*X_copy[2] + X_copy[5]*X_copy[5];

        // g'*X
        phi_X_gX  += g[j]*X_copy[0] + g[j+1]*X_copy[3];
    }
    for (; i < N*SMPC_NUM_VAR; i += SMPC_NUM_CONTROL_VAR)
    {
        // X'*H*X
        phi_X_jerk += X[i] * X[i] + X[i+1] * X[i+1];
    }

    return (Q[0]*phi_X_pos + Q[1]*phi_X_vel + Q[2]*phi_X_acc + P*phi_X_jerk + phi_X_gX);
}


/**
 * @brief Find initial value of alpha.
 *
 * @note sets alpha to 0, if it is too small.
 */
double qp_ip::init_alpha()
{
    double min_alpha = 1.0;
    double alpha = 1.0;

    for (int i = 0; i < 2*N; i++)
    {
        // lower bound may be violated
        if (dX[i*3] < 0)
        {
            const double tmp_alpha = (lb[i]-X[i*3])/dX[i*3];
            if (tmp_alpha < min_alpha)
            {
                min_alpha = tmp_alpha;
            }
        }
        // upper bound may be violated
        else if (dX[i*3] > 0)
        {
            const double tmp_alpha = (ub[i]-X[i*3])/dX[i*3];
            if (tmp_alpha < min_alpha)
            {
                min_alpha = tmp_alpha;
            }
        }
    }

    if (min_alpha > tol)
    {
        while (alpha > min_alpha)
        {
            alpha *= bs_beta;
        }
    }
    else
    {
        alpha = 0;
    }

    return (alpha);
}



/**
 * @brief Forms bs_alpha * grad' * dX.
 *
 * @return result of multiplication.
 */
double qp_ip::form_bs_alpha_grad_dX ()
{
    double res_pos = 0;
    double res_vel = 0;
    double res_acc = 0;
    double res_jerk = 0;
    
    for (int i = 0, j = 0; i < N*SMPC_NUM_STATE_VAR; i += SMPC_NUM_STATE_VAR, j += 2)
    {
        res_pos += grad[j]   * dX[i]
                 + grad[j+1] * dX[i+3];

        res_vel += X[i+1] * dX[i+1]  //grad[i+1] * dX[i+1]
                 + X[i+4] * dX[i+4]; //grad[i+4] * dX[i+4]

        res_acc += X[i+2] * dX[i+2]  //grad[i+2] * dX[i+2]
                 + X[i+5] * dX[i+5]; //grad[i+5] * dX[i+5];
    }
    for (int i = N*SMPC_NUM_STATE_VAR; i < N*SMPC_NUM_VAR; i += SMPC_NUM_CONTROL_VAR)
    {
        res_jerk += X[i] * dX[i] + X[i+1] * dX[i+1];
            // grad[i+6] * dX[i+6]
            // grad[i+7] * dX[i+7];
    }

    return ((res_pos + res_vel/i2Q[1] + res_acc/i2Q[2] + res_jerk/i2P)*bs_alpha);
}


/**
 * @brief Forms phi(X+alpha*dX)
 *
 * @param[in] kappa logarithmic barrier multiplicator.
 *
 * @return a value of phi.
 */
double qp_ip::form_phi_X_tmp (const double kappa, const double alpha)
{
    int i,j;


    // phi_X += X'*H*X
    double res_gX = 0;
    double res_pos = 0;
    double res_vel = 0;
    double res_acc = 0;
    double res_jerk = 0;
    double res_logbar = 1.0;
    for (i = 0,j = 0; i < 2*N; i+=2, j += SMPC_NUM_STATE_VAR)
    {
        double X_tmp[6] = {
            X[j]   + alpha * dX[j],
            X[j+1] + alpha * dX[j+1],
            X[j+2] + alpha * dX[j+2],
            X[j+3] + alpha * dX[j+3],
            X[j+4] + alpha * dX[j+4],
            X[j+5] + alpha * dX[j+5],
        };

        // logarithmic barrier
        res_logbar *= (-lb[i]   + X_tmp[0]) * (ub[i]   - X_tmp[0])
                    * (-lb[i+1] + X_tmp[3]) * (ub[i+1] - X_tmp[3]);

        // phi_X += g'*X
        res_gX  += g[i] * X_tmp[0] + g[i+1] * X_tmp[3];

        // phi_X += X'*H*X // states
        res_pos += X_tmp[0]*X_tmp[0] + X_tmp[3]*X_tmp[3];
        res_vel += X_tmp[1]*X_tmp[1] + X_tmp[4]*X_tmp[4];
        res_acc += X_tmp[2]*X_tmp[2] + X_tmp[5]*X_tmp[5];
    }
    // phi_X += X'*H*X // controls
    for (i = N*SMPC_NUM_STATE_VAR; i < N*SMPC_NUM_VAR; i += SMPC_NUM_CONTROL_VAR)
    {
        double X_tmp[2] = {
            X[i]   + alpha * dX[i],
            X[i+1] + alpha * dX[i+1]
        };

        res_jerk += X_tmp[0] * X_tmp[0] + X_tmp[1] * X_tmp[1];
    }

    return (res_gX + Q[0]*res_pos + Q[1]*res_vel + Q[2]*res_acc + P*res_jerk - kappa * log(res_logbar));
}



/**
 * @brief Set parameters of interior-point method.
 *
 * @param[in] t_ logarithmic barrier parameter
 * @param[in] mu_ multiplier of t, >1.
 * @param[in] bs_alpha_ backtracking search parameter alpha
 * @param[in] bs_beta_  backtracking search parameter beta
 * @param[in] max_iter_ maximum number of internal loop iterations
 * @param[in] tol_out_ tolerance of the outer loop
 */
void qp_ip::set_ip_parameters (
        const double t_, 
        const double mu_, 
        const double bs_alpha_, 
        const double bs_beta_, 
        const unsigned int max_iter_,
        const double tol_out_)
{
    t = t_;
    mu = mu_;
    bs_alpha = bs_alpha_;
    bs_beta = bs_beta_;
    max_iter = max_iter_;
    tol_out = tol_out_;
}



/**
 * @brief Solve QP using interior-point method.
 *
 * @param[in,out] obj_log a vector of objective function values
 *
 * @return 0 if ok, negative number otherwise.
 */
void qp_ip::solve(vector<double> &obj_log)
{
    if (obj_computation_enabled)
    {
        obj_log.clear();
        obj_log.push_back(compute_obj());
    }

    double kappa = 1/t;
    double duality_gap = 2*N*kappa;

    int_loop_counter = 0;
    ext_loop_counter = 0;
    bs_counter = 0;

    while (duality_gap > tol_out)
    {
        ++ext_loop_counter;
        while (int_loop_counter < max_iter)
        {
            ++int_loop_counter;
            if(!solve_onestep(kappa, obj_log))
            {
                break;
            }
        }
        if (int_loop_counter == max_iter)
        {
            break;
        }

        kappa /= mu;
        duality_gap = 2*N*kappa;
    }
}


double qp_ip::form_decrement()
{
    double decrement_pos = 0;
    double decrement_vel = 0;
    double decrement_acc = 0;
    double decrement_jerk = 0;
    for (int i = 0, j = 0; i < N*2; i += 2, j += SMPC_NUM_STATE_VAR)
    {
        decrement_pos += dX[j]   * dX[j]   / i2hess[i]
                       + dX[j+3] * dX[j+3] / i2hess[i+1];

        decrement_vel += dX[j+1] * dX[j+1]
                       + dX[j+4] * dX[j+4];

        decrement_acc += dX[j+2] * dX[j+2]
                       + dX[j+5] * dX[j+5];
    }
    for (int i = N*SMPC_NUM_STATE_VAR; i < N*SMPC_NUM_VAR; i += SMPC_NUM_CONTROL_VAR)
    {
        decrement_jerk += dX[i]   * dX[i]
                        + dX[i+1] * dX[i+1];
    }
    return (decrement_pos + decrement_vel/i2Q[1] + decrement_acc/i2Q[2] + decrement_jerk/i2P);
}



/**
 * @brief One step of interior point method.
 *
 * @param[in] kappa logarithmic barrier multiplier
 * @param[in,out] obj_log a vector of objective function values
 *
 * @return true if a step was made, false if alpha or dX are too small.
 */
bool qp_ip::solve_onestep (const double kappa, vector<double> &obj_log)
{
    /// Value of phi(X), where phi is the cost function + log barrier.
    double phi_X;
    phi_X = form_grad_i2hess_logbar (kappa);
    phi_X += form_phi_X ();

    chol.solve (*this, i2hess_grad, i2hess, X, dX);


    // stopping criterion (decrement)
    if (form_decrement () < tol)
    {
        return (false);
    }


    // A number from 0 to 1, which controls depth of descent #X = #X + #alpha*#dX.
    double alpha = init_alpha ();
    // stopping criterion (step size)
    if (alpha < tol)
    {
        return (false); // done
    }


    // backtracking search
    const double bs_alpha_grad_dX = form_bs_alpha_grad_dX ();
    for (;;)
    {
        ++bs_counter;
        if (form_phi_X_tmp (kappa, alpha) <= phi_X + alpha * bs_alpha_grad_dX)
        {
            break;
        }

        alpha = bs_beta * alpha;

        // stopping criterion (step size)
        if (alpha < tol)
        {
            return (false); // done
        }
    }


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

    return (true);
}



/**
 * @brief Generates an initial feasible point. 
 * First we perform a change of variable to @ref pX_tilde "X_tilde"
 * generate a feasible point, and then we go back to @ref pX_bar "X_bar".
 *
 * @param[in] x_coord x coordinates of points satisfying constraints
 * @param[in] y_coord y coordinates of points satisfying constraints
 * @param[in] init_state current state
 * @param[in] state_tilde if true the state is assumed to be in @ref pX_tilde "X_tilde" form
 * @param[in,out] X_ initial guess / solution of optimization problem
 */
void qp_ip::form_init_fp (
        const double *x_coord, 
        const double *y_coord, 
        const double *init_state,
        const bool tilde_state,
        double* X_)
{
    X = X_;
    form_init_fp_tilde<problem_parameters>(*this, x_coord, y_coord, init_state, tilde_state, X);

    // go back to bar states
    double *cur_state = X;
    for (int i=0; i<N; i++)
    {
        state_handling::tilde_to_bar (spar[i].sin, spar[i].cos, cur_state);
        cur_state = &cur_state[SMPC_NUM_STATE_VAR];
    }
}


/**
 * @brief Computes value of the objective function.
 *
 * @return value of the objective function.
 */
double qp_ip::compute_obj()
{
    int i,j;
    double obj_pos = 0;
    double obj_vel = 0;
    double obj_acc = 0;
    double obj_jerk = 0;
    double obj_gX = 0;
    double obj_ref = 0;

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

        // g'*X
        obj_gX += g[j]*X_copy[0] + g[j+1]*X_copy[3];

        const int k = i / SMPC_NUM_STATE_VAR;
        obj_ref += zref_x[k]*zref_x[k] + zref_y[k]*zref_y[k];
    }
    for (; i < N*SMPC_NUM_VAR; i += SMPC_NUM_CONTROL_VAR)
    {
        // X'*H*X
        obj_jerk += X[i] * X[i] + X[i+1] * X[i+1];
    }

    return (Q[0]*obj_pos 
            + Q[1]*obj_vel 
            + Q[2]*obj_acc 
            + P*obj_jerk 
            + obj_gX
            + Q[0]*obj_ref);
}
