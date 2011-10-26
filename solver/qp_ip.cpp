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


#include <cmath> // log, fabs

/****************************************
 * FUNCTIONS
 ****************************************/


//==============================================
// qp_ip

/** @brief Constructor: initialization of the constant parameters

    @param[in] N_ Number of sampling times in a preview window
    @param[in] Alpha Velocity gain
    @param[in] Beta Position gain
    @param[in] Gamma Jerk gain
    @param[in] regularization regularization
    @param[in] tol_ tolerance
*/
qp_ip::qp_ip(
        const int N_, 
        const double Alpha, 
        const double Beta, 
        const double Gamma, 
        const double regularization, 
        const double tol_) : 
    qp_solver (N_, Alpha, Beta, Gamma, regularization, tol_),
    chol (N_)
{
    g = new double[2*N];
    i2hess = new double[2*N];
    i2hess_grad = new double[N*NUM_VAR];
    grad = new double[N*NUM_VAR];

    Q2[0] = Beta;
    Q2[1] = Alpha;
    Q2[2] = regularization*2;
    P2 = Gamma;
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
}


/** @brief Initializes quadratic problem.

    @param[in] T Sampling time (for the moment it is assumed to be constant) [sec.]
    @param[in] h Height of the Center of Mass divided by gravity
    @param[in] angle Rotation angle for each state in the preview window
    @param[in] zref_x reference values of z_x
    @param[in] zref_y reference values of z_y
    @param[in] lb_ array of lower bounds for z_x and z_y
    @param[in] ub_ array of upper bounds for z_x and z_y
*/
void qp_ip::set_parameters(
        const double* T, 
        const double* h, 
        const double* angle,
        const double* zref_x,
        const double* zref_y,
        const double* lb_,
        const double* ub_)
{
    set_state_parameters (T, h, angle);

    lb = lb_;
    ub = ub_;

    form_g (zref_x, zref_y);
}



/**
 * @brief Forms vector @ref pg "g".
 *
 * @param[in] zref_x x coordinates of reference ZMP positions
 * @param[in] zref_y y coordinates of reference ZMP positions
 */
void qp_ip::form_g (const double *zref_x, const double *zref_y)
{
    double p0, p1;
    double cosA, sinA;

    for (int i = 0; i < N; i++)
    {
        cosA = spar[i].cos;
        sinA = spar[i].sin;

        // zref
        p0 = zref_x[i];
        p1 = zref_y[i];

        // inv (2*H) * R' * Cp' * zref
        g[i*2] = -(cosA*p0 + sinA*p1)*gain_beta;
        g[i*2 + 1] = -(-sinA*p0 + cosA*p1)*gain_beta; 
    }
}



/**
 * @brief Compute gradient of phi, varying elements of i2hess, 
 *  logarithmic barrier of phi.
 *
 * @param[in] kappa 1/t, a logarithmic barrier multiplicator.
 */
void qp_ip::form_grad_i2hess_logbar (const double kappa)
{
    phi_X = 0;
    int i,j;


    // grad = H*X + g + kappa * b;
    // initialize grad = H*X
    for (i = 0; i < N*NUM_STATE_VAR; i++)
    {
        grad[i] = Q2[i%3] * X[i];
    }
    for (; i < N*NUM_VAR; i++)
    {
        grad[i] = P2 * X[i];
    }

    // finish initalization of grad 
    // initialize inverted hessian
    // initialize logarithmic barrier in the function
    for (i = 0, j = 0; i < 2*N; i++, j += 3)
    {
        double lb_diff = -lb[i] + X[j];
        double ub_diff =  ub[i] - X[j];

        // logarithmic barrier
        phi_X -= log(lb_diff) + log(ub_diff);

        lb_diff = 1/lb_diff;
        ub_diff = 1/ub_diff;

        // grad += g + kappa * (ub_diff - lb_diff)
        grad[j] += g[i] + kappa * (ub_diff - lb_diff);

        // only elements 1:3:N*NUM_STATE_VAR on the diagonal of hessian 
        // can change
        // hess = H + kappa * (ub_diff^2 - lb_diff^2)
        i2hess[i] = 1/(Q2[0] + kappa * (ub_diff*ub_diff + lb_diff*lb_diff));
    }
    phi_X *= kappa;
}



/**
 * @brief Finish initialization of i2hess_grad = -i2hess*grad.
 */
void qp_ip::form_i2hess_grad ()
{
    int i,j;
    for (i = 0, j = 0; i < N*2; i++)
    {
        i2hess_grad[j] = - grad[j] * i2hess[i]; 
        j++;
        i2hess_grad[j] = - grad[j] * i2Q[1]; 
        j++;
        i2hess_grad[j] = - grad[j] * i2Q[2]; 
        j++;
    }
    for (i = N*NUM_STATE_VAR; i < N*NUM_VAR; i++)
    {
        i2hess_grad[i] = - grad[i] * i2P;
    }
}



/**
 * @brief Compute phi_X for initial point, phi_X must already store
 *      logarithmic barrier term.
 */
void qp_ip::form_phi_X ()
{
    int i;

    // phi_X += X'*H*X
    for (i = 0; i < N*NUM_STATE_VAR; i++)
    {
        phi_X += 0.5*Q2[i%3] * X[i] * X[i];
    }
    for (; i < N*NUM_VAR; i++)
    {
        phi_X += 0.5*P2 * X[i] * X[i];
    }

    // phi_X += g'*X
    for (i = 0; i < 2*N; i++)
    {
        phi_X += g[i] * X[i*3];
    }
}


/**
 * @brief Find initial value of alpha.
 *
 * @note sets alpha to 0, if it is too small.
 */
void qp_ip::init_alpha()
{
    double min_alpha = 1;
    alpha = 1;

    for (int i = 0; i < 2*N; i++)
    {
        // lower bound may be violated
        if (dX[i*3] < 0)
        {
            double tmp_alpha = (lb[i]-X[i*3])/dX[i*3];
            if (tmp_alpha < min_alpha)
            {
                min_alpha = tmp_alpha;
            }
        }
        // upper bound may be violated
        else if (dX[i*3] > 0)
        {
            double tmp_alpha = (ub[i]-X[i*3])/dX[i*3];
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
}



/**
 * @brief Forms bs_alpha * grad' * dX.
 *
 * @return result of multiplication.
 */
double qp_ip::form_bs_alpha_grad_dX ()
{
    double res = 0;
    
    for (int i = 0; i < N*NUM_VAR; i++)
    {
        res += grad[i]*dX[i];
    }

    return (res*bs_alpha);
}


/**
 * @brief Forms phi(X+alpha*dX)
 *
 * @param[in] kappa logarithmic barrier multiplicator.
 *
 * @return a value of phi.
 */
double qp_ip::form_phi_X_tmp (const double kappa)
{
    int i,j;
    double res = 0;
    double X_tmp;


    // phi_X += X'*H*X
    for (i = 0,j = 0; i < 2*N; i++)
    {
        X_tmp = X[j] + alpha * dX[j];
        j++;

        // logarithmic barrier
        res -= kappa * (log(-lb[i] + X_tmp) + log(ub[i] - X_tmp));

        // phi_X += g'*X
        res += g[i] * X_tmp;


        // phi_X += X'*H*X // states
        res += 0.5*Q2[0] * X_tmp*X_tmp;

        X_tmp = X[j] + alpha * dX[j];
        j++;
        res += 0.5*Q2[1] * X_tmp*X_tmp;

        X_tmp = X[j] + alpha * dX[j];
        j++;
        res += 0.5*Q2[2] * X_tmp*X_tmp;
    }
    // phi_X += X'*H*X // controls
    for (i = N*NUM_STATE_VAR; i < N*NUM_VAR; i++)
    {
        X_tmp = X[i] + alpha * dX[i];
        res += 0.5*P2 * X_tmp * X_tmp;
    }

    return (res);
}



/**
 * @brief Set parameters of interior-point method.
 *
 * @param[in] t_ logarithmic barrier parameter
 * @param[in] mu_ multiplier of t, >1.
 * @param[in] bs_alpha_ backtracking search parameter alpha
 * @param[in] bs_beta_  backtracking search parameter beta
 * @param[in] max_iter_ maximum number of internal loop iterations
 */
void qp_ip::set_ip_parameters (
        const double t_, 
        const double mu_, 
        const double bs_alpha_, 
        const double bs_beta_, 
        const int max_iter_)
{
    t = t_;
    mu = mu_;
    bs_alpha = bs_alpha_;
    bs_beta = bs_beta_;
    max_iter = max_iter_;
}



/**
 * @brief Solve QP using interior-point method.
 *
 * @return 0 if ok, negative number otherwise.
 */
int qp_ip::solve()
{
    double kappa = 1/t;
    double duality_gap = 2*N*kappa;
    int i;


    while (duality_gap > tol)
    {
        for (i = 0; (i < max_iter) && solve_onestep(kappa); i++);
        if (i == max_iter)
        {
            return (-1);
        }

        kappa /= mu;
        duality_gap = 2*N*kappa;
    }
    return (0);
}


/**
 * @brief One step of interior point method.
 *
 * @param[in] kappa logarithmic barrier multiplier
 *
 * @return true if a step was made, false if alpha or dX are too small.
 */
bool qp_ip::solve_onestep (const double kappa)
{
    form_grad_i2hess_logbar (kappa);
    form_phi_X ();
    form_i2hess_grad ();

    chol.solve (this, i2hess_grad, i2hess, X, dX);


    double max_dX = 0;
    for (int i = 0; i < N*NUM_VAR; i++)
    {
        if (max_dX < fabs(dX[i]))
        {
            max_dX = fabs(dX[i]);
        }
    }
    if (max_dX < tol)
    {
        return (false);
    }


    init_alpha ();
    if (alpha < tol)
    {
        return (false); // done
    }

    double bs_alpha_grad_dX = form_bs_alpha_grad_dX ();
    for (;;)
    {
        if (form_phi_X_tmp (kappa) <= phi_X + alpha * bs_alpha_grad_dX)
        {
            break;
        }

        alpha = bs_beta * alpha;

        if (alpha < tol)
        {
            //XXX hm?
            break;
        }
    }

    if (alpha < tol)
    {
        return (false); // done
    }

    // Move in the feasible descent direction
    for (int j = 0; j < N*NUM_VAR ; j++)
    {
        X[j] += alpha * dX[j];
    }

    return (true);
}
