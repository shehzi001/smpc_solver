/** 
 * @file
 * @brief  
 *
 * @author Alexander Sherikov
 * @date 19.07.2011 22:30:13 MSD
 * @todo add description
 */


/****************************************
 * INCLUDES 
 ****************************************/
#include "qp_as.h"

#include <string.h> // in order to use memcpy and memmove
#include <cmath> // cos, sin

/****************************************
 * FUNCTIONS
 ****************************************/

//==============================================
// bound
/**
 * @brief Set parameters of the bound
 *
 * @param[in] var_num_ index of variable in the vector of states
 * @param[in] lb_ lower bound
 * @param[in] ub_ upper bound
 * @param[in] active activity of the bound
 */
void bound::set(int var_num_, double lb_, double ub_, int active)
{
    var_num = var_num_;
    lb = lb_;
    ub = ub_;
    isActive = active;
}



//==============================================
// qp_as

/** @brief Constructor: initialization of the constant parameters

    @param[in] N_ Number of sampling times in a preview window
    @param[in] Alpha Velocity gain
    @param[in] Beta Position gain
    @param[in] Gamma Jerk gain
*/
qp_as::qp_as(int N_, double Alpha, double Beta, double Gamma) : chol (N_)
{
    N = N_;

    gain_alpha = Alpha;
    gain_beta  = Beta;
    gain_gamma = Gamma;

    chol_param.i2Q[0] = 1/(2*(Beta/2));
    chol_param.i2Q[1] = 1/(2*(Alpha/2));
    chol_param.i2Q[2] = 1/(2*REGULARIZATION);

    chol_param.i2P = 1/(2 * (Gamma/2));

    chol_param.iHg = new double[NUM_VAR*N]();

#ifdef QPAS_VARIABLE_T_h
    chol_param.T = NULL;
    chol_param.h = NULL;
    chol_param.dh = NULL;
#endif

    chol_param.angle_cos = NULL;
    chol_param.angle_sin = NULL;


    // there are no inequality constraints in the initial working set (no hot-starting for the moment)
    W = new int[2*N];
#ifdef QPAS_DOWNDATE
    W_sign = new int[2*N];
#endif
    alpha = 1;
    
    dX = new double[NUM_VAR*N]();

    Bounds.resize(2*N);
}


/** Destructor */
qp_as::~qp_as()
{
    if (W  != NULL)
        delete [] W;

    if (dX  != NULL)
        delete [] dX;

    if (chol_param.angle_cos != NULL)
        delete chol_param.angle_cos;

    if (chol_param.angle_sin != NULL)
        delete chol_param.angle_sin;

    if (chol_param.iHg != NULL)
        delete chol_param.iHg;

#ifdef QPAS_VARIABLE_T_h
    if (chol_param.dh != NULL)
        delete chol_param.dh;
#endif
#ifdef QPAS_DOWNDATE
    if (W_sign != NULL)
        delete W_sign;
#endif
}


/** @brief Initializes quadratic problem.

    @param[in] T_ Sampling time (for the moment it is assumed to be constant) [sec.]
    @param[in] h_ Height of the Center of Mass divided by gravity
    @param[in] angle Rotation angle for each state in the preview window
    @param[in] zref_x reference values of z_x
    @param[in] zref_y reference values of z_y
    @param[in] lb array of lower bounds for z
    @param[in] ub array of upper bounds for z
    @param[in,out] X_ initial guess / solution of optimization problem
*/
void qp_as::init(
#ifdef QPAS_VARIABLE_T_h
        double* T_, 
        double* h_, 
#else
        double T_, 
        double h_, 
#endif
        double* angle,
        double* zref_x,
        double* zref_y,
        double* lb,
        double* ub,
        double* X_)
{
    nW = 0;

    chol_param.T = T_;
    chol_param.h = h_;
#ifdef QPAS_VARIABLE_T_h
    chol_param.dh = new double[N-1]();
    for (int i = 0; i < N-1; i++)
    {
        chol_param.dh[i] = chol_param.h[i+1] - chol_param.h[i];
    }
#endif


    X = X_;

    form_bounds(lb, ub);

    chol_param.angle_cos = new double[N];
    chol_param.angle_sin = new double[N];

    for (int i = 0; i < N; i++)
    {
        chol_param.angle_cos[i] = cos(angle[i]);
        chol_param.angle_sin[i] = sin(angle[i]);
    }

    form_iHg(zref_x, zref_y);
#ifdef SMPCS_DEBUG
    zref_x_copy = zref_x;
    zref_y_copy = zref_y;
#endif /*SMPCS_DEBUG*/
}


/**
 * @brief Forms inv(2*H) * g.
 *
 * @param[in] zref_x reference values of z_x
 * @param[in] zref_y reference values of z_y
 */
void qp_as::form_iHg(double *zref_x, double *zref_y)
{
    double p0, p1;
    double cosA, sinA;

    for (int i = 0; i < N; i++)
    {
        cosA = chol_param.angle_cos[i];
        sinA = chol_param.angle_sin[i];

        // zref
        p0 = zref_x[i];
        p1 = zref_y[i];

        // inv (2*H) * R' * Cp' * zref
        chol_param.iHg[i*NUM_STATE_VAR + 0] = 
            -chol_param.i2Q[0] * (cosA*p0 + sinA*p1)*gain_beta;
        chol_param.iHg[i*NUM_STATE_VAR + 3] = 
            -chol_param.i2Q[0] * (-sinA*p0 + cosA*p1)*gain_beta; 
    }
}




/**
 * @brief Forms the upper and lower bounds.
 *
 * @param[in] lb array of lower bounds for z
 * @param[in] ub array of upper bounds for z
 */
void qp_as::form_bounds(double *lb, double *ub)
{
    for (int i=0; i < N; i++)
    {
        Bounds[i*2].set(NUM_STATE_VAR*i, lb[i*2], ub[i*2], 0);
        Bounds[i*2+1].set(NUM_STATE_VAR*i+3, lb[i*2+1], ub[i*2+1], 0);
    }
}


/**
 * @brief Check for blocking bounds.
 *
 * @return sequential number of constraint to be added, -1 if no constraints.
 */
int qp_as::check_blocking_bounds()
{
    alpha = 1;

    /** Index to include in the working set #W, -1 if no constraint have to be included. */
    int activated_var_num = -1;


    for (int i = 0; i < 2*N; i++)
    {
        // Check only inactive constraints for violation. 
        // The constraints in the working set will not be violated regardless of 
        // the depth of descent
        if (Bounds[i].isActive == 0)
        {
            int ind = Bounds[i].var_num;
#ifdef QPAS_DOWNDATE
            int sign = 1;
#endif
            double t = 1;

            if ( dX[ind] < -TOL )
            {
                t = (Bounds[i].lb - X[ind])/dX[ind];
#ifdef QPAS_DOWNDATE
                sign = -1;
#endif
            }
            else if ( dX[ind] > TOL ) 
            {
                t = (Bounds[i].ub - X[ind])/dX[ind];
            }
            else
            {
                // do nothing because dX[Bounds[i].var_num] = 0 (numerically speaking)
                t = 1;
            }

            if (t < alpha)
            {
                alpha = t;
                activated_var_num = i;
#ifdef QPAS_DOWNDATE
                W_sign[nW] = sign;
#endif
            }
        }
    }
    if (activated_var_num != -1)
    {
        W[nW] = activated_var_num;    
        nW++;
        Bounds[activated_var_num].isActive = 1;
    }

    return (activated_var_num);
}  


#ifdef QPAS_DOWNDATE
int qp_as::choose_excl_constr (double *lambda)
{
    double min_lambda = -TOL;
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
        Bounds[W[ind_exclude]].isActive = 0;
        for (int i = ind_exclude; i < nW-1; i++)
        {
            W[i] = W[i + 1];
            W_sign[i] = W_sign[i + 1];
        }
        nW--;
    }

    return (ind_exclude);
}
#endif /* QPAS_DOWNDATE */


/**
 * @brief Solve QP problem.
 */
int qp_as::solve ()
{
    // obtain dX
    chol.solve(chol_param, X, dX);

    for (;;)
    {
        int activated_var_num = check_blocking_bounds();

        // Move in the feasible descent direction
        for (int i = 0; i < N*NUM_VAR ; i++)
        {
            X[i] += alpha * dX[i];
        }

        // no new inequality constraints
        if (activated_var_num == -1)
        {
#ifdef QPAS_DOWNDATE
            int ind_exclude = choose_excl_constr (chol.get_lambda());
            if (ind_exclude != -1)
            {
                chol.down_resolve (chol_param, nW, W, ind_exclude, X, dX);
            }
            else
            {
                break;
            }
#else
            break;
#endif /*QPAS_DOWNDATE*/
        }
        else
        {
            // add row to the L matrix and find new dX
            chol.up_resolve (chol_param, nW, W, X, dX);
        }
    }

    return (nW);
}


#ifdef SMPCS_DEBUG
/**
 * @brief Prints value of  X'*H*X + X'*g
 */
void qp_as::print_objective()
{
    int i;
    double obj = 0;


    // H_c
    for (i = 0; i < N*NUM_STATE_VAR; i++)
    {
        obj += (1/(2*chol_param.i2Q[i%3])) * X[i] * X[i];
    }

    // H_u
    for (; i < N*NUM_VAR; i++)
    {
        obj += (1/(2*chol_param.i2P)) * X[i] * X[i];
    }


    // X'*g
    for (i = 0; i < N; i++)
    {
        // rotation angle
        double cosA = chol_param.angle_cos[i];
        double sinA = chol_param.angle_sin[i];

        // zref
        double p0 = zref_x_copy[i];
        double p1 = zref_y_copy[i];

        // g = -[q_1 ... q_N]'
        // q = beta * R' * Cp' * zref
        obj -= X[i*NUM_STATE_VAR + 0]*(cosA*p0 + sinA*p1)*gain_beta;
        obj -= X[i*NUM_STATE_VAR + 3]*(-sinA*p0 + cosA*p1)*gain_beta; 
    }

    printf ("DEBUG: objective function = % 8e\n", obj);
}
#endif /*SMPCS_DEBUG*/
