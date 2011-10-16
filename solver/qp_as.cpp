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
void bound::set(const int var_num_, const double lb_, const double ub_, const int active)
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
#ifdef QPAS_DOWNDATE
    W_sign = new int[2*N];
#endif

    Bounds.resize(2*N);
}


/** Destructor */
qp_as::~qp_as()
{
    if (W  != NULL)
        delete W;

    if (iHg != NULL)
        delete iHg;

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
    @param[in] lb array of lower bounds for z_x and z_y
    @param[in] ub array of upper bounds for z_x and z_y
*/
void qp_as::set_parameters(
        const double* T_, 
        const double* h_, 
        const double* angle,
        const double* zref_x,
        const double* zref_y,
        const double* lb,
        const double* ub)
{
    nW = 0;

    set_state_parameters (T_, h_, angle);
    form_iHg (zref_x, zref_y);
    form_bounds(lb, ub);
}



/**
 * @brief Forms inv(2*H) * g.
 *
 * @param[in] zref_x reference values of z_x
 * @param[in] zref_y reference values of z_y
 */
void qp_as::form_iHg(const double *zref_x, const double *zref_y)
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
        iHg[i*2] =     -i2Q[0] * (cosA*p0 + sinA*p1)*gain_beta;
        iHg[i*2 + 1] = -i2Q[0] * (-sinA*p0 + cosA*p1)*gain_beta; 
    }
}



/**
 * @brief Forms the upper and lower bounds.
 *
 * @param[in] lb array of lower bounds for z
 * @param[in] ub array of upper bounds for z
 */
void qp_as::form_bounds(const double *lb, const double *ub)
{
    for (int i=0; i < N; i++)
    {
        Bounds[i*2].set(NUM_STATE_VAR*i, lb[i*2], ub[i*2], 0);
        Bounds[i*2+1].set(NUM_STATE_VAR*i+3, lb[i*2+1], ub[i*2+1], 0);
    }
}



/**
 * @brief Checks for blocking bounds.
 *
 * @return sequential number of constraint to be added, -1 if no constraints.
 */
int qp_as::check_blocking_bounds()
{
    alpha = 1;

    /* Index to include in the working set, -1 if no constraint have to be included. */
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

            if ( dX[ind] < -tol )
            {
                t = (Bounds[i].lb - X[ind])/dX[ind];
#ifdef QPAS_DOWNDATE
                sign = -1;
#endif
            }
            else if ( dX[ind] > tol ) 
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
 *
 * @return number of activated constraints
 */
int qp_as::solve ()
{
    // obtain dX
    chol.solve(this, iHg, X, dX);

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
            int ind_exclude = choose_excl_constr (chol.get_lambda(this));
            if (ind_exclude != -1)
            {
                chol.down_resolve (this, iHg, nW, W, ind_exclude, X, dX);
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
            chol.up_resolve (this, iHg, nW, W, X, dX);
        }
    }

    return (nW);
}
