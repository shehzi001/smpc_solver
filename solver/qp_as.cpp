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
    @param[in] regularization regularization
    @param[in] tol_ tolerance
*/
qp_as::qp_as(int N_, double Alpha, double Beta, double Gamma, double regularization, double tol_) : chol (N_)
{
    N = N_;
    tol = tol_;

    gain_alpha = Alpha;
    gain_beta  = Beta;
    gain_gamma = Gamma;

    chol_param.i2Q[0] = 1/(2*(Beta/2));
    chol_param.i2Q[1] = 1/(2*(Alpha/2));
    chol_param.i2Q[2] = 1/(2*regularization);

    chol_param.i2P = 1/(2 * (Gamma/2));

    chol_param.iHg = new double[2*N];

    chol_param.T = NULL;
    chol_param.h = NULL;
#ifdef QPAS_VARIABLE_T_h
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
    @param[in] lb array of lower bounds for z_x and z_y
    @param[in] ub array of upper bounds for z_x and z_y
    @param[in] X_tilde current state
    @param[in,out] X_ initial guess / solution of optimization problem
*/
void qp_as::init(
        double* T_, 
        double* h_, 
        double* angle,
        double* zref_x,
        double* zref_y,
        double* lb,
        double* ub,
        double* X_tilde,
        double* X_)
{
    nW = 0;

    chol_param.T = T_;
    chol_param.h = h_;
#ifdef QPAS_VARIABLE_T_h
    chol_param.dh = new double[N-1];
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

    form_init_fp (zref_x, zref_y, X_tilde);
    form_iHg (zref_x, zref_y);
}



/** 
 * @brief Generates an initial feasible point. 
 * First we perform a change of variable to @ref pX_tilde "X_tilde"
 * generate a feasible point, and then we go back to @ref pX_bar "X_bar".
 */
void qp_as::form_init_fp(double *zref_x, double *zref_y, double *X_tilde)
{
    double *control = &X[NUM_STATE_VAR*N];
    double *cur_state = X;
    double *prev_state = X_tilde;

#ifndef QPAS_VARIABLE_T_h    
    //------------------------------------
    double T = chol_param.T[0];
    double T2 = T*T/2;
    double h = chol_param.h[0];


    /* Control matrix. */
    double B[3] = {T2*T/3 - h*T, T2, T};


    /* inv(Cp*B). This is a [2 x 2] diagonal matrix (which is invertible if T^3/6-h*T is
     * not equal to zero). The two elements one the main diagonal are equal, and only one of them 
     * is stored, which is equal to
        1/(T^3/6 - h*T)
     */
    double iCpB = 1/(B[0]);


    /* inv(Cp*B)*Cp*A. This is a [2 x 6] matrix with the following structure
        iCpB_CpA = [a b c 0 0 0;
                    0 0 0 a b c];

        a = iCpB
        b = iCpB*T
        c = iCpB*T^2/2
     * Only a,b and c are stored.
     */
    double iCpB_CpA[3] = {iCpB, iCpB*T, iCpB*T2};
    //------------------------------------
#endif /*QPAS_VARIABLE_T_h*/

    
    for (int i=0; i<N; i++)
    {
#ifdef QPAS_VARIABLE_T_h
        //------------------------------------
        double T = chol_param.T[i];
        double T2 = T*T/2;
        double h = chol_param.h[i];


        /* Control matrix. */
        double B[3] = {T2*T/3 - h*T, T2, T};


        /* inv(Cp*B). This is a [2 x 2] diagonal matrix (which is invertible if T^3/6-h*T is
         * not equal to zero). The two elements one the main diagonal are equal, and only one of them 
         * is stored, which is equal to
            1/(T^3/6 - h*T)
         */
        double iCpB = 1/(B[0]);


        /* inv(Cp*B)*Cp*A. This is a [2 x 6] matrix with the following structure
            iCpB_CpA = [a b c 0 0 0;
                        0 0 0 a b c];

            a = iCpB
            b = iCpB*T
            c = iCpB*T^2/2
         * Only a,b and c are stored.
         */
        double iCpB_CpA[3] = {iCpB, iCpB*T, iCpB*T2};
        //------------------------------------
#endif /*QPAS_VARIABLE_T_h*/


        control[0] = -iCpB_CpA[0]*prev_state[0] - iCpB_CpA[1]*prev_state[1] - iCpB_CpA[2]*prev_state[2] + iCpB*zref_x[i];
        control[1] = -iCpB_CpA[0]*prev_state[3] - iCpB_CpA[1]*prev_state[4] - iCpB_CpA[2]*prev_state[5] + iCpB*zref_y[i];

        cur_state[0] = prev_state[0] + T*prev_state[1] + T2*prev_state[2] + B[0]*control[0];
        cur_state[1] =                   prev_state[1] +  T*prev_state[2] + B[1]*control[0];
        cur_state[2] =                                      prev_state[2] + B[2]*control[0];
        cur_state[3] = prev_state[3] + T*prev_state[4] + T2*prev_state[5] + B[0]*control[1];
        cur_state[4] =                   prev_state[4] +  T*prev_state[5] + B[1]*control[1];
        cur_state[5] =                                      prev_state[5] + B[2]*control[1];


        prev_state = &X[NUM_STATE_VAR*i];
        cur_state = &X[NUM_STATE_VAR*(i+1)];
        control = &control[NUM_CONTROL_VAR];
    }


    // go back to bar states
    cur_state = X;
    for (int i=0; i<N; i++)
    {
        state_handling::tilde_to_bar (cur_state, chol_param.angle_sin[i], chol_param.angle_cos[i]);
        cur_state = &cur_state[NUM_STATE_VAR];
    }
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
        chol_param.iHg[i*2] = 
            -chol_param.i2Q[0] * (cosA*p0 + sinA*p1)*gain_beta;
        chol_param.iHg[i*2 + 1] = 
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
int qp_as::choose_excl_constr (double *lambda)
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
