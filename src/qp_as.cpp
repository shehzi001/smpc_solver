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
        double* T_, 
        double* h_, 
        double* angle,
        double* zref_x,
        double* zref_y,
        double* lb,
        double* ub,
        double* X_)
{
    nW = 0;

#ifdef QPAS_VARIABLE_T_h
    chol_param.T = T_;
    chol_param.h = h_;
    chol_param.dh = new double[N-1]();
    for (int i = 0; i < N-1; i++)
    {
        chol_param.dh[i] = chol_param.h[i+1] - chol_param.h[i];
    }
#else
    chol_param.T = T_[0];
    chol_param.h = h_[0];
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

    form_init_fp (zref_x, zref_y);
    form_iHg (zref_x, zref_y);
#ifdef SMPCS_DEBUG
    zref_x_copy = zref_x;
    zref_y_copy = zref_y;
#endif /*SMPCS_DEBUG*/
}



/** 
 * @brief Generates an initial feasible point. 
 * First we perform a change of variable to "tilde states",
 * generate a feasible point, and then we go back to "bar states".
 */
        /** \brief The state of the linear model defined as follows 
            
            \verbatim
            X[0] - x ZMP position [meter]
            X[1] - x CoM velocity [meter/s]
            X[2] - x CoM acceleration [meter/s^2]
            X[3] - y ZMP position [meter]
            X[4] - y CoM velocity [meter/s]
            X[5] - y CoM acceleration [meter/s^2]
            \endverbatim

            \note This is tilde{c}_k in the paper. It is used when generating an initial feasible point.
        */
        /** \brief The state of the linear model defined as follows 
            \verbatim
            X[0] - x (bar) ZMP position [meter]
            X[1] - x CoM velocity [meter/s]
            X[2] - x XoM acceleration [meter/s^2]
            X[3] - y (bar) ZMP position [meter]
            X[4] - y CoM velocity [meter/s]
            X[5] - y CoM acceleration [meter/s^2]
            \endverbatim

            \note This is bar{c}_k in the paper
        */
void qp_as::form_init_fp(double *zref_x, double *zref_y)
{
    double *control = &X[NUM_STATE_VAR*N];
    double *cur_state = X;
    double *prev_state = X;


    for (int i=0; i<N; i++)
    {
        //------------------------------------
    ///@todo constant T
        double T = chol_param.T[i];
        double T2 = T*T/2;
        double h = chol_param.h[i];


        /** \brief Control matrix. */
        double B[3] = {T2*T/3 - h*T, T2, T};


        /** \brief inv(Cp*B). This is a [2 x 2] diagonal matrix (which is invertible if #T^3/6-#h*#T is
         * not equal to zero). The two elements one the main diagonal are equal, and only one of them 
         * is stored, which is equal to
            \verbatim
            1/(T^3/6 - h*T)
            \endverbatim      
         */
        double iCpB = 1/(B[0]);


        /** \brief inv(Cp*B)*Cp*A. This is a [2 x 6] matrix with the following structure
            \verbatim
            iCpB_CpA = [a b c 0 0 0;
                        0 0 0 a b c];

            a = iCpB
            b = iCpB*T
            c = iCpB*T^2/2
            \endverbatim

         * Only a,b and c are stored.
         */
        double iCpB_CpA[3] = {iCpB, iCpB*T, iCpB*T2};
        //------------------------------------


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
/*
        double tmp   =  chol_param.angle_cos[i]*cur_state[0] + chol_param.angle_sin[i]*cur_state[3];
        cur_state[3] = -chol_param.angle_sin[i]*cur_state[0] + chol_param.angle_cos[i]*cur_state[3];
        cur_state[0] = tmp;
*/
        tilde_to_bar (cur_state, chol_param.angle_sin[i], chol_param.angle_cos[i]);
        cur_state = &cur_state[NUM_STATE_VAR];
    }
}

void qp_as::tilde_to_bar (double *state, double sinA, double cosA)
{
    double tmp =  cosA*state[0] + sinA*state[3];
    state[3]   = -sinA*state[0] + cosA*state[3];
    state[0]   = tmp;
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
 * @return numebr of activated constraints
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


/**
 * @brief Determines coordinates of ZMP and CoM.
 *
 * @param[out] ZMP_x x coordinate of ZMP.
 * @param[out] ZMP_y y coordinate of ZMP.
 * @param[out] CoM_x x coordinate of CoM.
 * @param[out] CoM_y y coordinate of CoM.
 */
void qp_as::get_ZMP_CoM (double *ZMP_x, double *ZMP_y, double *CoM_x, double *CoM_y)
{
    *ZMP_x = chol_param.angle_cos[0] * X[0] - chol_param.angle_sin[0] * X[3];
    *ZMP_y = chol_param.angle_sin[0] * X[0] + chol_param.angle_cos[0] * X[3];

#ifdef QPAS_VARIABLE_T_h
    *CoM_x = *ZMP_x + chol_param.h[0] * X[2];
    *CoM_y = *ZMP_y + chol_param.h[0] * X[5];
#else
    *CoM_x = *ZMP_x + chol_param.h * X[2];
    *CoM_y = *ZMP_y + chol_param.h * X[5];
#endif
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
