/**
 * @file
 * @brief 
 * Contains code that (approximately) solves the following 
 * 
 * ---------------------------------
 * | Quadratic Program             |
 * ---------------------------------
 * |  min f(x) = 0.5*x'*H*x + x'*g |
 * |   x                           |
 * |                               |
 * |  st. E*x = e                  |
 * |      lb <= x <= ub            |
 * |                               |
 * | Note:                         |
 * | ------                        |
 * | H, g, E, e have structure as  |
 * | defined in my IROS11 paper    | 
 * ---------------------------------
 *
 * The "solver" uses a Primal Active Set scheme (with range space linear algebra). Inequality
 * constraints are only included in the active set (they are not removed). For the purposes of MPC
 * this appears to be reasonable (of course in general it is not).
 *
 * @author Alexander Sherikov
 * @note Based on the code originally developed by Dimitar Dimitrov.
 * @date 19.07.2011 15:56:18 MSD
 */


#ifndef QPAS_H
#define QPAS_H

/****************************************
 * INCLUDES 
 ****************************************/
#include "smpc_common.h"
#include "chol_solve.h"

#include <vector>


/****************************************
 * Defines
 ****************************************/
#define TOL 1e-7 // tolerance

#define REGULARIZATION 0.01


using namespace std;


/** \brief Defines simple bounds associated on a variable. */
class bound
{
    public:
        void set (int, double, double, int);


        /** Variable number (on which to impose the bounds). */
        int var_num;

        /** Lower bound. */
        double lb;

        /** Upper bound. */
        double ub;

        /** 
         * If isActive = 1 then one of the bounds is in the 
         * working set (if isActive = 0 it is not). 
         */
        int isActive;
};



/** 
 * Approximately solve a quadratic program with a specific structure. 
 * qp_as = Quadratic Programming / Active Set
 */
class qp_as
{
    public:
        qp_as(int N_, double Alpha = 150.0, double Beta = 2000.0, double Gamma = 1.0);
        ~qp_as();

#ifdef QPAS_VARIABLE_T_h
        void init(double*, double*, double*, double*, double*, double*, double*, double*);
#else
        void init(double, double, double*, double*, double*, double*, double*, double*);
#endif
   
        void get_ZMP_CoM (double *, double *, double *, double *);

        ///@todo return W
        int solve ();


    private:

// functions        
        void form_iHg(double *, double *);
        void initialize_bounds();
        void form_bounds(double *, double *);
        int check_blocking_bounds();

#ifdef QPAS_DOWNDATE
        int choose_excl_constr (double *);
#endif

#ifdef SMPCS_DEBUG
        void print_objective();
#endif

// variables
        /** Number of iterations in a preview window. */
        int N;

        /// Parameters, which are fed to the methods of #chol_solve class.
        chol_solve_param chol_param;

        /// An instance of #chol_solve class.
        chol_solve chol;


    // gains        
        double gain_alpha;
        double gain_beta;
        double gain_gamma;

    // variables and descent direction
        /** Variables for the QP (contain the states + control variables).
            Initial feasible point with respect to the equality and inequality constraints. */
        double *X;
     
        /** Feasible descent direction (to be used for updating #X). */
        double *dX;

        /** A number from 0 to 1, which controls depth of descent #X = #X + #alpha*#dX. */
        double alpha;

    // active set        
        /** Working set (contains the indexes of only inequality constraints). It is assumed that
            the only inequality constraints are simple bounds. Since only 2*#N variables in #X are
            subjecto to bounds, #W is defined to have size 2*#N. */
        int *W;

#ifdef QPAS_DOWNDATE
        /** 
         * Since we do not distinguish lower/upper bounds of active constraints (<= and => 
         * inequlities are treated in the same way), we have to adjust signs of lagrange 
         * multipliers before downdate. 
         */
        int *W_sign;
#endif

        /** Number of indexes of inequality constraints already included in #W. #W(#nW-1) is the
            index of the last inequality constraint added to #W. */
        int nW;

        /// Vector of bounds.
        std::vector <bound> Bounds;

#ifdef SMPCS_DEBUG
        double *zref_x_copy;
        double *zref_y_copy;
#endif
};


#endif /*QPAS_H*/
