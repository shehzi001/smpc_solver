/**
 * @file
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


using namespace std;

/// @addtogroup gINTERNALS
/// @{

/** \brief Defines simple bounds associated with variables. */
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
 * @brief Solve a quadratic program with a specific structure. 
 * qp_as = Quadratic Programming / Active Set
 */
class qp_as
{
    public:
// functions        
        qp_as(
                int N_, 
                double Alpha = 150.0, 
                double Beta = 2000.0, 
                double Gamma = 1.0,
                double regularization = 0.01,
                double tol = 1e-7);
        ~qp_as();

        void init(double*, double*, double*, double*, double*, double*, double*, double*, double*);
   

        ///@todo Do we need to return active set (W)?
        int solve ();

        void get_next_state_tilde (double *);
        void get_next_state (double *);


// variables
        /// Parameters, which are fed to the methods of #chol_solve class.
        chol_solve_param chol_param;

        /** Variables for the QP (contain the states + control variables).
            Initial feasible point with respect to the equality and inequality 
            constraints. */
        double *X;


    private:

// functions        
        void form_init_fp(double *, double *, double *);
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


        /// An instance of #chol_solve class.
        chol_solve chol;


    ///@{
    /// Gains used in @ref pPDObj "objective function".
        double gain_alpha;
        double gain_beta;
        double gain_gamma;
    ///@}

        /// tolerance
        double tol;

    // variables and descent direction
     
        /** Feasible descent direction (to be used for updating #X). */
        double *dX;

        /** A number from 0 to 1, which controls depth of descent #X = #X + #alpha*#dX. */
        double alpha;

    // active set        
        /** Working set (contains the indexes of only inequality constraints). It is assumed that
            the only inequality constraints are simple bounds. See also '@ref pBounds'. */
        int *W;

#ifdef QPAS_DOWNDATE
        /** 
         * Since we do not distinguish lower/upper bounds of active constraints (<= and => 
         * inequlities are treated in the same way), we have to adjust signs of lagrange 
         * multipliers before downdate. See also '@ref pBounds'.
         */
        int *W_sign;
#endif

        /** Number of inequality constraints already included in #W. #W(#nW-1) is the
            index of the last inequality constraint added to #W. */
        int nW;

        /// Vector of bounds.
        std::vector <bound> Bounds;
};

///@}
#endif /*QPAS_H*/
