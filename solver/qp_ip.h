/**
 * @file
 * @author Alexander Sherikov
 * @note Based on the code originally developed by Dimitar Dimitrov.
 * @date 19.07.2011 15:56:18 MSD
 */


#ifndef QPIP_H
#define QPIP_H

/****************************************
 * INCLUDES 
 ****************************************/
#include "smpc_solver.h"
#include "smpc_common.h"
#include "ip_chol_solve.h"
#include "ip_problem_param.h"

#include <vector>


/****************************************
 * Defines
 ****************************************/


using namespace std;
using namespace smpc;

/// @addtogroup gIP
/// @{

/** 
 * @brief Solve a quadratic program with a specific structure. 
 * qp_ip = Quadratic Programming / Interior-point method
 */
class qp_ip : public IP::problem_parameters
{
    public:
// functions        
        qp_ip(
                const int N_, 
                const double, 
                const double, 
                const double,
                const double,
                const double,
                const bool,
                const backtrackingSearchType);
        ~qp_ip();

        void set_parameters(
                const double*, 
                const double*, 
                const double, 
                const double*, 
                const double*, 
                const double*, 
                const double*, 
                const double*);

        void form_init_fp (
                const double *, 
                const double *, 
                const double *, 
                const bool,
                double *);


        void set_ip_parameters (
                const double, 
                const double, 
                const double, 
                const double, 
                const unsigned int,
                const double);

        void solve(vector<double> &);

        /** Variables for the QP (contain the states + control variables).
            Initial feasible point with respect to the equality and inequality 
            constraints. */
        double *X;

        
        unsigned int int_loop_counter;
        unsigned int ext_loop_counter;
        unsigned int bs_counter;


    private:
    // parameters
        double gain_position;

        /// tolerance
        double tol;

        bool obj_computation_on;
        backtrackingSearchType bs_type;


    // variables and descent direction
     
        /** Feasible descent direction (to be used for updating #X). */
        double *dX;

        /// 2*#N non-zero elements of vector @ref pg "g".
        double *g;

        /// Inverted hessian: non-repeating diagonal elements
        /// 1:3:#N*#SMPC_NUM_STATE_VAR, 2*#N in total.
        double *i2hess;

        /// Inverted hessian * gradient (#N*#SMPC_NUM_VAR vector)
        double *i2hess_grad;

        /// 2*#N gradient vector, only the elements that correspond to the ZMP
        /// positions are computed, it is faster to compute the others on the fly.
        /// Hint: the computed terms are affected by the logarithmic barrier.
        double *grad;


        ///@{
        /// Diagonal elements of H.
        double Q[3];
        double P;
        ///@}
        

        /// An instance of #chol_solve_ip class.
        IP::chol_solve chol;


        ///@{
        /// lower and upper bounds
        const double *lb;
        const double *ub;
        ///@}

        const double *zref_x;
        const double *zref_y;


// IP parameters
        double t; /// logarithmic barrier parameter
        double mu; /// multiplier of t, >1.
        double bs_alpha; /// backtracking search parameter alpha
        double bs_beta; /// backtracking search parameter beta
        unsigned int max_iter; /// maximum number of internal loop iterations (in total)
        double tol_out; /// tolerance of the outer loop


// functions        
        double init_alpha();
        double form_bs_alpha_obj_dX ();
        double form_phi_X_tmp (const double, const double);
        bool solve_onestep (const double, vector<double> &);
        void form_g (const double *, const double *);
        double form_grad_i2hess_logbar (const double);
        double form_phi_X ();
        double form_decrement();
        double compute_obj(const bool);
};

///@}
#endif /*QPIP_H*/
