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
#include "smpc_common.h"
#include "qp_solver.h"
#include "chol_solve_ip.h"

#include <vector>


/****************************************
 * Defines
 ****************************************/


using namespace std;

/// @addtogroup gIP
/// @{

/** 
 * @brief Solve a quadratic program with a specific structure. 
 * qp_ip = Quadratic Programming / Interior-point method
 */
class qp_ip : public qp_solver
{
    public:
// functions        
        qp_ip(
                const int N_, 
                const double, 
                const double, 
                const double,
                const double,
                const double);
        ~qp_ip();

        void set_parameters(
                const double*, 
                const double*, 
                const double*, 
                const double*, 
                const double*, 
                const double*, 
                const double*);


        void set_ip_parameters (
                const double, 
                const double, 
                const double, 
                const double, 
                const int,
                const double);
        int solve();


    private:
        /// 2*#N non-zero elements of vector @ref pg "g".
        double *g;

        /// Inverted hessian: non-repeating diagonal elements
        /// 1:3:#N*#SMPC_NUM_STATE_VAR, 2*#N in total.
        double *i2hess;

        /// Inverted hessian * gradient (#N*#SMPC_NUM_VAR vector)
        double *i2hess_grad;

        /// #N*#SMPC_NUM_VAR gradient vector
        double *grad;

        /// Value of phi(X), where phi is the cost function + log barrier.
        double phi_X; 


        ///@{
        /// Diagonal elements of H.
        double Q[3];
        double P;
        ///@}
        

        /// An instance of #chol_solve_ip class.
        chol_solve_ip chol;


        ///@{
        /// lower and upper bounds
        const double *lb;
        const double *ub;
        ///@}
        
// IP parameters
        double t; /// logarithmic barrier parameter
        double mu; /// multiplier of t, >1.
        double bs_alpha; /// backtracking search parameter alpha
        double bs_beta; /// backtracking search parameter beta
        int max_iter; /// maximum number of internal loop iterations
        double tol_out; /// tolerance of the outer loop


// functions        
        void init_alpha();
        double form_bs_alpha_grad_dX ();
        double form_phi_X_tmp (const double);
        bool solve_onestep (const double);
        void form_g (const double *, const double *);
        void form_grad_i2hess_logbar (const double);
        void form_i2hess_grad ();
        void form_phi_X ();
};

///@}
#endif /*QPIP_H*/
