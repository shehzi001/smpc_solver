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

/// @addtogroup gINTERNALS
/// @{

/** 
 * @brief Solve a quadratic program with a specific structure. 
 * qp_as = Quadratic Programming / Active Set
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


        void solve (const double);


    private:
        void form_g (const double *, const double *);
        double *g;
        double *i2hess;
        double *i2hess_grad;
        double phi_X;

        void form_grad_hess_logbar (const double);
        void form_i2hess_grad ();
        void form_phi_X ();

        double Q2[3];
        double P2;

        /// An instance of #chol_solve_ip class.
        chol_solve_ip chol;

// functions        
        ///@{
        /// lower and upper bounds
        const double *lb;
        const double *ub;
        ///@}
};

///@}
#endif /*QPIP_H*/
