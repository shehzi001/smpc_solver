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
#include "chol_solve_as.h"
#include "qp_solver.h"
#include "constraint.h"

#include <vector>


/****************************************
 * Defines
 ****************************************/


using namespace std;

/// @addtogroup gAS
/// @{

/** 
 * @brief Solve a quadratic program with a specific structure. 
 * qp_as = Quadratic Programming / Active Set
 */
class qp_as : public qp_solver
{
    public:
// functions        
        qp_as(
                const int N_, 
                const double, 
                const double, 
                const double,
                const double,
                const double);
        ~qp_as();

        void set_parameters(
                const double*, 
                const double*, 
                const double, 
                const double*, 
                const double*, 
                const double*, 
                const double*, 
                const double*);


        int solve ();


    private:

// functions        
        void form_iHg(const double *, const double *);
        void form_constraints(const double *, const double *, const double *);
        int check_blocking_constraints();

        int choose_excl_constr (const double *);

        /// An instance of #chol_solve_as class.
        chol_solve_as chol;


        /// @ref piHg "inv(H) * g"
        double *iHg;


    // active set        
        /** Working set (contains the indexes of only inequality constraints). */
        int *W;

        /** 
         * Since we do not distinguish lower/upper bounds of active constraints (<= and => 
         * inequlities are treated in the same way), we have to adjust signs of lagrange 
         * multipliers before downdate. See also '@ref pBounds'.
         */
        int *W_sign;

        /** Number of inequality constraints already included in #W. #W(#nW-1) is the
            index of the last inequality constraint added to #W. */
        int nW;

        /// Vector of constraints.
        std::vector <constraint> constraints;
};

///@}
#endif /*QPAS_H*/
