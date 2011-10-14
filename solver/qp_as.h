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
        void set (const int, const double, const double, const int);


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
class qp_as : public qp_solver
{
    public:
// functions        
        qp_as(
                const int N_, 
                const double Alpha = 150.0, 
                const double Beta = 2000.0, 
                const double Gamma = 1.0,
                const double regularization = 0.01,
                const double tol = 1e-7);
        ~qp_as();

        void set_parameters(
                const double*, 
                const double*, 
                const double*, 
                const double*, 
                const double*, 
                const double*, 
                const double*);


        int solve ();


    private:

// functions        
        void form_iHg(const double *, const double *);
        void form_bounds(const double *, const double *);
        int check_blocking_bounds();

#ifdef QPAS_DOWNDATE
        int choose_excl_constr (const double *);
#endif

        /// An instance of #chol_solve class.
        chol_solve_as chol;


        /// @ref piHg "inv(H) * g"
        double *iHg;


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
