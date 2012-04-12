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
#include "constraint.h"
#include "problem_param.h"

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
class qp_as : public problem_parameters
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

        void form_init_fp(const double *, const double *, const double *, double *);


        /** Variables for the QP (contain the states + control variables).
            Initial feasible point with respect to the equality and inequality 
            constraints. */
        double *X;


    private:

// functions        
        int check_blocking_constraints();

        int choose_excl_constr (const double *);

        /// An instance of #chol_solve_as class.
        chol_solve_as chol;


        double *zref_x;
        double *zref_y;


    // active set        
        /// A set of active constraints.
        vector <constraint> active_set;

        /// Vector of constraints.
        vector <constraint> constraints;



        /// tolerance
        double tol;

    // variables and descent direction
     
        /** Feasible descent direction (to be used for updating #X). */
        double *dX;

        /** A number from 0 to 1, which controls depth of descent #X = #X + #alpha*#dX. */
        double alpha;

        /// Height of the CoM at initial state divided by the gravity, this initial state
        /// precede the first state in the preview window.
        double h_initial;
};

///@}
#endif /*QPAS_H*/
