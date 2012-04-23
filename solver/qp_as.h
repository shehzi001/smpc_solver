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
#include "as_chol_solve.h"
#include "as_constraint.h"
#include "as_problem_param.h"

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
class qp_as : public AS::problem_parameters
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


        void solve ();
        void formInitialFP (
                const double *, 
                const double *, 
                const double *, 
                const bool, 
                double *);


        /** Variables for the QP (contain the states + control variables).
            Initial feasible point with respect to the equality and inequality 
            constraints. */
        double *X;


    // counters
        unsigned int added_constraints_num;
        unsigned int removed_constraints_num;
        unsigned int active_set_size;
    // limits
        bool constraint_removal_enabled;
        unsigned int max_added_constraints_num;


    private:

// functions        
        int check_blocking_constraints();

        int choose_excl_constr (const double *);

        /// An instance of #as_chol_solve class.
        AS::chol_solve chol;


        const double *zref_x;
        const double *zref_y;

        /// tolerance
        double tol;

    // active set        
        /// A set of active constraints.
        vector <AS::constraint> active_set;

        /// Vector of constraints.
        vector <AS::constraint> constraints;


    // descent direction
        /** Feasible descent direction (to be used for updating #X). */
        double *dX;

        /** A number from 0 to 1, which controls depth of descent #X = #X + #alpha*#dX. */
        double alpha;
};

///@}
#endif /*QPAS_H*/
