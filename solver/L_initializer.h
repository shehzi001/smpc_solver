/**
 * @file
 * @author Alexander Sherikov
 * @date 28.08.2011 13:43:12 MSD
 */


#ifndef L_INITIALIZER_H
#define L_INITIALIZER_H


/****************************************
 * INCLUDES 
 ****************************************/

#include "smpc_common.h"


/****************************************
 * DEFINES
 ****************************************/


/****************************************
 * TYPEDEFS 
 ****************************************/

using namespace std;

/// @addtogroup gINTERNALS
/// @{

/**
 * @brief Initializes lower diagonal matrix @ref pCholesky "L".
 */
class L_initializer
{
    public:
        /*********** Constructors / Destructors ************/
        L_initializer();
        ~L_initializer();

        void form_L (const solver_parameters*, const int, double *);


    private:
        void chol_dec (double *);

        void form_iQBiPB (const double *, const double *, const double);
        void form_iQAT (const double, const double, const double *);
        void form_AiQATiQBiPB (const double, const double);

        void form_L_non_diag(const double *, double *);
        void form_L_diag(const double *, double *);


        // intermediate results used in computation of L
        double *iQBiPB;     /// inv(Q) + B * inv(P) * B'
        double *iQAT;       /// inv(Q) * A'
        double *AiQATiQBiPB;/// A * inv(Q) * A' + inv(Q) + B * inv(P) * B'
};
/// @}

#endif /*L_INITIALIZER_H*/
