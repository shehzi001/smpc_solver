/**
 * @file
 * @brief 
 *
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

class L_initializer
{
    public:
        /*********** Constructors / Destructors ************/
        L_initializer();
        ~L_initializer();

        void form_L(chol_solve_param, int, double *);


    private:
        void chol_dec (double *);

        void form_iQBiPB (double *, double *, double);
        void form_iQAT (double, double, double *);
        void form_AiQATiQBiPB (double, double);

        void form_L_non_diag(double *, double *);
        void form_L_diag(double *, double *);


        // intermediate results used in computation of L
        double *iQBiPB;     /// inv(Q) + B * inv(P) * B'
        double *iQAT;       /// inv(Q) * A'
        double *AiQATiQBiPB;/// A * inv(Q) * A' + inv(Q) + B * inv(P) * B'
};

#endif /*L_INITIALIZER_H*/
