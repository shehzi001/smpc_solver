/** 
 * @file
 * @author Alexander Sherikov
 * @date 27.09.2011 18:59:03 MSD
 */

#include "WMG.h"


smpc_parameters::smpc_parameters(
        const unsigned int N,
        const double hCoM_,
        const double gravity_)
{
    hCoM = hCoM_;
    gravity = gravity_;

    X = new double[SMPC_NUM_VAR*N];

    T = new double[N];
    h = new double[N];

    h0 = hCoM/gravity;
    for (unsigned int i = 0; i < N; i++)
    {
        h[i] = h0;
    }

    angle = new double[N];
    zref_x = new double[N];
    zref_y = new double[N];
    fp_x = new double[N];
    fp_y = new double[N];
    lb = new double[2*N];
    ub = new double[2*N];
}



smpc_parameters::~smpc_parameters()
{
    if (X != NULL)
    {
        delete [] X;
        X = NULL;
    }

    if (T != NULL)
    {
        delete T;
    }
    if (h != NULL)
    {
        delete h;
    }

    if (angle != NULL)
    {
        delete angle;
    }
    if (zref_x != NULL)
    {
        delete zref_x;
    }
    if (zref_y != NULL)
    {
        delete zref_y;
    }
    if (fp_x != NULL)
    {
        delete fp_x;
    }
    if (fp_y != NULL)
    {
        delete fp_y;
    }
    if (lb != NULL)
    {
        delete lb;
    }
    if (ub != NULL)
    {
        delete ub;
    }
}
