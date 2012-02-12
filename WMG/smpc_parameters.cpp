/** 
 * @file
 * @author Alexander Sherikov
 * @date 27.09.2011 18:59:03 MSD
 */

#include "WMG.h"


/** \brief Default constructor. */
smpc_parameters::smpc_parameters ()
{
    X = NULL;

    T = NULL;
    h = NULL;

    angle = NULL;
    zref_x = NULL;
    zref_y = NULL;
    fp_x = NULL;
    fp_y = NULL;
    lb = NULL;
    ub = NULL;
}



/** \brief Default destructor. */
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



/**
 * @brief Allocate memory and initialize some of the parameters.
 *
 * @param[in] N preview window length
 * @param[in] h_ height of center of mass divided by gravity.
 */
void smpc_parameters::init(
        const unsigned int N,
        const double h_)
{
    X = new double[SMPC_NUM_VAR*N];

    T = new double[N];
    h = new double[N];

    h0 = h_;
    for (unsigned int i = 0; i < N; i++)
    {
        h[i] = h_;
    }

    angle = new double[N];
    zref_x = new double[N];
    zref_y = new double[N];
    fp_x = new double[N];
    fp_y = new double[N];
    lb = new double[2*N];
    ub = new double[2*N];
}
