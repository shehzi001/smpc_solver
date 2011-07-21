#include <iostream>
#include <fstream>

#include <cmath> // abs
#include <cstdio>


#include "WMG.h"

// QPAS_VARIABLE_AB is taken from this header, the results
// must be the same
#include "smpc_common.h" 

#include "qp_as.h"


using namespace std;


int main()
{
    //-----------------------------------------------------------
    // initialize
    WMG wmg(15, 0.1, 0.261);

    double d[4] = {0.09 , 0.025, 0.03, 0.075};
    wmg.AddFootstep(0.0, 0.05, 0.0, 3, 3, 1, d);

    double z = 5.0*M_PI/180.0;
    double step_x = 0.035;
    double step_y = 0.1;

    d[3] = 0.025;
    wmg.AddFootstep(0.0   , -step_y, 0.0 , 4,  4, -1, d);
    wmg.AddFootstep(step_x,  step_y, z);
    wmg.AddFootstep(step_x, -step_y, z);
    wmg.AddFootstep(step_x,  step_y, z);
    wmg.AddFootstep(step_x, -step_y, z);
    wmg.AddFootstep(step_x,  step_y, z);
    wmg.AddFootstep(step_x, -step_y, z);
    wmg.AddFootstep(step_x,  step_y, z);
    wmg.AddFootstep(step_x, -step_y, z);
    wmg.AddFootstep(step_x,  step_y, 0.0, 30, 30);
    wmg.AddFootstep(0.0   , -step_y, 0.0);
    //-----------------------------------------------------------
  


    qp_as solver(wmg.N);

    double err = 0;
    double max_err = 0;
   
    double angle[wmg.N];
    double zref_x[wmg.N];
    double zref_y[wmg.N];
    double lb[2*wmg.N];
    double ub[2*wmg.N];
  
    // reference states generated using thr implementation of
    // the algorithm in Octave/MATLAB
    ifstream inFile;
    inFile.open ("./data/test_01_1_states.dat");


    printf ("\n################################\n test 01\n################################\n");
    for(;;)
    {
        //------------------------------------------------------
        wmg.FormPreviewWindow();    
        if (wmg.halt)
        {
            cout << "EXIT (halt = 1)" << endl;
            break;
        }
        wmg.form_FP_init(); 

        for (int i = 0; i < wmg.N; i++)
        {
            angle[i] = wmg.FS[wmg.ind[i]].angle;
            zref_x[i] = wmg.FS[wmg.ind[i]].p.x;
            zref_y[i] = wmg.FS[wmg.ind[i]].p.y;

            lb[i*2] = -wmg.FS[wmg.ind[i]].ctr.d[2];
            ub[i*2] = wmg.FS[wmg.ind[i]].ctr.d[0];

            lb[i*2 + 1] = -wmg.FS[wmg.ind[i]].ctr.d[3];
            ub[i*2 + 1] = wmg.FS[wmg.ind[i]].ctr.d[1];
        }
        //------------------------------------------------------


//**************************************************************************
// SOLVER IS USED HERE
//**************************************************************************
#ifdef QPAS_VARIABLE_AB
        solver.init(wmg.Tk, wmg.hk, angle, zref_x, zref_y, lb, ub, wmg.FP_init);
#else
        solver.init(wmg.T, wmg.h, angle, zref_x, zref_y, lb, ub, wmg.FP_init);
#endif
        solver.solve();
//**************************************************************************


        //------------------------------------------------------
        // compare with reference results
        for (int i = 0; i < 6; i++)
        {
            double dataref;

            wmg.X[i] = wmg.FP_init[i];

            inFile >> dataref;
            err = abs(wmg.X[i] - dataref);
            if (err > max_err)
            {
                max_err = err;
            }
            printf("value: % 8e   ref: % 8e   err: % 8e\n", wmg.X[i], dataref, err);
        }
        cout << "Max. error (over all steps): " << max_err << endl;
        //------------------------------------------------------

        wmg.slide();
    }
    inFile.close();
    printf ("################################\n");

    return 0;
}
