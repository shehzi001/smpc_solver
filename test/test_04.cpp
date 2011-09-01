#include <iostream>
#include <fstream>

#include <cmath> // abs
#include <cstdio>


#include "WMG.h"

// QPAS_VARIABLE_T_h is taken from this header, the results
// must be the same
#include "smpc_solver.h" 


#define PREVIEW_SIZE 15 // Size of the preview window

using namespace std;


int main(int argc, char **argv)
{
    //-----------------------------------------------------------
    // initialize
    WMG wmg(PREVIEW_SIZE, 0.1, 0.261);

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
  


    smpc_solver solver(PREVIEW_SIZE);

    double angle[PREVIEW_SIZE];
    double zref_x[PREVIEW_SIZE];
    double zref_y[PREVIEW_SIZE];
    double lb[2*PREVIEW_SIZE];
    double ub[2*PREVIEW_SIZE];

    int nW;
  

    // change the first period (all others are equal)

    printf ("\n################################\n %s \n################################\n", argv[0]);
    for(;;)
    {
        wmg.T[0] -= 0.01;
        if (abs(wmg.T[0]) < 1e-10)
        {
            wmg.T[0] = 0.1;
        }
        for (int i=0; i<PREVIEW_SIZE; i++)
        {
            cout << wmg.T[i] << "   ";
        }
        cout << endl;

        //------------------------------------------------------
        wmg.FormPreviewWindow();    
        if (wmg.halt)
        {
            cout << "EXIT (halt = 1)" << endl;
            break;
        }

        wmg.form_FP_init(); 

        for (int i = 0; i < PREVIEW_SIZE; i++)
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
        solver.init(wmg.T, wmg.h, angle, zref_x, zref_y, lb, ub, wmg.FP_init);
        nW = solver.solve();
//**************************************************************************


        //------------------------------------------------------
        // compare with reference results
        printf("Num. of activated constraints: %d\n", nW);
        for (int i = 0; i < 6; i++)
        {
            wmg.X[i] = wmg.FP_init[i];
            printf("value: % 8e\n", wmg.FP_init[i]);
        }
        //------------------------------------------------------

        wmg.slide();
    }
    printf ("################################\n");

    return 0;
}
