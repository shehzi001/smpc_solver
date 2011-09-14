#include <iostream>
#include <fstream>

#include <cmath> // abs
#include <cstdio>


#include "WMG.h"

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
        //------------------------------------------------------


//**************************************************************************
// SOLVER IS USED HERE
//**************************************************************************
        solver.init(wmg.T, wmg.h, wmg.angle, wmg.zref_x, wmg.zref_y, wmg.lb, wmg.ub, wmg.FP_init);
        nW = solver.solve();
//**************************************************************************


        //------------------------------------------------------
        // compare with reference results
        printf("Num. of activated constraints: %d\n", nW);
        for (int i = 0; i < 6; i++)
        {
            printf("value: % 8e\n", wmg.FP_init[i]);
        }
        //------------------------------------------------------

        wmg.slide();
    }
    printf ("################################\n");

    return 0;
}
