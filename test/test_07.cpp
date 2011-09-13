#include <iostream>
#include <fstream>

#include <cmath> // abs


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

    //wmg.FS2file(); // output results for later use in Matlab/Octave
    //-----------------------------------------------------------
  


    smpc_solver solver(PREVIEW_SIZE);

    printf ("\n################################\n %s \n################################\n", argv[0]);


    double ZMP_x, ZMP_y, CoM_x, CoM_y;

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
        //------------------------------------------------------


//**************************************************************************
// SOLVER IS USED HERE
//**************************************************************************
        solver.init(wmg.T, wmg.h, wmg.angle, wmg.zref_x, wmg.zref_y, wmg.lb, wmg.ub, wmg.FP_init);
        solver.solve();
//**************************************************************************


        //------------------------------------------------------
        // compare with reference results
        for (int i = 0; i < 6; i++)
        {
            wmg.X[i] = wmg.FP_init[i];
        }
        //------------------------------------------------------


        wmg.CoM.x = wmg.X_tilde[0] + wmg.h[0]*(wmg.X_tilde[2]);
        wmg.CoM.y = wmg.X_tilde[3] + wmg.h[0]*(wmg.X_tilde[5]);

        wmg.ZMP.x = wmg.X_tilde[0];
        wmg.ZMP.y = wmg.X_tilde[3];

        if (wmg.counter != 0)   // get_ZMP_CoM returns coordinates of ZMP and CoM
        {                       // from the next simulation step
            printf ("ZMP and CoM coordinates check: % 8e\n", 
                    abs(wmg.CoM.x - CoM_x) + 
                    abs(wmg.CoM.y - CoM_y) +
                    abs(wmg.ZMP.x - ZMP_x) + 
                    abs(wmg.ZMP.y - ZMP_y));
            /*
            printf ("CoM coord. check: % 6e  % 6e | % 6e  % 6e\n",
                    wmg.CoM.x, CoM_x, wmg.CoM.y, CoM_y);
                  
            printf ("ZMP coord. check: % 6e  % 6e | % 6e  % 6e\n", 
                    wmg.ZMP.x, ZMP_x, wmg.ZMP.y, ZMP_y);
            */
        }
        solver.get_ZMP_CoM (&ZMP_x, &ZMP_y, &CoM_x, &CoM_y);

        wmg.slide();
    }
    printf ("################################\n");

    return 0;
}
