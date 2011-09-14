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

    //wmg.FS2file(); // output results for later use in Matlab/Octave
    //-----------------------------------------------------------
  


    smpc_solver solver(PREVIEW_SIZE);

    double err = 0;
    double max_err = 0;
   
  
    // reference states generated using thr implementation of
    // the algorithm in Octave/MATLAB
    ifstream inFile;
    inFile.open ("./data/states_chol_downdate.dat");


    printf ("\n################################\n %s \n################################\n", argv[0]);


    for(;;)
    {
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
        solver.init(wmg.T, wmg.h, wmg.angle, wmg.zref_x, wmg.zref_y, wmg.lb, wmg.ub, wmg.X_tilde, wmg.X);
        solver.solve();
        solver.get_next_state_tilde (wmg.X_tilde);
//**************************************************************************


        //------------------------------------------------------
        // compare with reference results
        for (int i = 0; i < wmg.N*NUM_VAR; i++)
        {
            double dataref;

            inFile >> dataref;
            err = abs(wmg.X[i] - dataref);
            if (err > max_err)
            {
                max_err = err;
            }
            //printf("value: % 8e   ref: % 8e   err: % 8e\n", wmg.FP_init[i], dataref, err);
        }
        cout << "Max. error (over all steps): " << max_err << endl;
        //------------------------------------------------------

        wmg.slide();
    }
    inFile.close();
    printf ("################################\n");

    return 0;
}
