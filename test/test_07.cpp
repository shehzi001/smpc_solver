#include <iostream>
#include <fstream>

#include <cmath> // abs
#include <cstdio>


#include "WMG.h"
#include "init_WMG.h"

#include "smpc_solver.h" 


using namespace std;


int main(int argc, char **argv)
{
    //-----------------------------------------------------------
    // initialize
    WMG wmg;
    init_01(&wmg);
    //wmg.FS2file(); // output results for later use in Matlab/Octave
    //-----------------------------------------------------------
  


    smpc_solver solver(wmg.N);

    double err = 0;
    double max_err = 0;
    double max_err_first_state = 0;
   
    // reference states generated using thr implementation of
    // the algorithm in Octave/MATLAB
    ifstream inFile;
    inFile.open ("./data/states_inv_downdate.dat");


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
            if ((i < 6) && (err > max_err_first_state))
            {
                max_err_first_state = err;
            }
            if (err > max_err)
            {
                max_err = err;
            }
            //printf("value: % 8e   ref: % 8e   err: % 8e\n", wmg.FP_init[i], dataref, err);
        }
        cout << "Max. error (first state, all steps): " << max_err_first_state << endl;
        cout << "Max. error (all states, all steps): " << max_err << endl;
        //------------------------------------------------------

        wmg.slide();
    }
    inFile.close();
    printf ("################################\n");

    return 0;
}
