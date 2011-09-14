#include <iostream>
#include <fstream>

#include <cmath> // abs
#include <cstdio>

#include <sys/time.h>
#include <time.h>


#include "WMG.h"

#include "smpc_solver.h" 

#define PREVIEW_SIZE 15 // Size of the preview window

using namespace std;


int main(int argc, char **argv)
{
    struct timeval start, end;
    double CurrentCPUTime;

    WMG wmg(PREVIEW_SIZE, 0.1, 0.261);

    double d[4] = {0.09 , 0.025, 0.03, 0.075};
    wmg.AddFootstep(0.0, 0.05, 0.0, 3, 3, 1, d);
    
    double z = 5.0*M_PI/180.0;
    double step_x = 0.035;
    double step_y = 0.1;

    d[3] = 0.025;

    // use this for smaller feet (and we will have more active constraints)
    d[0] = 0.03; d[1] = 0.01; d[2] = 0.01; d[3] = 0.01;

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



        int NN = 1000;
        int nW;
        gettimeofday(&start,0);
        for(int kk=0; kk<NN ;kk++)
        {
//**************************************************************************
// SOLVER IS USED HERE
//**************************************************************************
            solver.init(wmg.T, wmg.h, wmg.angle, wmg.zref_x, wmg.zref_y, wmg.lb, wmg.ub, wmg.X_tilde, wmg.X);
            nW = solver.solve();
//**************************************************************************
        }
        solver.get_next_state_tilde (wmg.X_tilde);
        gettimeofday(&end,0);             
        CurrentCPUTime = end.tv_sec - start.tv_sec + 0.000001 * (end.tv_usec - start.tv_usec);
        double TT = CurrentCPUTime/NN;
        printf("(%i) time = % f (%i)\n", wmg.counter, TT, nW);
        //------------------------------------------------------


        wmg.slide();
    }
    printf ("################################\n");

    return 0;
}
