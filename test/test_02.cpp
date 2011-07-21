#include <iostream>
#include <fstream>

#include <cmath> // abs
#include <cstdio>

#include <sys/time.h>
#include <time.h>


#include "WMG.h"

// QPAS_VARIABLE_AB is taken from this header, the results
// must be the same
#include "smpc_common.h" 

#include "qp_as.h"


using namespace std;


int main()
{
    struct timeval start, end;
    double CurrentCPUTime;

    WMG wmg(15, 0.1, 0.261);

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
 

    qp_as solver(wmg.N);

    double angle[wmg.N];
    double zref_x[wmg.N];
    double zref_y[wmg.N];
    double lb[2*wmg.N];
    double ub[2*wmg.N];
  


    printf ("\n################################\n test 02\n################################\n");
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



        int NN = 1000;
        int nW;
        gettimeofday(&start,0);             
        for(int kk=0; kk<NN ;kk++)
        {
            wmg.form_FP_init(); // reset X
//**************************************************************************
// SOLVER IS USED HERE
//**************************************************************************
#ifdef QPAS_VARIABLE_AB
            solver.init(wmg.Tk, wmg.hk, angle, zref_x, zref_y, lb, ub, wmg.FP_init);
#else
            solver.init(wmg.T, wmg.h, angle, zref_x, zref_y, lb, ub, wmg.FP_init);
#endif
            nW = solver.solve();
//**************************************************************************
        }
        gettimeofday(&end,0);             
        CurrentCPUTime = end.tv_sec - start.tv_sec + 0.000001 * (end.tv_usec - start.tv_usec);
        double TT = CurrentCPUTime/NN;
        printf("(%i) time = % f (%i)\n", wmg.counter, TT, nW);
        //------------------------------------------------------


        for (int i = 0; i < 6; i++)
        {
            wmg.X[i] = wmg.FP_init[i];
        }

        wmg.slide();
    }
    printf ("################################\n");

    return 0;
}
