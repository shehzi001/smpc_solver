#include <iostream>
#include <fstream>

#include <cmath> // abs
#include <cstdio>

#include <sys/time.h>
#include <time.h>


#include "WMG.h"
#include "init_WMG.h"

#include "smpc_solver.h" 


using namespace std;


int main(int argc, char **argv)
{
    struct timeval start, end;
    double CurrentCPUTime;

    WMG wmg;
    init_02 (&wmg);
    //-----------------------------------------------------------
 

    smpc_solver solver(wmg.N);


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
