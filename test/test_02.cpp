/** 
 * @file
 * @author Alexander Sherikov
 * @brief Performs a full simulation and measures average time
 *  required to solve a QP.
 */


#include <sys/time.h>
#include <time.h>

#include "tests_common.h" 

///@addtogroup gTEST
///@{

int main(int argc, char **argv)
{
    struct timeval start, end;
    double CurrentCPUTime;

    WMG wmg;
    init_02 (&wmg);
    //-----------------------------------------------------------
 

    test_start(argv[0]);

    smpc::solver solver(wmg.N);
    for(int counter = 0 ;; counter++)
    {
        //------------------------------------------------------
        if (wmg.formPreviewWindow() == WMG_HALT)
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
            solver.set_parameters (wmg.T, wmg.h, wmg.h[0], wmg.angle, wmg.fp_x, wmg.fp_y, wmg.lb, wmg.ub);
            solver.form_init_fp (wmg.fp_x, wmg.fp_y, wmg.init_state, wmg.X);
            nW = solver.solve();
        }
        wmg.init_state.get_next_state (solver);
        gettimeofday(&end,0);             
        CurrentCPUTime = end.tv_sec - start.tv_sec + 0.000001 * (end.tv_usec - start.tv_usec);
        double TT = CurrentCPUTime/NN;
        printf("(%i) time = % f (%i)\n", counter, TT, nW);
        //------------------------------------------------------
    }


    test_end(argv[0]);
    return 0;
}
///@}
