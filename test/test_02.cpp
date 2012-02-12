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
    smpc_parameters par;
    init_02 (&wmg);
    par.init(wmg.N, wmg.hCoM/wmg.gravity);
    //-----------------------------------------------------------
 

    test_start(argv[0]);

    smpc::solver solver(wmg.N);
    for(int counter = 0 ;; counter++)
    {
        //------------------------------------------------------
        if (wmg.formPreviewWindow(par) == WMG_HALT)
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
            solver.set_parameters (par.T, par.h, par.h0, par.angle, par.fp_x, par.fp_y, par.lb, par.ub);
            solver.form_init_fp (par.fp_x, par.fp_y, par.init_state, par.X);
            nW = solver.solve();
        }
        par.init_state.get_next_state (solver);
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
