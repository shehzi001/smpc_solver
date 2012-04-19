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

    init_02 test_02("test_02");
    //-----------------------------------------------------------
 

    smpc::solver_as solver(test_02.wmg->N);
    for(int counter = 0 ;; counter++)
    {
        //------------------------------------------------------
        if (test_02.wmg->formPreviewWindow(*test_02.par) == WMG_HALT)
        {
            cout << "EXIT (halt = 1)" << endl;
            break;
        }
        //------------------------------------------------------



        int NN = 1000;
        gettimeofday(&start,0);
        for(int kk=0; kk<NN ;kk++)
        {
            solver.set_parameters (test_02.par->T, test_02.par->h, test_02.par->h0, test_02.par->angle, test_02.par->fp_x, test_02.par->fp_y, test_02.par->lb, test_02.par->ub);
            solver.form_init_fp (test_02.par->fp_x, test_02.par->fp_y, test_02.par->init_state, test_02.par->X);
            solver.solve();
        }
        solver.get_next_state(test_02.par->init_state);
        gettimeofday(&end,0);             
        CurrentCPUTime = end.tv_sec - start.tv_sec + 0.000001 * (end.tv_usec - start.tv_usec);
        double TT = CurrentCPUTime/NN;
        printf("(%i) time = % f (%i)\n", counter, TT, solver.active_set_size);
        //------------------------------------------------------
    }

    return 0;
}
///@}
