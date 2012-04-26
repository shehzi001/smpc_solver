/** 
 * @file
 * @author Alexander Sherikov
 * @brief Performs a full simulation and measures average time
 *  required to solve a QP using IP method.
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

    init_02 test_11("test_11");
    //-----------------------------------------------------------

    smpc::solver_ip solver(
            test_11.wmg->N,  // N
            10,             // max_iter
            2000,            // gain_position
            150,             // gain_velocity
            0.02,            // gain_acceleration
            1,               // gain_jerk
            1e-2,            // tol
            1e-2,            // tol_out
            1e-4,              // t
            1.0,              // mu
            0.01,            // bs_alpha
            0.9);            // bs_beta          
    
    for(int counter = 0; counter < 2; counter++)
    {
        //------------------------------------------------------
        if (test_11.wmg->formPreviewWindow(*test_11.par) == WMG_HALT)
        {
            cout << "EXIT (halt = 1)" << endl;
            break;
        }
        //------------------------------------------------------



        int NN = 1000;
        gettimeofday(&start,0);
        for(int kk=0; kk<NN ;kk++)
        {
            solver.set_parameters (test_11.par->T, test_11.par->h, test_11.par->h0, test_11.par->angle, test_11.par->fp_x, test_11.par->fp_y, test_11.par->lb, test_11.par->ub);
            solver.form_init_fp (test_11.par->fp_x, test_11.par->fp_y, test_11.par->init_state, test_11.par->X);
            solver.solve();
        }
        solver.get_next_state(test_11.par->init_state);
        gettimeofday(&end,0);             
        CurrentCPUTime = end.tv_sec - start.tv_sec + 0.000001 * (end.tv_usec - start.tv_usec);
        double TT = CurrentCPUTime/NN;
        printf("(%i) time = % f (ext = %d // int = %d // bs = %d)\n", 
                counter, 
                TT, 
                solver.ext_loop_iterations,
                solver.int_loop_iterations,
                solver.bt_search_iterations);
        //------------------------------------------------------
    }

    return 0;
}
///@}
