/** 
 * @file
 * @author Alexander Sherikov
 */


#include <sys/time.h>
#include <time.h>

#include "tests_common.h" 

///@addtogroup gTEST
///@{

int main(int argc, char **argv)
{
    const int control_sampling_time_ms = 20;
//    const int preview_sampling_time_ms = 40;
//    const int next_preview_len_ms = 0;

    struct timeval start, end;
    double CurrentCPUTime;

    ofstream IP_log;


    IP_log.open("test_17_ip.m", fstream::trunc);
    IP_log.precision (numeric_limits<double>::digits10);

    const int stop_iter = 140;
    const int NN = 1;
    const int ip_limit_num = 7;
    const int ip_limit[ip_limit_num] = {0, 20, 10, 7, 6, 5, 4};

    for(int j = 0; j < ip_limit_num; ++j)
    {
        // Circular
        /*
        init_11 IP_test("");
        */
        // Straight
        init_10 IP_test("");
        
        //-----------------------------------------------------------

        IP_test.par->init_state.set (0.019978839010709938, -6.490507362468014e-05);
        IP_test.X_tilde.set (0.019978839010709938, -6.490507362468014e-05);

        IP_test.wmg->T_ms[0] = control_sampling_time_ms;
        IP_test.wmg->T_ms[1] = control_sampling_time_ms;


        IP_log << endl << endl;
        IP_log << "OBJ" << j << "= [";


        smpc::solver_ip IP_solver(
                IP_test.wmg->N,          // N
                8000,                    // gain_position
                1.0,                     // gain_velocity
                0.02,                    // gain_acceleration
                1.0,                     // gain_jerk
                1e-2,                    // tol
                1e+4,                    // tol_out | if it is big only one iteration is made.
                1.0,                     // t
                1.0,                     // mu
                0.01,                    // bs_alpha
                0.95,                     // bs_beta
                ip_limit[j],              // max_iter
                smpc::SMPC_IP_BS_LOGBAR, // backtracking search
                true);                   // obj


        for(int counter = 0; /*counter < stop_iter*/; ++counter)
        {
            //------------------------------------------------------
            if (IP_test.wmg->formPreviewWindow(*IP_test.par) == WMG_HALT)
            {
                cout << "EXIT (halt = 1)" << endl;
                break;
            }
            //------------------------------------------------------

            if (counter == stop_iter - 1)
            {
                IP_test.par->init_state.x() += 0.01;
                IP_test.par->init_state.y() -= 0.03;
            }

            gettimeofday(&start,0);
            for(int kk=0; kk<NN ;kk++)
            {
                IP_solver.set_parameters (IP_test.par->T, IP_test.par->h, IP_test.par->h0, IP_test.par->angle, IP_test.par->zref_x, IP_test.par->zref_y, IP_test.par->lb, IP_test.par->ub);
                IP_solver.form_init_fp (IP_test.par->fp_x, IP_test.par->fp_y, IP_test.par->init_state, IP_test.par->X);
                IP_solver.solve();
            }
            IP_solver.get_next_state(IP_test.par->init_state);
            gettimeofday(&end,0);             
            CurrentCPUTime = end.tv_sec - start.tv_sec + 0.000001 * (end.tv_usec - start.tv_usec);
            printf("(%3i)  IP: time = % f (ext = %d // int = %d // bs = %d)\n", 
                    counter, 
                    CurrentCPUTime/NN, 
                    IP_solver.ext_loop_iterations,
                    IP_solver.int_loop_iterations,
                    IP_solver.bt_search_iterations);
            printf ("           OBJ = ");
            for (unsigned int i = 0; i < IP_solver.objective_log.size(); ++i)
            {
                printf ("% 8e ", IP_solver.objective_log[i]);
            }
            printf ("\n");
            IP_solver.get_next_state(IP_test.X_tilde);
            IP_log << endl << IP_solver.objective_log.back() << ";";
        }

        IP_log << "];" << endl;
        IP_log << "figure;" << "semilogy (OBJ" << j << "(:,1));" << endl;
        IP_log << "legend ('Limit = " << ip_limit[j] << "');"<< endl;
    }
    IP_log.close();

    return 0;
}
///@}
