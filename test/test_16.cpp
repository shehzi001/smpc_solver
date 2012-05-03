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

    ofstream AS_log;


    AS_log.open("test_16_as.m", fstream::trunc);
    AS_log.precision (numeric_limits<double>::digits10);

    const int stop_iter = 59;
    const int NN = 1;
    const int as_limit_num = 9;
    const int as_limit[as_limit_num] = {0, 30, 25, 24, 23, 22, 21, 20, 15};

    for(int j = 0; j < as_limit_num; ++j)
    {
        // Circular
        /*
        init_11 AS_test("");
        */
        // Straight
        init_10 AS_test("");
        
        //-----------------------------------------------------------

        AS_test.par->init_state.set (0.019978839010709938, -6.490507362468014e-05);
        AS_test.X_tilde.set (0.019978839010709938, -6.490507362468014e-05);

        AS_test.wmg->T_ms[0] = control_sampling_time_ms;
        AS_test.wmg->T_ms[1] = control_sampling_time_ms;


        AS_log << endl << endl;
        AS_log << "OBJ" << j << "= [";


        smpc::solver_as AS_solver(
                AS_test.wmg->N, // size of the preview window
                8000.0,         // gain_position
                1.0,            // gain_velocity
                0.02,           // gain_acceleration
                1.0,            // gain_jerk
                1e-7,           // tolerance
                as_limit[j],    // limit number of activated constraints
                true,           // remove constraints
                true);          // obj

        for(int counter = 0; /*counter < stop_iter*/; ++counter)
        {
            //------------------------------------------------------
            if (AS_test.wmg->formPreviewWindow(*AS_test.par) == WMG_HALT)
            {
                cout << "EXIT (halt = 1)" << endl;
                break;
            }
            //------------------------------------------------------

            if (counter == stop_iter - 1)
            {
                AS_test.par->init_state.x() += 0.01;
                AS_test.par->init_state.y() -= 0.03;
            }

            gettimeofday(&start,0);
            for(int kk=0; kk<NN ;kk++)
            {
                AS_solver.set_parameters (AS_test.par->T, AS_test.par->h, AS_test.par->h0, AS_test.par->angle, AS_test.par->zref_x, AS_test.par->zref_y, AS_test.par->lb, AS_test.par->ub);
                AS_solver.form_init_fp (AS_test.par->fp_x, AS_test.par->fp_y, AS_test.par->init_state, AS_test.par->X);
                AS_solver.solve();
            }
            AS_solver.get_next_state(AS_test.par->init_state);
            gettimeofday(&end,0);             
            CurrentCPUTime = end.tv_sec - start.tv_sec + 0.000001 * (end.tv_usec - start.tv_usec);
            printf("       AS: time = % f (AS size = %i // Added = %i // Removed = %i)\n", 
                    CurrentCPUTime/NN, 
                    AS_solver.active_set_size,
                    AS_solver.added_constraints_num,
                    AS_solver.removed_constraints_num);
            printf("           OBJ = ");
            for (unsigned int i = 0; i < AS_solver.objective_log.size(); ++i)
            {
                printf ("% 8e ", AS_solver.objective_log[i]);
            }
            printf ("\n ---\n");
            AS_solver.get_next_state(AS_test.X_tilde);
            AS_log << endl << AS_solver.objective_log.back() << ";";
            //------------------------------------------------------
        }

        AS_log << "];" << endl;
        AS_log << "figure;" << "semilogy (OBJ" << j << "(:,1));" << endl;
        AS_log << "legend ('Limit = " << as_limit[j] << "');"<< endl;
    }
    AS_log.close();

    return 0;
}
///@}
