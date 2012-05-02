/** 
 * @file
 * @author Alexander Sherikov
 * @brief Comparison of IP and AS.
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
    const int N = 40;


    smpc::solver_ip IP_solver(
            N,                       // N
            8000,                    // gain_position
            1.0,                     // gain_velocity
            0.02,                    // gain_acceleration
            1.0,                     // gain_jerk
            1e-1,                    // tol
            1e+4,                    // tol_out | if it is big only one iteration is made.
            1e-1,                    // t
            1.0,                     // mu
            0.01,                    // bs_alpha
            0.9,                     // bs_beta
            5,                       // max_iter
            smpc::SMPC_IP_BS_LOGBAR, // backtracking search
            true);                   // obj

    smpc::solver_as AS_solver(
            N,              // size of the preview window
            8000.0,         // gain_position
            1.0,            // gain_velocity
            0.02,           // gain_acceleration
            1.0,            // gain_jerk
            1e-7,           // tolerance
            15,             // limit number of activated constraints
            true,           // remove constraints
            true);          // obj



    for(int stop_iter = 1; stop_iter < 200; ++stop_iter)
    {
        // Circular
        /*
        init_11 AS_test("test_14_as");
        init_11 IP_test("test_14_ip");
        */
        // Straight
        init_10 AS_test("");
        init_10 IP_test("");
        
        //-----------------------------------------------------------

        AS_test.par->init_state.set (0.019978839010709938, -6.490507362468014e-05);
        IP_test.par->init_state.set (0.019978839010709938, -6.490507362468014e-05);
        AS_test.X_tilde.set (0.019978839010709938, -6.490507362468014e-05);
        IP_test.X_tilde.set (0.019978839010709938, -6.490507362468014e-05);

        AS_test.wmg->T_ms[0] = control_sampling_time_ms;
        AS_test.wmg->T_ms[1] = control_sampling_time_ms;
        IP_test.wmg->T_ms[0] = control_sampling_time_ms;
        IP_test.wmg->T_ms[1] = control_sampling_time_ms;


        double AS_obj = 0.0;
        double IP_obj = 0.0;

        for(int counter = 0; /*counter < stop_iter*/; ++counter)
        {
            //------------------------------------------------------
            if (IP_test.wmg->formPreviewWindow(*IP_test.par) == WMG_HALT)
            {
                break;
            }
            if (AS_test.wmg->formPreviewWindow(*AS_test.par) == WMG_HALT)
            {
                break;
            }
            //------------------------------------------------------

            if (counter == stop_iter - 1)
            {
                IP_test.par->init_state.x() += 0.01;
                IP_test.par->init_state.y() -= 0.03;

                AS_test.par->init_state.x() += 0.01;
                AS_test.par->init_state.y() -= 0.03;
                /*
                IP_test.par->init_state.x() -= 0.05;
                IP_test.par->init_state.y() += 0.01;

                AS_test.par->init_state.x() -= 0.05;
                AS_test.par->init_state.y() += 0.01;
                */
            }
    //        IP_test.par->init_state =  AS_test.par->init_state;

            IP_solver.set_parameters (IP_test.par->T, IP_test.par->h, IP_test.par->h0, IP_test.par->angle, IP_test.par->zref_x, IP_test.par->zref_y, IP_test.par->lb, IP_test.par->ub);
            IP_solver.form_init_fp (IP_test.par->fp_x, IP_test.par->fp_y, IP_test.par->init_state, IP_test.par->X);
            IP_solver.solve();
            IP_solver.get_next_state(IP_test.par->init_state);
            IP_obj = IP_solver.objective_log.back();


            AS_solver.set_parameters (AS_test.par->T, AS_test.par->h, AS_test.par->h0, AS_test.par->angle, AS_test.par->zref_x, AS_test.par->zref_y, AS_test.par->lb, AS_test.par->ub);
            AS_solver.form_init_fp (AS_test.par->fp_x, AS_test.par->fp_y, AS_test.par->init_state, AS_test.par->X);
            AS_solver.solve();
            AS_solver.get_next_state(AS_test.par->init_state);
            AS_obj = AS_solver.objective_log.back();
        }
        printf ("(%3i)    AS: % 8e     IP: % 8e\n", stop_iter, AS_obj, IP_obj);
    }

    return 0;
}
///@}
