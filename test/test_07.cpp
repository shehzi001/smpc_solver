/** 
 * @file
 * @author Alexander Sherikov
 */


#include <cstring> //strcmp
#include "tests_common.h"

///@addtogroup gTEST
///@{

int main(int argc, char **argv)
{
    //-----------------------------------------------------------
    // initialize
    WMG wmg;
    init_04 (&wmg);
    //-----------------------------------------------------------


    test_start(argv[0]);
    smpc_solver solver(wmg.N);


    for(;;)
    {
        int num_iter_in_ss;
        int num_iter_in_ss_passed;
        double prev_swing_pos[3];
        double next_swing_pos[3];



        //------------------------------------------------------
        if (wmg.FormPreviewWindow() == WMG_HALT)
        {
            cout << "EXIT (halt = 1)" << endl;
            break;
        }
        //------------------------------------------------------

        wmg.get_swing_foot_pos (
                prev_swing_pos,
                next_swing_pos,
                &num_iter_in_ss,
                &num_iter_in_ss_passed);

        //------------------------------------------------------
        solver.set_parameters (wmg.T, wmg.h, wmg.angle, wmg.zref_x, wmg.zref_y, wmg.lb, wmg.ub);
        solver.form_init_fp (wmg.zref_x, wmg.zref_y, wmg.X_tilde, wmg.X);
        solver.solve();
        solver.get_next_state_tilde (wmg.X_tilde);
        //------------------------------------------------------
    }
    test_end(argv[0]);

    return 0;
}
///@}
