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
        double swing_foot_position[3];
        double angle;



        //------------------------------------------------------
        if (wmg.FormPreviewWindow() == WMG_HALT)
        {
            cout << "EXIT (halt = 1)" << endl;
            break;
        }
        //------------------------------------------------------
        wmg.getSwingFootPosition (
                WMG_SWING_PARABOLA,
                10,
                0,
                swing_foot_position,
                &angle);

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
