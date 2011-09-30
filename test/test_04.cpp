/** 
 * @file
 * @author Alexander Sherikov
 * @brief Simulation with double support
 */


#include <string>
#include <cstring> //strcmp
#include "tests_common.h"

///@addtogroup gTEST
///@{

int main(int argc, char **argv)
{
    bool dump_to_stdout = false;

    //-----------------------------------------------------------
    // initialize
    WMG wmg;
    init_03 (&wmg);

    string out_file("test_04_fs.m");
    wmg.FS2file(out_file); // output results for later use in Matlab/Octave
    //-----------------------------------------------------------


    if ((argc == 2) && (strcmp (argv[1], "stdout") == 0))
    {
        dump_to_stdout = true;
        cout.precision (numeric_limits<double>::digits10);
    }


    if (!dump_to_stdout)
        test_start(argv[0]);


    smpc_solver solver(wmg.N);


    for(;;)
    {
        //------------------------------------------------------
        if (wmg.FormPreviewWindow() == WMG_HALT)
        {
            cout << "EXIT (halt = 1)" << endl;
            break;
        }
        //------------------------------------------------------


        //------------------------------------------------------
        solver.init(wmg.T, wmg.h, wmg.angle, wmg.zref_x, wmg.zref_y, wmg.lb, wmg.ub, wmg.X_tilde, wmg.X);
        solver.solve();
        solver.get_next_state_tilde (wmg.X_tilde);
        //------------------------------------------------------


        if (dump_to_stdout)
        {
            for (int i = 0; i < wmg.N*NUM_VAR; i++)
            {
                cout << wmg.X[i] << endl;
            }
        }
    }


    if (!dump_to_stdout)
        test_end(argv[0]);

    return 0;
}
///@}
