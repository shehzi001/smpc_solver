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
    ofstream fs_out;

    //-----------------------------------------------------------
    // initialize
    WMG wmg;
    smpc_parameters par;
    init_03 (&wmg);
    par.init(wmg.N, wmg.hCoM/wmg.gravity);

    string fs_out_filename("test_04_fs.m");
    wmg.FS2file(fs_out_filename); // output results for later use in Matlab/Octave
    //-----------------------------------------------------------


    if ((argc == 2) && (strcmp (argv[1], "stdout") == 0))
    {
        dump_to_stdout = true;
        cout.precision (numeric_limits<double>::digits10);
    }


    if (!dump_to_stdout)
        test_start(argv[0]);


    fs_out.open(fs_out_filename.c_str(), fstream::app);
    fs_out.precision (numeric_limits<double>::digits10);
    fs_out << endl << endl;
    fs_out << "CoM_ZMP = [";


    smpc::solver solver(wmg.N);


    for(;;)
    {
        //------------------------------------------------------
        if (wmg.formPreviewWindow(par) == WMG_HALT)
        {
            cout << "EXIT (halt = 1)" << endl;
            break;
        }
        //------------------------------------------------------


        //------------------------------------------------------
        solver.set_parameters (par.T, par.h, par.h0, par.angle, par.fp_x, par.fp_y, par.lb, par.ub);
        solver.form_init_fp (par.fp_x, par.fp_y, par.init_state, par.X);
        solver.solve();
        par.init_state.get_next_state (solver);
        //------------------------------------------------------

        wmg.X_tilde.get_next_state (solver);
        fs_out << endl << par.init_state.x() << " " << par.init_state.y() << " " << wmg.X_tilde.x() << " " << wmg.X_tilde.y() << ";";
    
        if (dump_to_stdout)
        {
            for (unsigned int i = 0; i < wmg.N*SMPC_NUM_VAR; i++)
            {
                cout << par.X[i] << endl;
            }
        }
    }

    fs_out << "];" << endl;
    fs_out << "plot (CoM_ZMP(:,1), CoM_ZMP(:,2), 'b');" << endl;
    fs_out << "plot (CoM_ZMP(:,3), CoM_ZMP(:,4), 'ks','MarkerSize',5);" << endl;
    fs_out.close();


    if (!dump_to_stdout)
        test_end(argv[0]);

    return 0;
}
///@}
