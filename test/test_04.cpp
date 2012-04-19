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
    string test_name = "";

    //-----------------------------------------------------------
    // initialize
    //-----------------------------------------------------------


    if ((argc == 2) && (strcmp (argv[1], "stdout") == 0))
    {
        dump_to_stdout = true;
        cout.precision (numeric_limits<double>::digits10);
    }
    else
    {
        test_name = "test_04";
    }


    init_03 test_04(test_name);



    fs_out.open(test_04.fs_out_filename.c_str(), fstream::app);
    fs_out.precision (numeric_limits<double>::digits10);
    fs_out << endl << endl;
    fs_out << "CoM_ZMP = [";


    smpc::solver_as solver(test_04.wmg->N);


    for(;;)
    {
        //------------------------------------------------------
        if (test_04.wmg->formPreviewWindow(*test_04.par) == WMG_HALT)
        {
            if (dump_to_stdout == false)
            {
                cout << "EXIT (halt = 1)" << endl;
            }
            break;
        }
        //------------------------------------------------------


        //------------------------------------------------------
        solver.set_parameters (test_04.par->T, test_04.par->h, test_04.par->h0, test_04.par->angle, test_04.par->fp_x, test_04.par->fp_y, test_04.par->lb, test_04.par->ub);
        solver.form_init_fp (test_04.par->fp_x, test_04.par->fp_y, test_04.par->init_state, test_04.par->X);
        solver.solve();
        solver.get_next_state(test_04.par->init_state);
        //------------------------------------------------------

        solver.get_next_state(test_04.X_tilde);
        fs_out << endl << test_04.par->init_state.x() << " " << test_04.par->init_state.y() << " " << test_04.X_tilde.x() << " " << test_04.X_tilde.y() << ";";
    
        if (dump_to_stdout)
        {
            for (unsigned int i = 0; i < test_04.wmg->N*SMPC_NUM_VAR; i++)
            {
                cout << test_04.par->X[i] << endl;
            }
        }
    }

    fs_out << "];" << endl;
    fs_out << "plot (CoM_ZMP(:,1), CoM_ZMP(:,2), 'b');" << endl;
    fs_out << "plot (CoM_ZMP(:,3), CoM_ZMP(:,4), 'ks','MarkerSize',5);" << endl;
    fs_out.close();


    return 0;
}
///@}
