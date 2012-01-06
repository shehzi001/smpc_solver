/** 
 * @file
 * @author Alexander Sherikov
 * @brief Performs a full simulation and compares results with 
 *  reference data produced by Octave.
 */


#include <cstring> //strcmp
#include "tests_common.h"

///@addtogroup gTEST
///@{

int main(int argc, char **argv)
{
    bool dump_to_stdout = false;
    ifstream inFile;
    ofstream fs_out;

    //-----------------------------------------------------------
    // initialize
    WMG wmg;
    init_01 (&wmg);

    std::string fs_out_filename("test_01_fs.m");
    wmg.FS2file(fs_out_filename); // output results for later use in Matlab/Octave
    //-----------------------------------------------------------


    if ((argc == 2) && (strcmp (argv[1], "stdout") == 0))
    {
        dump_to_stdout = true;
        cout.precision (numeric_limits<double>::digits10);
    }
    else
    {
        // reference states generated using thr implementation of
        // the algorithm in Octave/MATLAB
        inFile.open ("./data/as_states_inv_downdate.dat");
        //inFile.open ("./data/states_chol_downdate.dat");
    }


    if (!dump_to_stdout)
        test_start(argv[0]);


    smpc_solver solver(wmg.N);

    double err = 0;
    double max_err = 0;
    double max_err_first_state = 0;

    fs_out.open(fs_out_filename.c_str(), fstream::app);
    fs_out.precision (numeric_limits<double>::digits10);
    fs_out << endl << endl;
    fs_out << "CoM_ZMP = [";
   


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
        solver.set_parameters (wmg.T, wmg.h, wmg.h[0], wmg.angle, wmg.fp_x, wmg.fp_y, wmg.lb, wmg.ub);
        solver.form_init_fp (wmg.fp_x, wmg.fp_y, wmg.init_state, wmg.X);
        solver.solve();
        solver.get_next_state (wmg.init_state);
        //------------------------------------------------------

        solver.get_next_state_tilde (wmg.X_tilde);
        fs_out << endl << wmg.init_state[0] << " " << wmg.init_state[3] << " " << wmg.X_tilde[0] << " " << wmg.X_tilde[3] << ";";

        if (dump_to_stdout)
        {
            for (int i = 0; i < wmg.N*SMPC_NUM_VAR; i++)
            {
                cout << wmg.X[i] << endl;
            }
        }
        else
        {
            //------------------------------------------------------
            // compare with reference results
            for (int i = 0; i < wmg.N*SMPC_NUM_VAR; i++)
            {
                double dataref;

                inFile >> dataref;
                err = abs(wmg.X[i] - dataref);
                if ((i < 6) && (err > max_err_first_state))
                {
                    max_err_first_state = err;
                }
                if (err > max_err)
                {
                    max_err = err;
                }
                //printf("value: % 8e   ref: % 8e   err: % 8e\n", wmg.X[i], dataref, err);
            }
            cout << "Max. error (first state, all steps): " << max_err_first_state << endl;
            cout << "Max. error (all states, all steps): " << max_err << endl;
            //------------------------------------------------------
        }
    }
    inFile.close();

    fs_out << "];" << endl;
    fs_out << "plot (CoM_ZMP(:,1), CoM_ZMP(:,2), 'b');" << endl;
    fs_out << "plot (CoM_ZMP(:,3), CoM_ZMP(:,4), 'ks','MarkerSize',5);" << endl;
    fs_out.close();

    if (!dump_to_stdout)
        test_end(argv[0]);

    return 0;
}
///@}
