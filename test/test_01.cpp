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
    string test_name;

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
        // reference states generated using thr implementation of
        // the algorithm in Octave/MATLAB
        inFile.open ("./data/as_states_inv_downdate.dat");
        //inFile.open ("./data/states_chol_downdate.dat");
        test_name = "test_01";
    }

    init_01 test_01 (test_name);

    smpc::solver_as solver(test_01.wmg->N);
    smpc::enable_fexceptions();


    double err = 0;
    double max_err = 0;
    double max_err_first_state = 0;

    fs_out.open(test_01.fs_out_filename.c_str(), fstream::app);
    fs_out.precision (numeric_limits<double>::digits10);
    fs_out << endl << endl;
    fs_out << "CoM_ZMP = [";
   


    for(;;)
    {
        //------------------------------------------------------
        if (test_01.wmg->formPreviewWindow(*test_01.par) == WMG_HALT)
        {
            if (dump_to_stdout == false)
            {
                cout << "EXIT (halt = 1)" << endl;
            }
            break;
        }
        //------------------------------------------------------


        //------------------------------------------------------
        solver.set_parameters (test_01.par->T, test_01.par->h, test_01.par->h0, test_01.par->angle, test_01.par->fp_x, test_01.par->fp_y, test_01.par->lb, test_01.par->ub);
        solver.form_init_fp (test_01.par->fp_x, test_01.par->fp_y, test_01.par->init_state, test_01.par->X);
        solver.solve();
        test_01.par->init_state.get_next_state (solver);
        //------------------------------------------------------

        test_01.X_tilde.get_next_state (solver);
        fs_out << endl << test_01.par->init_state.x() << " " << test_01.par->init_state.y() << " " << test_01.X_tilde.x() << " " << test_01.X_tilde.y() << ";";

        if (dump_to_stdout)
        {
            for (unsigned int i = 0; i < test_01.wmg->N*SMPC_NUM_VAR; i++)
            {
                cout << test_01.par->X[i] << endl;
            }
        }
        else
        {
            //------------------------------------------------------
            // compare with reference results
            for (unsigned int i = 0; i < test_01.wmg->N*SMPC_NUM_VAR; i++)
            {
                double dataref;

                inFile >> dataref;
                if (i%3 == 0)
                {
                    // Test reference data was generated using a variant of the solver 
                    // with two variable substitutions, and consequently the coordinates
                    // of ZMP cannot be compared with the current implementation, which
                    // makes only one variable substitution.
                    printf("value: % 8e   ref: % 8e   err: -skip-\n", test_01.par->X[i], dataref);
                }
                else
                {
                    err = abs(test_01.par->X[i] - dataref);
                    if ((i < 6) && (err > max_err_first_state))
                    {
                        max_err_first_state = err;
                    }
                    if (err > max_err)
                    {
                        max_err = err;
                    }
                    printf("value: % 8e   ref: % 8e   err: % 8e\n", test_01.par->X[i], dataref, err);
                }
            }
            cout << "Max. error (first state, all preview windows): " << max_err_first_state << endl;
            cout << "Max. error (all states, all preview windows): " << max_err << endl;
            //------------------------------------------------------
        }
    }
    inFile.close();

    fs_out << "];" << endl;
    fs_out << "plot (CoM_ZMP(:,1), CoM_ZMP(:,2), 'b');" << endl;
    fs_out << "plot (CoM_ZMP(:,3), CoM_ZMP(:,4), 'ks','MarkerSize',5);" << endl;
    fs_out.close();

    return 0;
}
///@}
