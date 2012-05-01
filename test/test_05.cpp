/** 
 * @file
 * @author Alexander Sherikov
 * @brief Performs a full simulation using interior point method and compares 
 *  results with reference data produced by Octave.
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
        // reference states generated using thr implementation of
        // the algorithm in Octave/MATLAB
        inFile.open ("./data/ip_states_inv.dat");
        //inFile.open ("./data/states_chol_downdate.dat");
        test_name = "test_05";
    }
    init_01 test_05 (test_name);



    smpc::solver_ip solver(
            test_05.wmg->N, 
            2000, 
            150, 
            0.02, 
            1, 
            1e-3,
            1e-2,
            100, 
            15, 
            0.01, 
            0.5,
            0,
            true,
            true);


    double err = 0;
    double max_err = 0;
    double max_err_first_state = 0;


    fs_out.open(test_05.fs_out_filename.c_str(), fstream::app);
    fs_out.precision (numeric_limits<double>::digits10);
    fs_out << endl << endl;
    fs_out << "CoM_ZMP = [";
   


    for(;;)
    {
        //------------------------------------------------------
        if (test_05.wmg->formPreviewWindow(*test_05.par) == WMG_HALT)
        {
            if (dump_to_stdout == false)
            {
                cout << "EXIT (halt = 1)" << endl;
            }
            break;
        }
        //------------------------------------------------------

        //------------------------------------------------------
        solver.set_parameters (test_05.par->T, test_05.par->h, test_05.par->h0, test_05.par->angle, test_05.par->fp_x, test_05.par->fp_y, test_05.par->lb, test_05.par->ub);
        solver.form_init_fp (test_05.par->fp_x, test_05.par->fp_y, test_05.par->init_state, test_05.par->X);
        solver.solve();
        solver.get_next_state(test_05.par->init_state);
        //------------------------------------------------------

        solver.get_next_state(test_05.X_tilde);
        fs_out << endl << test_05.par->init_state.x() << " " << test_05.par->init_state.y() << " " << test_05.X_tilde.x() << " " << test_05.X_tilde.y() << ";";


        if (dump_to_stdout)
        {
            for (unsigned int i = 0; i < test_05.wmg->N*SMPC_NUM_VAR; i++)
            {
                cout << test_05.par->X[i] << endl;
            }
        }
        else
        {
            printf ("-------------------------------------\nOBJ: ");
            for (unsigned int i = 0; i < solver.objective_log.size(); ++i)
            {   
                printf ("% 8e ", solver.objective_log[i]);
            }
            printf ("\n-------------------------------------\n");


            //------------------------------------------------------
            // compare with reference results
            for (unsigned int i = 0; i < test_05.wmg->N*SMPC_NUM_VAR; i++)
            {
                double dataref;

                inFile >> dataref;
                err = abs(test_05.par->X[i] - dataref);
                if ((i < 6) && (err > max_err_first_state))
                {
                    max_err_first_state = err;
                }
                if (err > max_err)
                {
                    max_err = err;
                }
                printf("value: % 8e   ref: % 8e   err: % 8e\n", test_05.par->X[i], dataref, err);
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

    return 0;
}
///@}
