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

    //-----------------------------------------------------------
    // initialize
    WMG wmg;
    init_01 (&wmg);
    //wmg.FS2file(); // output results for later use in Matlab/Octave
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
        inFile.open ("./data/states_inv_downdate.dat");
        //inFile.open ("./data/states_chol_downdate.dat");
        //inFile.open ("./data/states_chol_nodowndate.dat");
    }


    if (!dump_to_stdout)
        test_start(argv[0]);


    smpc_solver solver(wmg.N);

    double err = 0;
    double max_err = 0;
    double max_err_first_state = 0;
   


    for(;;)
    {
        //------------------------------------------------------
        wmg.FormPreviewWindow();    
        if (wmg.halt)
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
        else
        {
            //------------------------------------------------------
            // compare with reference results
            for (int i = 0; i < wmg.N*NUM_VAR; i++)
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
                //printf("value: % 8e   ref: % 8e   err: % 8e\n", wmg.FP_init[i], dataref, err);
            }
            cout << "Max. error (first state, all steps): " << max_err_first_state << endl;
            cout << "Max. error (all states, all steps): " << max_err << endl;
            //------------------------------------------------------
        }
    }
    inFile.close();


    if (!dump_to_stdout)
        test_end(argv[0]);

    return 0;
}
///@}
