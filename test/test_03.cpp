/** 
 * @file
 * @author Alexander Sherikov
 * @brief Performs a full simulation with a dumb emulation of variable
 *  sampling time in the preview window: the differing period repeated
 *  in the preview window in a cycle.
 */


#include "tests_common.h" 

///@addtogroup gTEST
///@{
int main(int argc, char **argv)
{
    //-----------------------------------------------------------
    // initialize
    init_01 test_03 ("test_03");
    //-----------------------------------------------------------
  

    smpc::solver solver(test_03.wmg->N);
    int nW;

    unsigned int j=0;
    for(;;)
    {
        test_03.wmg->T_ms[(test_03.wmg->N-1) - j] = 50;
        if (j != 0)
        {
            test_03.wmg->T_ms[(test_03.wmg->N-1) - j + 1] = 100;
            if (j == test_03.wmg->N-1)
            {
                j = 0;
            }
            else
            {
                j++;
            }
        }
        else
        {
            test_03.wmg->T_ms[j] = 100;
            j++;
        }
        for (unsigned int i=0; i < test_03.wmg->N; i++)
        {
            cout << test_03.wmg->T_ms[i] << "   ";
        }
        cout << endl;

        //------------------------------------------------------
        if (test_03.wmg->formPreviewWindow(*test_03.par) == WMG_HALT)
        {
            cout << "EXIT (halt = 1)" << endl;
            break;
        }
        //------------------------------------------------------


        //------------------------------------------------------
        solver.set_parameters (test_03.par->T, test_03.par->h, test_03.par->h0, test_03.par->angle, test_03.par->fp_x, test_03.par->fp_y, test_03.par->lb, test_03.par->ub);
        solver.form_init_fp (test_03.par->fp_x, test_03.par->fp_y, test_03.par->init_state, test_03.par->X);
        nW = solver.solve();
        test_03.par->init_state.get_next_state (solver);
        //------------------------------------------------------


        //------------------------------------------------------
        printf("Num. of activated constraints: %d\n", nW);
        for (int i = 0; i < 6; i++)
        {
            printf("value: % 8e\n", test_03.par->X[i]);
        }
        //------------------------------------------------------
    }

    return 0;
}
///@}
