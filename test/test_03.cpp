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
    WMG wmg;
    smpc_parameters par;
    init_01 (&wmg);
    par.init(wmg.N, wmg.hCoM/wmg.gravity);
    //-----------------------------------------------------------
  

    test_start (argv[0]);

    smpc::solver solver(wmg.N);
    int nW;

    unsigned int j=0;
    for(;;)
    {
        wmg.T_ms[(wmg.N-1) - j] = 50;
        if (j != 0)
        {
            wmg.T_ms[(wmg.N-1) - j + 1] = 100;
            if (j == wmg.N-1)
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
            wmg.T_ms[j] = 100;
            j++;
        }
        for (unsigned int i=0; i < wmg.N; i++)
        {
            cout << wmg.T_ms[i] << "   ";
        }
        cout << endl;

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
        nW = solver.solve();
        par.init_state.get_next_state (solver);
        //------------------------------------------------------


        //------------------------------------------------------
        printf("Num. of activated constraints: %d\n", nW);
        for (int i = 0; i < 6; i++)
        {
            printf("value: % 8e\n", par.X[i]);
        }
        //------------------------------------------------------
    }

    test_end (argv[0]);
    return 0;
}
///@}
