/** 
 * @file
 * @author Alexander Sherikov
 * @brief Performs a full simulation with a dumb emulation of variable
 *  sampling time in the preview window: the first period in the preview
 *  window is decreased by a small fraction step by step.
 *
 * @todo Suspiciosly many constraints are activated.
 */


#include "tests_common.h" 


///@addtogroup gTEST
///@{
int main(int argc, char **argv)
{
    //-----------------------------------------------------------
    // initialize
    WMG wmg;
    init_01 (&wmg);
    //-----------------------------------------------------------
  


    test_start (argv[0]);
    smpc_solver solver(wmg.N);

    int nW;
  

    // change the first period (all others are equal)
    for(;;)
    {
        wmg.T[0] -= 0.01;
        if (abs(wmg.T[0]) < 1e-10)
        {
            wmg.T[0] = 0.1;
        }
        for (int i=0; i<wmg.N; i++)
        {
            cout << wmg.T[i] << "   ";
        }
        cout << endl;

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
        nW = solver.solve();
        solver.get_next_state_tilde (wmg.X_tilde);
        //------------------------------------------------------


        //------------------------------------------------------
        printf("Num. of activated constraints: %d\n", nW);
        for (int i = 0; i < 6; i++)
        {
            printf("value: % 8e\n", wmg.X[i]);
        }
        //------------------------------------------------------
    }

    test_end (argv[0]);
    return 0;
}
///@}
