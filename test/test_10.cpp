/** 
 * @file
 * @author Alexander Sherikov
 * @brief Performs a full simulation and measures average time
 *  required to solve a QP.
 */


#include <sys/time.h>
#include <time.h>

#include "tests_common.h" 

///@addtogroup gTEST
///@{

int main(int argc, char **argv)
{
    init_02 test_10("test_10");
    //-----------------------------------------------------------
 

    // N periods
    ofstream fT("test_10__T.out", ios::out);
    // N CoM heights
    ofstream fh("test_10__h.out", ios::out);
    // N rotation angles
    ofstream fa("test_10__angle.out", ios::out);
    // N x coordinates of feasible points
    ofstream ffpx("test_10__FPx.out", ios::out);
    // N y coordinates of feasible points
    ofstream ffpy("test_10__FPy.out", ios::out);
    // N lower bounds for ZMP position
    ofstream flb("test_10__lower_bounds.out", ios::out);
    // N upper bounds for ZMP position
    ofstream fub("test_10__upper_bounds.out", ios::out);


    fT.precision (numeric_limits<double>::digits10);
    fh.precision (numeric_limits<double>::digits10);
    fa.precision (numeric_limits<double>::digits10);
    ffpx.precision (numeric_limits<double>::digits10);
    ffpy.precision (numeric_limits<double>::digits10);
    flb.precision (numeric_limits<double>::digits10);
    fub.precision (numeric_limits<double>::digits10);


    smpc::solver_as solver(test_10.wmg->N);
    for(int counter = 0 ;; counter++)
    {
        //------------------------------------------------------
        if (test_10.wmg->formPreviewWindow(*test_10.par) == WMG_HALT)
        {
            cout << "EXIT (halt = 1)" << endl;
            break;
        }
        //------------------------------------------------------

        for (unsigned int i = 0; i < test_10.wmg->N; i++)
        {
            fT << test_10.par->T[i] << endl;
            fh << test_10.par->h[i] << endl;
            // all h are equal
            //test_10.par->h0;
            fa << test_10.par->angle[i] << endl;
            ffpx << test_10.par->fp_x[i] << endl;
            ffpy << test_10.par->fp_y[i] << endl;
            flb << test_10.par->lb[i] << endl;
            fub << test_10.par->ub[i] << endl;
        }

        break;
        solver.set_parameters (test_10.par->T, test_10.par->h, test_10.par->h0, test_10.par->angle, test_10.par->fp_x, test_10.par->fp_y, test_10.par->lb, test_10.par->ub);
        solver.form_init_fp (test_10.par->fp_x, test_10.par->fp_y, test_10.par->init_state, test_10.par->X);
        solver.solve();
        test_10.par->init_state.get_next_state (solver);
        //------------------------------------------------------
    }

    fT.close();
    fh.close();
    fa.close();
    ffpx.close();
    ffpy.close();
    flb.close();
    fub.close();

    return 0;
}
///@}
