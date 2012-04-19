/** 
 * @file
 * @author Alexander Sherikov
 * @brief Perform simulation using different footstep patterns.
 */


#include <sstream>
#include "tests_common.h"

///@addtogroup gTEST
///@{

int main(int argc, char **argv)
{
    test_init_base* test_06 = NULL;


    for (int j = 0; j < 6; j++)
    {
        ofstream fs_out;

        switch (j)
        {
            case 0:
                test_06 = new init_01 ("test_06_init_01");
                break;
            case 1:
                test_06 = new init_02 ("test_06_init_02");
                break;
            case 2:
                test_06 = new init_03 ("test_06_init_03");
                break;
            case 3:
                test_06 = new init_04 ("test_06_init_04");
                break;
            case 4:
                test_06 = new init_05 ("test_06_init_05");
                break;
            case 5:
                test_06 = new init_06 ("test_06_init_06");
                break;
        }


        smpc::solver_as solver(test_06->wmg->N);


        fs_out.open(test_06->fs_out_filename.c_str(), fstream::app);
        fs_out.precision (numeric_limits<double>::digits10);
        fs_out << endl << endl;
        fs_out << "CoM_ZMP = [";
       


        for(;;)
        {
            //------------------------------------------------------
            if (test_06->wmg->formPreviewWindow(*test_06->par) == WMG_HALT)
            {
                cout << "EXIT (halt = 1)" << endl;
                break;
            }
            //------------------------------------------------------

            //------------------------------------------------------
            solver.set_parameters (test_06->par->T, test_06->par->h, test_06->par->h0, test_06->par->angle, test_06->par->fp_x, test_06->par->fp_y, test_06->par->lb, test_06->par->ub);
            solver.form_init_fp (test_06->par->fp_x, test_06->par->fp_y, test_06->par->init_state, test_06->par->X);
            solver.solve();
            solver.get_next_state(test_06->par->init_state);
            //------------------------------------------------------

            solver.get_next_state(test_06->X_tilde);
            fs_out << endl << test_06->par->init_state.x() << " " << test_06->par->init_state.y() << " " << test_06->X_tilde.x() << " " << test_06->X_tilde.y() << ";";
        }

        fs_out << "];" << endl;
        fs_out << "plot (CoM_ZMP(:,1), CoM_ZMP(:,2), 'b');" << endl;
        fs_out << "plot (CoM_ZMP(:,3), CoM_ZMP(:,4), 'ks','MarkerSize',5);" << endl;
        fs_out.close();

        delete test_06;
    }
    return 0;
}
///@}
