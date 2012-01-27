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
    test_start(argv[0]);


    for (int j = 0; j < 6; j++)
    {
        ofstream fs_out;
        WMG wmg;
        std::string fs_out_filename("");

        switch (j)
        {
            case 0:
                init_01 (&wmg);
                fs_out_filename = "test_06_init_01_fs.m";
                break;
            case 1:
                init_02 (&wmg);
                fs_out_filename = "test_06_init_02_fs.m";
                break;
            case 2:
                init_03 (&wmg);
                fs_out_filename = "test_06_init_03_fs.m";
                break;
            case 3:
                init_04 (&wmg);
                fs_out_filename = "test_06_init_04_fs.m";
                break;
            case 4:
                init_05 (&wmg);
                fs_out_filename = "test_06_init_05_fs.m";
                break;
            case 5:
                init_06 (&wmg);
                fs_out_filename = "test_06_init_06_fs.m";
                break;
        }

        wmg.FS2file(fs_out_filename); // output results for later use in Matlab/Octave


        smpc::solver solver(wmg.N);


        fs_out.open(fs_out_filename.c_str(), fstream::app);
        fs_out.precision (numeric_limits<double>::digits10);
        fs_out << endl << endl;
        fs_out << "CoM_ZMP = [";
       


        for(;;)
        {
            //------------------------------------------------------
            if (wmg.formPreviewWindow() == WMG_HALT)
            {
                cout << "EXIT (halt = 1)" << endl;
                break;
            }
            //------------------------------------------------------

            //------------------------------------------------------
            solver.set_parameters (wmg.T, wmg.h, wmg.h[0], wmg.angle, wmg.fp_x, wmg.fp_y, wmg.lb, wmg.ub);
            solver.form_init_fp (wmg.fp_x, wmg.fp_y, wmg.init_state, wmg.X);
            solver.solve();
            wmg.init_state.get_next_state (solver);
            //------------------------------------------------------

            wmg.X_tilde.get_next_state (solver);
            fs_out << endl << wmg.init_state.x() << " " << wmg.init_state.y() << " " << wmg.X_tilde.x() << " " << wmg.X_tilde.y() << ";";
        }

        fs_out << "];" << endl;
        fs_out << "plot (CoM_ZMP(:,1), CoM_ZMP(:,2), 'b');" << endl;
        fs_out << "plot (CoM_ZMP(:,3), CoM_ZMP(:,4), 'ks','MarkerSize',5);" << endl;
        fs_out.close();
    }
    test_end(argv[0]);
    return 0;
}
///@}
