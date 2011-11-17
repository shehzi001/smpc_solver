/** 
 * @file
 * @author Alexander Sherikov
 * @brief Performs a full simulation and compares results with 
 *  reference data produced by Octave.
 */


#include <sstream>
#include "tests_common.h"

///@addtogroup gTEST
///@{

int main(int argc, char **argv)
{
    test_start(argv[0]);


    for (int j = 0; j < 3; j++)
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
        }

        wmg.FS2file(fs_out_filename); // output results for later use in Matlab/Octave


        smpc_solver solver(wmg.N);


        double X[6];
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
            solver.set_parameters (wmg.T, wmg.h, wmg.angle, wmg.zref_x, wmg.zref_y, wmg.lb, wmg.ub);
            solver.form_init_fp (wmg.zref_x, wmg.zref_y, wmg.X_tilde, wmg.X);
            solver.solve();
            solver.get_next_state_tilde (wmg.X_tilde);
            //------------------------------------------------------


            solver.get_next_state (X);
            fs_out << endl << X[0] << " " << X[3] << " " << wmg.X_tilde[0] << " " << wmg.X_tilde[3] << ";";
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
