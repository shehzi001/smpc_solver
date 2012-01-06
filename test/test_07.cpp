/** 
 * @file
 * @author Alexander Sherikov
 * @brief Simulate control loop, which is shorter than preview window iteration.
 */


#include <cstring> //strcmp
#include "tests_common.h"

///@addtogroup gTEST
///@{

int main(int argc, char **argv)
{
    //-----------------------------------------------------------
    // the numbers must correspond to the numbers in init_04()
    int control_sampling_time_ms = 10;
    int preview_sampling_time_ms = 100;
    int next_preview_len_ms = 0;

    // initialize
    WMG wmg;
    init_04 (&wmg);

    std::string fs_out_filename("test_07_fs.m");
    wmg.FS2file(fs_out_filename); // output results for later use in Matlab/Octave
    //-----------------------------------------------------------


    test_start(argv[0]);
    //-----------------------------------------------------------
    smpc_solver solver(
            wmg.N, // size of the preview window
            300.0,  // Alpha
            800.0,  // Beta
            1.0,    // Gamma
            0.01,   // regularization
            1e-7);  // tolerance
    //-----------------------------------------------------------



    //-----------------------------------------------------------
    wmg.initABMatrices ((double) control_sampling_time_ms / 1000);
    wmg.init_state[0] = 0.019978839010709938;
    wmg.init_state[1] = 0;
    wmg.init_state[2] = 0;
    wmg.init_state[3] = -6.490507362468014e-05;
    wmg.init_state[4] = 0;
    wmg.init_state[5] = 0;
    double cur_control[2];
    cur_control[0] = cur_control[1] = 0;
    //-----------------------------------------------------------


    FILE *file_op = fopen(fs_out_filename.c_str(), "a");
    fprintf(file_op,"hold on\n");

    vector<double> ZMP_ref_x;
    vector<double> ZMP_ref_y;
    vector<double> ZMP_x;
    vector<double> ZMP_y;
    vector<double> CoM_x;
    vector<double> CoM_y;

    vector<double> swing_foot_x;
    vector<double> swing_foot_y;
    vector<double> swing_foot_z;


    for(int i=0 ;; i++)
    {
        if (next_preview_len_ms == 0)
        {
            WMGret wmg_retval = wmg.FormPreviewWindow();

            if (wmg_retval == WMG_HALT)
            {
                cout << "EXIT (halt = 1)" << endl;
                break;
            }

            ZMP_ref_x.push_back(wmg.zref_x[0]);
            ZMP_ref_y.push_back(wmg.zref_y[0]);

            next_preview_len_ms = preview_sampling_time_ms;
        }   

        
        //------------------------------------------------------
        wmg.T[0] = (double) next_preview_len_ms / 1000; // get seconds
        solver.set_parameters (wmg.T, wmg.h, wmg.h[0], wmg.angle, wmg.zref_x, wmg.zref_y, wmg.lb, wmg.ub);
        solver.form_init_fp (wmg.fp_x, wmg.fp_y, wmg.init_state, wmg.X);
        solver.solve();
        //------------------------------------------------------
        // update state
        solver.get_first_controls (cur_control);
        wmg.calculateNextState(cur_control, wmg.init_state);
        //-----------------------------------------------------------


        solver.get_next_state_tilde (wmg.X_tilde);
        ZMP_x.push_back(wmg.X_tilde[0]);
        ZMP_y.push_back(wmg.X_tilde[3]);
        double X[SMPC_NUM_STATE_VAR];
        solver.get_next_state (X);
        CoM_x.push_back(X[0]);
        CoM_y.push_back(X[3]);
    

        // support foot and swing foot position/orientation
        double LegPos[3];
        double angle;
        /* wrong, but makes nice graph */
        wmg.getSwingFootPosition (
                WMG_SWING_2D_PARABOLA,
                preview_sampling_time_ms / control_sampling_time_ms,
                (preview_sampling_time_ms - next_preview_len_ms) / control_sampling_time_ms,
                LegPos,
                &angle);
                
        /* correct
        wmg.getSwingFootPosition (
                WMG_SWING_2D_PARABOLA,
                1,
                1,
                LegPos,
                &angle);
                */

        swing_foot_x.push_back(LegPos[0]);
        swing_foot_y.push_back(LegPos[1]);
        swing_foot_z.push_back(LegPos[2]);
        /*
        fprintf(file_op, "plot3([%f %f], [%f %f], [%f %f])\n",
                LegPos[0], LegPos[0] + cos(angle)*0.005,
                LegPos[1], LegPos[1] + sin(angle)*0.005,
                LegPos[2], LegPos[2]);
        */
        
        next_preview_len_ms -= control_sampling_time_ms;
    }

    /*
    fprintf(file_op,"SFP = [\n");
    for (unsigned int i=0; i < swing_foot_x.size(); i++)
    {
        fprintf(file_op, "%f %f %f;\n", swing_foot_x[i], swing_foot_y[i], swing_foot_z[i]);
    }
    fprintf(file_op, "];\n\n plot3(SFP(:,1), SFP(:,2), SFP(:,3), 'r')\n");
    */


    fprintf(file_op,"ZMP = [\n");
    for (unsigned int i=0; i < ZMP_x.size(); i++)
    {
        fprintf(file_op, "%f %f;\n", ZMP_x[i], ZMP_y[i]);
    }
    fprintf(file_op, "];\n\n plot(ZMP(:,1), ZMP(:,2), 'k')\n");

    fprintf(file_op,"ZMPref = [\n");
    for (unsigned int i=0; i < ZMP_ref_x.size(); i++)
    {
        fprintf(file_op, "%f %f;\n", ZMP_ref_x[i], ZMP_ref_y[i]);
    }
    fprintf(file_op, "];\n\n plot(ZMPref(:,1), ZMPref(:,2), 'ko')\n");

    fprintf(file_op,"CoM = [\n");
    for (unsigned int i=0; i < CoM_x.size(); i++)
    {
        fprintf(file_op, "%f %f;\n", CoM_x[i], CoM_y[i]);
    }
    fprintf(file_op, "];\n\n plot(CoM(:,1), CoM(:,2), 'b')\n");


    fprintf(file_op,"hold off\n");
    fclose(file_op);
    test_end(argv[0]);

    return 0;
}
///@}
