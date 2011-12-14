/** 
 * @file
 * @author Alexander Sherikov
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
    double control_sampling_time = (double) control_sampling_time_ms / 1000;
    double A[9];
    double B[3];
    double cur_control[2];
    
    A[0] = A[4] = A[8] = 1;
    A[1] = A[2] = A[5] = 0;
    A[3] = A[7] = control_sampling_time;
    A[6] = control_sampling_time * control_sampling_time/2 /*- delta_hCoM = 0*/;

    B[0] = control_sampling_time * control_sampling_time * control_sampling_time / 6
        - wmg.hCoM/wmg.gravity * control_sampling_time;
    B[1] = control_sampling_time * control_sampling_time/2;
    B[2] = control_sampling_time;

    for (int i = 0; i < 6; i++)
    {
        wmg.X_tilde[i] = 0;
    }
    cur_control[0] = cur_control[1] = 0;
    //-----------------------------------------------------------


    FILE *file_op = fopen(fs_out_filename.c_str(), "a");
    fprintf(file_op,"hold on\n");

    double X[SMPC_NUM_STATE_VAR];
    vector<double> ZMP_x;
    vector<double> ZMP_y;
    vector<double> CoM_x;
    vector<double> CoM_y;

    vector<double> swing_foot_x;
    vector<double> swing_foot_y;
    vector<double> swing_foot_z;
    for(int i=0 ;; i++)
    {
        //-----------------------------------------------------------
        // update state
        wmg.X_tilde[0] = wmg.X_tilde[0] * A[0] 
                         + wmg.X_tilde[1] * A[3]
                         + wmg.X_tilde[2] * A[6]
                         + cur_control[0] * B[0];

        wmg.X_tilde[1] = wmg.X_tilde[1] * A[4]
                         + wmg.X_tilde[2] * A[7]
                         + cur_control[0] * B[1];

        wmg.X_tilde[2] = wmg.X_tilde[2] * A[8]
                         + cur_control[0] * B[2];

        wmg.X_tilde[3] = wmg.X_tilde[3] * A[0] 
                         + wmg.X_tilde[4] * A[3]
                         + wmg.X_tilde[5] * A[6]
                         + cur_control[1] * B[0];

        wmg.X_tilde[4] = wmg.X_tilde[4] * A[4]
                         + wmg.X_tilde[5] * A[7]
                         + cur_control[1] * B[1];

        wmg.X_tilde[5] = wmg.X_tilde[5] * A[8]
                         + cur_control[1] * B[2];
        //-----------------------------------------------------------




        if (next_preview_len_ms == 0)
        {
            WMGret wmg_retval = wmg.FormPreviewWindow();

            if (wmg_retval == WMG_HALT)
            {
                cout << "EXIT (halt = 1)" << endl;
                break;
            }

            next_preview_len_ms = preview_sampling_time_ms;
        }   

       
        /// @attention wmg.X_tilde does not always satisfy the lower and upper bounds!
        /// (but wmg.X does)
        
        wmg.T[0] = (double) next_preview_len_ms / 1000; // get seconds
        //------------------------------------------------------
        solver.set_parameters (wmg.T, wmg.h, wmg.angle, wmg.zref_x, wmg.zref_y, wmg.lb, wmg.ub);
        solver.form_init_fp (wmg.zref_x, wmg.zref_y, wmg.X_tilde, wmg.X);
        solver.solve();
        solver.get_next_state (X);
        solver.get_first_controls (cur_control);
        //------------------------------------------------------
        ZMP_x.push_back(wmg.X_tilde[0]);
        ZMP_y.push_back(wmg.X_tilde[3]);
        CoM_x.push_back(X[0]);
        CoM_y.push_back(X[3]);
    

        // support foot and swing foot position/orientation
        double LegPos[3];
        double angle;
        wmg.getSwingFootPosition (
                WMG_SWING_PARABOLA,
                preview_sampling_time_ms / control_sampling_time_ms,
                (preview_sampling_time_ms - next_preview_len_ms) / control_sampling_time_ms,
                LegPos,
                &angle);

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

    fprintf(file_op,"SFP = [\n");
    for (unsigned int i=0; i < swing_foot_x.size(); i++)
    {
        fprintf(file_op, "%f %f %f;\n", swing_foot_x[i], swing_foot_y[i], swing_foot_z[i]);
    }
    fprintf(file_op, "];\n\n plot3(SFP(:,1), SFP(:,2), SFP(:,3), 'r')\n");


    fprintf(file_op,"ZMP = [\n");
    for (unsigned int i=0; i < ZMP_x.size(); i++)
    {
        fprintf(file_op, "%f %f %f;\n", ZMP_x[i], ZMP_y[i], 0.0);
    }
    fprintf(file_op, "];\n\n plot3(ZMP(:,1), ZMP(:,2), ZMP(:,3), 'k')\n");


    fprintf(file_op,"CoM = [\n");
    for (unsigned int i=0; i < CoM_x.size(); i++)
    {
        fprintf(file_op, "%f %f %f;\n", CoM_x[i], CoM_y[i], 0.0);
    }
    fprintf(file_op, "];\n\n plot3(CoM(:,1), CoM(:,2), CoM(:,3), 'b')\n");


    fprintf(file_op,"hold off\n");
    fclose(file_op);
    test_end(argv[0]);

    return 0;
}
///@}
