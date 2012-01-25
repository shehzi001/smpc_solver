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
    smpc::solver solver(
            wmg.N, // size of the preview window
            300.0,  // Alpha
            800.0,  // Beta
            1.0,    // Gamma
            0.01,   // regularization
            1e-7);  // tolerance
    //-----------------------------------------------------------



    //-----------------------------------------------------------
    wmg.initABMatrices ((double) control_sampling_time_ms / 1000);
    wmg.init_state.set (0.019978839010709938, -6.490507362468014e-05);
    // state_tilde = state_orig, when velocity = acceleration = 0
    wmg.X_tilde.set (0.019978839010709938, -6.490507362468014e-05);
    //-----------------------------------------------------------


    FILE *file_op = fopen(fs_out_filename.c_str(), "a");
    fprintf(file_op,"hold on\n");

    vector<double> ZMP_ref_x;
    vector<double> ZMP_ref_y;
    vector<double> ZMP_x;
    vector<double> ZMP_y;
    vector<double> CoM_x;
    vector<double> CoM_y;

    vector<double> left_foot_x;
    vector<double> left_foot_y;
    vector<double> left_foot_z;

    vector<double> right_foot_x;
    vector<double> right_foot_y;
    vector<double> right_foot_z;


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
        wmg.next_control.get_first_controls (solver);
        wmg.calculateNextState(wmg.next_control, wmg.init_state);
        //-----------------------------------------------------------


        if (next_preview_len_ms == preview_sampling_time_ms)
        {
            // if the values are saved on each iteration the plot becomes sawlike.
            // better solution - more frequent sampling.
            ZMP_x.push_back(wmg.X_tilde.x());
            ZMP_y.push_back(wmg.X_tilde.y());
            wmg.X_tilde.get_next_state (solver);
        }
        CoM_x.push_back(wmg.init_state.x());
        CoM_y.push_back(wmg.init_state.y());
    

        // feet position/orientation
        double left_foot_pos[3+1];
        double right_foot_pos[3+1];
        /* wrong, but makes nice graph */
        wmg.getFeetPositions (
                preview_sampling_time_ms / control_sampling_time_ms,
                (preview_sampling_time_ms - next_preview_len_ms) / control_sampling_time_ms + 1,
                left_foot_pos,
                right_foot_pos);
                
        left_foot_x.push_back(left_foot_pos[0]);
        left_foot_y.push_back(left_foot_pos[1]);
        left_foot_z.push_back(left_foot_pos[2]);

        right_foot_x.push_back(right_foot_pos[0]);
        right_foot_y.push_back(right_foot_pos[1]);
        right_foot_z.push_back(right_foot_pos[2]);
        
        next_preview_len_ms -= control_sampling_time_ms;
    }

    // feet positions    
    fprintf(file_op,"LFP = [\n");
    for (unsigned int i=0; i < left_foot_x.size(); i++)
    {
        fprintf(file_op, "%f %f %f;\n", left_foot_x[i], left_foot_y[i], left_foot_z[i]);
    }
    fprintf(file_op, "];\n\n plot3(LFP(:,1), LFP(:,2), LFP(:,3), 'r')\n");
    
    fprintf(file_op,"RFP = [\n");
    for (unsigned int i=0; i < right_foot_x.size(); i++)
    {
        fprintf(file_op, "%f %f %f;\n", right_foot_x[i], right_foot_y[i], right_foot_z[i]);
    }
    fprintf(file_op, "];\n\n plot3(RFP(:,1), RFP(:,2), RFP(:,3), 'r')\n");


    // ZMP
    fprintf(file_op,"ZMP = [\n");
    for (unsigned int i=0; i < ZMP_x.size(); i++)
    {
        fprintf(file_op, "%f %f;\n", ZMP_x[i], ZMP_y[i]);
    }
    fprintf(file_op, "];\n\n plot(ZMP(:,1), ZMP(:,2), 'k')\n");

    // reference ZMP points
    fprintf(file_op,"ZMPref = [\n");
    for (unsigned int i=0; i < ZMP_ref_x.size(); i++)
    {
        fprintf(file_op, "%f %f;\n", ZMP_ref_x[i], ZMP_ref_y[i]);
    }
    fprintf(file_op, "];\n\n plot(ZMPref(:,1), ZMPref(:,2), 'ko')\n");

    // CoM
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
