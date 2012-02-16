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
    int control_sampling_time_ms = 20;
    int preview_sampling_time_ms = 40;
    int next_preview_len_ms = 0;

    // initialize
    init_07 test_08 ("test_08", false);
    //-----------------------------------------------------------


    //-----------------------------------------------------------
    smpc::solver solver(
            test_08.wmg->N, // size of the preview window
            400.0,  // Alpha
            4000.0,  // Beta
            1.0,    // Gamma
            0.01,   // regularization
            1e-7);  // tolerance
    solver.enable_fexceptions();
    //-----------------------------------------------------------



    //-----------------------------------------------------------
    test_08.par->init_state.set (0.019978839010709938, -6.490507362468014e-05);
    // state_tilde = state_orig, when velocity = acceleration = 0
    test_08.X_tilde.set (0.019978839010709938, -6.490507362468014e-05);
    //-----------------------------------------------------------


    FILE *file_op = fopen(test_08.fs_out_filename.c_str(), "a");
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

    test_08.wmg->T_ms[0] = control_sampling_time_ms;
    test_08.wmg->T_ms[1] = control_sampling_time_ms;

    for(int i=0 ;; i++)
    {
        if (next_preview_len_ms == 0)
        {
            next_preview_len_ms = preview_sampling_time_ms;
        }   


        test_08.wmg->T_ms[2] = next_preview_len_ms;

        cout << test_08.wmg->isSupportSwitchNeeded() << endl;
        if (test_08.wmg->formPreviewWindow(*test_08.par) == WMG_HALT)
        {
            cout << "EXIT (halt = 1)" << endl;
            break;
        }

        ZMP_ref_x.push_back(test_08.par->zref_x[0]);
        ZMP_ref_y.push_back(test_08.par->zref_y[0]);


        
        //------------------------------------------------------
        solver.set_parameters (test_08.par->T, test_08.par->h, test_08.par->h0, test_08.par->angle, test_08.par->zref_x, test_08.par->zref_y, test_08.par->lb, test_08.par->ub);
        solver.form_init_fp (test_08.par->fp_x, test_08.par->fp_y, test_08.par->init_state, test_08.par->X);
        solver.solve();
        //------------------------------------------------------
        // update state
        test_08.par->init_state.get_next_state(solver);
        //-----------------------------------------------------------


        if (next_preview_len_ms == preview_sampling_time_ms)
        {
            // if the values are saved on each iteration the plot becomes sawlike.
            // better solution - more frequent sampling.
            ZMP_x.push_back(test_08.X_tilde.x());
            ZMP_y.push_back(test_08.X_tilde.y());
            test_08.X_tilde.get_next_state (solver);
        }
        CoM_x.push_back(test_08.par->init_state.x());
        CoM_y.push_back(test_08.par->init_state.y());
    

        // feet position/orientation
        double left_foot_pos[16];
        double right_foot_pos[16];
        test_08.wmg->getFeetPositions (control_sampling_time_ms, left_foot_pos, right_foot_pos);
                
        left_foot_x.push_back(left_foot_pos[12]);
        left_foot_y.push_back(left_foot_pos[13]);
        left_foot_z.push_back(left_foot_pos[14]);

        right_foot_x.push_back(right_foot_pos[12]);
        right_foot_y.push_back(right_foot_pos[13]);
        right_foot_z.push_back(right_foot_pos[14]);
        
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

    return 0;
}
///@}
