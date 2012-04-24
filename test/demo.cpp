#include <cmath> // M_PI


#include "WMG.h"    // a helper class
#include "smpc_solver.h"    // the interface of the library


///@addtogroup gTEST
///@{

int main(int argc, char **argv)
{
    //-----------------------------------------------------------
    // initialize

    // WMG is a helper class, which provides data for QP on each
    // iteration of simulation
    WMG wmg(15,     // size of the preview window
            100);   // sampling time [ms.]
    smpc_parameters param (
            wmg.N,  // size of the preview window
            0.261); // height of the center of mass [meter]

    // each step is defined relatively to the previous step
    double z = 5.0*M_PI/180.0;  // relative angle
    double step_x = 0.035;      // relative X position
    double step_y = wmg.def_constraints.support_distance_y; // relative Y position


    // Initial double support
    // ZMP, CoM are at [0;0]
    wmg.setFootstepParameters (2, 1, 1);
    wmg.addFootstep(0.0, 0.0, 0.0, FS_TYPE_DS);

    // all subsequent steps have normal feet size
    // 2 reference ZMP positions in single support 
    // 1 in double support
    wmg.addFootstep(0.0   , -step_y/2, 0.0);
    wmg.addFootstep(step_x,  step_y, z);
    wmg.addFootstep(step_x, -step_y, z);
    wmg.addFootstep(step_x,  step_y, z);
    wmg.addFootstep(step_x, -step_y, z);
    wmg.addFootstep(step_x,  step_y, z);
    wmg.addFootstep(step_x, -step_y, z);
    wmg.addFootstep(step_x,  step_y, 0.0);

    // here we give many reference points, since otherwise we 
    // would not have enough steps in preview window to reach 
    // the last footsteps
    wmg.setFootstepParameters (30, 0, 0);
    wmg.addFootstep(0.0   , -step_y/2, 0.0, FS_TYPE_DS);
    //-----------------------------------------------------------


    smpc::solver_as solver(
            wmg.N); // size of the preview window
                    // other parameters are optional
    //-----------------------------------------------------------

  
    for(;;)
    {
        if (wmg.formPreviewWindow(param) == WMG_HALT) // initialize input for QP
        {
            // not enough time steps left (<15)
            break;
        }
        //------------------------------------------------------


        // initialize quadratic problem
        solver.set_parameters(
                // sampling time for each time step [sec.]
                param.T,
                // height of the center of mass divided by gravity for each time step
                param.h, 
                // current height of the center of mass divided by gravity
                param.h0, 
                // rotation angle for each state relative to the world frame
                param.angle, 
                // reference values of x coordinate of ZMP
                param.fp_x, 
                // reference values of y coordinate of ZMP
                param.fp_y, 
                // array of lower bounds for coordinates of ZMP
                param.lb, 
                // array of upper bounds for coordinates of ZMP
                param.ub);

        solver.form_init_fp (
                // coordinates of points, that can be used to generate
                // an initial feasible point
                param.fp_x, 
                param.fp_y, 
                // initial state
                param.init_state, 
                // solution of the optimization problem
                param.X);

        //------------------------------------------------------

        // solve quadratic problem
        solver.solve();
        //------------------------------------------------------

        // obtain the next state
        solver.get_next_state(param.init_state);
    }

    return 0;
}
//---------------------------------------------------------------
///@}
