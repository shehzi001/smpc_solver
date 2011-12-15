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
    WMG wmg;

    wmg.init (
            15);    // size of the preview window
    wmg.init_param (
            0.1,    // sampling time [sec.]
            0.261); // height of the center of mass [meter]


    // Initial double support
    double d[4] = {0.09 , 0.075, 0.03, 0.025};
    wmg.AddFootstep(0.0, 0.0, 0.0, 2, 3, d, FS_TYPE_DS);
    // ZMP, CoM are at [0;0]

    // each step is defined relatively to the previous step
    double z = 5.0*M_PI/180.0;  // relative angle
    double step_x = 0.035;      // relative X position
    double step_y = 0.1;        // relative Y position

    // all subsequent steps have normal feet size
    d[0] = 0.09;
    d[1] = 0.025;
    d[2] = 0.03;
    d[3] = 0.025;
    // 2 reference ZMP positions in single support 
    // 1 in double support
    // 1 + 2 = 3
    wmg.AddFootstep(0.0   , -step_y/2, 0.0 , 2,  3, d);
    wmg.AddFootstep(step_x,  step_y, z);
    wmg.AddFootstep(step_x, -step_y, z);
    wmg.AddFootstep(step_x,  step_y, z);
    wmg.AddFootstep(step_x, -step_y, z);
    wmg.AddFootstep(step_x,  step_y, z);
    wmg.AddFootstep(step_x, -step_y, z);
    wmg.AddFootstep(step_x,  step_y, 0.0);

    // here we give many reference points, since otherwise we 
    // would not have enough steps in preview window to reach 
    // the last footsteps
    d[0] = 0.09;
    d[1] = 0.025;
    d[2] = 0.03;
    d[3] = 0.075;
    wmg.AddFootstep(0.0   , -step_y/2, 0.0, 30, 30, d, FS_TYPE_DS);
    //-----------------------------------------------------------


    smpc_solver solver(
            wmg.N); // size of the preview window
                    // other parameters are optional
    //-----------------------------------------------------------

  
    for(;;)
    {
        if (wmg.FormPreviewWindow() == WMG_HALT) // initialize input for QP
        {
            // not enough time steps left (<15)
            break;
        }
        //------------------------------------------------------


        // initialize quadratic problem
        solver.set_parameters(
                // sampling time for each time step [sec.]
                wmg.T,
                // height of the center of mass divided by gravity for each time step
                wmg.h, 
                // rotation angle for each state relative to the world frame
                wmg.angle, 
                // reference values of x coordinate of ZMP
                wmg.fp_x, 
                // reference values of y coordinate of ZMP
                wmg.fp_y, 
                // array of lower bounds for coordinates of ZMP
                wmg.lb, 
                // array of upper bounds for coordinates of ZMP
                wmg.ub);

        solver.form_init_fp (
                // coordinates of points, that can be used to generate
                // an initial feasible point
                wmg.fp_x, 
                wmg.fp_y, 
                // current state
                wmg.X_tilde, 
                // solution of the optimization problem
                wmg.X);

        //------------------------------------------------------

        // solve quadratic problem
        solver.solve();
        //------------------------------------------------------

        // obtain the next state as X_tilde
        solver.get_next_state_tilde (wmg.X_tilde);
    }

    return 0;
}
//---------------------------------------------------------------
///@}
