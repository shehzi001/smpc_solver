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
            15,     // size of the preview window
            0.1,    // sampling time [sec.]
            0.261); // height of the center of mass [meter]

    // Double supports are not yet supported, hence the initial
    // double support is simulated by a single support with one
    // large foot.
    double d[4] = {0.09 , 0.025, 0.03, 0.075};
    wmg.AddFootstep(0.0, 0.05, 0.0, 3, 3, 1, d);

    // each step is defined relatively to the previous step
    double z = 5.0*M_PI/180.0;  // relative angle
    double step_x = 0.035;      // relative X position
    double step_y = 0.1;        // relative Y position

    d[3] = 0.025;   // all subsequent steps have normal feet size
    wmg.AddFootstep(0.0   , -step_y, 0.0 , 4,  4, -1, d);
    wmg.AddFootstep(step_x,  step_y, z);
    wmg.AddFootstep(step_x, -step_y, z);
    wmg.AddFootstep(step_x,  step_y, z);
    wmg.AddFootstep(step_x, -step_y, z);
    wmg.AddFootstep(step_x,  step_y, z);
    wmg.AddFootstep(step_x, -step_y, z);
    wmg.AddFootstep(step_x,  step_y, z);
    wmg.AddFootstep(step_x, -step_y, z);
    wmg.AddFootstep(step_x,  step_y, 0.0, 30, 30);
    wmg.AddFootstep(0.0   , -step_y, 0.0);

    //wmg.FS2file(); // output results for later use in Matlab/Octave
    //-----------------------------------------------------------


    smpc_solver solver(
            wmg.N); // size of the preview window
                    // other parameters are optional
    //-----------------------------------------------------------

  
    for(;;)
    {
        wmg.FormPreviewWindow();    // initialize input for QP
        if (wmg.halt)               // not enough time steps left (<15)
        {
            break;
        }
        //------------------------------------------------------


        // initialize quadratic problem
        solver.init(
                // sampling time for each time step [sec.]
                wmg.T,
                // height of the center of mass divided by gravity for each time step
                wmg.h, 
                // rotation angle for each state relative to the world frame
                wmg.angle, 
                // reference values of x coordinate of ZMP
                wmg.zref_x, 
                // reference values of y coordinate of ZMP
                wmg.zref_y, 
                // array of lower bounds for coordinates of ZMP
                wmg.lb, 
                // array of upper bounds for coordinates of ZMP
                wmg.ub, 
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
