#include <iostream>
#include <fstream>
#include <cstdio>
#include <limits>
#include <cmath> // abs, M_PI


#include "WMG.h"
#include "smpc_solver.h" 

using namespace std;

///@addtogroup gTEST
///@{

void test_start (char *name)
{
    printf ("\n################################\n %s \n################################\n", name);
}

void test_end (char *name)
{
    printf ("################################\n");
}


//**********************************************************************

void init_01 (WMG *wmg)
{
    wmg->init(15, 0.1, 0.261);

    double d[4] = {0.09 , 0.025, 0.03, 0.075};
    wmg->AddFootstep(0.0, 0.05, 0.0, 3, 3, d, FS_TYPE_DS);

    double z = 5.0*M_PI/180.0;
    double step_x = 0.035;
    double step_y = 0.1;

    d[3] = 0.025;
    wmg->AddFootstep(0.0   , -step_y, 0.0 , 4,  4, d);
    wmg->AddFootstep(step_x,  step_y, z);
    wmg->AddFootstep(step_x, -step_y, z);
    wmg->AddFootstep(step_x,  step_y, z);
    wmg->AddFootstep(step_x, -step_y, z);
    wmg->AddFootstep(step_x,  step_y, z);
    wmg->AddFootstep(step_x, -step_y, z);
    wmg->AddFootstep(step_x,  step_y, z);
    wmg->AddFootstep(step_x, -step_y, z);
    wmg->AddFootstep(step_x,  step_y, 0.0);
    wmg->AddFootstep(0.0   , -step_y, 0.0, 30, 30);
}


void init_02 (WMG *wmg)
{
    wmg->init(15, 0.1, 0.261);

    double d[4] = {0.03 , 0.01, 0.01, 0.11};
    wmg->AddFootstep(0.0, 0.05, 0.0, 3, 3, d);
    
    double z = 5.0*M_PI/180.0;
    double step_x = 0.035;
    double step_y = 0.1;

    d[3] = 0.025;

    // use this for smaller feet (and we will have more active constraints)
    d[0] = 0.03; d[1] = 0.01; d[2] = 0.01; d[3] = 0.01;

    wmg->AddFootstep(0.0   , -step_y, 0.0 , 4,  4, d);
    wmg->AddFootstep(step_x,  step_y, z);
    wmg->AddFootstep(step_x, -step_y, z);
    wmg->AddFootstep(step_x,  step_y, z);
    wmg->AddFootstep(step_x, -step_y, z);
    wmg->AddFootstep(step_x,  step_y, z);
    wmg->AddFootstep(step_x, -step_y, z);
    wmg->AddFootstep(step_x,  step_y, z);
    wmg->AddFootstep(step_x, -step_y, z);
    wmg->AddFootstep(step_x,  step_y, 0.0, 30, 30);
    wmg->AddFootstep(0.0   , -step_y, 0.0);
}



void init_03 (WMG *wmg)
{
    wmg->init(15, 0.1, 0.261);

    // Initial double support
    double d[4] = {0.09 , 0.075, 0.03, 0.025};
    wmg->AddFootstep(0.0, 0.0, 0.0, 2, 3, d, FS_TYPE_DS);
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
    wmg->AddFootstep(0.0   , -step_y/2, 0.0 , 2,  3, d);
    wmg->AddFootstep(step_x,  step_y, z);
    wmg->AddFootstep(step_x, -step_y, z);
    wmg->AddFootstep(step_x,  step_y, z);
    wmg->AddFootstep(step_x, -step_y, z);
    wmg->AddFootstep(step_x,  step_y, z);
    wmg->AddFootstep(step_x, -step_y, z);
    wmg->AddFootstep(step_x,  step_y, 0.0);

    // here we give many reference points, since otherwise we 
    // would not have enough steps in preview window to reach 
    // the last footsteps
    d[0] = 0.09;
    d[1] = 0.025;
    d[2] = 0.03;
    d[3] = 0.075;
    wmg->AddFootstep(0.0   , -step_y/2, 0.0, 30, 30, d, FS_TYPE_DS);
}

///@}
