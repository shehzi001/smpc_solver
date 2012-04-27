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

class test_init_base
{
    public:
        test_init_base(const string& test_name, const bool plot_ds_)
        {
            name = test_name;
            plot_ds = plot_ds_;

            if (!name.empty())
            {
                cout << "################################" << endl;
                cout << name << endl;
                cout << "################################" << endl;
                fs_out_filename = name + "_fs.m";
            }
        }
        ~test_init_base()
        {
            if (!name.empty())
            {
                cout << "################################" << endl;
            }
        }



        smpc::state_tilde X_tilde;
        WMG* wmg;
        smpc_parameters* par;

        string name;
        string fs_out_filename;
        bool plot_ds;
};


//**********************************************************************

class init_01 : public test_init_base
{
    public:
        init_01 (const string & test_name, const bool plot_ds_ = true) : 
            test_init_base (test_name, plot_ds_)
        {
            wmg = new WMG (15, 100);
            par = new smpc_parameters (wmg->N, 0.261);


            double z = 5.0*M_PI/180.0;
            double step_x = 0.035;
            double step_y = wmg->def_constraints.support_distance_y;


            // the test data is generated for this rectangular 
            // constraint which is not properly chosen (though it
            // is ok for testing purposes)
            wmg->user_constraints[0] = 0.09;
            wmg->user_constraints[1] = 0.025;
            wmg->user_constraints[2] = 0.03;
            wmg->user_constraints[3] = 0.075;


            wmg->setFootstepParameters (3, 0, 0, true);
            wmg->addFootstep(0.0, 0.05, 0.0, FS_TYPE_DS);
            wmg->setFootstepParameters (4, 0, 0, false);
            wmg->addFootstep(0.0   , -step_y, 0.0);
            wmg->addFootstep(step_x,  step_y, z);
            wmg->addFootstep(step_x, -step_y, z);
            wmg->addFootstep(step_x,  step_y, z);
            wmg->addFootstep(step_x, -step_y, z);
            wmg->addFootstep(step_x,  step_y, z);
            wmg->addFootstep(step_x, -step_y, z);
            wmg->addFootstep(step_x,  step_y, z);
            wmg->addFootstep(step_x, -step_y, z);
            wmg->addFootstep(step_x,  step_y, 0.0);
            wmg->setFootstepParameters (30, 0, 0);
            wmg->addFootstep(0.0   , -step_y, 0.0);

            if (!name.empty())
            {
                wmg->FS2file(fs_out_filename, plot_ds); 
            }
        }
};



class init_02 : public test_init_base
{
    public:
        init_02 (const string & test_name, const bool plot_ds_ = true) : 
            test_init_base (test_name, plot_ds_)
        {
            wmg = new WMG (15, 100);
            par = new smpc_parameters (wmg->N, 0.261);


            double z = 5.0*M_PI/180.0;
            double step_x = 0.035;
            double step_y = wmg->def_constraints.support_distance_y;


            wmg->user_constraints[0] = 0.03;
            wmg->user_constraints[1] = 0.01;
            wmg->user_constraints[2] = 0.01;
            wmg->user_constraints[3] = 0.11;

            wmg->setFootstepParameters (3, 0, 0, true);
            wmg->addFootstep(0.0, 0.05, 0.0, FS_TYPE_DS);
            

            // use this for smaller feet (and we will have more active constraints)
            wmg->user_constraints[0] = 0.03;
            wmg->user_constraints[1] = 0.01;
            wmg->user_constraints[2] = 0.01;
            wmg->user_constraints[3] = 0.01;

            wmg->setFootstepParameters (4, 0, 0, true);
            wmg->addFootstep(0.0   , -step_y, 0.0);
            wmg->addFootstep(step_x,  step_y, z);
            wmg->addFootstep(step_x, -step_y, z);
            wmg->addFootstep(step_x,  step_y, z);
            wmg->addFootstep(step_x, -step_y, z);
            wmg->addFootstep(step_x,  step_y, z);
            wmg->addFootstep(step_x, -step_y, z);
            wmg->addFootstep(step_x,  step_y, z);
            wmg->addFootstep(step_x, -step_y, z);
            wmg->setFootstepParameters (30, 0, 0, true);
            wmg->addFootstep(step_x,  step_y, 0.0);
            wmg->addFootstep(0.0   , -step_y, 0.0);

            if (!name.empty())
            {
                wmg->FS2file(fs_out_filename, plot_ds); 
            }
        }
};



class init_03 : public test_init_base
{
    public:
        init_03 (const string & test_name, const bool plot_ds_ = true) : 
            test_init_base (test_name, plot_ds_)
        {
            wmg = new WMG (15, 100);
            par = new smpc_parameters (wmg->N, 0.261);

            // each step is defined relatively to the previous step
            double z = 5.0*M_PI/180.0;  // relative angle
            double step_x = 0.035;      // relative X position
            double step_y = wmg->def_constraints.support_distance_y;        // relative Y position


            // Initial double support
            wmg->user_constraints[0] = 0.09;
            wmg->user_constraints[1] = 0.075;
            wmg->user_constraints[2] = 0.03;
            wmg->user_constraints[3] = 0.025;
            wmg->setFootstepParameters (2, 1, 1, true);
            wmg->addFootstep(0.0, 0.0, 0.0, FS_TYPE_DS);
            // ZMP, CoM are at [0;0]

            // all subsequent steps have normal feet size
            // 2 reference ZMP positions in single support 
            // 1 in double support
            wmg->setFootstepParameters (2, 1, 1, false);
            wmg->addFootstep(0.0   , -step_y/2, 0.0);
            wmg->addFootstep(step_x,  step_y, z);
            wmg->addFootstep(step_x, -step_y, z);
            wmg->addFootstep(step_x,  step_y, z);
            wmg->addFootstep(step_x, -step_y, z);
            wmg->addFootstep(step_x,  step_y, z);
            wmg->addFootstep(step_x, -step_y, z);
            wmg->addFootstep(step_x,  step_y, 0.0);

            // here we give many reference points, since otherwise we 
            // would not have enough steps in preview window to reach 
            // the last footsteps
            wmg->user_constraints[0] = 0.09;
            wmg->user_constraints[1] = 0.025;
            wmg->user_constraints[2] = 0.03;
            wmg->user_constraints[3] = 0.075;
            wmg->setFootstepParameters (30, 0, 0, true);
            wmg->addFootstep(0.0   , -step_y/2, 0.0, FS_TYPE_DS);

            if (!name.empty())
            {
                wmg->FS2file(fs_out_filename, plot_ds); 
            }
        }
};



/**
 * @brief Walk straight
 */
class init_04 : public test_init_base
{
    public:
        init_04 (const string & test_name, const bool plot_ds_ = true) : 
            test_init_base (test_name, plot_ds_)
        {
            wmg = new WMG (15, 100);
            par = new smpc_parameters (wmg->N, 0.261);

            // each step is defined relatively to the previous step
            double step_x = 0.035;      // relative X position
            double step_y = wmg->def_constraints.support_distance_y;       // relative Y position

            wmg->setFootstepParameters (0, 0, 0);
            wmg->addFootstep(0.0, step_y/2, 0.0, FS_TYPE_SS_L);

            // Initial double support
            wmg->setFootstepParameters (2, 0, 0);
            wmg->addFootstep(0.0, -step_y/2, 0.0, FS_TYPE_DS);
            // ZMP, CoM are at [0;0]


            // all subsequent steps have normal feet size
            wmg->setFootstepParameters (6, 1, 2);
            wmg->addFootstep(0.0   , -step_y/2, 0.0);
            wmg->addFootstep(step_x,  step_y, 0.0);
            wmg->addFootstep(step_x, -step_y, 0.0);
            wmg->addFootstep(step_x,  step_y, 0.0);
            wmg->addFootstep(step_x, -step_y, 0.0);
            wmg->addFootstep(step_x,  step_y, 0.0);

            // here we give many reference points, since otherwise we 
            // would not have enough steps in preview window to reach 
            // the last footsteps
            wmg->setFootstepParameters (30, 0, 0);
            wmg->addFootstep(0.0   , -step_y/2, 0.0, FS_TYPE_DS);
            wmg->setFootstepParameters (0, 0, 0);
            wmg->addFootstep(0.0   , -step_y/2, 0.0, FS_TYPE_SS_R);

            if (!name.empty())
            {
                wmg->FS2file(fs_out_filename, plot_ds); 
            }
        }
};


/**
 * @brief 'Staircase' walking
 */
class init_05 : public test_init_base
{
    public:
        init_05 (const string & test_name, const bool plot_ds_ = true) : 
            test_init_base (test_name, plot_ds_)
        {
            wmg = new WMG (15, 100);
            par = new smpc_parameters (wmg->N, 0.261);


            // each step is defined relatively to the previous step
            double step_x = 0.03;      // relative X position
            double step_y = wmg->def_constraints.support_distance_y;       // relative Y position
            double shift = -0.02;


            // Initial double support
            wmg->user_constraints[0] = 0.09;
            wmg->user_constraints[1] = 0.075;
            wmg->user_constraints[2] = 0.03;
            wmg->user_constraints[3] = 0.025;
            wmg->setFootstepParameters (2, 1, 1, true);
            wmg->addFootstep(0.0, 0.0, 0.0, FS_TYPE_DS);

            // all subsequent steps have normal feet size
            // 2 reference ZMP positions in single support 
            // 1 in double support
            wmg->setFootstepParameters (2, 1, 1, false);
            wmg->addFootstep(0.0   , -step_y/2, 0.0);
            wmg->addFootstep(step_x,  step_y + shift, 0.0);
            wmg->addFootstep(step_x, -step_y + shift, 0.0);
            wmg->addFootstep(step_x,  step_y + shift, 0.0);
            wmg->addFootstep(step_x, -step_y + shift, 0.0);
            wmg->addFootstep(step_x,  step_y + shift, 0.0);
            wmg->addFootstep(step_x, -step_y + shift, 0.0);
            wmg->addFootstep(step_x,  step_y + shift, 0.0);
            wmg->addFootstep(step_x, -step_y + shift, 0.0);
            wmg->addFootstep(step_x,  step_y + shift, 0.0);

            // here we give many reference points, since otherwise we 
            // would not have enough steps in preview window to reach 
            // the last footsteps
            wmg->user_constraints[0] = 0.09;
            wmg->user_constraints[1] = 0.025;
            wmg->user_constraints[2] = 0.03;
            wmg->user_constraints[3] = 0.075;
            wmg->setFootstepParameters (30, 0, 0, true);
            wmg->addFootstep(0.0   , -step_y/2, 0.0, FS_TYPE_DS);

            if (!name.empty())
            {
                wmg->FS2file(fs_out_filename, plot_ds); 
            }
        }
};



/**
 * @brief Side walking
 */
class init_06 : public test_init_base
{
    public:
        init_06 (const string & test_name, const bool plot_ds_ = true) : 
            test_init_base (test_name, plot_ds_)
        {
            wmg = new WMG (15, 100);
            par = new smpc_parameters (wmg->N, 0.261);


            // each step is defined relatively to the previous step
            double step_x = 0.0;      // relative X position
            double step_y = wmg->def_constraints.support_distance_y;       // relative Y position
            double shift = -0.02;


            // Initial double support
            wmg->user_constraints[0] = 0.09;
            wmg->user_constraints[1] = 0.075;
            wmg->user_constraints[2] = 0.03;
            wmg->user_constraints[3] = 0.025;
            wmg->setFootstepParameters (3, 0, 0, true);
            wmg->addFootstep(0.0, 0.0, 0.0, FS_TYPE_DS);

            // all subsequent steps have normal feet size
            // 3 reference ZMP positions in single support 
            // 0 in double support
            wmg->setFootstepParameters (3, 0, 0, false);
            wmg->addFootstep(0.0   , -step_y/2,       0.0);
            wmg->addFootstep(step_x,  step_y + shift, 0.0);
            wmg->addFootstep(step_x, -step_y + shift, 0.0);
            wmg->addFootstep(step_x,  step_y + shift, 0.0);
            wmg->addFootstep(step_x, -step_y + shift, 0.0);
            wmg->addFootstep(step_x,  step_y + shift, 0.0);
            wmg->addFootstep(step_x, -step_y + shift, 0.0);
            wmg->addFootstep(step_x,  step_y + shift, 0.0);
            wmg->addFootstep(step_x, -step_y + shift, 0.0);
            wmg->addFootstep(step_x,  step_y + shift, 0.0);

            // here we give many reference points, since otherwise we 
            // would not have enough steps in preview window to reach 
            // the last footsteps
            wmg->user_constraints[0] = 0.09;
            wmg->user_constraints[1] = 0.025;
            wmg->user_constraints[2] = 0.03;
            wmg->user_constraints[3] = 0.075;
            wmg->setFootstepParameters (30, 0, 0, true);
            wmg->addFootstep(0.0   , -step_y/2, 0.0, FS_TYPE_DS);

            if (!name.empty())
            {
                wmg->FS2file(fs_out_filename, plot_ds); 
            }
        }
};



/**
 * @brief Walk straight
 */
class init_07 : public test_init_base
{
    public:
        init_07 (const string & test_name, const bool plot_ds_ = true) : 
            test_init_base (test_name, plot_ds_)
        {
            wmg = new WMG (40, 40, 0.015);
            par = new smpc_parameters (wmg->N, 0.261);

            // each step is defined relatively to the previous step
            double step_x = 0.035;      // relative X position
            double step_y = wmg->def_constraints.support_distance_y;       // relative Y position


            wmg->setFootstepParameters (0, 0, 0);
            wmg->addFootstep(0.0, step_y/2, 0.0, FS_TYPE_SS_L);

            // Initial double support
            wmg->setFootstepParameters (10, 0, 0);
            wmg->addFootstep(0.0, -step_y/2, 0.0, FS_TYPE_DS);
            // ZMP, CoM are at [0;0]


            // all subsequent steps have normal feet size
            wmg->setFootstepParameters (10, 1, 3);
            wmg->addFootstep(0.0   , -step_y/2, 0.0);
            wmg->addFootstep(step_x,  step_y, 0.0);
            wmg->addFootstep(step_x, -step_y, 0.0);
            wmg->addFootstep(step_x,  step_y, 0.0);
            wmg->addFootstep(step_x, -step_y, 0.0);
            wmg->addFootstep(step_x,  step_y, 0.0);
            wmg->addFootstep(step_x, -step_y, 0.0);
            wmg->addFootstep(step_x,  step_y, 0.0);
            wmg->addFootstep(step_x, -step_y, 0.0);
            wmg->addFootstep(step_x,  step_y, 0.0);

            // here we give many reference points, since otherwise we 
            // would not have enough steps in preview window to reach 
            // the last footsteps
            wmg->setFootstepParameters (60, 0, 0);
            wmg->addFootstep(0.0   , -step_y/2, 0.0, FS_TYPE_DS);
            wmg->setFootstepParameters (0, 0, 0);
            wmg->addFootstep(0.0   , -step_y/2, 0.0, FS_TYPE_SS_R);


            if (!name.empty())
            {
                wmg->FS2file(fs_out_filename, plot_ds); 
            }
        }
};



/**
 * @brief Walk straight
 */
class init_08 : public test_init_base
{
    public:
        init_08 (const string & test_name, const bool plot_ds_ = true) : 
            test_init_base (test_name, plot_ds_)
        {
            wmg = new WMG (25, 60, 0.02);
            par = new smpc_parameters (wmg->N, 0.261);
            int ss_time_ms = 420;
            int ds_time_ms = 60;
            int ds_number = 2;


            // each step is defined relatively to the previous step
            double step_x = 0.04;      // relative X position
            double step_y = wmg->def_constraints.support_distance_y;       // relative Y position


            wmg->setFootstepParametersMS (0, 0, 0);
            wmg->addFootstep(0.0, step_y/2, 0.0, FS_TYPE_SS_L);

            // Initial double support
            wmg->setFootstepParametersMS (3*ss_time_ms, 0, 0);
            wmg->addFootstep(0.0, -step_y/2, 0.0, FS_TYPE_DS);
            // ZMP, CoM are at [0;0]


            // all subsequent steps have normal feet size
            wmg->setFootstepParametersMS (ss_time_ms, ds_time_ms, ds_number);
//            wmg->setFootstepParameters (ss_time_ms, 0, 0);
            wmg->addFootstep(0.0   , -step_y/2, 0.0);
//            wmg->setFootstepParameters (ss_time_ms, ds_time_ms, ds_number);
            wmg->addFootstep(step_x,  step_y, 0.0);
            wmg->addFootstep(step_x, -step_y, 0.0);
            wmg->addFootstep(step_x,  step_y, 0.0);
            wmg->addFootstep(step_x, -step_y, 0.0);
            wmg->addFootstep(step_x,  step_y, 0.0);
            wmg->addFootstep(step_x, -step_y, 0.0);
            wmg->addFootstep(step_x,  step_y, 0.0);
            wmg->addFootstep(step_x, -step_y, 0.0);
            wmg->addFootstep(step_x,  step_y, 0.0);

            // here we give many reference points, since otherwise we 
            // would not have enough steps in preview window to reach 
            // the last footsteps
            wmg->setFootstepParametersMS (5*ss_time_ms, 0, 0);
            wmg->addFootstep(0.0   , -step_y/2, 0.0, FS_TYPE_DS);
            wmg->setFootstepParametersMS (0, 0, 0);
            wmg->addFootstep(0.0   , -step_y/2, 0.0, FS_TYPE_SS_R);

            if (!name.empty())
            {
                wmg->FS2file(fs_out_filename, plot_ds);
            }
        }
};


/**
 * @brief Walk straight
 */
class init_09 : public test_init_base
{
    public:
        init_09 (const string & test_name, const bool plot_ds_ = true) : 
            test_init_base (test_name, plot_ds_)
        {
            wmg = new WMG (
                    40, 
                    40, 
                    0.015,
                    1.0,
                    2.0,
                    0.01,
                    0.008,
                    true); // Use alternative constraints
            par = new smpc_parameters (wmg->N, 0.261);

            // each step is defined relatively to the previous step
            double step_x = 0.035;      // relative X position
            double step_y = wmg->def_constraints.support_distance_y;       // relative Y position
            double z = 5.0*M_PI/180.0;


            wmg->setFootstepParameters (0, 0, 0);
            wmg->addFootstep(0.0, step_y/2, 0.0, FS_TYPE_SS_L);

            // Initial double support
            wmg->setFootstepParameters (10, 0, 0);
            wmg->addFootstep(0.0, -step_y/2, 0.0, FS_TYPE_DS);
            // ZMP, CoM are at [0;0]


            // all subsequent steps have normal feet size
            wmg->setFootstepParameters (10, 1, 3);
            wmg->addFootstep(0.0   , -step_y/2, 0.0);
            wmg->addFootstep(step_x,  step_y, z);
            wmg->addFootstep(step_x, -step_y, z);
            wmg->addFootstep(step_x,  step_y, z);
            wmg->addFootstep(step_x, -step_y, z);
            wmg->addFootstep(step_x,  step_y, z);
            wmg->addFootstep(step_x, -step_y, z);
            wmg->addFootstep(step_x,  step_y, z);
            wmg->addFootstep(step_x, -step_y, z);
            wmg->addFootstep(step_x,  step_y, z);

            // here we give many reference points, since otherwise we 
            // would not have enough steps in preview window to reach 
            // the last footsteps
            wmg->setFootstepParameters (60, 0, 0);
            wmg->addFootstep(0.0   , -step_y/2, 0.0, FS_TYPE_DS);
            wmg->setFootstepParameters (0, 0, 0);
            wmg->addFootstep(0.0   , -step_y/2, 0.0, FS_TYPE_SS_R);


            if (!name.empty())
            {
                wmg->FS2file(fs_out_filename, plot_ds); 
            }
        }
};



/**
 * @brief Walk straight
 */
class init_10 : public test_init_base
{
    public:
        init_10 (const string & test_name, const bool plot_ds_ = true) : 
            test_init_base (test_name, plot_ds_)
        {
            int preview_sampling_time_ms = 40;
            wmg = new WMG (40, preview_sampling_time_ms, 0.02);
            par = new smpc_parameters (wmg->N, 0.252007);
            int ss_time_ms = 400;
            int ds_time_ms = 40;
            int ds_number = 3;


            // each step is defined relatively to the previous step
            double step_x = 0.04;      // relative X position
            double step_y = wmg->def_constraints.support_distance_y;       // relative Y position


            wmg->setFootstepParametersMS (0, 0, 0);
            wmg->addFootstep(0.0, -step_y/2, 0.0, FS_TYPE_SS_R);

            // Initial double support
            wmg->setFootstepParametersMS (3*ss_time_ms, 0, 0);
            wmg->addFootstep(0.0, step_y/2, 0.0, FS_TYPE_DS);


            // all subsequent steps have normal feet size
            wmg->setFootstepParametersMS (ss_time_ms, 0, 0);
            wmg->addFootstep(0.0   ,  step_y/2, 0.0);
            wmg->setFootstepParametersMS (ss_time_ms, ds_time_ms, ds_number);
            wmg->addFootstep(step_x, -step_y  , 0.0);


            for (int i = 0; i < 3; i++)
            {
                wmg->addFootstep(step_x,  step_y, 0.0);
                wmg->addFootstep(step_x, -step_y, 0.0);
            }

            // here we give many reference points, since otherwise we 
            // would not have enough steps in preview window to reach 
            // the last footsteps
            wmg->setFootstepParametersMS (5*ss_time_ms, 0, 0);
            wmg->addFootstep(0.0, step_y/2, 0.0, FS_TYPE_DS);
            wmg->setFootstepParametersMS (0, 0, 0);
            wmg->addFootstep(0.0, step_y/2, 0.0, FS_TYPE_SS_L);

            if (!name.empty())
            {
                wmg->FS2file(fs_out_filename, plot_ds);
            }
        }
};

///@}
