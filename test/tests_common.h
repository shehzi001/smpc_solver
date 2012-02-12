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
            double step_y = 0.1;       // relative Y position
            double d[4];

            d[0] = 0.09;
            d[1] = 0.025;
            d[2] = 0.03;
            d[3] = 0.025;
            wmg->AddFootstep(0.0, step_y/2, 0.0, 0, 0, d, FS_TYPE_SS_L);

            // Initial double support
            d[0] = 0.09;
            d[1] = 0.075;
            d[2] = 0.03;
            d[3] = 0.075;
            wmg->AddFootstep(0.0, -step_y/2, 0.0, 2, 2, d, FS_TYPE_DS);
            // ZMP, CoM are at [0;0]


            // all subsequent steps have normal feet size
            d[0] = 0.09;
            d[1] = 0.025;
            d[2] = 0.03;
            d[3] = 0.025;
            // 2 reference ZMP positions in single support 
            // 1 in double support
            // 1 + 2 = 3
            wmg->AddFootstep(0.0   , -step_y/2, 0.0 , 6,  8, d);
            wmg->AddFootstep(step_x,  step_y, 0.0);
            wmg->AddFootstep(step_x, -step_y, 0.0);
            wmg->AddFootstep(step_x,  step_y, 0.0);
            wmg->AddFootstep(step_x, -step_y, 0.0);
            wmg->AddFootstep(step_x,  step_y, 0.0);

            // here we give many reference points, since otherwise we 
            // would not have enough steps in preview window to reach 
            // the last footsteps
            d[0] = 0.09;
            d[1] = 0.075;
            d[2] = 0.03;
            d[3] = 0.075;
            wmg->AddFootstep(0.0   , -step_y/2, 0.0, 30, 30, d, FS_TYPE_DS);
            d[0] = 0.09;
            d[1] = 0.025;
            d[2] = 0.03;
            d[3] = 0.025;
            wmg->AddFootstep(0.0   , -step_y/2, 0.0 , 0,  0, d, FS_TYPE_SS_R);

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

            // Initial double support
            double d[4] = {0.09 , 0.075, 0.03, 0.025};
            wmg->AddFootstep(0.0, 0.0, 0.0, 2, 3, d, FS_TYPE_DS);
            // ZMP, CoM are at [0;0]

            // each step is defined relatively to the previous step
            double step_x = 0.03;      // relative X position
            double step_y = 0.1;       // relative Y position
            double shift = -0.02;

            // all subsequent steps have normal feet size
            d[0] = 0.09;
            d[1] = 0.025;
            d[2] = 0.03;
            d[3] = 0.025;
            // 2 reference ZMP positions in single support 
            // 1 in double support
            // 1 + 2 = 3
            wmg->AddFootstep(0.0   , -step_y/2, 0.0 , 2,  3, d);
            wmg->AddFootstep(step_x,  step_y + shift, 0.0);
            wmg->AddFootstep(step_x, -step_y + shift, 0.0);
            wmg->AddFootstep(step_x,  step_y + shift, 0.0);
            wmg->AddFootstep(step_x, -step_y + shift, 0.0);
            wmg->AddFootstep(step_x,  step_y + shift, 0.0);
            wmg->AddFootstep(step_x, -step_y + shift, 0.0);
            wmg->AddFootstep(step_x,  step_y + shift, 0.0);
            wmg->AddFootstep(step_x, -step_y + shift, 0.0);
            wmg->AddFootstep(step_x,  step_y + shift, 0.0);

            // here we give many reference points, since otherwise we 
            // would not have enough steps in preview window to reach 
            // the last footsteps
            d[0] = 0.09;
            d[1] = 0.025;
            d[2] = 0.03;
            d[3] = 0.075;
            wmg->AddFootstep(0.0   , -step_y/2, 0.0, 30, 30, d, FS_TYPE_DS);

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

            // Initial double support
            double d[4] = {0.09 , 0.075, 0.03, 0.025};
            wmg->AddFootstep(0.0, 0.0, 0.0, 3, 3, d, FS_TYPE_DS);
            // ZMP, CoM are at [0;0]

            // each step is defined relatively to the previous step
            double step_x = 0.0;      // relative X position
            double step_y = 0.1;       // relative Y position
            double shift = -0.02;

            // all subsequent steps have normal feet size
            d[0] = 0.09;
            d[1] = 0.025;
            d[2] = 0.03;
            d[3] = 0.025;
            // 3 reference ZMP positions in single support 
            // 0 in double support
            // 0 + 3 = 3
            wmg->AddFootstep(0.0   , -step_y/2, 0.0 , 3,  3, d);
            wmg->AddFootstep(step_x,  step_y + shift, 0.0);
            wmg->AddFootstep(step_x, -step_y + shift, 0.0);
            wmg->AddFootstep(step_x,  step_y + shift, 0.0);
            wmg->AddFootstep(step_x, -step_y + shift, 0.0);
            wmg->AddFootstep(step_x,  step_y + shift, 0.0);
            wmg->AddFootstep(step_x, -step_y + shift, 0.0);
            wmg->AddFootstep(step_x,  step_y + shift, 0.0);
            wmg->AddFootstep(step_x, -step_y + shift, 0.0);
            wmg->AddFootstep(step_x,  step_y + shift, 0.0);

            // here we give many reference points, since otherwise we 
            // would not have enough steps in preview window to reach 
            // the last footsteps
            d[0] = 0.09;
            d[1] = 0.025;
            d[2] = 0.03;
            d[3] = 0.075;
            wmg->AddFootstep(0.0   , -step_y/2, 0.0, 30, 30, d, FS_TYPE_DS);

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
            double step_y = 0.1;       // relative Y position

            double ds_constraint[4] = {
                wmg->def_ss_constraint[0],
                wmg->def_ss_constraint[1] + 0.5*step_y,
                wmg->def_ss_constraint[2],
                wmg->def_ss_constraint[3] + 0.5*step_y};


            wmg->AddFootstep(0.0, step_y/2, 0.0, 0, 0, wmg->def_ss_constraint, FS_TYPE_SS_L);

            // Initial double support
            wmg->AddFootstep(0.0, -step_y/2, 0.0, 10, 10, ds_constraint, FS_TYPE_DS);
            // ZMP, CoM are at [0;0]


            // all subsequent steps have normal feet size
            // 2 reference ZMP positions in single support 
            // 1 in double support
            // 1 + 2 = 3
            wmg->AddFootstep(0.0   , -step_y/2, 0.0 , 10,  13, wmg->def_ss_constraint);
            wmg->AddFootstep(step_x,  step_y, 0.0);
            wmg->AddFootstep(step_x, -step_y, 0.0);
            wmg->AddFootstep(step_x,  step_y, 0.0);
            wmg->AddFootstep(step_x, -step_y, 0.0);
            wmg->AddFootstep(step_x,  step_y, 0.0);
            wmg->AddFootstep(step_x, -step_y, 0.0);
            wmg->AddFootstep(step_x,  step_y, 0.0);
            wmg->AddFootstep(step_x, -step_y, 0.0);
            wmg->AddFootstep(step_x,  step_y, 0.0);

            // here we give many reference points, since otherwise we 
            // would not have enough steps in preview window to reach 
            // the last footsteps
            wmg->AddFootstep(0.0   , -step_y/2, 0.0, 60, 60, ds_constraint, FS_TYPE_DS);
            wmg->AddFootstep(0.0   , -step_y/2, 0.0 , 0,  0, wmg->def_ss_constraint, FS_TYPE_SS_R);


            if (!name.empty())
            {
                wmg->FS2file(fs_out_filename, plot_ds); 
            }
        }
};

///@}
