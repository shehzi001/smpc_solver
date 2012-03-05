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
            wmg->setFootstepDefaults (3, 0, d);
            wmg->addFootstep(0.0, 0.05, 0.0, FS_TYPE_DS);

            double z = 5.0*M_PI/180.0;
            double step_x = 0.035;
            double step_y = 0.1;

            wmg->setFootstepDefaults (4, 0, wmg->def_ss_constraint);
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
            wmg->setFootstepDefaults (30, 0);
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

            double d[4] = {0.03 , 0.01, 0.01, 0.11};
            wmg->setFootstepDefaults (3, 0, d);
            wmg->addFootstep(0.0, 0.05, 0.0);
            
            double z = 5.0*M_PI/180.0;
            double step_x = 0.035;
            double step_y = 0.1;

            // use this for smaller feet (and we will have more active constraints)
            d[0] = 0.03; d[1] = 0.01; d[2] = 0.01; d[3] = 0.01;

            wmg->setFootstepDefaults (4, 0, d);
            wmg->addFootstep(0.0   , -step_y, 0.0);
            wmg->addFootstep(step_x,  step_y, z);
            wmg->addFootstep(step_x, -step_y, z);
            wmg->addFootstep(step_x,  step_y, z);
            wmg->addFootstep(step_x, -step_y, z);
            wmg->addFootstep(step_x,  step_y, z);
            wmg->addFootstep(step_x, -step_y, z);
            wmg->addFootstep(step_x,  step_y, z);
            wmg->addFootstep(step_x, -step_y, z);
            wmg->setFootstepDefaults (30, 0);
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

            // Initial double support
            double d[4] = {0.09 , 0.075, 0.03, 0.025};
            wmg->setFootstepDefaults (2, 1, d);
            wmg->addFootstep(0.0, 0.0, 0.0, FS_TYPE_DS);
            // ZMP, CoM are at [0;0]

            // each step is defined relatively to the previous step
            double z = 5.0*M_PI/180.0;  // relative angle
            double step_x = 0.035;      // relative X position
            double step_y = 0.1;        // relative Y position

            // all subsequent steps have normal feet size
            // 2 reference ZMP positions in single support 
            // 1 in double support
            wmg->setFootstepDefaults (2, 1, wmg->def_ss_constraint);
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
            d[0] = 0.09;
            d[1] = 0.025;
            d[2] = 0.03;
            d[3] = 0.075;
            wmg->setFootstepDefaults (30, 0, d);
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
            double step_y = 0.1;       // relative Y position
            double d[4];

            wmg->setFootstepDefaults (0, 0, wmg->def_ss_constraint);
            wmg->addFootstep(0.0, step_y/2, 0.0, FS_TYPE_SS_L);

            // Initial double support
            d[0] = 0.09;
            d[1] = 0.075;
            d[2] = 0.03;
            d[3] = 0.075;
            wmg->setFootstepDefaults (2, 0, d);
            wmg->addFootstep(0.0, -step_y/2, 0.0, FS_TYPE_DS);
            // ZMP, CoM are at [0;0]


            // all subsequent steps have normal feet size
            wmg->setFootstepDefaults (6, 2, wmg->def_ss_constraint);
            wmg->addFootstep(0.0   , -step_y/2, 0.0);
            wmg->addFootstep(step_x,  step_y, 0.0);
            wmg->addFootstep(step_x, -step_y, 0.0);
            wmg->addFootstep(step_x,  step_y, 0.0);
            wmg->addFootstep(step_x, -step_y, 0.0);
            wmg->addFootstep(step_x,  step_y, 0.0);

            // here we give many reference points, since otherwise we 
            // would not have enough steps in preview window to reach 
            // the last footsteps
            d[0] = 0.09;
            d[1] = 0.075;
            d[2] = 0.03;
            d[3] = 0.075;
            wmg->setFootstepDefaults (30, 0, d);
            wmg->addFootstep(0.0   , -step_y/2, 0.0, FS_TYPE_DS);
            wmg->setFootstepDefaults (0, 0, wmg->def_ss_constraint);
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

            // Initial double support
            double d[4] = {0.09 , 0.075, 0.03, 0.025};
            wmg->setFootstepDefaults (2, 1, d);
            wmg->addFootstep(0.0, 0.0, 0.0, FS_TYPE_DS);
            // ZMP, CoM are at [0;0]

            // each step is defined relatively to the previous step
            double step_x = 0.03;      // relative X position
            double step_y = 0.1;       // relative Y position
            double shift = -0.02;

            // all subsequent steps have normal feet size
            // 2 reference ZMP positions in single support 
            // 1 in double support
            wmg->setFootstepDefaults (2, 1, wmg->def_ss_constraint);
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
            d[0] = 0.09;
            d[1] = 0.025;
            d[2] = 0.03;
            d[3] = 0.075;
            wmg->setFootstepDefaults (30, 0, d);
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

            // Initial double support
            double d[4] = {0.09 , 0.075, 0.03, 0.025};
            wmg->setFootstepDefaults (3, 0, d);
            wmg->addFootstep(0.0, 0.0, 0.0, FS_TYPE_DS);
            // ZMP, CoM are at [0;0]

            // each step is defined relatively to the previous step
            double step_x = 0.0;      // relative X position
            double step_y = 0.1;       // relative Y position
            double shift = -0.02;

            // all subsequent steps have normal feet size
            // 3 reference ZMP positions in single support 
            // 0 in double support
            wmg->setFootstepDefaults (3, 0, wmg->def_ss_constraint);
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
            d[0] = 0.09;
            d[1] = 0.025;
            d[2] = 0.03;
            d[3] = 0.075;
            wmg->setFootstepDefaults (30, 0, d);
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
            double step_y = 0.1;       // relative Y position

            double ds_constraint[4] = {
                wmg->def_ss_constraint[0],
                wmg->def_ss_constraint[1] + 0.5*step_y,
                wmg->def_ss_constraint[2],
                wmg->def_ss_constraint[3] + 0.5*step_y};


            wmg->setFootstepDefaults (0, 0, wmg->def_ss_constraint);
            wmg->addFootstep(0.0, step_y/2, 0.0, FS_TYPE_SS_L);

            // Initial double support
            wmg->setFootstepDefaults (10, 0, ds_constraint);
            wmg->addFootstep(0.0, -step_y/2, 0.0, FS_TYPE_DS);
            // ZMP, CoM are at [0;0]


            // all subsequent steps have normal feet size
            wmg->setFootstepDefaults (10, 3, wmg->def_ss_constraint);
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
            wmg->setFootstepDefaults (60, 0, ds_constraint);
            wmg->addFootstep(0.0   , -step_y/2, 0.0, FS_TYPE_DS);
            wmg->setFootstepDefaults (0, 0, wmg->def_ss_constraint);
            wmg->addFootstep(0.0   , -step_y/2, 0.0, FS_TYPE_SS_R);


            if (!name.empty())
            {
                wmg->FS2file(fs_out_filename, plot_ds); 
            }
        }
};

///@}
