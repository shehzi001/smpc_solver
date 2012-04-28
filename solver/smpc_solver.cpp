/** 
 * @file
 * @brief The interface class, essentially, a wrapper around #qp_as.
 *
 * @author Alexander Sherikov
 * @date 02.09.2011 00:49:30 MSD
 */



/****************************************
 * INCLUDES 
 ****************************************/

#include "solver_config.h"

#ifdef HAVE_FEENABLEEXCEPT
#include <fenv.h> // feenableexcept()
#endif

#include "qp_as.h"
#include "qp_ip.h"
#include "smpc_solver.h"
#include "state_handling.h"


/****************************************
 * FUNCTIONS 
 ****************************************/
namespace smpc
{
    void enable_fexceptions()
    {
#ifdef HAVE_FEENABLEEXCEPT
        feenableexcept(
                FE_ALL_EXCEPT &
                // ignore precision loss due to rounding
                ~FE_INEXACT);
#endif
    }


    solver_as::solver_as (
                    const int N,
                    const double gain_position, 
                    const double gain_velocity, 
                    const double gain_acceleration,
                    const double gain_jerk, 
                    const double tol,
                    const bool obj_computation_on,
                    const unsigned int max_added_constraints_num,
                    const bool constraint_removal_on)
    {
        qp_sol = new qp_as (
                N, 
                gain_position, gain_velocity, gain_acceleration, gain_jerk, 
                tol, 
                obj_computation_on,
                max_added_constraints_num, constraint_removal_on);
        added_constraints_num = 0;
        removed_constraints_num = 0;
        active_set_size = 0;
    }


    solver_as::~solver_as()
    {
        if (qp_sol != NULL)
        {
            delete qp_sol;
        }
    }




    void solver_as::set_parameters(
            const double* T, const double* h, const double h_initial,
            const double* angle,
            const double* zref_x, const double* zref_y,
            const double* lb, const double* ub)
    {
        if (qp_sol != NULL)
        {
            qp_sol->set_parameters(T, h, h_initial, angle, zref_x, zref_y, lb, ub);
        }
    }


    void solver_as::form_init_fp (
            const double *x_coord,
            const double *y_coord,
            const state_orig &init_state,
            double* X)
    {
        if (qp_sol != NULL)
        {
            qp_sol->form_init_fp (x_coord, y_coord, init_state.state_vector, false, X);
        }
    }


    void solver_as::form_init_fp (
            const double *x_coord,
            const double *y_coord,
            const state_tilde &init_state,
            double* X)
    {
        if (qp_sol != NULL)
        {
            qp_sol->form_init_fp (x_coord, y_coord, init_state.state_vector, true, X);
        }
    }


    void solver_as::solve()
    {
        if (qp_sol != NULL)
        {
            qp_sol->solve (objective_log);
            
            added_constraints_num   = qp_sol->added_constraints_num;
            removed_constraints_num = qp_sol->removed_constraints_num;
            active_set_size         = qp_sol->active_set_size;
        }
    }


    //************************************************************


    void solver_as::get_next_state (state_tilde &s)
    {
        get_state (s, 0);
    }


    void solver_as::get_state (state_tilde &s, const int ind)
    {
        if (qp_sol != NULL)
        {
            int index;
            if (ind >= qp_sol->N)
            {
                index = qp_sol->N - 1;
            }
            else
            {
                index = ind;
            }

            for (int i = 0; i < SMPC_NUM_STATE_VAR; i++)
            {
                s.state_vector[i] = qp_sol->X[index*SMPC_NUM_STATE_VAR + i];
            }
        }
    }


    //************************************************************


    void solver_as::get_next_state (state_orig &s)
    {
        get_state (s, 0);
    }


    void solver_as::get_state (state_orig &s, const int ind)
    {
        if (qp_sol != NULL)
        {
            int index;
            if (ind >= qp_sol->N)
            {
                index = qp_sol->N - 1;
            }
            else
            {
                index = ind;
            }

            for (int i = 0; i < SMPC_NUM_STATE_VAR; i++)
            {
                s.state_vector[i] = qp_sol->X[index*SMPC_NUM_STATE_VAR + i];
            }
            state_handling::tilde_to_orig (qp_sol->spar[index].h, s.state_vector);
        }
    }


    //************************************************************


    void solver_as::get_first_controls (control &c)
    {
        get_controls (c, 0);
    }


    void solver_as::get_controls (control &c, const int ind)
    {
        if (qp_sol != NULL)
        {
            state_handling::get_controls (
                    qp_sol->N, 
                    qp_sol->X, 
                    ind, 
                    c.control_vector);
        }
    }


//************************************************************
//************************************************************
//************************************************************


    solver_ip::solver_ip (
                    const int N,
                    const unsigned int max_iter,
                    const double gain_position, const double gain_velocity, const double gain_acceleration,
                    const double gain_jerk, 
                    const double tol, const double tol_out,
                    const double t,
                    const double mu,
                    const double bs_alpha, const double bs_beta,
                    const bool obj_computation_on,
                    const bool backtracking_search_on)
    {
        qp_sol = new qp_ip (
                N, 
                gain_position, gain_velocity, gain_acceleration, gain_jerk, 
                tol, 
                obj_computation_on, backtracking_search_on);
        qp_sol->set_ip_parameters (t, mu, bs_alpha, bs_beta, max_iter, tol_out);

        int_loop_iterations = 0;
        ext_loop_iterations = 0;
        bt_search_iterations = 0;
    }


    solver_ip::~solver_ip()
    {
        if (qp_sol != NULL)
        {
            delete qp_sol;
        }
    }


    void solver_ip::set_parameters(
            const double* T, const double* h, const double h_initial,
            const double* angle,
            const double* zref_x, const double* zref_y,
            const double* lb, const double* ub)
    {
        if (qp_sol != NULL)
        {
            qp_sol->set_parameters(T, h, h_initial, angle, zref_x, zref_y, lb, ub);
        }
    }



    void solver_ip::form_init_fp (
            const double *x_coord,
            const double *y_coord,
            const state_orig &init_state,
            double* X)
    {
        if (qp_sol != NULL)
        {
            qp_sol->form_init_fp (x_coord, y_coord, init_state.state_vector, false, X);
        }
    }

    void solver_ip::form_init_fp (
            const double *x_coord,
            const double *y_coord,
            const state_tilde &init_state,
            double* X)
    {
        if (qp_sol != NULL)
        {
            qp_sol->form_init_fp (x_coord, y_coord, init_state.state_vector, true, X);
        }
    }



    void solver_ip::solve()
    {
        if (qp_sol != NULL)
        {
            qp_sol->solve (objective_log);

            int_loop_iterations = qp_sol->int_loop_counter;
            ext_loop_iterations = qp_sol->ext_loop_counter;
            bt_search_iterations = qp_sol->bs_counter;
        }
    }


    //************************************************************


    void solver_ip::get_next_state (state_tilde &s)
    {
        get_state (s, 0);
    }


    void solver_ip::get_state (state_tilde &s, const int ind)
    {
        if (qp_sol != NULL)
        {
            int index;
            if (ind >= qp_sol->N)
            {
                index = qp_sol->N - 1;
            }
            else
            {
                index = ind;
            }

            for (int i = 0; i < SMPC_NUM_STATE_VAR; i++)
            {
                s.state_vector[i] = qp_sol->X[index*SMPC_NUM_STATE_VAR + i];
            }
            state_handling::bar_to_tilde (
                    qp_sol->spar[index].sin, 
                    qp_sol->spar[index].cos, 
                    s.state_vector);
        }
    }


    //************************************************************


    void solver_ip::get_next_state (state_orig &s)
    {
        get_state (s, 0);
    }


    void solver_ip::get_state (state_orig &s, const int ind)
    {
        if (qp_sol != NULL)
        {
            int index;
            if (ind >= qp_sol->N)
            {
                index = qp_sol->N - 1;
            }
            else
            {
                index = ind;
            }

            for (int i = 0; i < SMPC_NUM_STATE_VAR; i++)
            {
                s.state_vector[i] = qp_sol->X[index*SMPC_NUM_STATE_VAR + i];
            }
            state_handling::bar_to_tilde (
                    qp_sol->spar[index].sin, 
                    qp_sol->spar[index].cos, 
                    s.state_vector);
            state_handling::tilde_to_orig (qp_sol->spar[index].h, s.state_vector);
        }
    }


    //************************************************************


    void solver_ip::get_first_controls (control &c)
    {
        get_controls (c, 0);
    }


    void solver_ip::get_controls (control &c, const int ind)
    {
        if (qp_sol != NULL)
        {
            state_handling::get_controls (
                    qp_sol->N,
                    qp_sol->X,
                    ind,
                    c.control_vector);
        }
    }


//************************************************************
//************************************************************
//************************************************************


    control::control()
    {
        control_vector[0] = 0.0;
        control_vector[1] = 0.0;
    }


    state::state()
    {
        set (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    }


    void state::set (double x_, double y_)
    {
        set (x_, 0.0, 0.0, y_, 0.0, 0.0);
    }


    void state::set (
            double x_, double vx_, double ax_,
            double y_, double vy_, double ay_)
    {
        state_vector[0] = x_;
        state_vector[1] = vx_;
        state_vector[2] = ax_;
        state_vector[3] = y_;
        state_vector[4] = vy_;
        state_vector[5] = ay_;
    }
}
