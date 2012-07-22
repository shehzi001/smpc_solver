/** 
 * @file
 * @author Alexander Sherikov
 * @date 19.07.2011 22:30:13 MSD
 */



/****************************************
 * INCLUDES 
 ****************************************/

#include "as_chol_solve.h"

#include <cmath> // sqrt
#include <cstring> // memset, memmove


/****************************************
 * FUNCTIONS 
 ****************************************/
namespace AS
{
    //==============================================
    // constructors / destructors

    /**
     * @brief Constructor
     *
     * @param[in] N size of the preview window.
     */
    chol_solve::chol_solve (const int N) : ecL(N)
    {
        nu = new double[SMPC_NUM_VAR*N];
        z = new double[SMPC_NUM_VAR*N];

        icL = new double*[N*2];
        icL_mem = new double[SMPC_NUM_VAR*N*N*2];
        for(int i = 0; i < N*2; ++i)
        {
            icL[i] = &icL_mem[i * SMPC_NUM_VAR*N];
        }
    }


    /**
     * @brief Destructor
     */
    chol_solve::~chol_solve()
    {
        if (icL != NULL)
        {
            delete icL_mem;
            delete icL;
        }
        if (z != NULL)
        {
            delete z;
        }
        if (nu != NULL)
        {
            delete nu;
        }
    }
    //==============================================


    /**
     * @brief Forms row vector 's_a' (@ref pCholUp).
     *
     * @param[in] ppar parameters
     * @param[in] c activated constraint
     * @param[in] ic_len the length of new row in L
     * @param[out] row 's_a' row
     */
    void chol_solve::form_sa_row(
            const problem_parameters& ppar, 
            const constraint& c,
            const int ic_len, 
            double *row)
    {
        double i2H = ppar.i2Q[0]; // a'*inv(H) = a'*inv(H)*a
        int first_num = c.ind; // first !=0 element


        // reset memory
        memset(row, 0, ic_len * sizeof(double));


        // a * iH * -I
        row[first_num] = -i2H * c.coef_x;
        row[first_num + 3] = -i2H * c.coef_y;

        if (first_num/SMPC_NUM_STATE_VAR != ppar.N-1)
        {
            // a * iH * A'
            row[first_num + 6] = i2H * c.coef_x;
            row[first_num + 9] = i2H * c.coef_y;
        }

        // initialize the last element in the row
        // coef_x^2 + coef_y^2 = 1, since they are row in a rotation matrix
        row[ic_len] = i2H;
    }



    /**
     * @brief Determines feasible descent direction.
     *
     * @param[in] ppar   parameters.
     * @param[in] x    initial guess.
     * @param[out] dx   feasible descent direction, must be allocated.
     */
    void chol_solve::solve(
            const problem_parameters& ppar, 
            const double *x, 
            double *dx)
    {
        double *s_nu = nu;
        int i;


        // generate L
        ecL.form (ppar);

        // obtain s = E * x;
        E.form_Ex (ppar, x, s_nu);

        // obtain nu
        ecL.solve_forward(ppar.N, s_nu);
        // make copy of z - it is constant
        for (i = 0; i < ppar.N * SMPC_NUM_STATE_VAR; i += SMPC_NUM_STATE_VAR)
        {
            s_nu[i]   = -s_nu[i];
            s_nu[i+1] = -s_nu[i+1];
            s_nu[i+2] = -s_nu[i+2];
            s_nu[i+3] = -s_nu[i+3];
            s_nu[i+4] = -s_nu[i+4];
            s_nu[i+5] = -s_nu[i+5];
        }
        memmove(z, s_nu, sizeof(double) * ppar.N * SMPC_NUM_STATE_VAR);
        ecL.solve_backward(ppar.N, s_nu);

        // - i2H * E' * nu
        E.form_i2HETx (ppar, s_nu, dx);

        
        // dx = -(x + inv(H) * E' * nu)
        //            ~~~~~~~~~~~~~~~~
        // dx = -(x +       dx        ) 
        for (i = 0; i < ppar.N*SMPC_NUM_VAR; i += SMPC_NUM_VAR)
        {
            // dx for state variables
            dx[i]   -= x[i]  ;
            dx[i+1] -= x[i+1];
            dx[i+2] -= x[i+2];
            dx[i+3] -= x[i+3];
            dx[i+4] -= x[i+4];
            dx[i+5] -= x[i+5];
            dx[i+6] -= x[i+6];
            dx[i+7] -= x[i+7];
        }
    }


    /**
     * @brief A wrapper around private functions, which update Cholesky factor and 
     *  resolve the system.
     *
     * @param[in] ppar   parameters.
     * @param[in] active_set a vector of active constraints.
     * @param[in] x     initial guess.
     * @param[out] dx   feasible descent direction, must be allocated.
     */
    void chol_solve::up_resolve(
            const AS::problem_parameters& ppar, 
            const vector <AS::constraint>& active_set, 
            const double *x, 
            double *dx)
    {
        int ic_num = active_set.size()-1;
        constraint c = active_set.back();

        update (ppar, c, ic_num);
        update_z (ppar, c, ic_num, x);
        resolve (ppar, active_set, x, dx);
    }


    /**
     * @brief Adds a row corresponding to some inequality constraint to L, see
     * '@ref pCholUpAlg'.
     *
     * @param[in] ppar parameters.
     * @param[in] c activated constraint
     * @param[in] ic_num index of added constraint in W
     */
    void chol_solve::update (
            const problem_parameters& ppar, 
            const constraint& c,
            const int ic_num) 
    {
        int i, j, k;

        double *new_row = icL[ic_num]; // current row in icL
        // trailing elements of new_row corresponding to active constraints
        double *new_row_end = &new_row[ppar.N*SMPC_NUM_STATE_VAR]; 

        int last_num = ic_num + ppar.N*SMPC_NUM_STATE_VAR; // the last !=0 element


        // form row 'a' in the current row of icL
        form_sa_row(ppar, c, last_num, new_row);

        // Forward substitution using L for equality constraints
        ecL.solve_forward(ppar.N, new_row, c.cind/2);


        // update the trailing elements of new_row using the
        // elements computed using forward substitution above
        for(i = c.cind/2 * SMPC_NUM_STATE_VAR; 
            i < ppar.N * SMPC_NUM_STATE_VAR; 
            i += SMPC_NUM_STATE_VAR)
        {
            // make a copy for faster computations
            double tmp_copy_el[6] = {new_row[i], 
                                     new_row[i+1], 
                                     new_row[i+2], 
                                     new_row[i+3], 
                                     new_row[i+4], 
                                     new_row[i+5]};

            // update the last (diagonal) number in the row
            new_row[last_num] -= tmp_copy_el[0] * tmp_copy_el[0] 
                               + tmp_copy_el[1] * tmp_copy_el[1] 
                               + tmp_copy_el[2] * tmp_copy_el[2]
                               + tmp_copy_el[3] * tmp_copy_el[3] 
                               + tmp_copy_el[4] * tmp_copy_el[4] 
                               + tmp_copy_el[5] * tmp_copy_el[5];

            // update elements after N*SMPC_NUM_STATE_VAR using the previously added rows
            // in icL
            for (j = 0; j < ic_num; ++j)
            {
                new_row_end[j] -= tmp_copy_el[0] * icL[j][i]
                                + tmp_copy_el[1] * icL[j][i+1]
                                + tmp_copy_el[2] * icL[j][i+2] 
                                + tmp_copy_el[3] * icL[j][i+3]
                                + tmp_copy_el[4] * icL[j][i+4]
                                + tmp_copy_el[5] * icL[j][i+5];
            }
        }


        // update elements in the end of icL
        for(i = SMPC_NUM_STATE_VAR * ppar.N, k = 0; i < last_num; ++i, ++k)
        {
            new_row[i] /= icL[k][i];
            double tmp_copy_el = new_row[i];

            // determine number in the row of L

            // update the last (diagonal) number in the row
            new_row[last_num] -= tmp_copy_el * tmp_copy_el;

            for (j = k+1; j < ic_num; ++j)
            {
                new_row_end[j] -= tmp_copy_el * icL[j][i];
            }
        }

        // square root of the diagonal element
        new_row[last_num] = sqrt(new_row[last_num]);
    }



    /**
     * @brief Adjust vector '@ref pz "z"' after update.
     *
     * @param[in] ppar  parameters.
     * @param[in] c     activated constraint
     * @param[in] ic_num number of added constraint in the active set
     * @param[in] x     initial guess.
     */
    void chol_solve::update_z (
            const problem_parameters& ppar, 
            const constraint& c,
            const int ic_num, 
            const double *x)
    {
        // update lagrange multipliers
        const int zind = ppar.N*SMPC_NUM_STATE_VAR + ic_num;
        // sn
        const int first_num = c.ind; // first !=0 element

        double zn = -x[first_num] * c.coef_x
                    -x[first_num+3] * c.coef_y;

        // zn
        memmove (nu, z, zind * sizeof(double));
        for (int i = first_num; i < zind; ++i)
        {
            zn -= z[i] * icL[ic_num][i];
        }
        nu[zind] = z[zind] = zn/icL[ic_num][zind];
        return;
    }



    /**
     * @brief Determines feasible descent direction with respect to added
     *  inequality constraints.
     *
     * @param[in] ppar   parameters.
     * @param[in] active_set a vector of active constraints.
     * @param[in] x     initial guess.
     * @param[out] dx   feasible descent direction, must be allocated.
     */
    void chol_solve::resolve (
            const AS::problem_parameters& ppar, 
            const vector <AS::constraint>& active_set, 
            const double *x, 
            double *dx)
    {
        int i;
        const int nW = active_set.size();

        // backward substitution for icL
        for (i = nW-1; i >= 0; --i)
        {
            const int last_el_num = i + ppar.N*SMPC_NUM_STATE_VAR;
            nu[last_el_num] /= icL[i][last_el_num];

            for (int j = active_set[i].ind; j < last_el_num; ++j)
            {
                nu[j] -= nu[last_el_num] * icL[i][j];
            }
        }
        // backward substitution for ecL
        ecL.solve_backward(ppar.N, nu);


        // - i2H * E' * nu
        E.form_i2HETx (ppar, nu, dx);

        
        // dx = -(x + inv(H) * E' * nu)
        //            ~~~~~~~~~~~~~~~~
        // dx = -(x +       dx        ) 
        for (i = 0; i < ppar.N*SMPC_NUM_VAR; i += SMPC_NUM_VAR)
        {
            // dx for state variables
            dx[i]   -= x[i]  ;
            dx[i+1] -= x[i+1];
            dx[i+2] -= x[i+2];
            dx[i+3] -= x[i+3];
            dx[i+4] -= x[i+4];
            dx[i+5] -= x[i+5];
            dx[i+6] -= x[i+6];
            dx[i+7] -= x[i+7];
        }

        // -iH * A(W,:)' * lambda
        const double *lambda = get_lambda(ppar);
        const double i2Q0 = ppar.i2Q[0];
        for (i = 0; i < nW; ++i)
        {
            constraint c = active_set[i];
            dx[c.ind]   -= i2Q0 * c.coef_x * lambda[i];
            dx[c.ind+3] -= i2Q0 * c.coef_y * lambda[i];
        }
    }



    /**
     * @brief A wrapper around private functions, which downdate Cholesky factor and 
     *  resolve the system.
     *
     * @param[in] ppar   parameters.
     * @param[in] active_set a vector of active constraints.
     * @param[in] ind_exclude index of excluded constraint.
     * @param[in] x     initial guess.
     * @param[out] dx   feasible descent direction, must be allocated.
     *
     * @note Downdate of vector @ref pz 'z' is described on the page '@ref pRemoveICz'.
     */
    void chol_solve::down_resolve(
            const AS::problem_parameters& ppar, 
            const vector <AS::constraint>& active_set,
            const int ind_exclude, 
            const double *x, 
            double *dx)
    {
        const int nW = active_set.size();
        const int last_el_ind = ppar.N*SMPC_NUM_STATE_VAR + ind_exclude;

        // for each element of z affected by removed constraint
        // find a base that stays the same
        double z_tmp = 0;
        for (int i = nW; i > ind_exclude; --i)
        {
            const int zind = ppar.N*SMPC_NUM_STATE_VAR + i;
            double zn = z[zind] * icL[i][zind];
            z[zind] = z_tmp;

            for (int j = last_el_ind; j < zind; ++j)
            {
                zn += z[j] * icL[i][j];
            }
            z_tmp = zn;
        }
        z[last_el_ind] = z_tmp;


        // downdate L
        downdate (ppar, nW, ind_exclude, x);


        // recompute elements of z
        for (int i = ind_exclude; i < nW; i++)
        {
            int zind = ppar.N*SMPC_NUM_STATE_VAR + i;
            double zn = z[zind];

            // zn
            // start from the first !=0 element
            for (int j = last_el_ind; j < zind; j++)
            {
                zn -= z[j] * icL[i][j];
            }
            z[zind] = zn/icL[i][zind];
        }

        // copy z to nu
        memmove(nu, z, sizeof(double) * (ppar.N*SMPC_NUM_STATE_VAR + nW));

        resolve (ppar, active_set, x, dx);
    }



    /**
     * @return a pointer to the memory where current lambdas are stored.
     * @param[in] ppar parameters
     */
    double * chol_solve::get_lambda(const problem_parameters& ppar)
    {
        return(&nu[SMPC_NUM_STATE_VAR*ppar.N]);
    }




    /**
     * @brief Delete a line from icL, see page '@ref pCholDown'.
     *
     * @param[in] ppar   parameters.
     * @param[in] nW    number of added constrains.
     * @param[in] ind_exclude index of excluded constraint.
     * @param[in] x     initial guess.
     */
    void chol_solve::downdate(
            const problem_parameters& ppar, 
            const int nW, 
            const int ind_exclude, 
            const double *x)
    {
        // Shuffle memory pointers to avoid copying of the data.
        double * downdate_row = icL[ind_exclude];
        for (int i = ind_exclude + 1; i < nW + 1; i++)
        {
            icL[i-1] = icL[i];
        }
        icL[nW] = downdate_row;


        for (int i = ind_exclude; i < nW; i++)
        {
            int el_index = SMPC_NUM_STATE_VAR*ppar.N + i;
            double *cur_el = &icL[i][el_index];
            double x1 = cur_el[0];
            double x2 = cur_el[1];
            double cosT, sinT;


            // Givens rotation matrix
            if (abs(x2) >= abs(x1))
            {
                double t = x1/x2;
                sinT = 1/sqrt(1 + t*t);
                cosT = sinT*t;
            }
            else
            {
                double t = x2/x1;
                cosT = 1/sqrt(1 + t*t);
                sinT = cosT*t;
            }


            // update elements in the current line
            cur_el[0] = cosT*x1 + sinT*x2;
            cur_el[1] = 0;

            // change sign if needed (diagonal elements of Cholesky 
            // decomposition must be positive)
            double sign = copysign(1, cur_el[0]);
            cur_el[0] = fabs(cur_el[0]);

            // update the lines below the current one.
            for (int j = i + 1; j < nW; j++)
            {
                x1 = icL[j][el_index];
                x2 = icL[j][el_index + 1];

                icL[j][el_index] = sign * (cosT*x1 + sinT*x2);
                icL[j][el_index + 1] = -sinT*x1 + cosT*x2;
            }
        }
    }
}
