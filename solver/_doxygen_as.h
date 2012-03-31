/**
 * @file _doxygen_as.h
 * @brief This file contains only doxygen definitions and is not to be 
 *  included anywhere.
 *
 * @author Alexander Sherikov
 * @date 15.09.2011 13:24:49 MSD
 */


#ifndef DOXYGEN_AS_H
#define DOXYGEN_AS_H

/**
 * @page pAddIC Adding inequality constraints
 
    This page describes changes in @ref pKKT "KKT system" after addition of
    inequality constraint to the active set.

@section pCholUp Update of Cholesky factor

    Let 
    @f$\mbm{a}_i^T@f$ 
    be the normal to the i-th
    inequality constraint. It contains two nonzero elements from a row in
    the respective rotation matrix. Define the matrix 

    @f$
    \mbm{C} = \left[\begin{array}{c} \mbm{E} \\ \mbm{A}_{W}\end{array}\right]
    @f$

    where the rows of 
    @f$\mbm{A}_W@f$ 
    contain the normals to the inequality constraints in the working set.

    If a new constraint must be added to the active set, then the Schur 
    complement must be updated:

    @f$
    \frac{1}{2}
    \left[
    \begin{array}{c}
        \mbm{C}\\
        \mbm{a}^T\\
    \end{array}
    \right] \mbm{H}^{-1} \left[\mbm{C}^T \quad \mbm{a}\right] = 
    \frac{1}{2}
    \left[
    \begin{array}{cc}
        \mbm{C} \mbm{H}^{-1} \mbm{C}^T      & \mbm{C} \mbm{H}^{-1} \mbm{a}\\
        \mbm{a}^T \mbm{H}^{-1} \mbm{C}^T    & \mbm{a}^T \mbm{H}^{-1} \mbm{a}\\
    \end{array}
    \right]
    @f$

    In a general case the last line is

    @f$
    \mbm{s_a}^T = 
    \frac{1}{2}
    \left[
        \mbm{a}^T \mbm{H}^{-1} \mbm{E}^T \quad
        \mbm{a}^T \mbm{H}^{-1} \mbm{A}_W^T \quad
        \mbm{a}^T \mbm{H}^{-1} \mbm{a}
    \right]
    @f$

    Note that 
    @f$
    \frac{1}{2}
    \mbm{a}^T \mbm{H}^{-1} \mbm{A}_W^T
    @f$
    is a vector of zeros since all normals are orthogonal.
    While 
    @f$
    \frac{1}{2}
    \mbm{a}^T \mbm{H}^{-1} \mbm{a} = \frac{1}{2\beta}
    @f$
    is a scalar.

    The last part 
    @f$
    \frac{1}{2}
    \mbm{a}^T \mbm{H}^{-1} \mbm{E}^T 
    @f$
    has at most 4 non-zero elements.

    The total number of non-zero elements in the new row of Schur complement
    is 5 or 3 (for the last state in the preview window).
\n\n

@subsection pCholUpAlg Algorithm of Cholesky factor update
@verbatim
Input:
    m_e % the number of equality constraints
    m_a % the current cardinality of the active set
    L   % Cholesky factor
    s_a % a row added to the Schur complement

Output:
    l   % a new (the last) row of L

    l = s_a
    first   % the index of the first !=0 element of s_a
    end = m_e + m_a + 1 % the index of the last element of s_a

    for i = l_s:m_e
        l(i) = l(i) / L(i,i)
        l(end) = l(end) - l(i)^2

        % Since ecL is sparse, no more than three subsequent elements with
        % (known) indexes 'k' <= 'end' in 'l' must be updated: 
        l(k) = l(k) - l(i) * L(k,i)

        for j = m_e+1:end-1
            l(j) = l(j) - l(i) * L_(j,i)
        end
    end

    for i = m_e+1:end-1
        l(i) = l(i) / L(i,i)
        l(end) = l(end) - l(i)^2

        for j = i+1:end-1
            l(j) = l(j) - l(i) * L(j,i)
        end
    end
    l(end) = sqrt(l(end))
@endverbatim
\n\n


@section pAddICz Update of z
    @anchor pz
    After Cholesky decomposition we get
    @f$
    \mbm{S}\mbm{\nu} = \mbm{L} \mbm{L}^T \mbm{\nu} = \mbm{L} \mbm{z} = \mbm{s}
    @f$

    When a constraint is added to the active set, there is no need to 
    perform full forward substitution in order to form 
    @f$\mbm{\nu}@f$ (but full backward substitution is still required).

    @f$
      \mbm{s}^{+} = -\left[\begin{array}{c} \mbm{C} \\ \mbm{a}_i^T\end{array}\right]
      \left((\mbm{x}+\alpha\Delta\mbm{x})+\mbm{H}^{-1}\mbm{g}\right) 
      = -\left[\begin{array}{c} \mbm{C}\mbm{x} + \mbm{C}\mbm{H}^{-1}\mbm{g} \\ 
          \mbm{a}_i^T(\mbm{x}+\alpha\Delta\mbm{x}) + \mbm{a}_i^T\mbm{H}^{-1}\mbm{g}\end{array}\right]
      = \left[\begin{array}{c} \mbm{s} \\ s_n \end{array}\right]. 
    @f$

    Note that 
    @f$\alpha\mbm{C}\Delta\mbm{x} = \mbm{0}@f$, 
    because 
    @f$\Delta\mbm{x}@f$ 
    is in the null space of the normals to the active constraints (stored in 
    @f$\mbm{C}@f$), 
    @f$\mbm{H}^{-1}\mbm{g}@f$ 
    is constant, and 
    @f$\mbm{x}+\alpha\Delta\mbm{x}@f$ 
    is already formed).
    \n\n


    Now consider the forward substitution 

    @f$
      \underbrace{\left[\begin{array}{cc} \mbm{L} & \mbm{0} \\ \mbm{l}^T & \ell \end{array}\right]}_{\mbm{L}^{+}}
      \underbrace{\left[\begin{array}{c} \mbm{z} \\ z_n \end{array}\right]}_{\mbm{z}^{+}} = 
      \underbrace{\left[\begin{array}{c} \mbm{s} \\ s_n \end{array}\right]}_{\mbm{s}^{+}}, 
    @f$

    where 
    @f$\left[\begin{array}{cc} \mbm{l}^T  \ell \end{array}\right]@f$ 
    is an appended row.
    These are two equations

    @f$\\
    \mbm{L}\mbm{z} = \mbm{s}  \\
    \mbm{l}^T\mbm{z} + bz_n = s_n 
    @f$

    From the second one we can compute (note that @f$\ell\neq0@f$)

    @f$
    z_n = \frac{s_n - \mbm{l}^T\mbm{z}}{\ell}.
    @f$

    Hence, forming @f$\mbm{z}^{+}@f$ amounts to performing one dot product.
 */

/**
 * @page pRemoveIC Removing inequality constraints
 *
 *  This page describes changes in @ref pKKT "KKT system" after removal of
 *  inequality constraint from the active set.
 *  \n\n

@section pCholDown Downdate of Cholesky factor
    Imagine, that we have selected an inequality constraint for removal,
    then corresponding line and column must be removed from matrix
    @f$\mbm{S}@f$. We can represent this by moving these lines to the end 
    and to the right of the matrix using permutation matrix:

    @f$
    \mbm{E}_{perm} \mbm{S} \mbm{E}_{perm}^T = 
    \mbm{E}_{perm} \mbm{L} \mbm{L}^T \mbm{E}_{perm}^T = 
    (\mbm{E}_{perm} \mbm{L}) (\mbm{E}_{perm} \mbm{L})^T
    @f$

    But matrix @f$ \mbm{E}_{perm} \mbm{L} @f$
    is not lower triangular. We can transform it using Givens rotation 
    matrices (which, obviously, cancel out):
@verbatim
    .                 .                     .                 .
    ..                ..                    ..                ..
    ...  -> remove -> .... <  -> Givens  -> ...  -> Givens -> ...
    ....     line     .....<     rotation   .....   rotation  ....
    .....               ^^^
                      these elements 
     L                must be adjusted
@endverbatim
    The rotations are performed in the following way:
    -# from two columns, which must be updated, topmost 2x2 matrix is taken;
    -# the rotation matrix is formed in such a way, that this 2x2 matrix
    becomes diagonal;
    -# the selected two columns are multiplied by the rotation matrix.

    The Cholesky decomposition is unique: given a positive-definite matrix A,
    there is only one lower triangular matrix L with strictly positive
    diagonal entries such that A = L*L'. The algorithm presented above can 
    produce L with negative diagonal entries. If a negative diagonal element
    is found, the sign of the whole column must be altered after rotation.
\n\n


@section pRemoveICz Downdate of z
    All elements of @f$\mbm{z}@f$ starting from the position corresponding 
    to the last non-zero (diagonal) element of removed row must be updated.

    Consider the following situation (based on formulas derived in section 
    '@ref pAddICz'):

    @f$
      \left[
        \begin{array}{cccc} 
            \mbm{L} & \mbm{0} & \mbm{0} & \mbm{0} \\
            & \mbox{removed line} &&\\
            & \mbox{ignored lines} &&\\
            \mbm{l}^T & \ell_{rem} & \mbm{l}^T_u& \ell_n
        \end{array}
      \right]
      \left[\begin{array}{c} \mbm{z} \\ z_{rem} \\ \mbm{z}_u \\ z_n \end{array}\right] = 
      \mbm{s}
    @f$

    Some lines are not important right now and they are marked as 'ignored'.
    Elements of
    @f$\mbm{s}@f$
    are computed on demand and there is no need to alter it.
    @f$\mbm{l}^T_u, \mbm{z}_u@f$
    must be updated.
    Vector
    @f$\mbm{l}^T@f$
    is not affected by update (see section '@ref pCholDown'). Hence the following
    number stays constant:

    @f$
    z_n \ell_n + \ell_{rem} z_{rem} + \mbm{l}^T_u \mbm{z}_u = s_n - \mbm{l}^T\mbm{z} = z_{n,const}
    @f$

    
    After Cholesky factor was updated we get:

    @f$
      \left[
        \begin{array}{ccc} 
            \mbm{L} & \mbm{0} & \mbm{0} \\
            & \mbox{ignored lines} &\\
            \mbm{l}^T & \mbm{l}^T_{u,new}& \ell_{n,new}
        \end{array}
      \right]
      \left[\begin{array}{c} \mbm{z} \\ \mbm{z}_{u,new} \\ z_{n,new} \end{array}\right] = 
      \mbm{s}
    @f$
    
    The new element of vector z can be computed:

    @f[
    z_{n,new} = \frac{z_{n,const} - \mbm{l}^T_{u,new} \mbm{z}_{u,new} }{\ell_{new}}
    @f]

    Note, that computation of new elements @f$z_{n,new}@f$ must be started from the
    element corresponding to the first changed line of L and continued towards the
    end of z.
 */


/**
 * @page pDetails Notes on implementation

@section pDetMatrices Representation of matrices
    Note, that parts of state and control matrices (@ref pPDModel) corresponding to x 
    and y coordinates are identical. This property is preserved in all subsequent
    transformations, hence there is no need to compute and store all elements of 
    static matrices (Cholesky factor for example). Consequently, we mainly operate on
    3x3 matrices, which are stored in the following way:
@verbatim
  0   3   6
  1   4   7
  2   5   8
@endverbatim


@subsection pDetCholesky Cholesky factor
    The Cholesky factor consists of two parts:
    - A constant and well structured 'upper' part, which corresponds to
      equality constraints:
@verbatim
           full L                           |   compressed L
            a                               |
            0   a                           |           a
            b   0   c                   0   |           b   c
            0   b   0   c                   |   ===>        d   e
                    d   0   e               |   ===>            f   g
                    0   d   0   e           |                       ...
                0           f   0   g       |
                            0   f   0   g   |
                                    ....    |
@endverbatim
      Here a,c,e,g - lower diagonal 3x3 matrices; b,d,f - upper triangular
      3x3 matrices. This matrix is reffered as ecL. It is stored as a 
      sequence of vectors of 9 elements:
@verbatim
        0  3  6 
        1  4  7
        2  5  8

        9  12 15  18 21 24
        10 13 16  19 22 25
        11 14 17  20 23 26

                    ...         ...
@endverbatim

    - A less structured part corresponding to inequality constraints. It is 
      altered each time the active set is changed. This part is stored as a
      set of 2*N vectors. The length of each vector could not be longer than
      N*#SMPC_NUM_VAR.
\n\n

@section pBounds Implementing bounds

    Consider the variable

    @f$
    \mbm{x} = \left[\begin{array}{c} x_1 \\ x_2 \\ x_3 \\ x_4 \\ x_5 \\ x_6 \end{array}\right].
    @f$
     
    Suppose that variables @f$x_1@f$ and @f$x_4@f$ are constrained (and the rest of the variables are are
    not subject to inequality constraints), i.e., 

    @f$\\
      lb_1 \leq a_1 x_1 + b_1 x_2 \leq ub_1  \\
      lb_2 \leq a_2 x_1 + b_2 x_2 \leq ub_2  
    @f$

    Note that both 
    @f$a_1 x_1 + b_1 x_2 \leq ub_1@f$ 
    and 
    @f$a_1 x_1 + b_1 x_2 \geq lb_1@f$ 
    cannot be in the working set at the same time (because if we are on one of the bounds we 
    cannot be on the other one).

    @attention Note, that if we do not distinguish lower and upper bounds, Lagrange 
    multipliers can be negative.
 */

#endif /*DOXYGEN_AS_H*/
 
