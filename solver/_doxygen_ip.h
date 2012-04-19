/**
 * @file _doxygen_ip.h
 * @brief This file contains only doxygen definitions and is not to be 
 *  included anywhere.
 *
 * @author Alexander Sherikov
 * @date 15.09.2011 13:24:49 MSD
 */


#ifndef DOXYGEN_IP_H
#define DOXYGEN_IP_H

/**
 * @page pIP Derivations for a linearly constrained problem

    Objective function:

    @f$
        f(\mbm{x}) = 0.5 \mbm{x}^T \mbm{H} \mbm{x} + \mbm{g}^T \mbm{x}
    @f$

    Constraints:

    @f$\\
        \mbm{E}\mbm{x} = \mbm{e}\\
        \mbm{D}\mbm{x} - \mbm{d} \le \mbm{0}
    @f$


    Define the new objective function as the sum of the old one and logarithmic barrier function:

    @f$
        \phi(\mbm{x}) = 0.5 \mbm{x}^T \mbm{H} \mbm{x} + \mbm{g}^T \mbm{x} 
            -\frac{1}{t} \sum^m_{i=1} \ln {(-\mbm{D}_i \mbm{x} + \mbm{d}_i)}
    @f$


    Approximation of @f$\phi@f$ near some feasible point @f$\mbm{x}@f$:

    @f$
        \phi(\mbm{x} + \mbm{\Delta x}) = \phi(\mbm{x}) + \nabla\phi({\mbm{x}})\mbm{\Delta x}
            + 0.5 \mbm{\Delta x} \nabla^2 \phi(\mbm{x}) \mbm{\Delta x}
    @f$

    where

    @f$\\
        \nabla\phi(\mbm{x}) = \mbm{H}\mbm{x} + \mbm{g} + \frac{1}{t} \sum^m_{i=1} 
            \left( 
                \frac{1}{-\mbm{D}_i\mbm{x} + \mbm{d}_i} \mbm{D}_i^T 
            \right)\\
        \nabla^2\phi(\mbm{x}) = \mbm{H} + \frac{1}{t} \sum^m_{i=1} 
            \left( 
                \frac{1}{(-\mbm{D}_i\mbm{x} + \mbm{d}_i)^2} \mbm{D}_i^T \mbm{D}_i
            \right)
    @f$


    KKT system:

    @f$
        \left(
            \begin{array}{cc}
                \nabla^2\phi(\mbm{x})   &   \mbm{E}^T\\
                \mbm{E}                 &   \mbm{0}\\
            \end{array}
        \right)
        \left(
            \begin{array}{c}
                \mbm{\Delta x}\\
                \mbm{\omega}\\
            \end{array}
        \right)
        =
        \left(
            \begin{array}{c}
                -\nabla\phi(\mbm{x})\\
                \mbm{0}\\
            \end{array}
        \right)
    @f$


    Solution of the KKT system:

    @f$\\
        \mbm{\Delta x} = (\nabla^2\phi(\mbm{x}))^{-1}
            (-\nabla\phi(\mbm{x}) - \mbm{E}^T \mbm{\omega})\\
        \mbm{E} (\nabla^2\phi(\mbm{x}))^{-1}
            (-\nabla\phi(\mbm{x}) - \mbm{E}^T \mbm{\omega}) = 0\\
        \mbm{E} (\nabla^2\phi(\mbm{x}))^{-1} \mbm{E}^T \mbm{\omega} =
            - \mbm{E} (\nabla^2\phi(\mbm{x}))^{-1} \nabla\phi(\mbm{x})
    @f$

    We can find @f$\mbm{\omega}@f$ using Cholesky decomposition and then determine
    the descent direction. Since approximation of @f$\phi@f$ is used, we have to
    resolve the system until @f$\mbm{\Delta x}@f$ or feasible step size is small 
    enough.
 */


/**
 * @page pIPSMPC Derivations for the SMPC scheme

    Inequality constraints:

    @f$
        \mbm{l}_i \le \mbm{x}_{3(i-1) + 1} \le \mbm{u}_i \;\;\; i = 1:2N
    @f$

    Objective function + logarithmic barrier:

    @f$
        \phi(\mbm{x}) = 0.5 \mbm{x}^T \mbm{H} \mbm{x} + \mbm{g}^T \mbm{x} 
            -\frac{1}{t} \sum^{2N}_{i=1} 
            \left(
            \ln (-\mbm{x}_{3(i-1) + 1} + \mbm{u}_i)
            +
            \ln (\mbm{x}_{3(i-1) + 1} - \mbm{l}_i)
            \right)
    @f$


    @f$\nabla\phi(\mbm{x})@f$ and @f$\nabla^2\phi(\mbm{x})@f$ are defined as follows:

    @f$\\
     \nabla\phi(\mbm{x}) = \mbm{H}\mbm{x} + \mbm{g} + \frac{1}{t} \mbm{b}\\
     \mbm{b}_{3(i-1) + 1} = 
        \frac{1}{-\mbm{x}_{3(i-1) + 1} + \mbm{u}_i}
        - 
        \frac{1}{\mbm{x}_{3(i-1) + 1} - \mbm{l}_i}
    @f$

    @f$\\
     \nabla^2\phi(\mbm{x}) = \mbm{H} + \frac{1}{t} \mbm{B}\\
     \mbm{B}_{(3(i-1) + 1),(3(i-1) + 1)} = 
        \frac{1}{(-\mbm{x}_{3(i-1) + 1} + \mbm{u}_i)^2}
        + 
        \frac{1}{(\mbm{x}_{3(i-1) + 1} - \mbm{l}_i)^2}
    @f$

    A diagonal 6x6 matrix, which lies on the main diagonal of @f$\mbm{B}@f$
    and corresponds to state k, is denoted @f$\mbm{B}_k@f$. Each such matrix
    has only two non-zero elements: (1,1) and (4,4).

@section pIPSchur Schur complement (IP method)

    The Schur complement for the new objective function is defined in a similar way
    as for the original objective (see '@ref pSchurComplement'). The part of the hessian 
    corresponding to the vector of states is defined differently though.

    @f$\\
        \bar{\mbm{E}}_c\tilde{\mbm{H}}_c^{-1}\bar{\mbm{E}}_c^T =  \\
    = \left[
      \begin{array}{ccccc} 
        \mbm{M}_{11}    &  -\mbm{M}_{11}\mbm{A}^T    &  \mbm{0}  & \mbm{0}    & \mbm{0}  \\
        -\mbm{A}\mbm{M}_{11} & \mbm{A}\mbm{M}_{11}\mbm{A}^T + \mbm{M}_{22}   &  -\mbm{M}_{22}\mbm{A}^T  & \mbm{0}    & \mbm{0}  \\
        \mbm{0}    &  -\mbm{A}\mbm{M}_{22} & \mbm{A}\mbm{M}_{22}\mbm{A}^T + \mbm{M}_{33}  & -\mbm{M}_{33}\mbm{A}^T    & \mbm{0}  \\
        \mbm{0}    &  \mbm{0}    &  -\mbm{A}\mbm{M}_{33}  & \mbm{A}\mbm{M}_{33}\mbm{A}^T + \mbm{M}_{44}  & -\mbm{M}_{44}\mbm{A}^T  \\
        \mbm{0}    &  \mbm{0}    &  \mbm{0}  & -\mbm{A}\mbm{M}_{44} & \mbm{A}\mbm{M}_{44}\mbm{A}^T + \mbm{M}_{55}
      \end{array}
      \right] 
    @f$

    where @f$\mbm{M}_{ii} = \bar{\mbm{R}}_i\hat{\mbm{Q}}^{-1}_i\bar{\mbm{R}}_i^T@f$ and
    @f$\hat{\mbm{Q}}_i = \tilde{\mbm{Q}_i} + \mbm{B}_i@f$.
 */

/**
 * @page pIPImplementation Notes on implementation
 
@section pIPChol Cholesky decomposition

    The Cholesky factor of the Schur complement must be formed on each iteration.
    Due to the differences in the Schur complement (see '@ref pIPSchur') caused 
    by addition of the logarithic barrier, the Cholesky factor is not so well 
    structured as in the case of the active set method.

    Each 6x6 @f$\mbm{L}_{ij}@f$ matrix of the factor (see '@ref pCholesky') must
    be stored. The matrices are stored sequentially from the top left to the bottom
    right.
 */

#endif /*DOXYGEN_IP_H*/
 
