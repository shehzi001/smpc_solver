/**
 * @file
 * @brief This file contains only doxygen definitions and is not to be 
 *  included anywhere.
 *
 * @author Alexander Sherikov
 * @date 15.09.2011 13:24:49 MSD
 *
 *
 * @todo (low priority) Interface with Matlab/Octave.
 * @todo Error processing. Exceptions?
 * @todo Downdate complexity
 */


#ifndef DOXYGEN_H
#define DOXYGEN_H

/**
 * @mainpage A sparse MPC solver for walking motion generation.
 *
 * @par Contents & links
 * - @ref MainIntro
 * - @ref MainLicense
 * - @ref MainRef
 * - @ref MainOverview
 * - @ref MainHowTo
 * - @ref MainFormulas
 * - @ref MainSrcDocs
 * - @ref MainNotes
 * @par
 * - <a href="http://github.com/asherikov/smpc_solver/">Sources on GitHub</a>
 * - <a href="./v1/index.html">Old version of the solver</a>
 *
 * \n
 *
 *
 * @section MainIntro Introduction
 * @verbinclude "README.md"
 * \n
 *
 *
 * @section MainLicense License
 * @verbinclude "LICENSE"
 * \n
 *
 *
 * @section MainRef References
 *
 * Dimitar Nikolaev Dimitrov, Alexander Sherikov, and Pierre-Brice Wieber\n
 * <a href="http://www.aass.oru.se/Research/Learning/drdv_dir/publications/iros11/iros11_1.html">
 * A sparse model predictive control formulation for walking motion generation</a>\n
 * IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS),\n
 * September 25-30, 2011, San Francisco, California
 * \n\n
 *
 *
 * @section MainOverview Purpose of the library
 * Walk of a robot can be controlled using the following scheme:
 * - (1) determine a desired position of the center of mass (CoM);
 * - (2) determine desired positions of the end-effectors (feet);
 * - (3) using the model of the robot compute necessary changes of joint positions.
 *
 * The smpc_solver library addresses the first step. The position of the CoM is
 * determined using MPC scheme, where it is modeled using inverted pendulum in 3 
 * dimensions. MPC implies, that we have to solve an optimization problem. The 
 * objective function is quadratic, additional constraints are imposed by position 
 * of the support foot/feet: we need to satisfy certain requrements to
 * prevent robot from falling.
 *
 * The WMG library contains functions, which are necessary to define footsteps, 
 * prepare parameters for the solver library, determine positions of the feet.
 *
 * Terms and abbreviations:
 * - CoM -- center of mass.
 * - ZMP -- zero moment point.
 * - MPC -- model predictive control.
 * - SMPC -- sparse model predictive control, this term is introduced by us.
 * - Support foot -- the foot, on which a robot is standing.
 * - Single support -- a situation, when a robot stands on only one foot.
 * - Double support -- a situation, when a robot stands on both feet.
 * - Reference foot -- even when a robot is in double support, we use one foot 
 *   as the reference.
 * - Reference ZMP coordinates -- the objective function contains term, that
 *   tries to minimize the difference between the solution an reference points.
 *   This is allows tuning of the solution.
 * 
 * For more information refer to the papers listed in '@ref MainRef'. Also, the
 * '@ref pProblemDef' section contains more detailed explanations.
 * \n\n
 *
 *
 * @section MainHowTo API & examples
 * @par
 * - @ref gAPI 
 * - @ref gWMG_API
 * - @ref pDemo "A simple demo"
 *
 * \n
 *
 *
 * @section MainFormulas Derivations and algorithms
 * @par Conventions & definitions
 * - All letters referencing matrices or vectors in formulas are typed in bold.
 * - Matrices and vectors are usually denoted by the upper case and lower case
 *   letters respectively.
 * - 'N' denotes the size of the preview window.
 * @par Further information can be found on the following pages
 * - @ref pProblemDef
 *      - @ref pPDModel
 *      - @ref pPDVarSub
 *      - @ref pPDObj
 *      - @ref pPD_EC
 *      - @ref pPD_IC
 * - @ref pKKT
 * - @ref pInitGuess
 * - @ref pSchurComplement
 * - @ref pCholesky
 * - Active set method
 *      - @ref pAddIC
 *          - @ref pCholUp
 *          - @ref pAddICz
 *      - @ref pRemoveIC
 *          - @ref pCholDown
 *          - @ref pRemoveICz
 *      - @ref pDetails 
 *          - @ref pDetMatrices
 *          - @ref pBounds
 *
 * \n
 *
 * @section MainSrcDocs Internal implementation of the libraries
 * @par
 * - @ref gINTERNALS
 *      - @ref gAS
 * - @ref gWMG_INTERNALS
 *
 * \n
 *
 * @section MainNotes Important notes
 * @todo The library is not thoroughly tested with variable height of CoM.
 */


/**
 * @defgroup gAPI API of the solver library
 *  
 * @defgroup gINTERNALS Internal classes, functions and definitions of the solver library
 *
 * @defgroup gAS Active set method
 * @ingroup gINTERNALS
 *
 * @defgroup gWMG_API API of the simulation support library
 *
 * @defgroup gWMG_INTERNALS Internal classes, functions and definitions of the simulation support library.
 */


/**
 * @page pDemo A simple demo
 * @dontinclude "test/demo.cpp"
 *
 * @section pDemo1 Include headers:
 * @until -------
 *
 * @section pDemo2 Define footsteps:
 * @until -------
 *
 * @section pDemo3 Setup solver:
 * @until -------
 * The gains are defined @ref pGains "here". The default values of optional
 * parameters can be found here: smpc_solver#smpc_solver.
 *
 * @section pDemo4 Start simulation:
 * @until -------
 *
 * @section pDemo5 Solve QP on each iteration of simulation:
 *
 * @subsection pDemo51 Initialize QP
 * @until -------
 * See also smpc_solver#init. Note, that the current state must be given
 * after @ref pX_tilde "first variable substitution".
 *
 * @subsection pDemo52 Solve QP
 * @until -------
 * See also smpc_solver#solve.
 *
 * @subsection pDemo53 Obtain the next state
 * @until -------
 * See also smpc_solver#get_next_state_tilde.
 */

 


/**
@page pProblemDef Problem definition

@section pPDModel Model of the system
    3D linear inverted pendulum is used as an approximate model
    of a humanoid robot.

    @f$
    \tilde{\mbm{c}}_{k+1} = \tilde{\mbm{A}}_k\tilde{\mbm{c}}_{k}+\tilde{\mbm{B}}_k\dddot{\mbm{c}}_{k}, \quad
    \dddot{\mbm{c}}_{k} = (\dddot{c}_k^x, \dddot{c}_k^y)
    @f$

    @f$
    \mbm{A} = \left[\hspace{-0.1cm}
      \begin{array}{cccccc} 
        1 & T_k & T_k^{2}/2 & 0 & 0 & 0 \vspace{0.05cm}\\ 
        0 & 1 & T_k & 0 & 0 & 0 \vspace{0.05cm}\\ 
        0 & 0 & 1 & 0 & 0 & 0 \vspace{0.05cm}\\
        0 & 0 & 0 & 1 & T_k & T_k^{2}/2 \vspace{0.05cm}\\
        0 & 0 & 0 & 0 & 1 & T_k \vspace{0.05cm}\\
        0 & 0 & 0 & 0 & 0 & 1
      \end{array}
      \hspace{-0.1cm}\right], \quad 

    \mbm{B}_k = \left[\hspace{-0.1cm}
      \begin{array}{cc}
        T_k^{3}/6 & 0 \vspace{0.05cm} \\ 
        T_k^{2}/2 & 0 \vspace{0.05cm} \\
        T_k      & 0\\
        0      & T_k^{3}/6 \vspace{0.05cm} \\
        0      & T_k^{2}/2 \vspace{0.05cm} \\
        0      & T_k
      \end{array}
    \right]
    @f$

    Where @f$T_k@f$ is a time sampling period in the preview window.

    Originally the state vector is defined as
    @f$
    \hat{\mbm{c}}_{k} = (c_k^x,\dot{c}_k^x,\ddot{c}_k^x,c_k^y,\dot{c}_k^y,\ddot{c}_k^y)
    @f$
    where @f$c^x_k, c^y_k@f$ are coordintes of the center of mass.
\n\n


@section pPDVarSub Variable substitution
    @anchor pX_tilde

    After the variable substitution we get 
    @f$
    \tilde{\mbm{c}}_{k} = (z_k^x,\dot{c}_k^x,\ddot{c}_k^x,z_k^y,\dot{c}_k^y,\ddot{c}_k^y)
    @f$
    where @f$z^x_k, z^y_k@f$ are coordintes of the ZMP.

    The state and control input matrices are changed accordingly:

    @f$
    \tilde{\mbm{A}} =
    \left[
      \begin{array}{cccccc} 
        1 & T_k & T_k^{2}/2-\Delta h_k & 0 & 0 & 0 \vspace{0.05cm}\\ 
        0 & 1 & T_k & 0 & 0 & 0 \vspace{0.05cm}\\ 
        0 & 0 & 1 & 0 & 0 & 0 \vspace{0.05cm}\\
        0 & 0 & 0 & 1 & T_k & T_k^{2}/2-\Delta h_k \vspace{0.05cm}\\
        0 & 0 & 0 & 0 & 1 & T_k \vspace{0.05cm}\\
        0 & 0 & 0 & 0 & 0 & 1
      \end{array}
      \right], \quad
    \tilde{\mbm{B}} = 
    \left[
      \begin{array}{cc}
        T^{3}/6-hT & 0 \\ 
        T^{2}/2 & 0\\ 
        T      & 0\\
        0      & T^{3}/6-hT \\
        0      & T^{2}/2\\
        0      & T
      \end{array}
    \right]
    @f$

    @anchor ph
    Here @f$h = c_k^z/g@f$, i.e. the height of center of mass divided by the 
    norm of gravitational acceleration.
\n\n


@section pPDObj Objective function
    Output matrices for position and velocity:

    @f$
    \mbm{C}_p =
    \left[
      \begin{array}{cccccc} 
        1 & 0 & 0 & 0 & 0 & 0 \\
        0 & 0 & 0 & 1 & 0 & 0 \\
      \end{array}\right], \quad  

    \mbm{C}_v =
    \left[
      \begin{array}{cccccc} 
        0 & 1 & 0 & 0 & 0 & 0 \\
        0 & 0 & 0 & 0 & 1 & 0 \\
      \end{array}\right], 
    @f$


    @f$
    \tilde{f}(\tilde{\mbm{v}}) =
    \left[\hspace{-0.1cm}\begin{array}{c} \tilde{\mbm{v}}_c \\ \mbm{v}_u \end{array}\hspace{-0.1cm}\right]^T
    \left[\hspace{-0.1cm}\begin{array}{cc} \tilde{\mbm{H}}_c & \mbm{0} \\ \mbm{0} & \mbm{H}_u \end{array}\hspace{-0.1cm}\hspace{-0.1cm}\right]
    \left[\hspace{-0.1cm}\begin{array}{c} \tilde{\mbm{v}}_c \\ \mbm{v}_u \end{array}\hspace{-0.1cm}\right] + 
    \left[\hspace{-0.1cm}\begin{array}{c} \tilde{\mbm{v}}_c \\ \mbm{v}_u \end{array}\hspace{-0.1cm}\right]^T
    \left[\hspace{-0.1cm}\begin{array}{c} \tilde{\mbm{g}}_c \\ \mbm{0} \end{array}\hspace{-0.1cm}\right] 
    @f$

    where 

    @f$\tilde{\mbm{v}}_c @f$ is a column vector containing state vectors and
    @f${\mbm{v}}_u @f$ is a column vector containing control inputs.

    or

    @f$
    f(\mbm{v}) = \frac{\gamma}{2}\sum_{k=0}^{N-1}\left(\dddot{\mbm{c}}_k^T\dddot{\mbm{c}}_k\right) + 
    \frac{\alpha}{2}\sum_{k=1}^{N}\left(\dot{\mbm{c}}_k^T\dot{\mbm{c}}_k\right) + 
    \frac{\beta}{2}\sum_{k=1}^{N}\left(\mbm{z}_k^T\mbm{z}_k - 2\mbm{z}_k^T\mbm{z}^{\mbox{ref}}_k\right),
    @f$

    @anchor pGains
    where @f$\alpha, \beta, \gamma > 0@f$ are gains.

    @anchor pg
    @f$
    \frac{\beta}{2}\mbm{z}_k^T\mbm{z}_k - \beta\mbm{z}_k^T\mbm{z}^{\mbox{ref}}_k =

    \tilde{\mbm{c}}_k^T\frac{\beta}{2}\mbm{C}_p^T\mbm{C}_p\tilde{\mbm{c}}_k -
    \tilde{\mbm{c}}_k^T\underbrace{\beta\mbm{C}_p^T\mbm{z}^{\mbox{ref}}_k}_{\tilde{\mbm{q}}_k},  \\

    \tilde{\mbm{v}}_c = \left[\begin{array}{c} \tilde{\mbm{c}}_1 \\ \vdots \\ \tilde{\mbm{c}}_N \end{array} \right], \quad

    \tilde{\mbm{g}}_c = \left[\begin{array}{c} -\tilde{\mbm{q}}_1 \\ \vdots \\ -\tilde{\mbm{q}}_N \end{array} \right], \quad

    \tilde{\mbm{v}} = \left[\begin{array}{c} \tilde{\mbm{v}}_c \\ \mbm{v}_u \end{array} \right]
    @f$


    @f$\\
      \mbm{H}_u = 
      \left[
        \begin{array}{ccc}
          \mbm{P} & \dots  & \mbm{0} \\ 
          \vdots  & \ddots & \vdots  \\ 
          \mbm{0} & \dots  & \mbm{P} \\ 
        \end{array}
      \right],\quad

    \mbm{P} = 
      \left[
        \begin{array}{cc}
          \frac{\gamma}{2}  & 0  \\ 
          0                 &  \frac{\gamma}{2} \\ 
        \end{array}
      \right];\\

    \tilde{\mbm{H}}_c =
    \left[
      \begin{array}{ccc}
        \tilde{\mbm{Q}} & \dots  & \mbm{0} \\ 
        \vdots  & \ddots & \vdots  \\ 
        \mbm{0} & \dots  & \tilde{\mbm{Q}} \\ 
      \end{array}
    \right], \quad 

    \tilde{\mbm{Q}} =
    \left[
      \begin{array}{cccccc}
        \frac{\beta}{2} & 0                 & 0 & 0 & 0 & 0\\ 
        0               & \frac{\alpha}{2}  & 0 & 0 & 0 & 0\\ 
        0               & 0                 & r & 0 & 0 & 0\\ 
        0 & 0 & 0 & \frac{\beta}{2} & 0                 & 0 \\ 
        0 & 0 & 0 & 0               & \frac{\alpha}{2}  & 0 \\ 
        0 & 0 & 0 & 0               & 0                 & r \\ 
      \end{array}
    \right]
    @f$

    @anchor RegFactor
    Here @f$r@f$ is a regularization factor, which makes the matrix nonsingular.
\n\n


@section pPD_EC Equality constraints

    The states of the system can be found using equations

    @f$\\
    \tilde{\mbm{c}}_1 = \tilde{\mbm{A}}_0\tilde{\mbm{c}}_0 + \tilde{\mbm{B}}_0\mbm{u}_0  \\
    \tilde{\mbm{c}}_2 = \tilde{\mbm{A}}_1\tilde{\mbm{c}}_1 + \tilde{\mbm{B}}_1\mbm{u}_1  \\
    \dots
    @f$

    From these equations we can build equality constraints

    @f$
    \tilde{\mbm{E}}_c\tilde{\mbm{v}}_c + \tilde{\mbm{E}}_u\mbm{v}_u = \tilde{\mbm{e}}, 
    @f$

    Where

    @f$\tilde{\mbm{e}} = (-\tilde{\mbm{A}}_0\tilde{\mbm{c}}_0, \mbm{0}, \dots, \mbm{0})@f$,

    @f$
      \tilde{\mbm{E}}_c =
      \left[
        \begin{array}{cccccc} 
          -\mbm{I}           &  \mbm{0} &  \mbm{0} & \dots  & \mbm{0}               & \mbm{0}  \\
           \tilde{\mbm{A}}_0 & -\mbm{I} &  \mbm{0} & \dots  & \mbm{0}               & \mbm{0}  \\
           \vdots            &  \vdots  &  \vdots  & \ddots & \vdots                & \vdots   \\
           \mbm{0}           &  \mbm{0} &  \mbm{0} & \dots  & \tilde{\mbm{A}}_{N-1} & -\mbm{I} \\
        \end{array}
      \right], \quad

      \tilde{\mbm{E}}_u =
      \left[
        \begin{array}{cccc} 
          \tilde{\mbm{B}}_0 & \dots  & \mbm{0} \\
          \vdots     & \ddots & \vdots  \\
          \mbm{0}    & \dots  & \tilde{\mbm{B}}_{N-1} \\
        \end{array}
      \right]. 
    @f$
\n\n


@section pPD_IC Inequality constraints
    @f$
    \left[
      \begin{array}{cc} 
        -1 & 0  \\
        0  & -1 \\
        1  & 0  \\
        0  & 1  
      \end{array}
    \right]\mbm{R}^T_k \mbm{C}_p \tilde{\mbm{c}}_k + \mbm{d}_{k} \geq \mbm{0},
    @f$

    Where 
    @f$
    \mbm{R}_k =
    \left[
      \begin{array}{cc} 
        \cos\theta_k & -\sin\theta_k \\
        \sin\theta_k & \cos\theta_k
      \end{array}
    \right] 
    @f$
    is a rotation matrix for the corresponding rectangular support.
 */


/**
@page pKKT KKT system
    @f$
      \left[
        \begin{array}{cc} 
            2\mbm{H} & \mbm{E}^T \\ 
            \mbm{E} & \mbm{0} 
        \end{array}
      \right]
      \left[
        \begin{array}{c} 
            \mbm{x}_{init} + \Delta\mbm{x}\\ 
            \mbm{\nu}
        \end{array}
      \right] =
      \left[
        \begin{array}{c} 
            -\mbm{g} \\ 
            \mbm{\bar{e}} 
        \end{array}
      \right]. 
    @f$

    Section '@ref pProblemDef' discusses formation of matrices
    @f$
    \mbm{H}, \mbm{E}^T, \mbm{g}, \mbm{\bar{e}} 
    @f$.

    We assume, that an initial guess satisfying all constraints is given
    (instructions on how to generate it are given in section '@ref pInitGuess')
    Our goal is to find delta between initial guess and optimal point.
    From the system presented above we can derive:

    @f$\\
    \frac{1}{2} \mbm{E} \mbm{H}^{-1} \mbm{E}^T \mbm{\nu} = 
        \mbm{S} \mbm{\nu} = 
        \mbm{E} (-\frac{1}{2} \mbm{H}^{-1} \mbm{g} - \mbm{x}_{init}) = \mbm{s}\\
    \Delta\mbm{x} = -\frac{1}{2} \mbm{H}^{-1} \mbm{g} - \mbm{x}_{init} - 
        \frac{1}{2} \mbm{H}^{-1} \mbm{E}^T \mbm{\nu}
    @f$

    Here
    @f$
    \mbm{S}
    @f$
    is a Schur complement, its structure is described in section '@ref pSchurComplement'.

    @anchor piHg
    Note, that
    @f$
    \frac{1}{2} \mbm{H}^{-1} \mbm{g}
    @f$
    is constant and due to @ref pPDObj "the structure of the matrix and vector"
    has only 2*N non-zero elements. Also, since Hessian is diagonal its 
    invertion is trivial and multiplication of inverted Hessian by any vector 
    is O(N).
 */

/**
 * @page pInitGuess Generation of an initial feasible point
    Consider @ref pPDModel "model of the system".
 
    We assume, that initial state is given (in @ref pX_tilde "X_tilde form").
    A feasible ZMP profile can be build by selecting reference points, which
    satisfy inequality constraints. Then based on this profile we can find
    all states and control inputs:

    -# since the coordinates of the current and the next ZMP positions are 
    known, we can compute control inputs necessary to change position;
    -# given control inputs and current state we can find the next state;
    -# if there are more ZMP positions in ZMP profile go to step 1.
 */


/**
 * @page pSchurComplement Schur complement

    In order to solve @ref pKKT we have to form Schur complement:

    @f$\\
    \mbm{S} = \frac{1}{2}\mbm{E}\mbm{H}^{-1}\mbm{E}^T = \frac{1}{2}\left[\begin{array}{cc}\bar{\mbm{E}}_c  \tilde{\mbm{E}}_u\end{array}\right]
    \left[\begin{array}{cc}\tilde{\mbm{H}}_c & \mbm{0} \\ \mbm{0} & \mbm{H}_u\end{array}\right]
    \left[\begin{array}{c}\bar{\mbm{E}}_c^T \\ \tilde{\mbm{E}}_u^T \end{array}\right] 

    = \frac{1}{2}\bar{\mbm{E}}_c\tilde{\mbm{H}}_c^{-1}\bar{\mbm{E}}_c^T + \frac{1}{2}\tilde{\mbm{E}}_u\mbm{H}_u^{-1}\tilde{\mbm{E}}_u^T. 
    @f$

    For @f$N = 4@f$ we have

    @f$\\
    \tilde{\mbm{H}}_c^{-1}\bar{\mbm{E}}_c^T = 
      \left[
        \begin{array}{ccccc} 
          \tilde{\mbm{Q}}^{-1} & \mbm{0}    & \mbm{0}    & \mbm{0}   & \mbm{0} \\
          \mbm{0}    & \tilde{\mbm{Q}}^{-1} & \mbm{0}    & \mbm{0}   & \mbm{0} \\
          \mbm{0}    & \mbm{0}    & \tilde{\mbm{Q}}^{-1} & \mbm{0}   & \mbm{0} \\
          \mbm{0}    & \mbm{0}    & \mbm{0}    & \tilde{\mbm{Q}}^{-1} & \mbm{0} \\
          \mbm{0}    & \mbm{0}    & \mbm{0}    & \mbm{0}  & \tilde{\mbm{Q}}^{-1}
        \end{array}
        \right]
      
    \left[
        \begin{array}{ccccc} 
          -\mbm{I} &  \tilde{\mbm{A}}^T_1 &  \mbm{0}             & \mbm{0}             & \mbm{0}  \\
           \mbm{0} & -\mbm{I}             &  \tilde{\mbm{A}}^T_2 & \mbm{0}             & \mbm{0}  \\
           \mbm{0} &  \mbm{0}             & -\mbm{I}             & \tilde{\mbm{A}}^T_3 & \mbm{0}  \\
           \mbm{0} &  \mbm{0}             &  \mbm{0}             & -\mbm{I}            & \tilde{\mbm{A}}^T_4  \\
           \mbm{0} &  \mbm{0}             &  \mbm{0}             & \mbm{0}             & -\mbm{I} \\
        \end{array}
      \right]  \\

    = \left[
      \begin{array}{ccccc} 
        -\tilde{\mbm{Q}}^{-1} & \tilde{\mbm{Q}}^{-1}\tilde{\mbm{A}}^T_1  & \mbm{0}    & \mbm{0}   & \mbm{0} \\
        \mbm{0}    & -\tilde{\mbm{Q}}^{-1} & \tilde{\mbm{Q}}^{-1}\tilde{\mbm{A}}^T_2    & \mbm{0}   & \mbm{0} \\
        \mbm{0}    & \mbm{0}    & -\tilde{\mbm{Q}}^{-1} & \tilde{\mbm{Q}}^{-1}\tilde{\mbm{A}}^T_3   & \mbm{0} \\
        \mbm{0}    & \mbm{0}    & \mbm{0}    & -\tilde{\mbm{Q}}^{-1} & \tilde{\mbm{Q}}^{-1}\tilde{\mbm{A}}^T_4 \\
        \mbm{0}    & \mbm{0}    & \mbm{0}    & \mbm{0}  & -\tilde{\mbm{Q}}^{-1}
      \end{array}
      \right]
    @f$

    @f$
    \bar{\mbm{E}}_c\tilde{\mbm{H}}_c^{-1}\bar{\mbm{E}}_c^T =  \\
    = \left[
        \begin{array}{ccccc} 
          -\mbm{I}    &  \mbm{0}    &  \mbm{0}  & \mbm{0}    & \mbm{0}  \\
           \tilde{\mbm{A}}_1 & -\mbm{I}    &  \mbm{0}  & \mbm{0}    & \mbm{0}  \\
           \mbm{0}    &  \tilde{\mbm{A}}_2 & -\mbm{I}  & \mbm{0}    & \mbm{0}  \\
           \mbm{0}    &  \mbm{0}    &  \tilde{\mbm{A}}_3  & -\mbm{I}    & \mbm{0}  \\
           \mbm{0}    &  \mbm{0}    &  \mbm{0}  & \tilde{\mbm{A}}_4 & -\mbm{I} \\
        \end{array}
      \right]

    \left[
      \begin{array}{ccccc} 
        -\tilde{\mbm{Q}}^{-1} & \tilde{\mbm{Q}}^{-1}\tilde{\mbm{A}}^T_1  & \mbm{0}    & \mbm{0}   & \mbm{0} \\
        \mbm{0}    & -\tilde{\mbm{Q}}^{-1} & \tilde{\mbm{Q}}^{-1}\tilde{\mbm{A}}^T_2    & \mbm{0}   & \mbm{0} \\
        \mbm{0}    & \mbm{0}    & -\tilde{\mbm{Q}}^{-1} & \tilde{\mbm{Q}}^{-1}\tilde{\mbm{A}}^T_3   & \mbm{0} \\
        \mbm{0}    & \mbm{0}    & \mbm{0}    & -\tilde{\mbm{Q}}^{-1} & \tilde{\mbm{Q}}^{-1}\tilde{\mbm{A}}^T_4 \\
        \mbm{0}    & \mbm{0}    & \mbm{0}    & \mbm{0}  & -\tilde{\mbm{Q}}^{-1}
      \end{array}
      \right]

    = \left[
      \begin{array}{ccccc} 
        \tilde{\mbm{Q}}^{-1}    &  -\tilde{\mbm{Q}}^{-1}\tilde{\mbm{A}}^T_1    &  \mbm{0}  & \mbm{0}    & \mbm{0}  \\
        -\tilde{\mbm{A}}_1\tilde{\mbm{Q}}^{-1} & \tilde{\mbm{A}}_1\tilde{\mbm{Q}}^{-1}\tilde{\mbm{A}}^T_1 + \tilde{\mbm{Q}}^{-1}   &  -\tilde{\mbm{Q}}^{-1}\tilde{\mbm{A}}^T_2  & \mbm{0}    & \mbm{0}  \\
        \mbm{0}    &  -\tilde{\mbm{A}}_2\tilde{\mbm{Q}}^{-1} & \tilde{\mbm{A}}_2\tilde{\mbm{Q}}^{-1}\tilde{\mbm{A}}^T_2 + \tilde{\mbm{Q}}^{-1}  & -\tilde{\mbm{Q}}^{-1}\tilde{\mbm{A}}^T_3    & \mbm{0}  \\
        \mbm{0}    &  \mbm{0}    &  -\tilde{\mbm{A}}_3\tilde{\mbm{Q}}^{-1}  & \tilde{\mbm{A}}_3\tilde{\mbm{Q}}^{-1}\tilde{\mbm{A}}^T_3 + \tilde{\mbm{Q}}^{-1}  & -\tilde{\mbm{Q}}^{-1}\tilde{\mbm{A}}^T_4  \\
        \mbm{0}    &  \mbm{0}    &  \mbm{0}  & -\tilde{\mbm{A}}_4\tilde{\mbm{Q}}^{-1} & \tilde{\mbm{A}}_4\tilde{\mbm{Q}}^{-1}\tilde{\mbm{A}}^T_4 + \tilde{\mbm{Q}}^{-1}
      \end{array}
      \right]
    @f$


    @f$\\
      \tilde{\mbm{E}}_u\mbm{H}_u^{-1}\tilde{\mbm{E}}_u^T =
      \left[
        \begin{array}{ccccc} 
          \tilde{\mbm{B}}_0 & \mbm{0}    & \mbm{0}    & \mbm{0}   & \mbm{0} \\
          \mbm{0}    & \tilde{\mbm{B}}_1 & \mbm{0}    & \mbm{0}   & \mbm{0} \\
          \mbm{0}    & \mbm{0}    & \tilde{\mbm{B}}_2 & \mbm{0}   & \mbm{0} \\
          \mbm{0}    & \mbm{0}    & \mbm{0}    & \tilde{\mbm{B}}_3 & \mbm{0} \\
          \mbm{0}    & \mbm{0}    & \mbm{0}    & \mbm{0}  & \tilde{\mbm{B}}_4\\
        \end{array}
      \right]
      
      \left[
      \begin{array}{ccccc} 
        \mbm{P}^{-1} & \mbm{0}    & \mbm{0}    & \mbm{0}   & \mbm{0} \\
        \mbm{0}    & \mbm{P}^{-1} & \mbm{0}    & \mbm{0}   & \mbm{0} \\
        \mbm{0}    & \mbm{0}    & \mbm{P}^{-1} & \mbm{0}   & \mbm{0} \\
        \mbm{0}    & \mbm{0}    & \mbm{0}    & \mbm{P}^{-1} & \mbm{0} \\
        \mbm{0}    & \mbm{0}    & \mbm{0}    & \mbm{0}  & \mbm{P}^{-1}
      \end{array}
      \right]
        
        \left[
          \begin{array}{ccccc} 
            \tilde{\mbm{B}}^T_0 & \mbm{0}    & \mbm{0}    & \mbm{0}   & \mbm{0} \\
            \mbm{0}    & \tilde{\mbm{B}}^T_1 & \mbm{0}    & \mbm{0}   & \mbm{0} \\
            \mbm{0}    & \mbm{0}    & \tilde{\mbm{B}}^T_2 & \mbm{0}   & \mbm{0} \\
            \mbm{0}    & \mbm{0}    & \mbm{0}    & \tilde{\mbm{B}}^T_3 & \mbm{0} \\
            \mbm{0}    & \mbm{0}    & \mbm{0}    & \mbm{0}  & \tilde{\mbm{B}}^T_4\\
          \end{array}
          \right]  \\
        
        =\left[
          \begin{array}{ccccc} 
            \tilde{\mbm{P}}_0 & \mbm{0}    & \mbm{0}    & \mbm{0}   & \mbm{0} \\
            \mbm{0}    & \tilde{\mbm{P}}_1 & \mbm{0}    & \mbm{0}   & \mbm{0} \\
            \mbm{0}    & \mbm{0}    & \tilde{\mbm{P}}_2  & \mbm{0}   & \mbm{0} \\
            \mbm{0}    & \mbm{0}    & \mbm{0}    & \tilde{\mbm{P}}_3  & \mbm{0} \\
            \mbm{0}    & \mbm{0}    & \mbm{0}    & \mbm{0}  & \tilde{\mbm{P}}_4\\
          \end{array}
          \right], 
    @f$

    where 
    @f$\tilde{\mbm{P}}_k = \tilde{\mbm{B}}_k\mbm{P}^{-1}\tilde{\mbm{B}}^T_k@f$.

    @f$\\
      2\mbm{S}_{11} = \tilde{\mbm{Q}}^{-1} + \tilde{\mbm{P}}_0,  \\
      2\mbm{S}_{kk} = \mbm{A}_{k-1}\tilde{\mbm{Q}}^{-1}\mbm{A}^T_{k-1} + \tilde{\mbm{Q}}^{-1} + \tilde{\mbm{P}}_{k-1},  \\
      2\mbm{S}_{k,k+1} = \mbm{S}_{k+1,k}^T = -\tilde{\mbm{Q}}^{-1}\mbm{A}^T_{k}. 
    @f$
*/


/**
 * @page pCholesky Cholesky decomposition of Schur complement
 
    Once Schur complement is formed we can use Cholesky decomposition
    @f$\mbm{S} = \mbm{L}\mbm{L}^T@f$.
    to obtain Langrange multipliers.

    @f$
    \mbm{L} = 
    \left[
      \begin{array}{cccccc} 
        \mbm{L}_{11}  & \mbm{0}  & \mbm{0}     & \dots  & \mbm{0} & \mbm{0}        \\
        \mbm{L}_{21}  & \mbm{L}_{22}  & \mbm{0} & \dots  & \mbm{0} & \mbm{0}        \\
        \mbm{0}      & \mbm{L}_{32}  & \mbm{L}_{33} & \dots  & \mbm{0} & \mbm{0}        \\
        \vdots       & \vdots       & \vdots      & \ddots & \vdots  & \vdots         \\
        \mbm{0}      & \mbm{0}      & \mbm{0}     & \dots  & \mbm{L}_{N-1,N-1} & \mbm{0} \\
        \mbm{0}      & \mbm{0}      & \mbm{0}     & \dots  &  \mbm{L}_{N,N-1} & \mbm{L}_{NN} 
      \end{array}
    \right], 
    @f$

    Directly from observation we have

    @f$
    \\
    \mbm{S}_{11} = \mbm{L}_{11}\mbm{L}_{11}^T,  \\
    \mbm{S}_{12} = \mbm{S}_{21}^T = \mbm{L}_{11}\mbm{L}_{21}^T, \quad \mbm{L}_{21}^T = \mbm{L}_{11}^{-1}\mbm{S}_{12},   \\
    \mbm{S}_{22} = \mbm{L}_{21}\mbm{L}_{21}^T + \mbm{L}_{22}\mbm{L}_{22}^T, \quad \dots 
    @f$
     
    In the second step 
    @f$\mbm{L}_{21}^T@f$ 
    is computed by forward substitution, and in the third step, forming 
    @f$\mbm{L}_{22}@f$ 
    requires the computation of the Cholesky factors of 
    @f$\mbm{S}_{22} - \mbm{L}_{21}\mbm{L}_{21}^T@f$. 
\n\n
 */

#endif /*DOXYGEN_H*/
 
