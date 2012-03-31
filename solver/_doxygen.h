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
 * @mainpage A sparse MPC solver for walking motion generation (old version).
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
 * - <a href="https://github.com/asherikov/smpc_solver/tree/Version_1.x">Sources on GitHub</a>
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
 *      - @ref pPD_EC (example: @ref pExampleEC)
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
 * - Interior-point method
 *      - @ref pIP
 *      - @ref pIPSMPC
 *          - @ref pIPSchur
 *      - @ref pIPImplementation
 *          - @ref pIPChol
 *
 * \n
 *
 * @section MainSrcDocs Internal implementation of the libraries
 * @par
 * - @ref gINTERNALS
 *      - @ref gAS
 *      - @ref gIP
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
 * @defgroup gIP Interior-point method
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


@section pPDVarSub Variable substitutions

@subsection pX_tilde The first substitution
    After the first variable substitution we get 
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
    norm of gravitational acceleration;

@subsection pX_bar The second substitution
    The last substitution rotates the state vector using matrix

    @f$
    \bar{\mbm{R}}_k =
    \left[
      \begin{array}{cccccc} 
        \cos\theta_k & 0 & 0 & -\sin\theta_k & 0 & 0 \\
        0 & 1 & 0 & 0 & 0 & 0 \\
        0 & 0 & 1 & 0 & 0 & 0 \\
        \sin\theta_k & 0 & 0 & \cos\theta_k & 0 & 0 \\
        0 & 0 & 0 & 0 & 1 & 0 \\
        0 & 0 & 0 & 0 & 0 & 1
      \end{array}
    \right]. 
    @f$
     
    where @f$\theta_k@f$ is an angle with respect to the world frame.

    @f$
    \bar{\mbm{c}}_{k} = \bar{\mbm{R}}_k^T \tilde{\mbm{c}}_{k} = 
    (\bar{z}_k^x,\dot{c}_k^x,\ddot{c}_k^x,\bar{z}_k^y,\dot{c}_k^y,\ddot{c}_k^y)
    @f$
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
    \bar{f}(\bar{\mbm{v}}) =
    \left[\hspace{-0.1cm}\begin{array}{c} \bar{\mbm{v}}_c \\ \mbm{v}_u \end{array}\hspace{-0.1cm}\right]^T
    \left[\hspace{-0.1cm}\begin{array}{cc} \tilde{\mbm{H}}_c & \mbm{0} \\ \mbm{0} & \mbm{H}_u \end{array}\hspace{-0.1cm}\hspace{-0.1cm}\right]
    \left[\hspace{-0.1cm}\begin{array}{c} \bar{\mbm{v}}_c \\ \mbm{v}_u \end{array}\hspace{-0.1cm}\right] + 
    \left[\hspace{-0.1cm}\begin{array}{c} \bar{\mbm{v}}_c \\ \mbm{v}_u \end{array}\hspace{-0.1cm}\right]^T
    \left[\hspace{-0.1cm}\begin{array}{c} \bar{\mbm{g}}_c \\ \mbm{0} \end{array}\hspace{-0.1cm}\right] 
    @f$

    where 

    @f$\bar{\mbm{v}}_c @f$ is a column vector containing state vectors and
    @f$\bar{\mbm{v}}_u @f$ is a column vector containing control inputs.

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

    \bar{\mbm{c}}_k^T\frac{\beta}{2}\mbm{C}_p^T\mbm{C}_p\bar{\mbm{c}}_k -
    \bar{\mbm{c}}_k^T\underbrace{\bar{\mbm{R}}^T_k\beta\mbm{C}_p^T\mbm{z}^{\mbox{ref}}_k}_{\bar{\mbm{q}}_k},  \\

    \bar{\mbm{v}}_c = \left[\begin{array}{c} \bar{\mbm{c}}_1 \\ \vdots \\ \bar{\mbm{c}}_N \end{array} \right], \quad

    \bar{\mbm{g}}_c = \left[\begin{array}{c} -\bar{\mbm{q}}_1 \\ \vdots \\ -\bar{\mbm{q}}_N \end{array} \right], \quad

    \bar{\mbm{v}} = \left[\begin{array}{c} \bar{\mbm{v}}_c \\ \mbm{v}_u \end{array} \right]
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
    @f$
    \bar{\mbm{E}}_c\bar{\mbm{v}}_c + \tilde{\mbm{E}}_u\mbm{v}_u = \bar{\mbm{e}}, 
    @f$

    where @f$\bar{\mbm{e}} = (-\mbm{A}\bar{\mbm{R}}_0\bar{\mbm{c}}_0, \mbm{0}, \dots, \mbm{0})@f$,

    @f$
      \bar{\mbm{E}}_c =
      \left[
        \begin{array}{cccccc} 
          -\bar{\mbm{R}}_1    &  \mbm{0}            &  \mbm{0}         & \dots  & \mbm{0}               & \mbm{0}  \\
           \mbm{A}\bar{\mbm{R}}_1 & -\bar{\mbm{R}}_2    &  \mbm{0}         & \dots  & \mbm{0}               & \mbm{0}  \\
           \mbm{0}            &  \mbm{A}\bar{\mbm{R}}_2 & -\bar{\mbm{R}}_3 & \dots  & \mbm{0}               & \mbm{0}  \\
           \vdots             &  \vdots             &  \vdots          & \ddots & \vdots                & \vdots   \\
           \mbm{0}            &  \mbm{0}            &  \mbm{0}         & \dots  & \mbm{A}\bar{\mbm{R}}_{N-1} & -\bar{\mbm{R}}_N \\
        \end{array}
      \right], \quad

      \tilde{\mbm{E}}_u =
      \left[
        \begin{array}{cccc} 
          \tilde{\mbm{B}} & \dots  & \mbm{0} \\
          \vdots     & \ddots & \vdots  \\
          \mbm{0}    & \dots  & \tilde{\mbm{B}} \\
        \end{array}
      \right]. 
    @f$

    See page '@ref pExampleEC' for example.
\n\n


@section pPD_IC Inequality constraints
    @f$
    \left[
      \begin{array}{cccccc} 
        -1 & 0 & 0 & 0 & 0 & 0 \\
        0 & 0 & 0 & -1 & 0 & 0 \\
        1 & 0 & 0 & 0 & 0 & 0 \\
        0 & 0 & 0 & 1 & 0 & 0
      \end{array}
    \right]\bar{\mbm{c}}_k + \mbm{d}_{k} \geq \mbm{0},
    @f$
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
    -# if there are more ZMP positions in ZMP profile go to step 1;
    -# perform variable substitution (rotation) to convert the feasible point
    to @ref pX_bar "X_bar form".
 */



/**
 * @page pExampleEC Derivation of the matrix of equality constraints
    @f$N = 4@f$

    @f$\\
    \tilde{\mbm{c}}_1 = \mbm{A}\tilde{\mbm{c}}_0 + \tilde{\mbm{B}}\mbm{u}_0,  \\
    \tilde{\mbm{c}}_2 = \mbm{A}\tilde{\mbm{c}}_1 + \tilde{\mbm{B}}\mbm{u}_1,  \\
    \tilde{\mbm{c}}_3 = \mbm{A}\tilde{\mbm{c}}_2 + \tilde{\mbm{B}}\mbm{u}_2,  \\
    \tilde{\mbm{c}}_4 = \mbm{A}\tilde{\mbm{c}}_3 + \tilde{\mbm{B}}\mbm{u}_3,  \\
    \tilde{\mbm{c}}_5 = \mbm{A}\tilde{\mbm{c}}_4 + \tilde{\mbm{B}}\mbm{u}_4.  
    @f$

    @f$
      \mbm{E}_c = 
      \left[
        \begin{array}{ccccc} 
          -\mbm{I}    &  \mbm{0}    &  \mbm{0}  & \mbm{0}    & \mbm{0}  \\
           \mbm{A} & -\mbm{I}    &  \mbm{0}  & \mbm{0}    & \mbm{0}  \\
           \mbm{0}    &  \mbm{A} & -\mbm{I}  & \mbm{0}    & \mbm{0}  \\
           \mbm{0}    &  \mbm{0}    &  \mbm{A}  & -\mbm{I}    & \mbm{0}  \\
           \mbm{0}    &  \mbm{0}    &  \mbm{0}  & \mbm{A} & -\mbm{I} \\
        \end{array}
      \right], \quad
        
      \tilde{\mbm{E}}_u = 
      \left[
        \begin{array}{ccccc} 
          \tilde{\mbm{B}} & \mbm{0}    & \mbm{0}    & \mbm{0}   & \mbm{0} \\
          \mbm{0}    & \tilde{\mbm{B}} & \mbm{0}    & \mbm{0}   & \mbm{0} \\
          \mbm{0}    & \mbm{0}    & \tilde{\mbm{B}} & \mbm{0}   & \mbm{0} \\
          \mbm{0}    & \mbm{0}    & \mbm{0}    & \tilde{\mbm{B}} & \mbm{0} \\
          \mbm{0}    & \mbm{0}    & \mbm{0}    & \mbm{0}  & \tilde{\mbm{B}}\\
        \end{array}
      \right], \quad 
    @f$

    @f$
      \mbm{e} = 
      \left[
        \begin{array}{c}
          -\mbm{A}\tilde{\mbm{c}}_0 \\ \mbm{0} \\ \mbm{0} \\ \vdots \\ \mbm{0} 
        \end{array}
      \right]. 
    @f$

    @f$\\
    \bar{\mbm{R}}_1\bar{\mbm{c}}_1 = \mbm{A}\bar{\mbm{R}}_0\bar{\mbm{c}}_0 + \tilde{\mbm{B}}\mbm{u}_0,  \\
    \bar{\mbm{R}}_2\bar{\mbm{c}}_2 = \mbm{A}\bar{\mbm{R}}_1\bar{\mbm{c}}_1 + \tilde{\mbm{B}}\mbm{u}_1,  \\
    \bar{\mbm{R}}_3\bar{\mbm{c}}_3 = \mbm{A}\bar{\mbm{R}}_2\bar{\mbm{c}}_2 + \tilde{\mbm{B}}\mbm{u}_2,  \\
    \bar{\mbm{R}}_4\bar{\mbm{c}}_4 = \mbm{A}\bar{\mbm{R}}_3\bar{\mbm{c}}_3 + \tilde{\mbm{B}}\mbm{u}_3,  \\
    \bar{\mbm{R}}_5\bar{\mbm{c}}_5 = \mbm{A}\bar{\mbm{R}}_4\bar{\mbm{c}}_4 + \tilde{\mbm{B}}\mbm{u}_4. 
    @f$

    @f$
      \bar{\mbm{E}}_c = 
      \left[
        \begin{array}{ccccc} 
          -\bar{\mbm{R}}_1    &  \mbm{0}    &  \mbm{0}  & \mbm{0}    & \mbm{0}  \\
           \mbm{A}\bar{\mbm{R}}_1 & -\bar{\mbm{R}}_2    &  \mbm{0}  & \mbm{0}    & \mbm{0}  \\
           \mbm{0}    &  \mbm{A}\bar{\mbm{R}}_2 & -\bar{\mbm{R}}_3  & \mbm{0}    & \mbm{0}  \\
           \mbm{0}    &  \mbm{0}    &  \mbm{A}\bar{\mbm{R}}_3  & -\bar{\mbm{R}}_4    & \mbm{0}  \\
           \mbm{0}    &  \mbm{0}    &  \mbm{0}  & \mbm{A}\bar{\mbm{R}}_4 & -\bar{\mbm{R}}_5 \\
        \end{array}
      \right], \quad
        
      \tilde{\mbm{E}}_u = 
      \left[
        \begin{array}{ccccc} 
          \tilde{\mbm{B}} & \mbm{0}    & \mbm{0}    & \mbm{0}   & \mbm{0} \\
          \mbm{0}    & \tilde{\mbm{B}} & \mbm{0}    & \mbm{0}   & \mbm{0} \\
          \mbm{0}    & \mbm{0}    & \tilde{\mbm{B}} & \mbm{0}   & \mbm{0} \\
          \mbm{0}    & \mbm{0}    & \mbm{0}    & \tilde{\mbm{B}} & \mbm{0} \\
          \mbm{0}    & \mbm{0}    & \mbm{0}    & \mbm{0}  & \tilde{\mbm{B}}\\
        \end{array}
      \right], \quad 
    @f$

    @f$
      \mbm{e} = 
      \left[
        \begin{array}{c}
          -\mbm{A}\bar{\mbm{R}}_0\bar{\mbm{c}}_0 \\ \mbm{0} \\ \mbm{0} \\ \vdots \\ \mbm{0} 
        \end{array}
      \right]. 
    @f$
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

    For @f$N = 4@f$ we have.

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
          -\bar{\mbm{R}}_1^T    &  \bar{\mbm{R}}_1^T\mbm{A}^T    &  \mbm{0}  & \mbm{0}    & \mbm{0}  \\
           \mbm{0} & -\bar{\mbm{R}}_2^T    &  \bar{\mbm{R}}_2^T\mbm{A}^T  & \mbm{0}    & \mbm{0}  \\
           \mbm{0}    &  \mbm{0} & -\bar{\mbm{R}}_3^T  & \bar{\mbm{R}}_3^T\mbm{A}^T    & \mbm{0}  \\
           \mbm{0}    &  \mbm{0}    &  \mbm{0} & -\bar{\mbm{R}}_4^T    & \bar{\mbm{R}}_4^T\mbm{A}^T  \\
           \mbm{0}    &  \mbm{0}    &  \mbm{0}  & \mbm{0} & -\bar{\mbm{R}}_5^T \\
        \end{array}
      \right]  \\

    = \left[
      \begin{array}{ccccc} 
        -\tilde{\mbm{Q}}^{-1}\bar{\mbm{R}}_1^T & \tilde{\mbm{Q}}^{-1}\bar{\mbm{R}}_1^T\mbm{A}^T  & \mbm{0}    & \mbm{0}   & \mbm{0} \\
        \mbm{0}    & -\tilde{\mbm{Q}}^{-1}\bar{\mbm{R}}_2^T & \tilde{\mbm{Q}}^{-1}\bar{\mbm{R}}_2^T\mbm{A}^T    & \mbm{0}   & \mbm{0} \\
        \mbm{0}    & \mbm{0}    & -\tilde{\mbm{Q}}^{-1}\bar{\mbm{R}}_3^T & \tilde{\mbm{Q}}^{-1}\bar{\mbm{R}}_3^T\mbm{A}^T   & \mbm{0} \\
        \mbm{0}    & \mbm{0}    & \mbm{0}    & -\tilde{\mbm{Q}}^{-1}\bar{\mbm{R}}_4^T & \tilde{\mbm{Q}}^{-1}\bar{\mbm{R}}_4^T\mbm{A}^T \\
        \mbm{0}    & \mbm{0}    & \mbm{0}    & \mbm{0}  & -\tilde{\mbm{Q}}^{-1}\bar{\mbm{R}}_5^T
      \end{array}
      \right]
    @f$

    @f$
    \bar{\mbm{E}}_c\tilde{\mbm{H}}_c^{-1}\bar{\mbm{E}}_c^T =  \\
    = \left[
        \begin{array}{ccccc} 
          -\bar{\mbm{R}}_1    &  \mbm{0}    &  \mbm{0}  & \mbm{0}    & \mbm{0}  \\
           \mbm{A}\bar{\mbm{R}}_1 & -\bar{\mbm{R}}_2    &  \mbm{0}  & \mbm{0}    & \mbm{0}  \\
           \mbm{0}    &  \mbm{A}\bar{\mbm{R}}_2 & -\bar{\mbm{R}}_3  & \mbm{0}    & \mbm{0}  \\
           \mbm{0}    &  \mbm{0}    &  \mbm{A}\bar{\mbm{R}}_3  & -\bar{\mbm{R}}_4    & \mbm{0}  \\
           \mbm{0}    &  \mbm{0}    &  \mbm{0}  & \mbm{A}\bar{\mbm{R}}_4 & -\bar{\mbm{R}}_5 \\
        \end{array}
      \right]

    \left[
      \begin{array}{ccccc} 
        -\tilde{\mbm{Q}}^{-1}\bar{\mbm{R}}_1^T & \tilde{\mbm{Q}}^{-1}\bar{\mbm{R}}_1^T\mbm{A}^T  & \mbm{0}    & \mbm{0}   & \mbm{0} \\
        \mbm{0}    & -\tilde{\mbm{Q}}^{-1}\bar{\mbm{R}}_2^T & \tilde{\mbm{Q}}^{-1}\bar{\mbm{R}}_2^T\mbm{A}^T    & \mbm{0}   & \mbm{0} \\
        \mbm{0}    & \mbm{0}    & -\tilde{\mbm{Q}}^{-1}\bar{\mbm{R}}_3^T & \tilde{\mbm{Q}}^{-1}\bar{\mbm{R}}_3^T\mbm{A}^T   & \mbm{0} \\
        \mbm{0}    & \mbm{0}    & \mbm{0}    & -\tilde{\mbm{Q}}^{-1}\bar{\mbm{R}}_4^T & \tilde{\mbm{Q}}^{-1}\bar{\mbm{R}}_4^T\mbm{A}^T \\
        \mbm{0}    & \mbm{0}    & \mbm{0}    & \mbm{0}  & -\tilde{\mbm{Q}}^{-1}\bar{\mbm{R}}_5^T
      \end{array}
      \right]\\

    = \left[
      \begin{array}{ccccc} 
        \bar{\mbm{R}}_1\tilde{\mbm{Q}}^{-1}\bar{\mbm{R}}_1^T    &  -\bar{\mbm{R}}_1\tilde{\mbm{Q}}^{-1}\bar{\mbm{R}}_1^T\mbm{A}^T    &  \mbm{0}  & \mbm{0}    & \mbm{0}  \\
        -\mbm{A}\bar{\mbm{R}}_1\tilde{\mbm{Q}}^{-1}\bar{\mbm{R}}_1^T & \mbm{A}\bar{\mbm{R}}_1\tilde{\mbm{Q}}^{-1}\bar{\mbm{R}}_1^T\mbm{A}^T + \bar{\mbm{R}}_2\tilde{\mbm{Q}}^{-1}\bar{\mbm{R}}_2^T   &  -\bar{\mbm{R}}_2\tilde{\mbm{Q}}^{-1}\bar{\mbm{R}}_2^T\mbm{A}^T  & \mbm{0}    & \mbm{0}  \\
        \mbm{0}    &  -\mbm{A}\bar{\mbm{R}}_2\tilde{\mbm{Q}}^{-1}\bar{\mbm{R}}_2^T & \mbm{A}\bar{\mbm{R}}_2\tilde{\mbm{Q}}^{-1}\bar{\mbm{R}}_2^T\mbm{A}^T + \bar{\mbm{R}}_3\tilde{\mbm{Q}}^{-1}\bar{\mbm{R}}_3^T  & -\bar{\mbm{R}}_3\tilde{\mbm{Q}}^{-1}\bar{\mbm{R}}_3^T\mbm{A}^T    & \mbm{0}  \\
        \mbm{0}    &  \mbm{0}    &  -\mbm{A}\bar{\mbm{R}}_3\tilde{\mbm{Q}}^{-1}\bar{\mbm{R}}_3^T  & \mbm{A}\bar{\mbm{R}}_3\tilde{\mbm{Q}}^{-1}\bar{\mbm{R}}_3^T\mbm{A}^T + \bar{\mbm{R}}_4\tilde{\mbm{Q}}^{-1}\bar{\mbm{R}}_4^T  & -\bar{\mbm{R}}_4\tilde{\mbm{Q}}^{-1}\bar{\mbm{R}}_4^T\mbm{A}^T  \\
        \mbm{0}    &  \mbm{0}    &  \mbm{0}  & -\mbm{A}\bar{\mbm{R}}_4\tilde{\mbm{Q}}^{-1}\bar{\mbm{R}}_4^T & \mbm{A}\bar{\mbm{R}}_4\tilde{\mbm{Q}}^{-1}\bar{\mbm{R}}_4^T\mbm{A}^T + \bar{\mbm{R}}_5\tilde{\mbm{Q}}^{-1}\bar{\mbm{R}}_5^T
      \end{array}
      \right]  \\

    = \left[
      \begin{array}{ccccc} 
        \mbm{M}_{11}    &  -\mbm{M}_{11}\mbm{A}^T    &  \mbm{0}  & \mbm{0}    & \mbm{0}  \\
        -\mbm{A}\mbm{M}_{11} & \mbm{A}\mbm{M}_{11}\mbm{A}^T + \mbm{M}_{22}   &  -\mbm{M}_{22}\mbm{A}^T  & \mbm{0}    & \mbm{0}  \\
        \mbm{0}    &  -\mbm{A}\mbm{M}_{22} & \mbm{A}\mbm{M}_{22}\mbm{A}^T + \mbm{M}_{33}  & -\mbm{M}_{33}\mbm{A}^T    & \mbm{0}  \\
        \mbm{0}    &  \mbm{0}    &  -\mbm{A}\mbm{M}_{33}  & \mbm{A}\mbm{M}_{33}\mbm{A}^T + \mbm{M}_{44}  & -\mbm{M}_{44}\mbm{A}^T  \\
        \mbm{0}    &  \mbm{0}    &  \mbm{0}  & -\mbm{A}\mbm{M}_{44} & \mbm{A}\mbm{M}_{44}\mbm{A}^T + \mbm{M}_{55}
      \end{array}
      \right], 
    @f$

    where 
    @f$ \mbm{M}_{ii} = 
    \bar{\mbm{R}}_i\tilde{\mbm{Q}}^{-1}\bar{\mbm{R}}_i^T = 
    \tilde{\mbm{Q}}^{-1}@f$ 
    (due to the special structure of 
    @f$ \tilde{\mbm{Q}}^{-1} @f$ and 
    @f$ \bar{\mbm{R}}_i @f$), 
    this is not true if logarithmic barrier is added to the objective see '@ref pIPSchur'.

    @f$\\
      \tilde{\mbm{E}}_u\mbm{H}_u^{-1}\tilde{\mbm{E}}_u^T =
      \left[
        \begin{array}{ccccc} 
          \tilde{\mbm{B}} & \mbm{0}    & \mbm{0}    & \mbm{0}   & \mbm{0} \\
          \mbm{0}    & \tilde{\mbm{B}} & \mbm{0}    & \mbm{0}   & \mbm{0} \\
          \mbm{0}    & \mbm{0}    & \tilde{\mbm{B}} & \mbm{0}   & \mbm{0} \\
          \mbm{0}    & \mbm{0}    & \mbm{0}    & \tilde{\mbm{B}} & \mbm{0} \\
          \mbm{0}    & \mbm{0}    & \mbm{0}    & \mbm{0}  & \tilde{\mbm{B}}\\
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
            \tilde{\mbm{B}}^T & \mbm{0}    & \mbm{0}    & \mbm{0}   & \mbm{0} \\
            \mbm{0}    & \tilde{\mbm{B}}^T & \mbm{0}    & \mbm{0}   & \mbm{0} \\
            \mbm{0}    & \mbm{0}    & \tilde{\mbm{B}}^T & \mbm{0}   & \mbm{0} \\
            \mbm{0}    & \mbm{0}    & \mbm{0}    & \tilde{\mbm{B}}^T & \mbm{0} \\
            \mbm{0}    & \mbm{0}    & \mbm{0}    & \mbm{0}  & \tilde{\mbm{B}}^T\\
          \end{array}
          \right]  \\
        
        =\left[
          \begin{array}{ccccc} 
            \tilde{\mbm{P}} & \mbm{0}    & \mbm{0}    & \mbm{0}   & \mbm{0} \\
            \mbm{0}    & \tilde{\mbm{P}} & \mbm{0}    & \mbm{0}   & \mbm{0} \\
            \mbm{0}    & \mbm{0}    & \tilde{\mbm{P}}  & \mbm{0}   & \mbm{0} \\
            \mbm{0}    & \mbm{0}    & \mbm{0}    & \tilde{\mbm{P}}  & \mbm{0} \\
            \mbm{0}    & \mbm{0}    & \mbm{0}    & \mbm{0}  & \tilde{\mbm{P}}\\
          \end{array}
          \right], 
    @f$

    where 
    @f$\tilde{\mbm{P}} = \tilde{\mbm{B}}\mbm{P}^{-1}\tilde{\mbm{B}}^T@f$.

    @f$\\
      2\mbm{S}_{11} = \mbm{M}_{11} + \tilde{\mbm{P}},  \\
      2\mbm{S}_{kk} = \mbm{A}\mbm{M}_{k-1,k-1}\mbm{A}^T + \mbm{M}_{kk} + \tilde{\mbm{P}},  \\
      2\mbm{S}_{k,k+1} = \mbm{S}_{k+1,k}^T = -\mbm{M}_{kk}\mbm{A}^T. 
    @f$

    @f$\\
      2\mbm{S}_{11} = \tilde{\mbm{Q}}^{-1} + \tilde{\mbm{P}},  \\
      2\mbm{S}_{kk} = \mbm{A}\tilde{\mbm{Q}}^{-1}\mbm{A}^T + \tilde{\mbm{Q}}^{-1} + \tilde{\mbm{P}},  \\
      2\mbm{S}_{k,k+1} = \mbm{S}_{k+1,k}^T = -\tilde{\mbm{Q}}^{-1}\mbm{A}^T. 
    @f$

    Hence, the matrix 
    @f$\mbm{S}@f$ is constant (if 
    @f$\mbm{A}@f$ and 
    @f$\mbm{B}@f$ 
    do not change). 
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
 
