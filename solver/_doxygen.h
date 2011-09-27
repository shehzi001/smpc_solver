/**
 * @file _doxygen.h
 * @brief This file contains only doxygen definitions and is not to be 
 *  included anywhere.
 *
 * @author Alexander Sherikov
 * @date 15.09.2011 13:24:49 MSD
 */


#ifndef DOXYGEN_H
#define DOXYGEN_H
/**
 * @mainpage A sparse MPC solver for walking motion generation.
 *
 * @section MainIntro Introduction
 * @verbinclude "README"
 *
 * For more information on the sparse MPC refer to the papers listed 
 * in the '@ref MainRef' section. Section '@ref MainHowTo' contains 
 * instructions for application programmer. Implementaion details are 
 * given in section '@ref MainInternals'.
 * \n\n
 *
 *
 * @section MainRef References
 *
 * Dimitar Nikolaev Dimitrov, Alexander Sherikov, and Pierre-Brice Wieber\n
 * <a href ="http://www.aass.oru.se/Research/Learning/drdv_dir/publications/iros11/iros11_1.html">
 * A sparse model predictive control formulation for walking motion generation</a>\n
 * IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS),\n
 * September 25-30, 2011, San Francisco, California
 * \n\n
 * 
 *
 * @section MainLicense License
 * @verbinclude "LICENSE"
 * \n
 *
 *
 * @section MainHowTo API & examples
 * - @ref gAPI 
 * - @ref pDemo "A simple demo"
 * \n\n
 *
 * @section MainInternals Implementation details
 * - @ref gINTERNALS
 * - @ref pFormulas 
 */


/**
 * @defgroup gAPI API of the library
 *  
 * @defgroup gINTERNALS Internal classes, functions and definitions
 *
 * @defgroup gTEST Tests and benchmarks
 * @todo Currently tests are not included in the doxygen documentation,
 * doxygen does not handle multiple main functions well.
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
 * @page pFormulas Derivations and algorithms
 *
 * @section pFormConv Conventions & definitions
 * - All letters referencing matrices or vectors in formulas are typed in bold.
 * - 'N' denotes the size of the preview window.
 *   
 * @section pFormTOC Further information can be found on the following pages
 * - @ref pProblemDef
 *      - @ref pPDModel
 *      - @ref pPDVarSub
 *      - @ref pPDObj
 *      - @ref pPD_EC (example: @ref pExampleEC)
 *      - @ref pPD_IC
 * - @ref pKKT
 * - @ref pInitGuess
 * - @ref pProjectedHessian
 * - @ref pCholesky
 * - @ref pAddIC
 *      - @ref pCholUp
 *      - @ref pAddICz
 * - @ref pRemoveIC
 *      - @ref pCholDown
 *      - @ref pRemoveICz
 * - @ref pDetails 
 *      - @ref pDetMatrices
 *      - @ref pBounds
 */


/**
@page pProblemDef Problem definition

@section pPDModel Model of the system
    3D linear inverted pendulum is used as an approximate model
    of a humanoid robot.

    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
    \newcommand{\dddot}[1]{{\mathop{#1}\limits^{\vbox to-1.4ex{\kern-2ex \hbox{\normalfont ...}\vss}}}}
    \tilde{\mbm{c}}_{k+1} = \tilde{\mbm{A}}_k\tilde{\mbm{c}}_{k}+\tilde{\mbm{B}}_k\dddot{\mbm{c}}_{k}, \quad
    \dddot{\mbm{c}}_{k} = (\dddot{c}_k^x, \dddot{c}_k^y)
    @f$

    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
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
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
    \hat{\mbm{c}}_{k} = (c_k^x,\dot{c}_k^x,\ddot{c}_k^x,c_k^y,\dot{c}_k^y,\ddot{c}_k^y)
    @f$
    where @f$c^x_k, c^y_k@f$ are coordintes of the center of mass.
\n\n


@section pPDVarSub Variable substitutions

@subsection pX_tilde The first substitution
    After the first variable substitution we get 
    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
    \tilde{\mbm{c}}_{k} = (z_k^x,\dot{c}_k^x,\ddot{c}_k^x,z_k^y,\dot{c}_k^y,\ddot{c}_k^y)
    @f$
    where @f$z^x_k, z^y_k@f$ are coordintes of the ZMP.

    The state and control input matrices are changed accordingly:

    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
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
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
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
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
    \bar{\mbm{c}}_{k} = \bar{\mbm{R}}_k^T \tilde{\mbm{c}}_{k} = 
    (\bar{z}_k^x,\dot{c}_k^x,\ddot{c}_k^x,\bar{z}_k^y,\dot{c}_k^y,\ddot{c}_k^y)
    @f$
\n\n


@section pPDObj Objective function
    Output matrices for position and velocity:

    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
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
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
    \bar{f}(\bar{\mbm{v}}) =
    \left[\hspace{-0.1cm}\begin{array}{c} \bar{\mbm{v}}_c \\ \mbm{v}_u \end{array}\hspace{-0.1cm}\right]^T
    \left[\hspace{-0.1cm}\begin{array}{cc} \tilde{\mbm{H}}_c & \mbm{0} \\ \mbm{0} & \mbm{H}_u \end{array}\hspace{-0.1cm}\hspace{-0.1cm}\right]
    \left[\hspace{-0.1cm}\begin{array}{c} \bar{\mbm{v}}_c \\ \mbm{v}_u \end{array}\hspace{-0.1cm}\right] + 
    \left[\hspace{-0.1cm}\begin{array}{c} \bar{\mbm{v}}_c \\ \mbm{v}_u \end{array}\hspace{-0.1cm}\right]^T
    \left[\hspace{-0.1cm}\begin{array}{c} \bar{\mbm{g}}_c \\ \mbm{0} \end{array}\hspace{-0.1cm}\right] 
    @f$

    where 

    @f$\newcommand{\mbm}[1]{\mbox{\boldmath $#1$}} \bar{\mbm{v}}_c @f$
    is a column vector containing state vectors and
    @f$\newcommand{\mbm}[1]{\mbox{\boldmath $#1$}} \bar{\mbm{v}}_u @f$
    is a column vector containing control inputs.

    or

    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
    \newcommand{\dddot}[1]{{\mathop{#1}\limits^{\vbox to-1.4ex{\kern-2ex \hbox{\normalfont ...}\vss}}}}

    f(\mbm{v}) = \frac{\gamma}{2}\sum_{k=0}^{N-1}\left(\dddot{\mbm{c}}_k^T\dddot{\mbm{c}}_k\right) + 
    \frac{\alpha}{2}\sum_{k=1}^{N}\left(\dot{\mbm{c}}_k^T\dot{\mbm{c}}_k\right) + 
    \frac{\beta}{2}\sum_{k=1}^{N}\left(\mbm{z}_k^T\mbm{z}_k - 2\mbm{z}_k^T\mbm{z}^{\mbox{ref}}_k\right),
    @f$

    @anchor pGains
    where @f$\alpha, \beta, \gamma > 0@f$ are gains.


    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
    \frac{\beta}{2}\mbm{z}_k^T\mbm{z}_k - \beta\mbm{z}_k^T\mbm{z}^{\mbox{ref}}_k =

    \bar{\mbm{c}}_k^T\frac{\beta}{2}\mbm{C}_p^T\mbm{C}_p\bar{\mbm{c}}_k -
    \bar{\mbm{c}}_k^T\underbrace{\bar{\mbm{R}}^T_k\beta\mbm{C}_p^T\mbm{z}^{\mbox{ref}}_k}_{\bar{\mbm{q}}_k},  \\

    \bar{\mbm{v}}_c = \left[\begin{array}{c} \bar{\mbm{c}}_1 \\ \vdots \\ \bar{\mbm{c}}_N \end{array} \right], \quad

    \bar{\mbm{g}}_c = \left[\begin{array}{c} -\bar{\mbm{q}}_1 \\ \vdots \\ -\bar{\mbm{q}}_N \end{array} \right], \quad

    \bar{\mbm{v}} = \left[\begin{array}{c} \bar{\mbm{v}}_c \\ \mbm{v}_u \end{array} \right]
    @f$


    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
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
      \begin{array}{ccc}
        \frac{\beta}{2} & 0                 & 0 \\ 
        0               & \frac{\alpha}{2}  & 0 \\ 
        0               & 0                 & r \\ 
      \end{array}
    \right]
    @f$

    @anchor RegFactor
    Here @f$r@f$ is a regularization factor, which makes the matrix nonsingular.
\n\n


@section pPD_EC Equality constraints
    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
    \bar{\mbm{E}}_c\bar{\mbm{v}}_c + \tilde{\mbm{E}}_u\mbm{v}_u = \bar{\mbm{e}}, 
    @f$

    where @f$\newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}
    \bar{\mbm{e}} = (-\mbm{A}\bar{\mbm{R}}_0\bar{\mbm{c}}_0, \mbm{0}, \dots, \mbm{0})@f$,

    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
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
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
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
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
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
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
    \mbm{H}, \mbm{E}^T, \mbm{g}, \mbm{\bar{e}} 
    @f$.

    We assume, that an initial guess satisfying all constraints is given
    (instructions on how to generate it are given in section '@ref pInitGuess')
    Our goal is to find delta between initial guess and optimal point.
    From the system presented above we can derive:

    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
    \frac{1}{2} \mbm{E} \mbm{H}^{-1} \mbm{E}^T \mbm{\nu} = 
        \mbm{S} \mbm{\nu} = 
        \mbm{E} (-\frac{1}{2} \mbm{H}^{-1} \mbm{g} - \mbm{x}_{init}) = \mbm{s}\\
    \Delta\mbm{x} = -\frac{1}{2} \mbm{H}^{-1} \mbm{g} - \mbm{x}_{init} - 
        \frac{1}{2} \mbm{H}^{-1} \mbm{E}^T \mbm{\nu}
    @f$

    Here
    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
    \mbm{S}
    @f$
    is a projected Hessian, its structure is described in section '@ref pProjectedHessian'.

    Note, that
    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
    \frac{1}{2} \mbm{H}^{-1} \mbm{g}
    @f$
    is constant. Also, since Hessian is diagonal its invertion is trivial and
    multiplication of inverted Hessian by any vector is O(N).
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

    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
    \tilde{\mbm{c}}_1 = \mbm{A}\tilde{\mbm{c}}_0 + \tilde{\mbm{B}}\mbm{u}_0,  \\
    \tilde{\mbm{c}}_2 = \mbm{A}\tilde{\mbm{c}}_1 + \tilde{\mbm{B}}\mbm{u}_1,  \\
    \tilde{\mbm{c}}_3 = \mbm{A}\tilde{\mbm{c}}_2 + \tilde{\mbm{B}}\mbm{u}_2,  \\
    \tilde{\mbm{c}}_4 = \mbm{A}\tilde{\mbm{c}}_3 + \tilde{\mbm{B}}\mbm{u}_3,  \\
    \tilde{\mbm{c}}_5 = \mbm{A}\tilde{\mbm{c}}_4 + \tilde{\mbm{B}}\mbm{u}_4.  
    @f$

    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
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
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
      \mbm{e} = 
      \left[
        \begin{array}{c}
          -\mbm{A}\tilde{\mbm{c}}_0 \\ \mbm{0} \\ \mbm{0} \\ \vdots \\ \mbm{0} 
        \end{array}
      \right]. 
    @f$

    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
    \bar{\mbm{R}}_1\bar{\mbm{c}}_1 = \mbm{A}\bar{\mbm{R}}_0\bar{\mbm{c}}_0 + \tilde{\mbm{B}}\mbm{u}_0,  \\
    \bar{\mbm{R}}_2\bar{\mbm{c}}_2 = \mbm{A}\bar{\mbm{R}}_1\bar{\mbm{c}}_1 + \tilde{\mbm{B}}\mbm{u}_1,  \\
    \bar{\mbm{R}}_3\bar{\mbm{c}}_3 = \mbm{A}\bar{\mbm{R}}_2\bar{\mbm{c}}_2 + \tilde{\mbm{B}}\mbm{u}_2,  \\
    \bar{\mbm{R}}_4\bar{\mbm{c}}_4 = \mbm{A}\bar{\mbm{R}}_3\bar{\mbm{c}}_3 + \tilde{\mbm{B}}\mbm{u}_3,  \\
    \bar{\mbm{R}}_5\bar{\mbm{c}}_5 = \mbm{A}\bar{\mbm{R}}_4\bar{\mbm{c}}_4 + \tilde{\mbm{B}}\mbm{u}_4. 
    @f$

    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
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
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
      \mbm{e} = 
      \left[
        \begin{array}{c}
          -\mbm{A}\bar{\mbm{R}}_0\bar{\mbm{c}}_0 \\ \mbm{0} \\ \mbm{0} \\ \vdots \\ \mbm{0} 
        \end{array}
      \right]. 
    @f$
 */


/**
 * @page pProjectedHessian Projected Hessian matrix

    In order to solve @ref pKKT we have to form Schur complement (projected Hessian):

    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
    \mbm{S} = \frac{1}{2}\mbm{E}\mbm{H}^{-1}\mbm{E}^T = \frac{1}{2}\left[\begin{array}{cc}\bar{\mbm{E}}_c  \tilde{\mbm{E}}_u\end{array}\right]
    \left[\begin{array}{cc}\tilde{\mbm{H}}_c & \mbm{0} \\ \mbm{0} & \mbm{H}_u\end{array}\right]
    \left[\begin{array}{c}\bar{\mbm{E}}_c^T \\ \tilde{\mbm{E}}_u^T \end{array}\right] 

    = \frac{1}{2}\bar{\mbm{E}}_c\tilde{\mbm{H}}_c^{-1}\bar{\mbm{E}}_c^T + \frac{1}{2}\tilde{\mbm{E}}_u\mbm{H}_u^{-1}\tilde{\mbm{E}}_u^T. 
    @f$

    For @f$N = 4@f$ we have.

    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
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
      \right]. 
    @f$

    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
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
      \right]  \\

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

    where @f$\newcommand{\mbm}[1]{\mbox{\boldmath $#1$}} 
    \mbm{M}_{ii} = 
    \bar{\mbm{R}}_i\tilde{\mbm{Q}}^{-1}\bar{\mbm{R}}_i^T = 
    \tilde{\mbm{Q}}^{-1}@f$ 
    (due to the special structure of 
    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}
    \tilde{\mbm{Q}}^{-1}
    @f$ and 
    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}
    \bar{\mbm{R}}_i
    @f$).

    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
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
    @f$\newcommand{\mbm}[1]{\mbox{\boldmath $#1$}} 
    \tilde{\mbm{P}} = \tilde{\mbm{B}}\mbm{P}^{-1}\tilde{\mbm{B}}^T@f$.

    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
      2\mbm{S}_{11} = \mbm{M}_{11} + \tilde{\mbm{P}},  \\
      2\mbm{S}_{kk} = \mbm{A}\mbm{M}_{k-1,k-1}\mbm{A}^T + \mbm{M}_{kk} + \tilde{\mbm{P}},  \\
      2\mbm{S}_{k,k+1} = \mbm{S}_{k+1,k}^T = -\mbm{M}_{kk}\mbm{A}^T. 
    @f$

    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
      2\mbm{S}_{11} = \tilde{\mbm{Q}}^{-1} + \tilde{\mbm{P}},  \\
      2\mbm{S}_{kk} = \mbm{A}\tilde{\mbm{Q}}^{-1}\mbm{A}^T + \tilde{\mbm{Q}}^{-1} + \tilde{\mbm{P}},  \\
      2\mbm{S}_{k,k+1} = \mbm{S}_{k+1,k}^T = -\tilde{\mbm{Q}}^{-1}\mbm{A}^T. 
    @f$

    Hence, the matrix 
    @f$\newcommand{\mbm}[1]{\mbox{\boldmath $#1$}} \mbm{S}@f$ 
    is constant (if 
    @f$\newcommand{\mbm}[1]{\mbox{\boldmath $#1$}} \mbm{A}@f$ and 
    @f$\newcommand{\mbm}[1]{\mbox{\boldmath $#1$}} \mbm{B}@f$ 
    do not change). 
*/


/**
 * @page pCholesky Cholesky decomposition of projected Hessian
    Once projected Hessian is formed we can use Cholesky decomposition
    @f$\newcommand{\mbm}[1]{\mbox{\boldmath $#1$}} \mbm{S} = \mbm{L}\mbm{L}^T@f$.
    to obtain Langrange multipliers.

    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
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
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
    \mbm{S}_{11} = \mbm{L}_{11}\mbm{L}_{11}^T,  \\
    \mbm{S}_{12} = \mbm{S}_{21}^T = \mbm{L}_{11}\mbm{L}_{21}^T, \quad \mbm{L}_{21}^T = \mbm{L}_{11}^{-1}\mbm{S}_{12},   \\
    \mbm{S}_{22} = \mbm{L}_{21}\mbm{L}_{21}^T + \mbm{L}_{22}\mbm{L}_{22}^T, \quad \dots 
    @f$
     
    In the second step 
    @f$\newcommand{\mbm}[1]{\mbox{\boldmath $#1$}} \mbm{L}_{21}^T@f$ 
    is computed by forward substitution, and in the third step, forming 
    @f$\newcommand{\mbm}[1]{\mbox{\boldmath $#1$}} \mbm{L}_{22}@f$ 
    requires the computation of the Cholesky factors of 
    @f$\newcommand{\mbm}[1]{\mbox{\boldmath $#1$}} \mbm{S}_{22} - \mbm{L}_{21}\mbm{L}_{21}^T@f$. 
\n\n
 */



/**
 * @page pAddIC Adding inequality constraints
    This page describes changes in @ref pKKT "KKT system" after addition of
    inequality constraint to the active set.

@section pCholUp Update of Cholesky factor

    Let 
    @f$\newcommand{\mbm}[1]{\mbox{\boldmath $#1$}} \mbm{a}_i^T@f$ 
    be the normal to the i-th
    inequality constraint (assumed to be a simple bound). Define the matrix 

    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}
    \mbm{C} = \left[\begin{array}{c} \mbm{E} \\ \mbm{A}_{W}\end{array}\right]
    @f$

    where the rows of 
    @f$\newcommand{\mbm}[1]{\mbox{\boldmath $#1$}} \mbm{A}_W@f$ 
    contain the normals to the inequality constraints in the working set.

    If a new constraint must be added to the active set, then the projected 
    Hessian matrix must be updated:

    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
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

    In general the last line is

    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
    \mbm{s_a}^T = 
    \frac{1}{2}
    \left[
        \mbm{a}^T \mbm{H}^{-1} \mbm{C}^T \quad
        \mbm{a}^T \mbm{H}^{-1} \mbm{A}_W^T \quad
        \mbm{a}^T \mbm{H}^{-1} \mbm{a}
    \right]
    @f$

    Note that 
    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}
    \frac{1}{2}
    \mbm{a}^T \mbm{H}^{-1} \mbm{A}_W^T
    @f$
    is a vector of zeros.
    While 
    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}
    \frac{1}{2}
    \mbm{a}^T \mbm{H}^{-1} \mbm{a} = \frac{1}{2\beta}
    @f$
    is a number.

    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}
    \frac{1}{2}
    \mbm{a}^T \mbm{H}^{-1} \mbm{C}^T 
    @f$
    selects and scales one column of E, this column corresponds to
    ZMP coordinates and can have at most 4 non-zero elements.

    The total number of non-zero elements in the new row of projected Hessian
    is 5 or 3 (for the last state in the preview window).
\n\n

@subsection pCholUpAlg Algorithm of Cholesky factor update
@verbatim
Input:
    m_e % the number of equality constraints
    m_a % the current cardinality of the active set
    L   % Cholesky factor
    s_a % a row added to the projected Hessian

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
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
    \mbm{S}\mbm{\nu} = \mbm{L} \mbm{L}^T \mbm{\nu} = \mbm{L} \mbm{z} = \mbm{s}
    @f$

    When a constraint is added to the active set, there is no need to 
    perform full forward substitution in order to form 
    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\mbm{\nu}
    @f$
    (but full backward substitution is still required).

    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
      \mbm{s}^{+} = -\left[\begin{array}{c} \mbm{C} \\ \mbm{a}_i^T\end{array}\right]
      \left((\mbm{x}+\alpha\Delta\mbm{x})+\mbm{H}^{-1}\mbm{g}\right) 
      = -\left[\begin{array}{c} \mbm{C}\mbm{x} + \mbm{C}\mbm{H}^{-1}\mbm{g} \\ 
          \mbm{a}_i^T(\mbm{x}+\alpha\Delta\mbm{x}) + \mbm{a}_i^T\mbm{H}^{-1}\mbm{g}\end{array}\right]
      = \left[\begin{array}{c} \mbm{s} \\ s_n \end{array}\right]. 
    @f$

    Note that 
    @f$\newcommand{\mbm}[1]{\mbox{\boldmath $#1$}} \alpha\mbm{C}\Delta\mbm{x} = \mbm{0}@f$, 
    because 
    @f$\newcommand{\mbm}[1]{\mbox{\boldmath $#1$}} \Delta\mbm{x}@f$ 
    is in the null space of the normals to the active constraints (stored in 
    @f$\newcommand{\mbm}[1]{\mbox{\boldmath $#1$}} \mbm{C}@f$). 
    Hence, given 
    @f$\newcommand{\mbm}[1]{\mbox{\boldmath $#1$}} \mbm{s}@f$, 
    computing
    @f$\newcommand{\mbm}[1]{\mbox{\boldmath $#1$}} \mbm{s}^{+}@f$ 
    amounts to performing two multiplications plus one addition (note that
    @f$\newcommand{\mbm}[1]{\mbox{\boldmath $#1$}} \mbm{H}^{-1}\mbm{g}@f$ 
    is constant, and 
    @f$\newcommand{\mbm}[1]{\mbox{\boldmath $#1$}} \mbm{x}+\alpha\Delta\mbm{x}@f$ 
    is already formed).
    \n\n


    Now consider the forward substitution 

    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
      \underbrace{\left[\begin{array}{cc} \mbm{L} & \mbm{0} \\ \mbm{l}^T & \ell \end{array}\right]}_{\mbm{L}^{+}}
      \underbrace{\left[\begin{array}{c} \mbm{z} \\ z_n \end{array}\right]}_{\mbm{z}^{+}} = 
      \underbrace{\left[\begin{array}{c} \mbm{s} \\ s_n \end{array}\right]}_{\mbm{s}^{+}}, 
    @f$

    where 
    @f$\newcommand{\mbm}[1]{\mbox{\boldmath $#1$}} \\
    \left[\begin{array}{cc} \mbm{l}^T  \ell \end{array}\right]@f$ 
    is an appended row.
    These are two equations

    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
    \mbm{L}\mbm{z} = \mbm{s}  \\
    \mbm{l}^T\mbm{z} + bz_n = s_n 
    @f$

    From the second one we can compute (note that @f$\ell\neq0@f$)

    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
    z_n = \frac{s_n - \mbm{l}^T\mbm{z}}{\ell}.
    @f$

    Hence, forming 
    @f$\newcommand{\mbm}[1]{\mbox{\boldmath $#1$}} \mbm{z}^{+}@f$ 
    amounts to performing one dot product.
 */

/**
 * @page pRemoveIC Removing inequality constraints
    This page describes changes in @ref pKKT "KKT system" after removal of
    inequality constraint from the active set.
    \n\n

@section pCholDown Downdate of Cholesky factor
    Imagine, that we have selected an inequality constraint for removal,
    then corresponding line and column must be removed from matrix
    @f$\newcommand{\mbm}[1]{\mbox{\boldmath $#1$}} \mbm{S}@f$. We can
    represent this by moving these lines to the end and to the right of
    the matrix using permutation matrix:

    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}} 
    \mbm{E}_{perm} \mbm{S} \mbm{E}_{perm}^T = 
    \mbm{E}_{perm} \mbm{L} \mbm{L}^T \mbm{E}_{perm}^T = 
    (\mbm{E}_{perm} \mbm{L}) (\mbm{E}_{perm} \mbm{L})^T
    @f$

    But matrix 
    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}} 
    \mbm{E}_{perm} \mbm{L}
    @f$
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
    All elements of @f$\newcommand{\mbm}[1]{\mbox{\boldmath $#1$}} \mbm{z}@f$
    starting from the position corresponding to the last non-zero (diagonal)
    element of removed row must be updated.

    Consider the following situation (based on formulas derived in section 
    '@ref pAddICz'):

    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
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
    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
    \mbm{s}
    @f$
    are computed on demand and there is no need to alter it.
    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
    \mbm{l}^T_u, \mbm{z}_u
    @f$
    must be updated.
    Vector
    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
    \mbm{l}^T
    @f$
    is not affected by update (see section '@ref pCholDown'). Hence the following
    number stays constant:

    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
    z_n \ell_n + \ell_{rem} z_{rem} + \mbm{l}^T_u \mbm{z}_u = s_n - \mbm{l}^T\mbm{z} = z_{n,const}
    @f$

    
    After Cholesky factor was updated we get:

    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
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
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
    z_{n,new} = \frac{z_{n,const} - \mbm{l}^T_{u,new} \mbm{z}_{u,new} }{\ell_{new}}
    @f]

    Note, that computation of new elements @f$z_{n,new}@f$ must be started from the
    element corresponding to the first changed line of L and continued towards the
    end of z.
 */


/**
 * @page pDetails Minor implemetation details
@section pDetNotes Notes
    It is implicitly supposed that we have 6 state variables and 2 control variables.

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
      N*#NUM_VAR.
\n\n

@section pBounds Implementing bounds

    Consider the variable

    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
    \mbm{x} = \left[\begin{array}{c} x_1 \\ x_2 \\ x_3 \\ x_4 \\ x_5 \\ x_6 \end{array}\right].
    @f$
     
    Suppose that variables @f$x_1@f$ and @f$x_4@f$ have simple bounds (and the rest of the variables are are
    not subject to inequality constraints), i.e., 

    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
      lb_1 \leq x_1 \leq ub_1  \\
      lb_2 \leq x_2 \leq ub_2  
    @f$

    The above four inequality constraints can be written as

    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
    \mbm{a}_1^T\mbm{x} \leq ub_1  \\
    \mbm{a}_1^T\mbm{x} \geq lb_1  \\
    \mbm{a}_2^T\mbm{x} \leq ub_2  \\
    \mbm{a}_2^T\mbm{x} \geq lb_1, 
    @f$

    where

    @f$
    \newcommand{\mbm}[1]{\mbox{\boldmath $#1$}}\\
      \mbm{a}_1^T = \left[\begin{array}{cccccc} 1 & 0 & 0 & 0 & 0 & 0 \end{array}\right]  \\
      \mbm{a}_2^T = \left[\begin{array}{cccccc} 0 & 0 & 0 & 1 & 0 & 0 \end{array}\right]. 
    @f$

    Note that both 
    @f$\newcommand{\mbm}[1]{\mbox{\boldmath $#1$}} \mbm{a}_1^T\mbm{x} \leq ub_1@f$ 
    and 
    @f$\newcommand{\mbm}[1]{\mbox{\boldmath $#1$}} \mbm{a}_1^T\mbm{x} \geq lb_1@f$ 
    can not be in the working set at the same time (because if we are on one of the bounds we 
    can not be on the other one).

    @attention Note, that if we do not distinguish lower and upper bounds, Lagrange 
    multipliers can be negative.
 */

#endif /*DOXYGEN_H*/
 
