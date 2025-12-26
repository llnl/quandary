# Quandary: Optimal Control for Open and Closed Quantum Systems

# Introduction
Quandary numerically simulates and optimizes the time evolution of closed and open quantum systems. The
underlying dynamics are modelled by either Schroedinger's equation (for closed systems), or Lindblad's master equation (for open systems that interact with the environment). Quandary solves the respective ordinary differential equation (ODE) numerically by applying a time-stepping integration scheme, and applies a gradient-based optimization
scheme to determine optimal control pulses that drive the quantum system to a desired target.
The target can be a unitary, i.e. optimizing for pulses that
realize a logical quantum operation, or state preparation that aims to drive the quantum system from one (or multiple) initial state to a desired target state, such as for example the ground state of zero energy level, or for the creation of entangled states.

Quandary is designed to solve optimal control problems in larger (potentially open) quantum systems, targeting modern high performance computing (HPC) platforms. Quandary utilizes distributed memory computations using the message passing paradigm that enables scalability to large number of compute cores. Implemented in C++, Quandary is portable and its object-oriented implementation allows developers to extend the predefined setup to suit their particular simulation and optimization requirements. For example, customized gates for Hamiltonian simulations can easily be added to supplement Quandary’s predefined gate set.
The Python interface allows for greater flexibility where custom Hamiltonian models can be used.

This document outlines the mathematical background and underlying equations, and summarizes their
implementation and usage in Quandary. Also refer to publications [@guenther2021quandary] [@guenther2021quantum].

# Model equation {#sec:model}
Quandary models composite quantum systems consisting of $Q$ subsystems ("oscillators", "qubits"), with $n_k$ energy levels for the
$k$-th subsystem, $k=0,\dots,Q-1$. The Hilbert space dimension is $N = \prod_{k=0}^{Q-1} n_k$.

For a **closed quantum system** (no environmental interactions), the quantum state is described by a complex-valued vector $\psi\in\C^N$, with $\|\psi\| = 1$. For a given initial state $\psi(t=0)$, the evolution of the state vector is modelled through **Schroedinger's equation**

\begin{align} \label{eq:schroedinger}
  \dot{\psi}(t)  = -i H(t) \psi(t)
\end{align}

For **open quantum systems** (taking interactions with the environment into account), the quantum state is described by the density matrix $\rho\in \C^{N\times N}$ (hermitian, trace one, positive semi definite). In that case, Quandary models the the quantum dynamics via **Lindblad's master equation**:

\begin{align}\label{mastereq}
  \dot{\rho}(t) = &-i(H(t)\rho(t) - \rho(t)H(t)) + \Ell(\rho(t)),
\end{align}


*Note:* In the remainder, the quantum state will mostly be denoted by $\rho$, which, depending on the context, either relates to the density matrix solving Lindblad's equation, or to the state vector $\psi\in \C^N$ solving Schroedinger's equation. Distinction will only be made explicit where necessary.

The **default Hamiltonian** in Quandary models superconducting (transmon) qubits, decomposing the Hamiltonian matrix into a constant system part and a time-dependent control part that models the action of control fields applied to each subsystem:

\begin{align}
  H(t) = H_d + H_c(t) \quad \text{where} \quad
  H_d &:= \sum_{k=0}^{Q-1} \omega_k a_k^{\dagger}a_k- \frac{\xi_k}{2} a_k^{\dagger}a_k^{\dagger}a_k a_k  + \sum_{l> k} J_{kl} \left( a_k^\dagger a_l + a_k a_l^\dagger \right) - \sum_{l> k}\xi_{kl} a_{k}^{\dagger}a_{k}   a_{l}^{\dagger} a_{l} \\
  H_c(t) &:= \sum_{k=0}^{Q-1} f^k(t) \left(a_k + a_k^\dagger \right)
\end{align}

where $\omega_k\geq 0$ denotes $0 \rightarrow 1$ transition frequencies for each oscillator $k$, $\xi_k\geq 0$ are the self-Kerr coefficients. Couplings can be specified through the cross resonance coefficients $J_{kl}\geq 0$ ("dipole-dipole interaction") or through $\xi_{kl}\geq 0$ ("zz-coupling"). 
Here, $a_k\in \C^{N\times N}$ denotes the lowering operator acting on subsystem $k$.
The control pulses $f^k(t)$ can be either specified or optimized for, compare section [Control pulse parameterization](#sec:controlpulses). **Custom system and control Hamiltonian operators** can be specified through Quandary's python interface.

The **default Lindbladian** operator $\Ell(\rho(t))$ is of the form

\begin{align} \label{eq:collapseop}
  \Ell(\rho(t)) = \sum_{k=0}^{Q-1} \sum_{l=1}^2 \Ell_{lk} \rho(t)
  \Ell_{lk}^{\dagger} - \frac 1 2 \left( \Ell_{lk}^{\dagger}\Ell_{lk}
  \rho(t) + \rho(t)\Ell_{lk}^{\dagger} \Ell_{lk}\right)
\end{align}

where the collapse operators $\Ell_{lk}$ model decay and dephasing in each subsystem $k$ with

- Decay  ("$T_1$"): $\Ell_{1k} = \frac{1}{\sqrt{T_1^k}} a_k$
- Dephasing  ("$T_2$"): $\Ell_{2k} = \frac{1}{\sqrt{T_2^k}} a_k^{\dagger}a_k$

<!-- Note that the main choice here is which equation should be solved for and which representation of the quantum state will be used (either Schroedinger with a state vector $\psi \in \C^N$, or Lindblad's equation for a density matrix $\rho \in \C^{N\times N}$). In the C++ configuration file, this choice is determined through the option `collapse_type`, where `none` will result in Schroedinger's equation and any other choice will result in Lindblad's equation being solved for. Further note, that choosing `collapse_type` $\neq$ `none`, together with a collapse time $T_{l}^k = 0.0$ will omit the evaluation of the corresponding term in the Lindblad operator $\eqref{eq:collapseop}$ (but will still solve Lindblad's equation for the density matrix). In the python interface, Lindblad's solver is enabled by passing decay and decoherence times `T1` and `T2` per oscillator to the Quandary object. -->

## Rotational frame 
Quandary uses the rotating wave approximation to slow down the time scale of the quantum dynamics. The user can specify the rotation frequencies $\omega_k^r$ for each oscillator. Under the rotating frame wave approximation, the Hamiltonians are transformed to

\begin{align}
  \tilde{H}_d(t) &:= \sum_{k=0}^{Q-1} \left(\omega_k - \omega_k^{r}\right)a_k^{\dagger}a_k- \frac{\xi_k}{2}
  a_k^{\dagger}a_k^{\dagger}a_k a_k
   - \sum_{l> k} \xi_{kl} a_{k}^{\dagger}a_{k}   a_{l}^{\dagger} a_{l} \notag \\
   & + \sum_{l>k} J_{kl} \cos(\eta_{kl}t) \left(a_k^\dagger a_l + a_k a_l^\dagger\right) + iJ_{kl} \sin(\eta_{kl}t)\left(a_k^\dagger a_l - a_k a_l^\dagger\right) \label{eq:Hd_rotating} \\
   %
   \tilde{H}_c(t) &:= \sum_{k=0}^{Q-1} p^k(t) (a_k +
   a_k^{\dagger}) + i q^k(t)(a_k - a_k^{\dagger})
    \label{eq:Hc_rotating}
\end{align}

where $\eta_{kl} := \omega_k^{r} - \omega_l^{r}$ are the differences in rotational frequencies between subsystems. Note that the dipole-dipole coupling is time-dependent if $\eta_{kl} \neq 0$. Using the rotating wave approximation, the rotating-frame control pulses $p^k(t)$ and $q^k(t)$ relate to the laboratory frame control pulse through 

\begin{align}
  f^k(t) = 2\mbox{Re}\left(d^k(t)e^{i\omega_k^r t}\right), \quad d^k(t) = p^k(t) + i q^k(t)
\end{align}


## Essential and non-essential energy levels {#sec:essential}
It is recommended to model the system with more energy levels than those that will be occupied during the dynamical evolution, in order to prevent leakage out of the computational subspace (modelling the infinite dimensional system with more accuracy by including more levels) and/or to allow the system to transition through higher energy levels in order to achieve a final-time target faster. In that case, *essential* levels, $n_k^e$, denote the computational subspace, e.g. $n_k^e = 2$ for qubits, whereas the total number of energy levels $n_k \geq n_k^e$ could be larger, e.g. $n_k=3$. Any non-essential level is considered a *guard* level.

## Control pulse parameterization {#sec:controlpulses}
The time-dependent rotating-frame control pulses $d^k(t) = p^k(t) + iq^k(t)$ are parameterized using Bsplines with $N_s^k$ basis functions $B_s(t)$, that act as an envelope for $N_f^k$ carrier waves:

\begin{align}\label{eq:spline-ctrl}
  d^k(\vec{\alpha}^k,t) = \sum_{f=1}^{N_f^k} e^{i\Omega^k_ft} \sum_{s=1}^{N_s^k} \alpha_{s,f}^k B_s(t) , \quad \alpha_{s,f}^k = \alpha_{s,f}^{k(1)} + i \alpha_{s,f}^{k(2)} \in \C
\end{align}

Using trigonometric identities, the real and imaginary part of the rotating-frame control $d^k(\vec{\alpha}^k,t) = p^k(\vec{\alpha}^k,t) + iq^k(\vec{\alpha}^k,t)$ can be written as

\begin{align}
  p^k(\vec{\alpha}^k,t) &= \sum_{f=1}^{N_f^k} \cos(\Omega_f^k t) B^{(1)}(t) 
    - \sin(\Omega_f^k t) B^{(2)}(t) \\
  q^k(\vec{\alpha}^k,t) &= \sum_{f=1}^{N_f^k} \sin(\Omega_f^k t) B^{(1)}(t) + \cos(\Omega_f^k t)B^{(2)}(t) 
\end{align}

where $B^{(1)}(t) = \sum_{s=1}^{N_s^k} \alpha^{k(1)}_{s,f} B_s(t)$ and $B^{(2)}(t) = \sum_{s=1}^{N_s^k} \alpha^{k(2)}_{s,f} B_s(t)$ evaluate the splines using the control coefficients $\alpha_{f,s}^{k(1)}, \alpha_{f,s}^{k(2)}\in \R$. 
By default, the basis functions are piecewise quadratic B-spline polynomials with compact support, centered on an equally spaced grid in time. To instead use a piecewise constant (0th order) Bspline basis, see Section [0-th order Bspline basis functions](#sec:bspline-0).

The control parameter vector $\boldsymbol{\alpha} = (\alpha_{f,s}^{k(i)})$ (*design* variables) can be either specified or can be optimized for in order to realize a desired system behavior (Section [The Optimal Control Problem](#sec:optim)).  

### Carrier wave frequencies
The rotating-frame carrier wave frequencies $\Omega^k_f \in \R$ should be chosen to trigger intrinsic system resonance frequencies. For example, when $\xi_{kl} << \xi_k$, the intrinsic qubit transition frequencies are $\omega_k - n\xi_k$. Thus by choosing $\Omega^k_f = \omega_k-\omega_k^r - n \xi_k$ in the rotating frame, one triggers transition between energy levels $n$ and $n+1$ in subsystem $k$. Choosing effective carrier wave frequencies is quite important for optimization performance, particulary when qubit interactions are desired, such as when optimizing for a CNOT gate. Using the python interface for Quandary, the carrier wave frequencies $\Omega^k_f$ are automatically computed based on an eigenvalue decomposition of the system Hamiltonian. For the C++ code, it is recommended to follow [@petersson2021optimal] for details on how to choose them effectively.  

### Alternative control parameterization based on B-spline amplitudes and time-constant phases
As an alternative parameterization, the user can choose to parameterize only the *amplitudes* of the control pulse with 2nd order B-splines, adding a time-constant phase per carrierwave:

\begin{align}
  d(t) = \sum_f e^{i\Omega_f t} a_f(t)e^{ib_f} \quad \text{where} \quad a_f(t) = \sum_s \alpha_{f,s} B_s(t) \\
  \Rightarrow d(t)= \sum_f\sum_s \alpha_{f,s}B_s(t)e^{i\Omega_ft + b_f}
\end{align}

where the control parameters are $b_f\in [-\pi, \pi]$ (phases for each carrier wave) and the amplitudes $\alpha_{f,s}\in \R$ for $s=1,\dots, N_s$, $f=1,\dots, N_f$. Hence for $Q$ oscillators, we have a total of $\sum_q (N_s^q + 1) N_f^q$ control parameters.

### Piecewise constant control parameterization (0-th order Bspline basis functions) {#sec:bspline-0}
Piecewise constant control pulses can be generated by using 0-th order Bspline basis functions. In this case, it is recommended to set the carrier wave frequencies to zero in the rotating frame. When optimizing with 0-th order B-spline basis functions, strong variations between consecutive control amplitudes can be suppressed by enabling the total variation penalty term, compare Section [Regularization, penalty terms, and leakage prevention](#sec:penalty).


# The Optimal Control Problem {#sec:optim}
In the most general form, Quandary can solve the following optimization problem:

\begin{align}\label{eq:minproblem}
  \min_{\boldsymbol{\alpha}} J\left(\{\rho^{target}_i, \rho_i(T)\}\right) +  \mbox{Regularization} + \mbox{Penalty}
\end{align}

where the (single or multiple) final-time states $\rho_i(T)$ solve either Lindblad's master equation $\eqref{mastereq}$ or Schroedinger's equation $\eqref{eq:schroedinger}$ in the rotating frame for (one or multiple) initial conditions $\rho_i(0)$, as specified in Section [Initial conditions](#sec:initcond), $i=1,\dots, n_{init}$. The first term in $\eqref{eq:minproblem}$ minimizes an objective function $J$ (see Section [Objective function](#sec:objectivefunctionals)) that quantifies the discrepancy between the realized states $\rho_i(T)$ at final time $T$ driven by the current control $\boldsymbol{\alpha}$ and the desired target $\rho^{target}_i$, see Section [Optimization targets](#sec:targets).
The remaining terms are regularization and penalty terms that can be added to stabilize convergence, or prevent leakage, compare Section [Regularization, penalty terms, and leakage prevention](#sec:penalty)

## Objective function {#sec:objectivefunctionals}
The following objective functions can be used for optimization in Quandary (config option `optim_objective`):

\begin{align}
 J_{Frobenius} &= \sum_{i=1}^{n_{init}} \frac{\beta_i}{2} \left\| \rho^{target}_i - \rho_i(T)\right\|^2_F \\
 J_{trace} &=
\begin{cases}
 1 - \sum_{i=1}^{n_{init}} \frac{\beta_i}{w_i} \mbox{Tr}\left((\rho^{target}_i)^\dagger\rho_i(T)\right) & \text{if Lindblad}\\
 1 - \left|\sum_{i=1}^{n_{init}} \beta_i (\psi^{target}_i)^\dagger\psi_i(T)\right|^2 & \text{if Schroedinger}
\end{cases}\\
 J_{measure} &= \sum_{i=1}^{n_{init}} \beta_i \mbox{Tr} \left( N_m \rho(T) \right) \label{eq:Jmeasure}
\end{align}

for default weights default $\beta_i = 1/n_{init}$ that can be chosen to scale different contribution of each initial-to-target state.
$J_{Frobenius}$ measures (weighted average of) the Frobenius norm between target and final states. $J_{Trace}$ measures the (weighted) infidelity in terms of the Hilbert-Schmidt overlap. Here, $w_i = \mbox{Tr}\left(\rho_i(0)^2\right)$ is the purity of the initial state. Both measures are common for optimization towards a unitary gate transformation, for example. $J_{measure}$ is (only) useful when considering unconditional pure-state preparation, see Section [Optimization targets](#sec:targets). Here, $m\in\N$ is a given integer, and $N_m$ is a diagonal matrix with diagonal elements being $|k-m|, k=0,\dots N-1$

<!-- The distinction for the Lindblad vs. Schroedinger solver is made explicit for $J_{trace}$ above. The other two measures apply naturally to either the density matrix version solving Lindblad's equation, or the state vector version solving Schroedinger's equation. -->

### Fidelity
As a measure of optimization success, Quandary reports on the **fidelity** computed from

\begin{align}\label{eq:fidelity}
  F = \begin{cases}
    \frac{1}{n_{init}} \sum_{i=1}^{n_{init}} \mbox{Tr}\left(\left(\rho^{target}_i\right)^\dagger\rho_i(T) \right) & \text{if Lindblad} \\
    \left|\frac{1}{n_{init}} \sum_{i=1}^{n_{init}} (\psi^{target}_i)^\dagger \psi_i(T) \right|^2 & \text{if Schroedinger}
  \end{cases}
\end{align}

The fidelity is an average of Hilbert-Schmidt overlaps of the target states and the evolved states: for the density matrix, the Hilbert-Schmidt overlap is $\langle \rho^{target}, \rho(t)\rangle = \mbox{Tr}\left(\left(\rho^{target}\right)^\dagger\rho(T)\right)$, which is *real* if both states are density matrices (which is always the case in Quandary, see definition of basis matrices). For the state vector (and the Schroedinger solver), the Hilbert-Schmidt overlap is $\langle \psi^{target}, \psi(T)\rangle = (\psi^{target})^{\dagger}\psi$, which is complex. Note that in the fidelity above (and also in the corresponding objective function $J_{trace}$, the absolute value is taken *outside* of the sum, hence relative phases are taken into account.

Further note that this fidelity is averaged over the chosen initial conditions, so the user should be careful how to interpret this number. E.g. if one optimizes for a logical gate while choosing the three initial condition as in Section [Initial conditions](#sec:initcond), the fidelity that is reported during optimization will be averaged over those three initial states, which is not sufficient to estimate the actual average fidelity over the entire space of potential initial states. It is advised to recompute the average fidelity **after** optimization has finished by propagating all basis states $B_{kj}$ to final time $T$ using the optimized control parameter, or by propagating only $N+1$ initial states to get an estimate thereof.

## Optimization targets {#sec:targets}

### Gate optimization
Quandary can be used to design control pulses that realize logical gate operations. Let $V\in \C^{N\times N}$ be the unitary matrix (gate), optimized control pulses drive any initial state $\rho(0)$ to the unitary transformation $\rho^{target} = V\rho(0)V^{\dagger}$ (Lindblad), or, in the Schroedinger case, drive any $\psi(0)$ to $\psi(T) =  V\psi(0)$.
Some default target gates that are readily available, or can be specified from file or through the Python interface. (File format: column-wise vectorization, first all real parts then all imaginary parts.)

Since *any* initial quantum state should be transformed by the control pulses, the corresponding initial conditions must span a basis with $n_{init} = N$ for Schroedinger solver, and $n_{init}=N^2$ for Lindblad solver, see Section [Initial conditions](#sec:initcond). 

Target gates will by default be rotated into the computational frame (Section [Model equation](#sec:model)). Alternatively, the user can specify the rotation of the target gate through the configuration option `gate_rot_freq`.


If guard levels are used ($n_k > n_k^e$, compare Section [Essential and guard levels](#sec:essential)), the gate should be defined in the essential-level dimensions only. Internally, the gate is projected upwards to the full dimensions by inserting identity blocks for rows/columns that correspond to a non-essential levels of the subsystems. Hence, a realization of the gate $\tilde{V}$ will not alter the occupation of higher (non-essential) energy level compared to their initial occupation at $t=0$.

### State preparation {#sec:stateprep}
Quandary can be used to optimize for pulses that drive (one or multiple) initial states to a fixed target state $\rho^{target}$. Depending on the choice of the [Initial conditions](#sec:initcond), this enables pulses for either direct **state-to-state transfer** (by choosing one specific initial condition, $n_{init}=1$), and one specific target state), or **unconditional state preparation** such as qubit reset (by spanning a basis of initial conditions, $n_{init}=N$ or $N^2, and one specific target state). Driving a basis of initial state to a common target will require to couple to a dissipative bath, which should be accounted for in the model setup. For unconditional *pure*-state preparation, it is shown in [@guenther2021quantum] that if one chooses the objective function $J_{measure}$ with corresponding measurement operator $N_m$ (see eq. $\eqref{eq:Jmeasure}$), one can reduce the number of initial conditions to only *one* being an ensemble of all basis states, and hence $n_{init}=1$ independent of $N$. Compare [@guenther2021quantum] for details.

If the target state is *pure*, internal computations are simplified and it is recommended to pass the specific identifier ``pure, i0, i1, i2, ...`` to the Quandary configuration for the optimization target, denoting a pure target state of the form $\psi = |i_0i_1i_2...\rangle$, or $\rho = \psi\psi^\dagger$


## Initial conditions {#sec:initcond}
The initial states $\rho_i(0)$ which are accounted for in the objective function eq. $\eqref{eq:minproblem}$ can be specified with the configuration option `initialcondition`. 


* **Basis states for gate optimization**: $n_{init}=N$ (Schroedinger case), or $n_{init}=N^2$ Lindblad case. For the Schroedinger case, the basis states are the unit vectors $\psi_i(0)=\boldsymbol{e}_i \in \R^N, i=0,\dots N-1$. For the Lindblad's case, the $N^2$ basis density matrices defined in [@guenther2021quantum] are used as initial states. 
In order to uniquely identify the different initial conditions in the Quandary code and in the output files, a
unique index $i \in \{0,\dots, N^2-1\}$ is assigned to each basis state with $B^i := B^{k(i), j(i)}$ with $k(i) := i \,\mbox{mod}\, N$ and $j(i) := \left\lfloor \frac{i}{N} \right\rfloor$
(column-wise vectorization of the matrix of basis matrices).
<br>
For composite systems of multiple subsystems, the user can specify a consecutive list of integer ID's to determine in which of the subsystems the basis states should be spanned. Other subsystems will then be initialized in the ground state only.
<br>
*Note:* The basis states are spanned in the *essential dimensions* of the system, if applicable.

* **Only diagonal density basis matrices**: $n_{init}=N$. For the Lindblad solver, one can choose to propagate only the *diagonal* basis matrices $\boldsymbol{e}_k\boldsymbol{e}_k^\dagger$. For the Schroedinger solver, this option is equivalent to all basis states.
<br>
*Note:* the diagonal states are spanned in the *essential dimensions* of the system, if applicable.

* **Three initial states for gate optimization**: $n_{init}=3$. *Only valid for the Lindblad solver.*
When considering gate optimization with Lindblad's solver, it is shown in [@goerz2014optimal] that it is enough to consider only three specific initial states during optimization ($n_{init}=3$), independent of the Hilbert space dimension. They are readily implemented in Quandary. Note that it is important to choose the weights $\beta_i, i=1,2,3$ in the objective function appropriately to achieve fast convergence.
<br>
Note: The three initial states are spanned in the *full* dimension of the system, including non-essential levels. The theory for gate optimization with three initial states had been developed for considering *only* essential levels (the gate is defined in the same dimension as the system state evolution), and at this point we are not certain if the theory generalizes to the case when non-essential levels are present. It is advised to optimize on the full basis if non-essential levels are present (or work on the theory, and let us know what you find.). The same holds for $N+1$ initial states below.

* **$N+1$ initial states for gate optimization**: $n_{init}=N+1$. *Only valid for the density matrix version, solving Lindblad's master equation.*
The three initial states from above do not suffice to estimate the fidelity of the realized gate (compare [@goerz2014optimal]). Instead, it is suggested in that same paper to choose $N+1$ initial states to compute the fidelity. Those $N+1$ initial states consist of the $N$ diagonal states $B^{kk}$ in the Hilbert space of dimension $N$, as well as the totally rotated state $\rho(0)_2$ from above. 
<br>
Note: The $N+1$ initial states are spanned in the *full* dimension of the system, including non-essential levels, see above for 3-state initialization.

* **Pure initial state for state-to-state transfer**: $n_{init} = 1$. The user can choose a pure initial state of the form $\psi(0) = |i_0, i_1, i_2, ...\rangle$, or $\rho(0) = \psi(0)\psi(0)^\dagger$, through the configuration option ``pure, i0, i1, i2, ...``

* **Arbitrary initial state for state-to-state transformation**: $n_{init}=1$. An arbitrary (non-pure) initial state can be passed to Quandary directly through the Python interface, or can be read from a file in the C++ code. File format: column-wise vectorized density matrix or the state vector, first all real parts, then all imaginary parts. 

* **Ensemble state for unconditional pure-state preparation**: $n_{init}=1$. *Only valid for Lindblad's solver.* When choosing the objective function $J_{measure}$ $\eqref{eq:Jmeasure}$, one can use the ensemble state $\rho_s(0) = \frac{1}{N^2}\sum_{i,j=0}^{N-1} B^{kj}$ as the only initial condition for optimizing for pulses that realize unconditional pure-state preparation, compare [@guenther2021quantum]). To specify the ensemble state in Quandary (C++), one can provide a list of consecutive integer ID's that determine in which of the subsystems the ensemble state should be spanned. Other subsystems will be initialized in the ground state.
<br>
*Note*: The ensemble state will be spanned in the *essential* levels of the (sub)system, if applicable, and will then be lifted up to the full dimension by inserting zero rows and columns.


## Regularization, penalty terms, and leakage prevention {#sec:penalty}
In order to regularize the optimization problem (stabilize optimization convergence), it is advised to add a Tikhonov regularization term to the objective function, by choosing a small $\gamma_1 > 0$

\begin{align}
  \mbox{Tikhonov} = \frac{\gamma_1}{2} \| \bfa \|^2_2
\end{align}

In addition, the following penalty terms can be added to the objective function, if desired:

\begin{align*}
  Penalty &= \frac{\gamma_2}{T} \int_0^T \mbox{Leakage}\left(\rho(t)\right) \, \mathrm{d} t   \hspace{3.5cm} \rightarrow \text{Leakage prevention}\\
         &+  \frac{\gamma_3}{T} \int_0^T \, \| \partial_{tt} \mbox{Pop}(\rho(t)) \|^2 \mathrm{d}t \hspace{3.2cm} \rightarrow \text{State variation penalty} \\
        &+\frac{\gamma_4}{T} \int_0^T \, \sum_k |p^k(\alpha^k,t)|^2 + |q^k(\alpha^k,t)|^2\, dt  \hspace{0.7cm}\rightarrow  \text{Control energy penalty}\\
        &+ \frac{\gamma_5}{2} Var(\vec{\alpha}) \hspace{6.2cm}\rightarrow  \text{Control variation penalty}\\
        &+\frac{\gamma_6}{T} \int_0^T \, w(t)J(t)\, dt  \hspace{4.5cm}\rightarrow  \text{Weighted running cost penalty}\\
\end{align*}

* **Leakage prevention:** Choose a small $\gamma_2 > 0$ to penalize (suppress) leakage into non-essential energy levels (if $n_k^e < n_k$ for at least $k$, compare Sec. [Essential and non-essential energy levels](#sec:essential)). This term penalizes the occupation of all *guard levels* with $\mbox{Leakage}(\rho(t)) = \sum_{r} \| \rho(t)_{rr} \|^2_2$, where $r$ iterates over all indices that correspond to a guard level (i.e., the final (highest) non-essential energy level) of at least one of the subsystems.

* **State variation penalty**: Choose a small $\gamma_3 > 0$ to encourage state evolutions whose populations vary slowly in time by penalizing the second derivative of the populations of the state.

* **Control energy penalty**: Choose a small $\gamma_4 > 0$ to encourage small control pulse amplitudes by penalizing the control pulse energy. This term can be useful if hardware bounds are given for the control pulse amplitudes: Rather than specifying control bounds for the optimization directly, which can lead to convergence deterioration, one can utilize this penalty term to favor control pulses with smaller amplitudes. Compare also [@gunther2023practical] for its usage to determine minimal gate durations.

* **Control variation penalty**: Choose a small $\gamma_5>0$ to penalize variations in control strength between consecutive B-spline coefficients. It is currently only implemented for piecewise zeroth order spline functions, see Section [Zeroth order B-spline basis functions](#sec:bspline-0), where it is useful to prevent noisy control pulses. Referring to the control function representation in $\eqref{eq:spline-ctrl}$, this penalty function takes the form:
$Var(\vec{\alpha}) = \sum_{k=1}^Q Var_k(\vec{\alpha})$ with $Var_k(\vec{\alpha}) = \sum_{f,s}|\alpha_{s,f}^k - \alpha_{s-1,f}^k|^2$.

* **Weighted running cost penalty**: Choose a small $\gamma_6>0$ to penalize a weighted objective function throughout the time domain. The weight is a gaussian centered at final time $T$: $w(t) =
  \frac{1}{a} e^{ -\left(\frac{t-T}{a} \right)^2}$ for a width parameter $0 \leq a \leq 1$. Note, that as $a\to 0$, the weighting function $w(t)$ converges to the Dirac delta distribution with peak at final time $T$, hence reducing $a$ leads to more emphasis on the final time $T$ while larger $a$ penalize the objective function at earlier times $t\leq T$.

Note: All regularization and penalty coefficients $\gamma_i$ should be chosen small enough so that they do not dominate the final-time objective function $J$. This might require some fine-tuning. It is recommended to always add $\gamma_1>0$, e.g. $\gamma_1 = 10^{-4}$, and add other penalties only if needed.

## Optimization algorithm
Quandary utilized Petsc's Toolkit for Advanced Optimization (TAO) package to solve the optimal control problem. In the current setup, Quasi-Newton updates are applied to the control parameters using L-BFGS Hessian approximations. A projected line-search is used to incorporate box constraints for control pulse amplitude bounds $|p^k(t)| \leq c^k_{max}$, $|q^k(t)| \leq c^k_{max}$ via

\begin{align}
  | \alpha_{s,f}^{k(1)}| \leq \frac{c^k_{max}}{\sqrt{2}N_f^k} \quad \text{and} \quad |
  \alpha_{s,f}^{k(2)} | \leq \frac{c^k_{max}}{\sqrt{2}N_f^k}.
\end{align}



# Implementation

## Real-valued and vectorized formulation
Quandary solves the quantum dynamical system in real-valued variables with $q(t) = u(t) + iv(t)$, evolving the real-valued
states $u(t), v(t)\in \R^{M}$ for $M=N$ (if Schroedinger's eq.) or $M=N^2$ (if Lindblad's eq., see below). Particularly, considering Schroedinger's equation first, the real-valued dynamical system is given by 

\begin{align}
  \dot{q}(t) = -iH(t) q(t) \quad \Leftrightarrow \quad &\begin{bmatrix} \dot{u}(t) \\ \dot{v}(t) \end{bmatrix} =
\begin{bmatrix} A(t) & -B(t) \\ B(t) & A(t) \end{bmatrix}
\begin{pmatrix} u(t) \\ v(t) \end{pmatrix}  \\
\notag\\
\text{where} \quad A(t) &= Re(-iH(t)) = Im(H)\notag \\
     B(t) &= Im(-iH(t)) = -Re(H(t))\notag
\label{realvaluedODE}
\end{align}

When solving Lindblad's master equation $\eqref{mastereq}$, Quandary uses the column-wise vectorization relations

\begin{align}
  \text{vec}(AB) &= (I_N\otimes A)\text{vec}(B) = (B^T\otimes I_N)\text{vec}(A)
  \\
  \text{vec}(ABC) &= (C^T\otimes A)\text{vec}(B)
\end{align}

to derive the vectorized master equation for $q(t) := \text{vec}(\rho(t)) \in \C^{N^2}$, and its real-valued formulation as

\begin{align}\label{mastereq_vectorized}
  &\dot{q}(t) = \left(I_N\otimes \left(-i H(t)\right) - \left(-iH(t)\right)^T \otimes I_N  + \vec L \right) q(t) \quad  \text{for} \quad
  \vec L := \sum_{k=0}^{Q-1}\sum_{l=1}^2 \gamma_{lk}
  \left( \Ell_{lk}\otimes \Ell_{lk} - \frac 1 2 \left( I_N\otimes
  \Ell^T_{lk}\Ell_{lk} + \Ell^T_{lk}\Ell_{lk} \otimes I_N \right) \right) \\
  &\Leftrightarrow \quad \begin{bmatrix} \dot{u}(t) \\ \dot{v}(t) \end{bmatrix} =
\begin{bmatrix} A(t) & -B(t) \\ B(t) & A(t) \end{bmatrix}
\begin{pmatrix} u(t) \\ v(t) \end{pmatrix}  \quad \text{where}\quad A(t) = I_N\otimes Re(-i H(t)) - Re(-iH(t))^T \otimes I_N + \vec L \notag \\
& \qquad \qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\quad  B(t) = I_N\otimes Im(-i H(t)) - Im(-iH(t))^T \otimes I_N)  \notag 
\end{align}


The real and imaginary parts of $q(t)$ are stored in blocked manner: For
  $q = u+iv$ with $u,v\in\R^{M}$, a vector of size $2M$ as $q=\begin{bmatrix} u\\v \end{bmatrix}.

## Time-stepping
To solve the resulting real-valued differential equation 

\begin{align*}
\dot q(t) = M(t) q(t) \forall t\in (0,T), \quad \text{with}\quad M(t) =\begin{bmatrix} A(t) & -B(t) \\ B(t) & A(t) \end{bmatrix} 
\end{align*}

Quandary applies a time-stepping integration scheme on a uniform time discretization grid $0=t_0 < \dots t_{N} = T$, with
$t_n = n \delta t$ and $\delta t = \frac{T}{N}$ to approximate the
solution at each discrete time-step $q^{n} \approx q(t_n)$. The default and recommended time-stepping scheme is the Implicit Midpoint Rule `IMR`. The implicit midpoint rule is a second-order accurate, symplectic time-stepping algorithm with Runge-Kutta scheme tableau
$\begin{array}{c|c}
1/2 & 1/2 \\
\hline
& 1
\end{array}$
Given a state $q^n$ at time $t_n$, the update formula to compute $q^{n+1}$
is given by

\begin{align}
  q^{n+1} = q^n + \delta t k_1 \quad \text{where} \, k_1 \, \text{solves}
  \quad \left( I-\frac{\delta t}{2} M^{n+1/2} \right) k_1 = M^{n+1/2}  q^n
\end{align}

where $M^{n+1/2} := M(t_n + \frac{\delta t}{2})$. In each time-step,
a linear equation is solved using GMRES to compute the stage variable $k_1$, which is then used it
to update $q^{n+1}$.

In addition to the IMR, two higher-order time-stepping schemes are available in Quandary, particularly a 4-th order and a 8-th order scheme which both are compsitional versions of the IMR. Those methods perform multiple composed IMR steps in each time-step interval to achieve higher order accuracy. Particularly, the 4-th order scheme (`IMR4`) performs 3 sub-steps per time interval, and the 8-th order (`IMR8`) performs 15 sub-steps per time time interval. Compared to the standard IMR, the higher-order methods allow for much larger time-steps to be taken to reach a certain accuracy tolerance, however, more work is done per time-step, creating a tradeoff at which the compositional methods can outperform the standard IMR scheme.

### Choice of the time-step size
<!-- The python interface to Quandary automatically computes a time-step size based on the fastest period of the system Hamiltonian. For the C++ code, it needs to be set by the user.  -->

In order to estimate a time-step size $\delta t$ that yields stable and accurate timestepping, Quandary's Python interface performs an eigenvalue decomposition of the system Hamiltonian $H_d$ to determine the fastest period of intrinsic dynamics. The user can adjust the desired number of timesteps ($P_{min}$) that are used to discretize per fastest period, and the time-step size is then computed from 

\begin{align} \label{eq:timestepsize}
  \delta t = \frac{\tau_{min}}{P_{min}} = \frac{2\pi}{P_{min} \lambda_{max}}
\end{align}

where $\lambda_{max}$ is the largest eigenvalue of the Hamiltonian. For the 2nd order IMR scheme, we recommend at least $P_{min}=80$. 

In order to test time-stepping accuracy, a standard $\Delta t$-test is recommended performed to ensure that the chosen time-step is small enough. For example, one can compute the *Richardson error estimator* of the current approximation error in $J^{\Delta t}$ compared to the true quantity $J^*$ with 

\begin{align}
  J^* - J^{\Delta t} = \frac{J^{\Delta t} - J^{\Delta t m}}{1-m^p} + O(\Delta t^{p+1})
\end{align}

where $p$ is the order of the time-stepping scheme (i.e. $p=2$ for the IMR and $p=8$ for the compositional IMR8), and $J^{\Delta t}, J^{\Delta tm}$ are approximations of the target quantity when using time-step sizes $\Delta t$ and $\Delta t m$, for a factor $m$.

## Sparse-matrix vs. matrix-free solver
Solving the differential equation with a time-stepping scheme requires efficient application of the right-hand-side (RHS) system matrix $M(t)$ to a (vectorized) state $q$. In Quandary, two versions to evaluate the matrix product $M(t)q(t)$ are available:

1. The *sparse-matrix solver* uses PETSc's sparse matrix format (sparse AIJ) to set up (and store) all building blocks inside $A(t)$ and $B(t)$, compare the appendix. Sparse matrix-vector products are then applied at each time-step to evaluate the products $A(t)u(t) - B(t) v(t)$ and $B(t)u(t) + A(t)v(t)$. For developers, the appendix provides details on each term within $A(t)$ and $B(t)$ which can be matched to the implementation in the code (class `MasterEq`).

2. The *matrix-free solver* considers the quantum state to be a tensor of rank $Q$ (Schroedinger) or $2Q$ (Lindblad). Instead of storing the building block matrices inside $M(t)$, the matrix-free solver applies tensor contractions to realize the action of $A(t)$ and $B(t)$ on the state vector. The matrix-free solver is much faster than the sparse-matrix solver (about 10x), no surprise. However the matrix-free solver is currently only implemented for systems of **2, 3, 4, or 5** oscillators.

<!-- **The matrix-free solver currently does not parallelize across the system dimension $N$**, hence the state vector is **not** distributed (i.e. no parallel Petsc!). The reason why we did not implement that yet is that $Q$ can often be large while each axis can be very short (e.g. modelling $Q=12$ qubits with $n_k=2$ energy levels per qubit), which yields a very high-dimensional tensor with very short axes. In that case, the standard (?) approach of parallelizing the tensor along its axes will likely lead to very poor scalability due to high communication overhead. We have not found a satisfying solution yet - if you have ideas, please reach out, we are happy to collaborate! -->



## Gradient computation via discrete adjoint back-propagation
Quandary computes the gradients of the objective function with respect to the design variables $\boldsymbol{\alpha}$ using the discrete adjoint method. The discrete adjoint approach yields exact and consistent gradients on the algorithmic level, at costs that are independent of the number of design variables.
To that end, the adjoint approach propagates local sensitivities backwards through the time-domain while concatenating contributions to the gradient using the chain-rule.

For the IMR timestepper, the discrete adjoint time-integration step for
adjoint variables denoted by $\bar q^{n}$ is given by

\begin{align}
  \bar q^{n} = \bar q^{n+1} + \delta t \left(M^{n+1/2}\right)^T \bar k_1
  \quad \text{where} \, \bar k_1 \, \text{solves} \quad \left(
  I-\frac{\delta t}{2} M^{n+1/2}\right)^T  \bar k_1 = \bar q^{n+1}
\end{align}

The contribution to the gradient $\nabla J$ for each time-step is

\begin{align}\label{eq:gradient}
  \nabla J += \delta t \left( \frac{\partial M^{n+1/2}}{\partial z}
  \left(q^n + \frac{\delta t}{2} k_1\right) \right)^T\bar k_1
\end{align}

Each evaluation of the gradient $\nabla J$ involves a forward solve of $n_{init}$ initial quantum states to evaluate the objective function at final time $T$, as well as $n_{init}$ backward solves to compute the adjoint states and the contributions to the gradient. Note that the gradient computation $\eqref{eq:gradient}$ requires the states and adjoint states at each time-step. For the Schroedinger solver, the primal states are recomputed by integrating Schroedinger's equation backwards in time, alongside the adjoint computation. For the Lindblad solver, the states $q^n$ are stored during forward propagation, and taken from storage during adjoint back-propagation (since we can't recompute it in case of Lindblad solver, due to dissipation).

For developers, a Central Finite Difference (CFD) test can be enabled by setting the compiler directive `TEST_FD_GRAD = 1` at the beginning of the `src/main.cpp` file. Quandary will then iterate over all elements in $\alpha$ and report the *relative* error of the implemented gradient with respect to the "true" gradient computed from CFD:


\begin{align*}
    \left(\nabla J(\boldsymbol{\alpha}) \right)_i \approx \frac{J(\bfa + \epsilon\bs{e}_i) - J(\bfa - \epsilon\bs{e}_i)}{2\epsilon} 
\end{align*}

for unit vectors $\bs{e}_i\in \R^d$, and $d$ being the dimension of $\bfa$.


## Storage of the control parameters
The control parameters $\bs{\alpha}$ are stored in the Quandary code in the following order: List oscillators first $(\vec{\alpha}^0, \dots, \vec{\alpha}^{Q-1})$, then for each $\vec{\alpha}^k \in
\R^{2 N_s^k N_f^k}$, iterate over all carrierwaves $\vec{\alpha}^k =
(\alpha^k_1,\dots, \alpha^k_{N_f})$ with $\alpha^k_f \in \R^{2 N_s^k}$, then each
$\alpha^k_f$ iterates over spline basis functions listing first all real then all imaginary
components: $\alpha^k_f = \alpha^{k(1)}_{1,f}, \dots, \alpha^{k(1)}_{N_s^k,f}, \alpha^{k(2)}_{1,f}, \dots, \alpha^{k(2)}_{N_s^k,f}$. Hence there are a total of $2\sum_k N_s^k N_f^k$ real-valued optimization parameters, which are stored in the following order:

\begin{align}
  \boldsymbol{\alpha} &:= \left( \vec{\alpha}^0, \dots, \vec{\alpha}^{Q-1} \right), \in
  \R^{2\sum_k N_s^k N_f^k} \quad \text{where}\\
  \vec{\alpha}^k = &\left( \alpha_{1,1}^{k(1)},\dots, \alpha_{N_s^k,1}^{k(1)}, \dots, \alpha_{1,N_f^k}^{k(1)}, \dots, \alpha_{N_s^k,N_f^k}^{k(1)} \right.\\
                 &  \left. \alpha_{1,1}^{k(2)},\dots, \alpha_{N_s^k,1}^{k(2)}, \dots, \alpha_{1,N_f^k}^{k(2)}, \dots, \alpha_{N_s^k,N_f^k}^{k(2)} \right)
\end{align}

iterating over $Q$ subsystems first, then $N_f^k$ carrier wave frequencies, then $N_s^k$ splines, listing first all real parts then all imaginary parts. To access an element $\alpha_{s,f}^{k(i)}$, $i=0,1$, from storage $\bfa$:

\begin{align}
  \alpha_{s,f}^{k(i)} = \bfa \left[ \left(\sum_{j=0}^{k-1} 2 N_s^j N_f^j\right) + f*2 N_s^k + s + i*N_s^k N_f^k \right],
\end{align}

# Parallelization
Quandary offers two levels of parallelization using MPI.

1. Parallelization over initial conditions: The $n_{init}$ initial conditions $\rho_i(0)$ can be distributed over `np_init` compute units. Since initial condition are propagated through the time-domain for solving Lindblad's or Schroedinger's equation independently from each other, speedup from distributed initial conditions is ideal.
2. Parallel linear algebra with Petsc (sparse-matrix solver only): For the sparse-matrix solver, Quandary utilizes Petsc's parallel sparse matrix and vector storage to distribute the state vector onto `np_petsc` compute units (spatial parallelization). To perform scaling results, make sure to disable code output (or reduce the output frequency to print only the last time-step), because writing the data files invokes additional MPI calls to gather data on the master node. Strong and weak scaling studies for parallel linear algebra are presented in [@guenther2021quantum].

Since those two levels of parallelism are orthogonal, Quandary splits the global communicator (MPI\_COMM\_WORLD) into
two sub-communicator such that the total number of executing MPI
processes ($np_{total}$) is split as

\begin{align*}
  np_{init} * np_{petsc} = np_{total}.
\end{align*}

Since parallelization over different initial conditions is perfect, Quandary automatically sets $np_{init} = n_{init}$, i.e. the total number of cores for distributing initial conditions is the total number of initial conditions that are considered in this run, as specified by the configuration option `intialcondition`. The number of cores for distributed linear algebra with Petsc is then computed from above.

It is currently required that the number of total cores for executing quandary is an integer divisor of multiplier of the number of initial conditions, such that each processor group owns the same number of initial conditions.

It is further required that the system dimension is an integer multiple of the number of cores used for distributed linear algebra from Petsc, i.e. it is required that $\frac{M}{np_{petsc}} \in \mathbb{N}$ where $M=N^2$ in the Lindblad solver case and $M=N$ in the Schroedinger case. This requirement is a little
  annoying, however the current implementation requires this due to the
  storage of the real and imaginary parts of the vectorized
  state.

# Output and plotting the results
Quandary generates various output files for system evolution of the current (optimized) controls as well as the optimization progress. All data files will be dumped into a user-specified folder through the config option `datadir`.

### Output options with regard to state evolution
For each subsystem $k$, the user can specify the desired state evolution output through the config option `output<k>`:

- `expectedEnergy`: This option prints the time evolution of the expected energy level of subsystem $k$ into files with naming convention `expected<k>.iinit<i>.dat`, where $i=1,\dots,n_{init}$ denotes the unique identifier for each initial condition $\rho_i(0)$ that was propagated through (see Section [Initial conditions](#sec:initcond)). This file contains two columns, the first row being the time values, the second one being the expectation value of the energy level of subsystem $k$ at that time point, computed from

    \begin{align}
      \langle N^{(n_k)}\rangle = \mbox{Tr}\left(N^{(n_k)} \rho^k(t)\right)
    \end{align}

    where $N^{(n_k)} = \left(a^{(n_k)}\right)^\dagger \left(a^{(n_k)}\right)$ denotes the number operator in subsystem $k$ and $\rho^k$ denotes the reduced density matrix or state for subsystem $k$. 
- `expectedEnergyComposite` Prints the time evolution of the expected energy level of the entire (full-dimensional) system state into files (one for each initial condition, as above): $mbox{Tr}\left(N \rho(t)\right)$ for the number operator $N$ in the full dimensions.
- `population`: This option prints the time evolution of the state's occupation in each energy level into files named `population<k>.iinit<i>.dat`, for each initial condition $i=1,\dots, n_{init}$ and each subsystem $k$. The files contain $n_k+1$ columns, the first one being the time values, the remaining ones correspond to the population of each level $l=0,\dots,n_k-1$ of the reduced density matrix or state vector at that time point. For Lindblad's solver, these are the diagonal elements of the reduced density matrix ($\rho_{ll}^k(t), l=0,\dots n_k-1$), for Schroedinger's solver it's the absolute values of the reduced state vector elements $|\psi^k_l(t)|^2, l=0,\dots n_k-1$. 
- `populationComposite`: Prints the time evolution of the state populations of the entire (full-dimensional) system into files (one for each initial condition, as above).
- `fullstate`: For smaller systems, one can choose to print out the full state $\rho(t)$ or $\psi(t)$ for each time point into the files `rho_Re.iinit<m>.dat` and `rho_Im.iinit<m>.dat`, for the real and imaginary parts. These files contain $N^2+1$ (Lindblad) or $N+1$ (Schroedinger) columns the first one being the time point value and the remaining ones contain the vectorized density matrix or the state vector for that time point. Note that these file become very big very quickly -- use with care!

The user can change the frequency of output in time (printing only every $j$-th time point) through the option `output_frequency`. This is particularly important when doing performance tests, as computing the reduced states for output requires extra computation and communication that might skew performance tests.

### Output with regard to simulation and optimization
- `config_log.toml` contains all configuration options that had been used for the current run of the C++ code.
- `params.dat` contains the control parameters $\bfa$ that had been used to determine the current control pulses. This file contains one column containing all parameters, ordered as stored, see Section [Control pulses](#sec:controlpulses).
- `control<k>.dat` contain the resulting control pulses applied to subsystem $k$ over time. It contains four columns, the first one being the time, second and third being $p^k(t)$ and $q^k(t)$ (rotating frame controls), and the last one is the corresponding lab-frame pulse $f^k(t)$. Note that the units of the control pulses are in frequency domain (divided by $2\pi)$. The unit matches the unit specified with the system parameters such as the qubit ground frequencies $\omega_k$.
- `optim_history.dat` contains information about the optimization progress in terms of the overall objective function and contribution from each term (cost at final time $T$ and contribution from the tikhonov regularization and the penalty term), as well the norm of the gradient and the fidelity, for each iteration of the optimization. If only a forward simulation is performed, this file still prints out the objective function and fidelity for the forward simulation.
Quandary always prints the current parameters and control pulses at the beginning of a simulation or optimization, and in addition at every $l$-th optimization iteration determined from the `optim_monitor_frequency` configuration option.

### Plotting
The format of all output files are very well suited for plotting with [Gnuplot](http://www.gnuplot.info), which is a command-line based plotting program that can output directly to screen, or into many other formats such as png, eps, or even tex. As an example, from within a Gnuplot session, you can plot e.g. the expected energy level of subsystem $k=0$ for initial condition $m=0$ by the simple command

``` console
gnuplot> plot 'expected0.iinit0000.dat' using 1:2 with lines title 'expected energy subsystem 0'
```

which plots the first against the second column of the file 'expected0.iinit0000.dat' to screen, connecting each point with a line. Additional lines (and files) can be added to the same plot by extending the above command with another file separated by comma. There are many example scripts for plotting with gnuplot online, and as a starting point I recommend looking into some scripts in the 'quandary/util/' folder.


# Testing
Quandary has a set of regression tests. Please take a look at the `tests/regression/README.md` document for instructions on how to run the regression tests.

# Acknowledgments
This work was performed under the auspices of the U.S. Department of Energy by Lawrence
Livermore National Laboratory under Contract DE-AC52-07NA27344. LLNL-SM-818073.

<!--
This document was prepared as an account of work sponsored by an agency of the United States
government. Neither the United States government nor Lawrence Livermore National Security, LLC,
nor any of their employees makes any warranty, expressed or implied, or assumes any legal
liability or responsibility for the accuracy, completeness, or usefulness of any information,
apparatus, product, or process disclosed, or represents that its use would not infringe
privately owned rights. Reference herein to any specific commercial product, process, or service
by trade name, trademark, manufacturer, or otherwise does not necessarily constitute or imply
its endorsement, recommendation, or favoring by the United States government or Lawrence
Livermore National Security, LLC. The views and opinions of authors expressed herein do not
necessarily state or reflect those of the United States government or Lawrence Livermore
National Security, LLC, and shall not be used for advertising or product endorsement purposes.
-->

# Appendix: Details for the real-valued system matrices (standard model)
The RHS system matrices $A(t) = Re(-iH(t))=Im(H(t))$ and $B(t) = Im(-iH(t)) = -Re(H(t))$ for the standard Hamiltonian model are assembled as follows. 

\begin{align}
A(t) &= A_d + \sum_k  q^k(\vec{\alpha}^k,t) A_c^k + \sum_{l>k} J_{kl} \sin(\eta_{kl}t)  A_d^{kl} \\
B(t) &=  B_d + \sum_k p^k(\vec{\alpha}^k,t) B_c^k + \sum_{kl} J_{kl} \cos(\eta_{kl}t)B_d^{kl}\\
\end{align}

**Schroedinger solver**:

\begin{align}
  A_d &:= 0\\
  A_c^k &:=  a_k - a_k^\dagger  \\
  A_d^{kl} &:=  a_k^\dagger a_l + a_k a_l^\dagger \\
\notag\\
  B_d &:= \sum_k -(\omega_k - \omega_k^{\text{rot}}) a_k^\dagger a_k  + \frac{\xi_k}{2}\left( a_k^\dagger a_k^\dagger a_k a_k \right)  + \sum_{l>k}  \xi_{kl}\left(a_k^\dagger a_k a_l^\dagger a_l  \right)\\
    B_c^k &:=  -(a_k + a_k^\dagger) \\
    B_d^{kl} &:=  - \left(a_k^\dagger a_l + a_k a_l^\dagger\right)  
\end{align}

**Lindblad solver**

\begin{align}
  A_d &:= \vec L\\
  A_c^k &:=  I_N \otimes \left(a_k - a_k^\dagger\right) - \left(a_k - a_k^\dagger\right)^T\otimes I_N \\
  A_d^{kl} &:=  I_N\otimes \left(a_k^\dagger a_l - a_k a_l^\dagger\right) - \left(a_k^\dagger a_l - a_k a_l^\dagger\right)^T\otimes I_N \\
\notag\\
  B_d &:= \sum_k -(\omega_k - \omega_k^{\text{rot}}) \left(I_N \otimes a_k^\dagger a_k - (a_k^\dagger a_k)^T \otimes I_N \right) + \frac{\xi_k}{2}\left( I_N \otimes a_k^\dagger a_k^\dagger a_k a_k - (a_k^\dagger a_k^\dagger a_k a_k )^T\otimes I_N\right)  \\
    &\quad + \sum_{l>k}  \xi_{kl} \left(I_N \otimes a_k^\dagger a_k a_l^\dagger a_l - (a_k^\dagger a_k a_l^\dagger a_l)^T \otimes I_N \right)\\
    B_c^k &:=  -\left( I_N \otimes \left(a_k + a_k^\dagger\right) + \left(a_k + a_k^\dagger\right)^T\otimes I_N \right)\\
    B_d^{kl} &:=  - \left(I_N\otimes \left(a_k^\dagger a_l + a_k a_l^\dagger\right) + \left(a_k^\dagger a_l + a_k a_l^\dagger\right)^T\otimes I_N \right)\\
\end{align}

The sparse-matrix solver initializes and stores the building blocks 
$A_d, A_d^{kl}, A_c^k, B_d, B_d^{kl}, B_c^k$ using Petsc's sparse-matrix MPIAIJ format. The matrix-free solver (2-5 subsystems) applies the actions of those building blocks to state vectors.

