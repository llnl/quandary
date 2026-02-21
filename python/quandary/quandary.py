import os, os.path, copy
import numpy as np
from subprocess import run, PIPE, Popen, call
import matplotlib.pyplot as plt
from dataclasses import dataclass, field, replace
from typing import List, Dict
## For some Matplotlib installations to work, you might need the below...
# import PyQt6.QtCore

@dataclass
class Quandary:
    """
    This class collects configuration options to run quandary and sets all defaults. Each parameter can be overwritten within the constructor. The number of time-steps required to resolve the time-domain, as well as the resonant carrier wave frequencies are computed within the constructor. If you attempt to change options of the configuration *after* construction by accessing them directly (e.g. myconfig.fieldname = mynewsetting), it is therefore advised to call myconfig.update() afterwards to recompute the number of time-steps and carrier waves.

    Environment Variables:
    --------------------
    QUANDARY_BASE_DATADIR  : Base directory for output files. If set, all output will be written to this directory
                             instead of the current working directory. Relative paths provided to methods will be
                             considered relative to this base directory.

    Parameters
    ----------
    # Quantum system specifications
    Ne           # Number of essential energy levels per qubit. Default: [3]        
    Ng           # Number of extra guard levels per qubit. Default: [0] (no guard levels)
    freq01       # 01-transition frequencies [GHz] per qubit. Default: [4.10595]
    selfkerr     # Anharmonicities [GHz] per qubit. Default: [0.2198]
    rotfreq      # Frequency of rotations for computational frame [GHz] per qubit. Default: freq01
    Jkl          # Dipole-dipole coupling strength [GHz]. Formated list [J01, J02, ..., J12, J13, ...] Default: 0
    crosskerr    # ZZ coupling strength [GHz]. Formated list [g01, g02, ..., g12, g13, ...] Default: 0
    T1           # Optional: T1-Decay time [ns] per qubit (invokes Lindblad solver). Default: 0
    T2           # Optional: T2-Dephasing time [ns] per qubit (invokes Lindblad solver). Default: 0

    # Optional: User-defined system and control Hamiltonian operators. Default: Superconducting Hamiltonian model
    Hsys                # Optional: User specified system Hamiltonian model. Array. 
    Hc_re               # Optional: User specified control Hamiltonian operators for each qubit (real-parts). List of Arrays
    Hc_im               # Optional: User specified control Hamiltonian operators for each qubit (real-parts) List of Arrays
    standardmodel       # Internal: Bool to use standard Hamiltonian model for superconduction qubits. Default: True

    # Time duration and discretization options
    T            # Pulse duration (simulation time). Default: 100ns
    Pmin         # Number of discretization points to resolve the shortest period of the dynamics (determines <nsteps>). Default: 150
    nsteps       # Number of time-discretization points (will be computed internally based on Pmin, or can be set here)
    timestepper  # Time-discretization scheme. Default: "IMR"

    # Optimization targets and initial states options
    targetgate          # Complex target unitary in the essential level dimensions for gate optimization. Default: none
    targetstate         # Complex target state vector for state-to-state optimization. Default: none
    initialcondition    # Choose from provided initial states at time t=0.0: "basis" (all basis states, default), "pure, 0,0,1,..." (one pure initial state |001...>), or pass a vector as initial state. Default: "basis" 
    gate_rot_freq       # Specify frequencies to rotate a target gate (one per oscillator). Default: no rotation (0.0 for each oscillator)

    # Control pulse options
    pcof0               # Optional: Pass an initial vector of control parameters. Default: none
    pcof0_filename      # Optional: Load initial control parameter vector from a file. Default: none
    randomize_init_ctrl # Randomize the initial control parameters (will be ignored if pcof0 or pcof0_filename are given). Default: True
    initctrl_MHz        # Amplitude [MHz] of initial control parameters. Float or List[float]. Default: 10 MHz.
    maxctrl_MHz         # Amplitude bounds for the control pulses [MHz]. Float or List[float]. Default: none
    control_enforce_BC  # Bool to let control pulses start and end at zero. Default: False
    spline_knot_spacing # Spacing of Bspline basis functions [ns]. The smaller this is, the larger the number of splines. Default: 3ns
    nsplines            # Number of Bspline basis functions. Default: T/spline_knot_spacing + 2
    spline_order        # Order of the B-spline basis (0 or 2). Default: 2
    carrier_frequency   # Carrier frequencies for each oscillator. List[List[float]]. Default will be computed based on Hsys.
    cw_amp_thres        # Threshold to ignore carrier wave frequencies whose growth rate is below this value. Default: 1e-7
    cw_prox_thres       # Threshold to distinguish different carrier wave frequencies from each other. Default: 1e-2

    # Optimization options
    maxiter             # Maximum number of optimization iterations. Default 200
    tol_infidelity      # Optimization stopping criterion based on the infidelity. Default 1e-5
    tol_costfunc        # Optimization stopping criterion based on the objective function value. Default 1e-4
    tol_gnorm_abs       # Optimization stopping criterion based on absolute gradient norm drop. Default 1e-4
    tol_gnorm_rel       # Optimization stopping criterion based on relative gradient norm drop. Default 1e-4
    costfunction        # Cost function measure: "Jtrace" or "Jfrobenius". Default: "Jtrace"
    optim_target        # Optional: Set other optimization target string, if not specified through the targetgate or targetstate. 
    gamma_tik0          # Parameter for Tikhonov regularization ||alpha||^2. Default 1e-4
    gamma_tik0_interpolate # Switch to use ||alpha-alpha_0||^2 instead, where alpha_0 is the initial guess. Default: False
    gamma_leakage       # Parameter for leakage prevention. Default: 0.1
    gamma_energy        # Parameter for integral penality term on the control pulse energy. Default: 0.1
    gamma_dpdm          # Parameter for integral penality term on second state derivative. Default: 0.01
    gamma_variation     # Parameter for penality term on variations in the control parameters: Default: 0.01

    # General options
    rand_seed            # Set a fixed random number generator seed. Default: None (non-reproducable)
    print_frequency_iter # Output frequency for optimization iterations. (Print every <x> iterations). Default: 1
    usematfree           # Switch to use matrix-free (rather than sparse-matrix) solver. Default: True
    verbose              # Switch to turn on more screen output for debugging. Default: False

    Internal variables. 
    -------------------
    _ninit                : int  # number of initial conditions that are propagated
    _lindblad_solver      : bool # Flag to determine whether lindblad solver vs schroedinger solver 
    _hamiltonian_filename_Hsys : str 
    _hamiltonian_filename_Hc : str 
    _gatefilename         : str 
    _initstatefilename    : str 
    _initialstate         : List[complex] = field(default_factory=list)


    Output parameters, available after Quandary has been executed (simulate or optimze)
    -----------------------------------------------------
    popt         # Optimized control palamters (Bspline coefficients). List[float]
    time         # Vector of discretized time points. List[float]
    optim_hist   # Optimization history: all fields as in Quandary's output file optim_history.dat. Dictionary.
    uT           # Evolved states at final time T. This is the (unitary) solution operator, if the initial conditions span the full basis. 
    """ 

    # Quantum system specifications
    Ne        : List[int]   = field(default_factory=lambda: [3])
    Ng        : List[int]   = field(default_factory=lambda: [0])
    freq01    : List[float] = field(default_factory=lambda: [4.10595])
    selfkerr  : List[float] = field(default_factory=lambda: [0.2198])
    rotfreq   : List[float] = field(default_factory=list)
    Jkl       : List[float] = field(default_factory=list)
    crosskerr : List[float] = field(default_factory=list)
    T1        : List[float] = field(default_factory=list)
    T2        : List[float] = field(default_factory=list)
    # Optiona: User-defined Hamiltonian model (Default = superconducting qubits model)
    Hsys                : List[complex]     = field(default_factory=list)
    Hc_re               : List[List[float]] = field(default_factory=list)
    Hc_im               : List[List[float]] = field(default_factory=list)
    standardmodel       : bool              = True
    # Time duration and discretization options
    T            : float = 100.0
    Pmin         : int   = 150
    nsteps       : int   = -1
    dT           : float = -1.0
    timestepper  : str   = "IMR"
    # Optimization targets and initial states options
    targetgate             : List[List[complex]] = field(default_factory=list) 
    targetstate            : List[complex] = field(default_factory=list) 
    initialcondition       : str = "basis"
    gate_rot_freq          : List[float] = field(default_factory=list)
    # Control pulse options
    pcof0               : List[float] = field(default_factory=list)   
    pcof0_filename      : str         = ""                            
    randomize_init_ctrl : bool        = True                          
    initctrl_MHz        : List[float] = field(default_factory=list)   
    maxctrl_MHz         : List[float] = field(default_factory=list)   
    control_enforce_BC  : bool        = False                         
    spline_knot_spacing : float       = 3.0
    nsplines            : int         = -1 
    spline_order        : int         = 2                           
    carrier_frequency   : List[List[float]] = field(default_factory=list) 
    cw_amp_thres        : float       = 1e-7
    cw_prox_thres       : float       = 1e-2                        
    # Optimization options
    maxiter                : int   = 200         
    tol_infidelity         : float = 1e-5        
    tol_costfunc           : float = 1e-4        
    tol_gnorm_abs          : float = 1e-4
    tol_gnorm_rel          : float = 1e-4
    costfunction           : str   = "Jtrace"                      
    optim_target           : str   = "gate, none"
    gamma_tik0             : float = 1e-4 
    gamma_tik0_interpolate : float = 0.0
    gamma_leakage          : float = 0.1 	       
    gamma_energy           : float = 0.1
    gamma_dpdm             : float = 0.01        
    gamma_variation        : float = 0.01        
    # General options
    rand_seed              : int  = None
    print_frequency_iter   : int  = 1
    usematfree             : bool = True 
    verbose                : bool = False
    # Internal configuration. Should not be changed by user.
    _ninit                : int         = -1
    _lindblad_solver      : bool        = False
    _hamiltonian_filename_Hsys : str         = ""
    _hamiltonian_filename_Hc : str         = ""
    _gatefilename         : str         = ""
    _initstatefilename    : str         = ""
    _initialstate         : List[complex] = field(default_factory=list)
    
    # Output parameters available after Quandary has been run
    popt        : List[float]   = field(default_factory=list)
    time        : List[float]   = field(default_factory=list)
    optim_hist  : Dict          = field(default_factory=dict)
    uT          : List[float]   = field(default_factory=list)


    def __post_init__(self):
        """
        This function sets all default options and overwrites those that specified by the user. It performs sanity checks on the user-specified input, and it computes
          - <nsteps>            : the number of time steps based on Hamiltonian eigenvalues and Pmin
          - <carrier_frequency> : carrier wave frequencies bases on system resonances
        """
        
        if self.spline_order == 0:
            minspline = 2
        elif self.spline_order == 2:
            minspline = 5 if self.control_enforce_BC else 3
        else:
            print("Error: spline order = ", self.spline_order, " is currently not available. Choose 0 or 2.")
            return -1

        # Set some defaults, if not set by the user
        if len(self.freq01) != len(self.Ne) and len(self.Hsys)<=0:
            self.Ne = [2 for _ in range(len(self.freq01))]
        if len(self.Ng) != len(self.Ne):
            self.Ng = [0 for _ in range(len(self.Ne))]
        if len(self.selfkerr) != len(self.Ne):
            self.selfkerr= np.zeros(len(self.Ne))
        if len(self.rotfreq) == 0:
            self.rotfreq = self.freq01
        if len(self.gate_rot_freq) == 0:
            self.gate_rot_freq = np.zeros(len(self.rotfreq))
        
        if isinstance(self.initctrl_MHz, float) or isinstance(self.initctrl_MHz, int):
            max_alloscillators = self.initctrl_MHz
            self.initctrl_MHz = [max_alloscillators for _ in range(len(self.Ne))]
        if len(self.initctrl_MHz) == 0:
            self.initctrl_MHz = [10.0 for _ in range(len(self.Ne))]
        if len(self.Hsys) > 0 and not self.standardmodel: # User-provided Hamiltonian operators 
            self.standardmodel=False   
            self.usematfree=False
        else: # Using standard Hamiltonian model. Set it up only if needed for computing dT or the carrier wave frequencies later
            self.standardmodel=True
        if len(self.targetstate) > 0:
            self.optim_target = "file"
        if len(self.targetgate) > 0:
            self.optim_target = "gate, file"
        if not isinstance(self.initialcondition, str):
            self._initialstate=self.initialcondition.copy()
            self.initialcondition = "file" 
        # Convert maxctrl_MHz to a list for each oscillator, if not so already
        if isinstance(self.maxctrl_MHz, float) or isinstance(self.maxctrl_MHz, int):
            max_alloscillators = self.maxctrl_MHz
            self.maxctrl_MHz = [max_alloscillators for _ in range(len(self.Ne))]
        
        # Store the number of initial conditions and solver flag
        self._lindblad_solver = True if (len(self.T1)>0) or (len(self.T2)>0) else False
        if self.initialcondition[0:4] == "file" or self.initialcondition[0:4] == "pure":
            self._ninit = 1
        else:
            self._ninit = np.prod(self.Ne)
        if self._lindblad_solver:
            self._ninit = self._ninit**2
        
        # Estimate the number of required time steps
        if self.dT < 0:
            if self.standardmodel==True: # set up the standard Hamiltonian first
                Ntot = [sum(x) for x in zip(self.Ne, self.Ng)]
                self.Hsys, self.Hc_re, self.Hc_im = hamiltonians(N=Ntot, freq01=self.freq01, selfkerr=self.selfkerr, crosskerr=self.crosskerr, Jkl=self.Jkl, rotfreq=self.rotfreq, verbose=self.verbose)
            self.nsteps = estimate_timesteps(T=self.T, Hsys=self.Hsys, Hc_re=self.Hc_re, Hc_im=self.Hc_im, maxctrl_MHz=self.maxctrl_MHz, Pmin=self.Pmin)
            self.dT = self.T/self.nsteps
        else:
            self.nsteps = int(np.ceil(self.T / self.dT))
            self.T = self.nsteps*self.dT
        if self.verbose:
            print("Final time: ",self.T,"ns, Number of timesteps: ", self.nsteps,", dt=", self.T/self.nsteps, "ns")
            print("Maximum control amplitudes: ", self.maxctrl_MHz, "MHz")

        # Get number of splines right
        if self.nsplines < 0:
            if self.spline_order == 0:
                self.nsplines = int(np.max([np.rint(self.nsteps*self.dT/self.spline_knot_spacing+1), minspline])) 
            else: 
                self.nsplines = int(np.max([np.ceil(self.T/self.spline_knot_spacing+ 2), minspline]))

            self.spline_knot_spacing = self.nsteps*self.dT / (self.nsplines-1) if self.spline_order == 0 else self.nsteps*self.dT / (self.nsplines-2)
        else:
            self.spline_knot_spacing= self.nsteps*self.dT/(self.nsplines-1) if self.spline_order == 0 else self.T/(self.nsplines - 2)

        # Estimate carrier wave frequencies
        if self.spline_order == 0 and len(self.carrier_frequency) == 0:
            self.carrier_frequency = [[0.0] for _ in range(len(self.freq01))]
        if len(self.carrier_frequency) == 0: 
            # set up the standard Hamiltonian first, if needed and if not done so already
            if self.standardmodel==True and len(self.Hsys)<=0:
                Ntot = [sum(x) for x in zip(self.Ne, self.Ng)]
                self.Hsys, self.Hc_re, self.Hc_im = hamiltonians(N=Ntot, freq01=self.freq01, selfkerr=self.selfkerr, crosskerr=self.crosskerr, Jkl=self.Jkl, rotfreq=self.rotfreq, verbose=self.verbose)
            self.carrier_frequency, _ = get_resonances(Ne=self.Ne, Ng=self.Ng, Hsys=self.Hsys, Hc_re=self.Hc_re, Hc_im=self.Hc_im, rotfreq=self.rotfreq, verbose=self.verbose, cw_amp_thres=self.cw_amp_thres, cw_prox_thres=self.cw_prox_thres, stdmodel=self.standardmodel)

        if self.verbose: 
            print("\n")
            print("Carrier frequencies (rot. frame): ", self.carrier_frequency)
            print("\n")

    def copy(self):
        """ Return a new instance of the class, copying all member variables. """
        return replace(self)

    def update(self):
        """
        Call this function if you have changed a configuration option outside of the Quandary constructor with "quandary.variablename = new_variable". This will ensure that the number of time steps and carrier waves are re-computed, given the new setting.
        """

        # Output parameters available after __run().
        popt_org = self.popt.copy()
        time_org = self.time.copy()
        opti_org = self.optim_hist.copy() 
        uT_org   = self.uT.copy()

        self.__post_init__()

        self.popt       = popt_org.copy()
        self.time       = time_org.copy()
        self.optim_hist = opti_org.copy()
        self.uT         = uT_org.copy()


    def simulate(self, *, pcof0=[], pt0=[], qt0=[], maxcores=-1, datadir="./run_dir", quandary_exec="", cygwinbash="", mpi_exec="mpirun -np ", batchargs=[]):
        """ 
        Simulate the quantm dynamics using the current settings. 

        Optional arguments:
        ------------------- 
        pcof0         : List of control parameters to be simulated. Default: Use use initial guess from the quandary object (pcof0, or pcof0_filename, or randomized initial guess)
        pt0           : List of ndarrays for the real part of the control function [MHz] for each oscillator, ndarray size = nsteps+1. Assumes spline_order == 0 and ignores the pcof0 argument. Default: []
        qt0           : Same as pt0, but for the imaginary part.
        maxcores      : Maximum number of processing cores. Default: number of initial conditions
        datadir       : Data directory for storing output files. Default: "./run_dir".
                        If $QUANDARY_BASE_DATADIR is set, this will be relative to that directory, otherwise relative to current working directory
        quandary_exec : Location of Quandary's C++ executable, if not in $PATH
        cygwinbash    : To run on Windows through Cygwin, set the path to Cygwin/bash.exe. Default: None.
        mpi_exec      : String for MPI launcher prefix, e.g. "mpirun -np" or "srun -n". The string should include the flag for core counts, but not the number of cores itself which will be appended automatically
        batchargs     : [str(time), str(accountname), int(nodes)] If given, submits a batch job rather than local execution. Specify the max. runtime (string), the account name (string) and the number of requested nodes (int). Note, the number of executing *cores* is defined through 'maxcores'. 

        Returns:
        --------
        time            :  List of time-points at which the controls are evaluated  (List)
        pt, qt          :  p,q-control pulses [MHz] at each time point for each oscillator (List of list)
        infidelity      :  Infidelity of the pulses for the specified target operation (Float)
        expectedEnergy  :  Evolution of the expected energy of each oscillator and each initial condition. Acces: expectedEnergy[oscillator][initialcondition]
        population      :  Evolution of the population of each oscillator, of each initial condition. (expectedEnergy[oscillator][initialcondition])
        """
        

        if len(pt0) > 0 and len(qt0) > 0:
            nsplines_org = self.nsplines
            spline_knot_spacing_org = self.spline_knot_spacing
            cw_org = self.carrier_frequency.copy()
            self.carrier_frequency = [[0.0] for _ in range(len(self.Ne))]
            spline_order_org = self.spline_order
            self.spline_order = 0
            self.spline_knot_spacing = self.dT 
            self.nsplines = np.max([2,int(np.ceil(self.nsteps*self.dT/self.spline_knot_spacing + 1))])

            pcof0 = self.downsample_pulses(pt0=pt0, qt0=qt0)

        result = self.__run(pcof0=pcof0, runtype="simulation", overwrite_popt=False, maxcores=maxcores, datadir=datadir, quandary_exec=quandary_exec, cygwinbash=cygwinbash,mpi_exec=mpi_exec, batchargs=batchargs)

        if len(pt0) > 0 and len(qt0) > 0:
            self.nsplines = nsplines_org
            self.spline_knot_spacing = spline_knot_spacing_org
            self.spline_order = spline_order_org
            self.carrier_frequency = cw_org.copy()
        
        return result


    def optimize(self, *, pcof0=[], pt0=[], qt0=[], maxcores=-1, datadir="./run_dir", quandary_exec="", cygwinbash="", mpi_exec="mpirun -np ", batchargs=[]):
        """ 
        Optimize the quantm dynamics using the current settings. 

        Optional arguments:
        ------------------- 
        pcof0          : List of control parameters to start the optimization from. Default: Use initial guess from the Quandary (pcof0, or pcof0_filename, or randomized initial guess)
        maxcores      : Maximum number of processing cores. Default: number of initial conditions
        pt, qt          :  p,q-control pulses [MHz] at each time point for each oscillator (List of list)
        datadir       : Data directory for storing output files. Default: "./run_dir".
                        If $QUANDARY_BASE_DATADIR is set, this will be relative to that directory, otherwise relative to current working directory
        quandary_exec : Location of Quandary's C++ executable, if not in $PATH
        mpi_exec      : String for MPI launcher prefix, e.g. "mpirun -np" or "srun -n". The string should include the flag for core counts, but not the number of cores itself which will be appended automatically
        cygwinbash    : To run on Windows through Cygwin, set the path to Cygwin/bash.exe. Default: None.
        batchargs     : [str(time), str(accountname), int(nodes)] If given, submits a batch job rather than local execution. Specify the max. runtime (string), the account name (string) and the number of requested nodes (int). Note, the number of executing *cores* is defined through 'maxcores'. 

        Returns:
        --------
        time            :  List of time-points at which the controls are evaluated  (List)
        pt, qt          :  p,q-control pulses at each time point for each oscillator (List of list)
        infidelity      :  Infidelity of the pulses for the specified target operation (Float)
        expectedEnergy  :  Evolution of the expected energy of each oscillator and each initial condition. Acces: expectedEnergy[oscillator][initialcondition]
        population      :  Evolution of the population of each oscillator, of each initial condition. (expectedEnergy[oscillator][initialcondition])
        """
        
        if len(pt0) > 0 and len(qt0) > 0:
            nsplines_org = self.nsplines
            spline_knot_spacing_org = self.spline_knot_spacing
            spline_order_org = self.spline_order
            cw_org = self.carrier_frequency.copy()
            self.carrier_frequency = [[0.0] for _ in range(len(self.Ne))]
            self.spline_order = 0
            self.spline_knot_spacing = self.dT 
            self.nsplines = np.max([2,int(np.ceil(self.nsteps*self.dT/self.spline_knot_spacing + 1))])
            pcof0 = self.downsample_pulses(pt0=pt0, qt0=qt0)

        result = self.__run(pcof0=pcof0, runtype="optimization", overwrite_popt=True, maxcores=maxcores, datadir=datadir, quandary_exec=quandary_exec, cygwinbash=cygwinbash, mpi_exec=mpi_exec, batchargs=batchargs)

        if len(pt0) > 0 and len(qt0) > 0:
            self.nsplines = nsplines_org
            self.spline_knot_spacing = spline_knot_spacing_org       
            self.spline_order = spline_order_org
            self.carrier_frequency = cw_org.copy()
        
        return result
    

    def evalControls(self, *, pcof0=[], points_per_ns=1,datadir="./run_dir", quandary_exec="", mpi_exec="mpirun -np ", cygwinbash=""):
        """
        Evaluate control pulses on a specific sample rate.       
        
        Optional arguments:
        --------------------
        pcof0         :  List of control parameters (bspline coefficients) that determine the controls pulse. If not given, the initial guess from Quandary class will be used (pcof0, or filename, or random initial control...)
        points_per_ns :  sample rate of the resulting controls. Default: 1ns 
        datadir       :  Directory for output files. Default: "./run_dir".
                         If $QUANDARY_BASE_DATADIR is set, this will be relative to that directory, otherwise relative to current working directory
        quandary_exec :  Path to Quandary's C++ executable if not in $PATH
        mpi_exec      : String for MPI launcher prefix, e.g. "mpirun -np" or "srun -n". The string should include the flag for core counts, but not the number of cores itself which will be appended automatically
        cygwinbash    : To run on Windows through Cygwin, set the path to Cygwin/bash.exe. Default: None.
    
        Returns:
        ---------
        time    :  List of time-points
        pt, qt  :  Control pulses for each oscillator at each time point
        """

        # Copy original setting and overwrite number of time steps for simulation
        nsteps_org = self.nsteps
        dT_org = self.dT
        self.nsteps = int(np.floor(self.T * points_per_ns))
        self.dT = self.T/self.nsteps
    
        datadir = resolve_datadir(datadir)

        # Execute quandary in 'evalcontrols' mode
        datadir_controls = datadir +"_ppns"+str(points_per_ns)
        os.makedirs(datadir_controls, exist_ok=True)
        runtype = 'evalcontrols'
        configfile_eval= self.__dump(pcof0=pcof0, runtype=runtype, datadir=datadir_controls)
        err = execute(runtype=runtype, ncores=1, config_filename=configfile_eval, datadir=datadir_controls, quandary_exec=quandary_exec, verbose=False, mpi_exec=mpi_exec, cygwinbash=cygwinbash)
        time, pt, qt, _, _, _, pcof, _, _ = self.get_results(datadir=datadir_controls, ignore_failure=True)

        # Save pcof to config.popt
        self.popt = pcof[:]
    
        # Restore original setting
        self.nsteps = nsteps_org
        self.dT = dT_org 

        return time, pt, qt


    def downsample_pulses(self, *, pt0=[], qt0=[]):
        if self.spline_order == 0: #specifying (pt, qt) only makes sense for piecewise constant B-splines
            Nsys = len(self.Ne)
            if len(pt0) == Nsys and len(qt0) == Nsys:
                sizes_ok = True
                for iosc in range(Nsys):
                    if sizes_ok and len(pt0[iosc]) >= 2 and len(pt0[iosc]) == len(qt0[iosc]):
                        sizes_ok = True
                    else:
                        sizes_ok = False
                # print("simulate(): sizes_ok = ", sizes_ok)
                if sizes_ok:
                    # do the downsampling and construct pcof0
                    pcof0 = np.zeros(0) # to hold the downsampled numpy array for the control vector
                    fact = 2e-3*np.pi # conversion factor from MHz to rad/ns
                    
                    for iosc in range(Nsys):
                        Nelem = np.size(pt0[iosc])
                        dt = (self.nsteps*self.dT)/(Nelem-1) # time step corresponding to (pt0, qt0)
                        p_seg = pt0[iosc]
                        q_seg = qt0[iosc]

                        seg_re = np.zeros(self.nsplines) # to hold downsampled amplitudes
                        seg_im = np.zeros(self.nsplines)
                        # downsample p_seg, q_seg
                        for i_spl in range(self.nsplines):
                            # the B-spline0 coefficients correspond to the time levels
                            t_spl = (i_spl)*self.spline_knot_spacing
                            # i = max(0, np.rint(t_spl/dt).astype(int))# given t_spl, find the closest time step index
                            i = int(np.rint(t_spl / dt))
                            i = min(i, Nelem-1) # make sure i is in range
                            seg_re[i_spl] = fact * p_seg[i]
                            seg_im[i_spl] = fact * q_seg[i]

                        pcof0 = np.append(pcof0, seg_re) # append segment to the global control vector
                        pcof0 = np.append(pcof0, seg_im)
                    # print("simulation(): downsampling of (pt0, qt0) completed")
                else:
                    print("simulation(): detected a mismatch in the sizes of (pt0, qt0)")
            elif len(pt0) > 0 or len(qt0) > 0:
                print("simulation(): the length of pt0 or qt0 != Nsys = ", Nsys)
        else:
            print("Downsampling (pt,qt) is only implemented for spline order 0, not ", self.spline_order)
        
        return pcof0


    def __run(self, *, pcof0=[], runtype="optimization", overwrite_popt=False, maxcores=-1, datadir="./run_dir", quandary_exec="", cygwinbash="", mpi_exec="mpirun -np ", batchargs=[]):
        """
        Internal helper function to launch processes to execute the C++ Quandary code:
          1. Writes quandary config files to file system
          2. Envokes subprocesses to run quandary (on multiple cores), or submits the jobs if 'batchargs' is given.
          3. Gathers results from Quandays output directory into python
          4. Evaluate controls on the input sample rate, if given
        """

        datadir = resolve_datadir(datadir)

        # Create quandary data directory and dump configuration file
        os.makedirs(datadir, exist_ok=True)
        config_filename = self.__dump(pcof0=pcof0, runtype=runtype, datadir=datadir)

        # Set default number of cores to the number of initial conditions, unless otherwise specified. Make sure ncores is an integer divisible of ninit.
        ncores_init = self._ninit
        if maxcores > -1:
            ncores_init = min(self._ninit, maxcores)
        for i in range(self._ninit, 0, -1):
            if self._ninit % i == 0:  # i is a factor of ninit
                if i <= ncores_init:
                    ncores_init = i
                    break
        # Set remaining number of cores for petsc
        ncores_petsc = 1
        if maxcores > ncores_init and maxcores % ncores_init == 0:
            ncores_petsc = int(maxcores / ncores_init)
        ncores = ncores_init * ncores_petsc

        # Execute subprocess to run Quandary
        err = execute(runtype=runtype, ncores=ncores, config_filename=config_filename, datadir=datadir, quandary_exec=quandary_exec, verbose=self.verbose, cygwinbash=cygwinbash, mpi_exec=mpi_exec, batchargs=batchargs)
        if self.verbose:
            print("Quandary data dir: ", datadir, "\n")

        # Get results from quandary output files
        if not err:
            # Get results form quandary's output folder and store some
            time, pt, qt, uT, expectedEnergy, population, popt, infidelity, optim_hist = self.get_results(datadir=datadir)
            if (overwrite_popt):
                self.popt = popt[:]
            self.optim_hist = optim_hist
            self.time = time[:]
            self.uT   = uT.copy()
            # if len(pcof0) == 0:
                # self.pcof0 = popt[:]
        else:
            time = []
            pt = []
            qt = []
            infidelity = 1.0
            expectedEnergy = []
            population = []
        
        if len(pcof0)>0:
            self.pcof0 = pcof0[:]

        return time, pt, qt, infidelity, expectedEnergy, population


    def __dump(self, *, pcof0=[], runtype="simulation", datadir="./run_dir"):
        """
        Internal helper function that dumps all configuration options (and target gate, pcof0, Hamiltonian operators) into files for Quandary C++ runs. Returns the name of the configuration file needed for executing Quandary. 
        """

        # If given, write the target gate to file
        if len(self.targetgate) > 0:
            gate_vectorized = np.concatenate((np.real(self.targetgate).ravel(order='F'), np.imag(self.targetgate).ravel(order='F')))
            self._gatefilename = "targetgate.dat"
            with open(os.path.join(datadir, self._gatefilename), "w", newline='\n') as f:
                for value in gate_vectorized:
                    f.write("{:20.13e}\n".format(value))
            if self.verbose:
                print("Target gate written to ", os.path.join(datadir, self._gatefilename))

        # If given, write the target state to file
        if len(self.targetstate) > 0:
            if self._lindblad_solver: # If Lindblad solver, make it a density matrix
                state = np.outer(self.targetstate, np.array(self.targetstate).conj())
            else:
                state = self.targetstate
            vectorized = np.concatenate((np.real(state).ravel(order='F'), np.imag(state).ravel(order='F')))
            self._gatefilename = "targetstate.dat"
            with open(os.path.join(datadir, self._gatefilename), "w", newline='\n') as f:
                for value in vectorized:
                    f.write("{:20.13e}\n".format(value))
            if self.verbose:
                print("Target state written to ", os.path.join(datadir, self._gatefilename))

        # If given, write the initial state to file
        if self.initialcondition[0:4]=="file":
            if self._lindblad_solver: # If Lindblad solver, make it a density matrix
                state = np.outer(self._initialstate, np.array(self._initialstate).conj())
            else:
                state = self._initialstate
            vectorized = np.concatenate((np.real(state).ravel(order='F'), np.imag(state).ravel(order='F')))
            self._initstatefilename = "initialstate.dat"
            with open(os.path.join(datadir, self._initstatefilename), "w", newline='\n') as f:
                for value in vectorized:
                    f.write("{:20.13e}\n".format(value))
            if self.verbose:
                print("Initial state written to ", os.path.join(datadir, self._initstatefilename))


        # If not standard Hamiltonian model, write provided Hamiltonians to a file
        if not self.standardmodel:
            # Write non-standard Hamiltonians to file  
            # Write system Hamiltonian (complex)
            self._hamiltonian_filename_Hsys= "hamiltonian_Hsys.dat"
            H = self.Hsys
            with open(os.path.join(datadir, self._hamiltonian_filename_Hsys), "w", newline='\n') as f:
                f.write("# row col Hsys_real Hsys_imag \n")
                nz = np.nonzero(H)
                for i, j in zip(*nz):
                    v = H[i, j]
                    f.write(f"{i} {j} {v.real:.13e} {v.imag:.13e}\n")
            # Write control Hamiltonians, if given
            if len(self.Hc_re)>0 or len(self.Hc_im)>0:
                self._hamiltonian_filename_Hc = "hamiltonian_Hc.dat"
                with open(os.path.join(datadir, self._hamiltonian_filename_Hc), "w", newline='\n') as f:
                    for iosc, (Hc_re, Hc_im) in enumerate(zip(self.Hc_re, self.Hc_im)):
                        Hc = np.array(Hc_re) + 1j * np.array(Hc_im)
                        nz = np.nonzero(Hc)
                        f.write(f"# oscillator row col Hc_real Hc_imag \n")
                        for i, j in zip(*nz):
                            v = Hc[i, j]
                            f.write(f"{iosc} {i} {j} {v.real:.13e} {v.imag:.13e}\n")
            if self.verbose:
                print("Hamiltonian operators written to ", os.path.join(datadir, self._hamiltonian_filename_Hsys), os.path.join(datadir, self._hamiltonian_filename_Hc))

        # Initializing the control parameter vector 'pcof0'
        # 1. If the initial parameter vector (list) is given with the 'pcof0' argument, the list will be dumped to a file with name self.pcof0_filename := "pcof0.dat". 
        # 2. If pcof0 is empty but self.pcof_filename is given, use that filename in Quandary 

        read_pcof0_from_file = False
        if len(self.pcof0) > 0 or len(pcof0) > 0:
            if len(self.pcof0) > 0:
                writeme = self.pcof0
            if len(pcof0) > 0: # pcof0 is an argument to __dump(), while self.pcof0 is stored in the object
                writeme = pcof0
            if len(writeme)>0:
                self.pcof0_filename = "pcof0.dat"
                with open(os.path.join(datadir, self.pcof0_filename), "w", newline='\n') as f:
                    for value in writeme:
                        f.write("{:20.13e}\n".format(value))
                if self.verbose:
                    print("Initial control parameters written to ", os.path.join(datadir, self.pcof0_filename))
                read_pcof0_from_file = True
        elif len(self.pcof0_filename) > 0:
            print("Using the provided filename '", self.pcof0_filename, "' in the control_initialization command")
            read_pcof0_from_file = True

        # Set up string for Quandary's config file
        Nt = [self.Ne[i] + self.Ng[i] for i in range(len(self.Ng))]
        mystring = "nlevels = " + ",".join([str(i) for i in Nt]) + "\n"
        mystring += "nessential= " + ",".join([str(i) for i in self.Ne]) + "\n"
        mystring += "ntime = " + str(self.nsteps) + "\n"
        # mystring += "dt = " + str(self.T / self.nsteps) + "\n"
        mystring += "dt = " + str(self.dT) + "\n"
        mystring += "transfreq = " + ",".join([str(i) for i in self.freq01]) + "\n"
        mystring += "rotfreq= " + ",".join([str(i) for i in self.rotfreq]) + "\n"
        mystring += "selfkerr = " + ",".join([str(i) for i in self.selfkerr]) + "\n"
        if len(self.crosskerr)>0:
            mystring += "crosskerr= " + ",".join([str(i) for i in self.crosskerr]) + "\n"
        else:
            mystring += "crosskerr= 0.0\n"
        if len(self.Jkl)>0:
            mystring += "Jkl= " + ",".join([str(i) for i in self.Jkl]) + "\n"
        else:
            mystring += "Jkl= 0.0\n"
        decay = dephase = False
        if len(self.T1) > 0: 
            decay = True
            mystring += "decay_time = " + ",".join([str(i) for i in self.T1]) + "\n"
        if len(self.T2) > 0:
            dephase = True
            mystring += "dephase_time = " + ",".join([str(i) for i in self.T2]) + "\n"
        if decay and dephase:
            mystring += "collapse_type = both\n"
        elif decay:
            mystring += "collapse_type = decay\n"
        elif dephase:
            mystring += "collapse_type = dephase\n"
        else:
            mystring += "collapse_type = none\n"
        if self.initialcondition[0:4] == "file":
            mystring += "initialcondition = " + str(self.initialcondition) + ", " + self._initstatefilename + "\n"
        else:
            mystring += "initialcondition = " + str(self.initialcondition) + "\n"
        for iosc in range(len(self.Ne)):
            if self.spline_order == 0:
                mystring += "control_segments" + str(iosc) + " = spline0, " + str(self.nsplines) + "\n"
            elif self.spline_order == 2:
                mystring += "control_segments" + str(iosc) + " = spline, " + str(self.nsplines) + "\n"
            else:
                print("Error: spline order = ", self.spline_order, " is currently not available. Choose 0 or 2.")
                return -1
            # if len(self.pcof0_filename)>0:
            if read_pcof0_from_file:
                initstring = "file, "+str(self.pcof0_filename) + "\n"
            else:
                # Scale initial control amplitudes by the number of carrier waves and convert to ns
                initamp = self.initctrl_MHz[iosc] /1000.0 / np.sqrt(2) / len(self.carrier_frequency[iosc])
                initstring = ("random, " if self.randomize_init_ctrl else "constant, ") + str(initamp) + "\n"
            mystring += "control_initialization" + str(iosc) + " = " + initstring 
            if len(self.maxctrl_MHz) == 0: # Disable bounds, if not specified
                boundval = 1e+12
            else:
                boundval = self.maxctrl_MHz[iosc]/1000.0  # Scale to ns
            mystring += "control_bounds" + str(iosc) + " = " + str(boundval) + "\n"
            mystring += "carrier_frequency" + str(iosc) + " = "
            omi = self.carrier_frequency[iosc]
            for j in range(len(omi)):
                mystring += str(omi[j]) + ", " 
            mystring += "\n"
        mystring += "control_enforceBC = " + str(self.control_enforce_BC)+ "\n"
        if len(self._gatefilename) > 0:
            mystring += "optim_target = " + self.optim_target + ", " + self._gatefilename + "\n"
        else: 
            mystring += "optim_target = " + str(self.optim_target) + "\n"
        mystring += "optim_objective = " + str(self.costfunction) + "\n"
        mystring += "gate_rot_freq = "
        for val in self.gate_rot_freq:
            mystring += str(val) + ", " 
        mystring += "\n"
        mystring += "optim_weights= 1.0\n"
        mystring += "optim_atol= " + str(self.tol_gnorm_abs) + "\n"
        mystring += "optim_rtol= " + str(self.tol_gnorm_rel) + "\n"
        mystring += "optim_ftol= " + str(self.tol_costfunc) + "\n"
        mystring += "optim_inftol= " + str(self.tol_infidelity) + "\n"
        mystring += "optim_maxiter= " + str(self.maxiter) + "\n"
        if self.gamma_tik0_interpolate > 0.0:
            mystring += "optim_regul= " + str(self.gamma_tik0_interpolate) + "\n"
            mystring += "optim_regul_tik0 = true\n"
        else:
            mystring += "optim_regul= " + str(self.gamma_tik0) + "\n"
            mystring += "optim_regul_tik0=false\n"
        mystring += "optim_penalty= " + str(self.gamma_leakage) + "\n"
        mystring += "optim_penalty_param= 0.0\n"
        mystring += "optim_penalty_dpdm= " + str(self.gamma_dpdm) + "\n"
        mystring += "optim_penalty_variation= " + str(self.gamma_variation) + "\n"
        mystring += "optim_penalty_energy= " + str(self.gamma_energy) + "\n"
        mystring += "datadir= ./\n"
        for iosc in range(len(self.Ne)):
            mystring += "output" + str(iosc) + "=expectedEnergy, population, fullstate\n"
        mystring += "output_frequency = 1\n"
        mystring += "optim_monitor_frequency = " + str(self.print_frequency_iter) + "\n"
        mystring += "runtype = " + runtype + "\n"
        if len(self.Ne) < 6:
            mystring += "usematfree = " + str(self.usematfree) + "\n"
        else:
            mystring += "usematfree = false\n"
        mystring += "linearsolver_type = gmres\n"
        mystring += "linearsolver_maxiter = 20\n"
        if not self.standardmodel:
            if len(self._hamiltonian_filename_Hsys) > 0:
                mystring += "hamiltonian_file_Hsys= "+str(self._hamiltonian_filename_Hsys)+"\n"
            if len(self._hamiltonian_filename_Hc) > 0:
                mystring += "hamiltonian_file_Hc= "+str(self._hamiltonian_filename_Hc)+"\n"
        mystring += "timestepper = "+str(self.timestepper)+ "\n"
        if self.rand_seed is not None and self.rand_seed >= 0:
            mystring += "rand_seed = "+str(int(self.rand_seed))+ "\n"

        # Write the file
        outpath = os.path.join(datadir, "config.cfg")
        with open(outpath, "w", newline='\n') as file:
            file.write(mystring)

        if self.verbose:
            print("Quandary config file written to:", outpath)

        return "./config.cfg"


    def get_results(self, *, datadir="./", ignore_failure=False):
        """ 
        Loads Quandary's output into python.

        Parameters:
        -----------
        datadir        (string) : Directory containing Quandary's output files.
                                  If $QUANDARY_BASE_DATADIR is set, this will be relative to that directory, otherwise relative to current working directory
        ignore_failure (bool)   : Flag to ignore warning when an expected file can't be found

        Returns:
        ---------
        time            :  List of time-points at which the controls are evaluated  (List)
        pt, qt          :  p,q-control pulses at each time point for each oscillator (List of list)
        infidelity      :  Infidelity of the pulses for the specified target operation (Float)
        expectedEnergy  :  Evolution of the expected energy of each oscillator and each initial condition. Acces: expectedEnergy[oscillator][initialcondition]
        population      :  Evolution of the population of each oscillator, of each initial condition. (expectedEnergy[oscillator][initialcondition])
        """

        datadir = resolve_datadir(datadir)

        # Get control parameters
        filename = os.path.join(datadir, "params.dat")
        try:
            pcof = np.loadtxt(filename).astype(float)
        except:
            if not ignore_failure:
                print("Can't read control coefficients from ", filename)
            pcof=[]
    
        # Get optimization history information
        filename = os.path.join(datadir, "optim_history.dat")
        try:
            optim_hist_tmp = np.loadtxt(filename)
        except:
            if not ignore_failure:
                print("Can't read optimization history from ", filename)
            optim_hist_tmp = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    
        if optim_hist_tmp.ndim == 2:
            optim_last = optim_hist_tmp[-1]
        else:
            optim_last = optim_hist_tmp
            optim_hist_tmp = np.array([optim_hist_tmp])
        infid_last = 1.0 - optim_last[4]
    
        # Put optimization history into a dictionary: 
        optim_hist = {
            "Iters"                  : optim_hist_tmp[:,0],
            "Gradient"               : optim_hist_tmp[:,2],
            "Fidelity"               : optim_hist_tmp[:,4],
            "Cost"                   : optim_hist_tmp[:,5],
            "Tikhonov"               : optim_hist_tmp[:,6],
            "Penalty-Leakage"        : optim_hist_tmp[:,7],
            "Penalty-StateVariation" : optim_hist_tmp[:,8],
            "Penalty-TotalEnergy"    : optim_hist_tmp[:,9],
        }


        # Number of colums to be returned. If Schroedinger: all. If Lindblad: Only diagonals. 
        ninits = self._ninit if not self._lindblad_solver else int(np.sqrt(self._ninit))
    
        # Get the time-evolution of the expected energy for each qubit, for each (diagonal) initial condition
        expectedEnergy = [[] for _ in range(len(self.Ne))] 
        for iosc in range(len(self.Ne)):
            for iinit in range(ninits):
                iid = iinit if not self._lindblad_solver else iinit*ninits + iinit
                filename = os.path.join(datadir, f"expected{iosc}.iinit{str(iid).zfill(4)}.dat")
                try:
                    x = np.loadtxt(filename)
                    expectedEnergy[iosc].append(x[:,1])    # 0th column is time, second column is expected energy
                except:
                    if not ignore_failure:
                        print("Can't read expected energy from ", filename)
    
        # Get population for each qubit, for each initial condition
        population = [[] for _ in range(len(self.Ne))]
        for iosc in range(len(self.Ne)):
            for iinit in range(ninits):
                iid = iinit if not self._lindblad_solver else iinit*ninits + iinit
                filename = os.path.join(datadir, f"population{iosc}.iinit{str(iid).zfill(4)}.dat")
                try:
                    x = np.loadtxt(filename)
                    population[iosc].append(x[:,1:].transpose())    # first column is time
                except:
                    if not ignore_failure:
                        print("Can't read population from ", filename)
    
        # Get last time-step unitary
        Ntot = [i+j for i,j in zip(self.Ne,self.Ng)]
        ndim = np.prod(Ntot) if not self._lindblad_solver else np.prod(Ntot)**2
        uT = np.zeros((ndim, self._ninit), dtype=complex)
        for iinit in range(self._ninit):
            file_index = str(iinit).zfill(4)
            try:
                xre = np.loadtxt(os.path.join(datadir, f"rho_Re.iinit{file_index}.dat"), skiprows=1, usecols=range(1, ndim+1))[-1]
                uT[:, iinit] = xre 
            except:
                name = os.path.join(datadir, f"rho_Re.iinit{file_index}.dat")
                if not ignore_failure:
                    print("Can't read from ", name)
            try:
                xim = np.loadtxt(os.path.join(datadir, f"rho_Im.iinit{file_index}.dat"), skiprows=1, usecols=range(1, ndim+1))[-1]
                uT[:, iinit] += 1j * xim
            except:
                name = os.path.join(datadir, f"rho_Im.iinit{file_index}.dat")
                if not ignore_failure:
                    print("Can't read from ", name)
            # uT[:, iinit] = xre + 1j * xim
    
        # Get the control pulses for each qubit
        pt = []
        qt = []
        ft = []
        for iosc in range(len(self.Ne)):
            # Read the control pulse file
            filename = os.path.join(datadir, f"control{iosc}.dat")
            try:
                x = np.loadtxt(filename)
            except:
                print("Can't read control pulses from ", filename)
                x = np.zeros((1,4))
            # Extract the pulses 
            time = x[:,0]   # Time domain
            pt.append([x[n,1]*1e+3 for n in range(len(x[:,0]))])     # Rot frame p(t), MHz
            qt.append([x[n,2]*1e+3 for n in range(len(x[:,0]))])     # Rot frame q(t), MHz
            ft.append([x[n,3]*1e+3 for n in range(len(x[:,0]))])     # Lab frame f(t)
    
        return time, pt, qt, uT, expectedEnergy, population, pcof, infid_last, optim_hist


def estimate_timesteps(*, T=1.0, Hsys=[], Hc_re=[], Hc_im=[], maxctrl_MHz=[], Pmin=40):
    """
    Helper function to estimate the number of time steps based on eigenvalues of the system Hamiltonian and maximum control Hamiltonians. Note: The estimate does not account for quickly varying signals or a large number of splines. Double check that at least 2-3 points per spline are present to resolve control function. #TODO: Automate this.
    """

    # Get estimated control pulse amplitude
    est_ctrl_MHz = maxctrl_MHz[:]
    if len(maxctrl_MHz) == 0:
        est_ctrl_MHz = [10.0 for _ in range(max(len(Hc_re), len(Hc_im)))] 

    # Set up Hsys +  maxctrl*Hcontrol
    K1 = np.copy(Hsys) 

    for i in range(len(Hc_re)):
        est_radns = est_ctrl_MHz[i]*2.0*np.pi/1e+3
        if len(Hc_re[i])>0:
            K1 += est_radns * Hc_re[i] 
    for i in range(len(Hc_im)):
        est_radns = est_ctrl_MHz[i]*2.0*np.pi/1e+3
        if len(Hc_im[i])>0:
            K1 = K1 + 1j * est_radns * Hc_im[i] # can't use += due to type!
    
    # Estimate time step
    eigenvalues = np.linalg.eigvals(K1)
    maxeig = np.max(np.abs(eigenvalues))
    # ctrl_fac = 1.2  # Heuristic, assuming that the total Hamiltonian is dominated by the system part.
    ctrl_fac = 1.0
    samplerate = ctrl_fac * maxeig * Pmin / (2 * np.pi)
#     print(f"{samplerate=}")
    nsteps = int(np.ceil(T * samplerate))

    return nsteps


def eigen_and_reorder(H0, verbose=False):
    """ Internal function that computes eigen decomposition and re-orders it to make the eigenvector matrix as close to the identity as possible """

    # Get eigenvalues and vectors and sort them in ascending order
    Ntot = H0.shape[0]
    evals, evects = np.linalg.eig(H0)
    reord = np.argsort(evals)
    evals = evals[reord]
    evects = evects[:,reord]

    # Find the column index corresponding to the largest element in each row of evects 
    max_col = np.zeros(Ntot, dtype=np.int32)
    for row in range(Ntot):
        max_col[row] = np.argmax(np.abs(evects[row,:]))

    # test the error detection
    # max_col[1] = max_col[0]

    # loop over all columns and check max_col for duplicates
    Ndup_col = 0 
    for row in range(Ntot-1): 
        for k in range(row+1, Ntot):
            if max_col[row] == max_col[k]:
                Ndup_col += 1
                print("Error: detected identical max_col =", max_col[row], "for rows", row, "and", k)


    if Ndup_col > 0:
        print("Found", Ndup_col, "duplicate column indices in max_col array")
        raise ValueError('Permutation of eigen-vector matrix failed')

    evects = evects[:,max_col]
    evals = evals[max_col]
    
    # Make sure all diagonal elements are positive
    for j in range(Ntot):
        if evects[j,j]<0.0:
            evects[:,j] = - evects[:,j]

    return evals, evects


def get_resonances(*, Ne, Ng, Hsys, Hc_re=[], Hc_im=[], rotfreq=[], cw_amp_thres=1e-7, cw_prox_thres=1e-2,verbose=True, stdmodel=True):
    """ 
    Computes system resonances, to be used as carrier wave frequencie.
    Returns resonance frequencies in GHz and corresponding growth rates.
    """ 

    if verbose:
        print("\nComputing carrier frequencies, ignoring growth rate slower than:", cw_amp_thres, "and frequencies closer than:", cw_prox_thres, "[GHz])")

    nqubits = len(Ne)
    n = Hsys.shape[0]
    
    # Get eigenvalues of system Hamiltonian (GHz)
    Hsys_evals, Utrans = eigen_and_reorder(Hsys, verbose)
    Hsys_evals = Hsys_evals.real  # Eigenvalues may have a small imaginary part due to numerical imprecision
    Hsys_evals = Hsys_evals / (2 * np.pi) 
            
    # Look for resonances in the symmetric and anti-symmetric control Hamiltonians for each qubit
    resonances = []
    speed = []
    for q in range(nqubits):
       
        # Transform symmetric and anti-symmetric control Hamiltonians using eigenvectors (reordered)
        Hsym_trans = Utrans.conj().T @ Hc_re[q] @ Utrans
        Hanti_trans = Utrans.conj().T @ Hc_im[q] @ Utrans

        resonances_a = []
        speed_a = []
        if verbose:
            print("  Resonances in oscillator #", q)
        
        for Hc_trans in (Hsym_trans, Hanti_trans):

            # Iterate over non-zero elements in transformed control
            for i in range(n):
                # Only consider transitions from lower to higher levels
                for j in range(i):

                    # Look for non-zero elements (skip otherwise)
                    if abs(Hc_trans[i,j]) < 1e-14: 
                        continue 

                    # Get the resonance frequency
                    delta_f = Hsys_evals[i] - Hsys_evals[j]
                    if abs(delta_f) < 1e-10:
                        delta_f = 0.0

                    # Get involved oscillator levels
                    ids_i = map_to_oscillators(i, Ne, Ng)
                    ids_j = map_to_oscillators(j, Ne, Ng)

                    # make sure both indices correspond to essential energy levels
                    is_ess_i = all(ids_i[k] < Ne[k] for k in range(len(Ne)))
                    is_ess_j = all(ids_j[k] < Ne[k] for k in range(len(Ne)))

                    if (is_ess_i and is_ess_j):
                        # Ignore resonances that are too close by comparing to all previous resonances
                        if any(abs(delta_f - f) < cw_prox_thres for f in resonances_a):
                            if verbose:
                                print("    Ignoring resonance from ", ids_j, "to ", ids_i, ", freq", delta_f, ", growth rate=", abs(Hc_trans[i, j]), "being too close to one that already exists.")
                        # Ignore resonances with growth rate smaller than user-defined threshold
                        elif abs(Hc_trans[i,j]) < cw_amp_thres:
                            if verbose:
                                print("    Ignoring resonance from ", ids_j, "to ", ids_i, ", freq", delta_f, ", growth rate=", abs(Hc_trans[i, j]), "growth rate is too slow.")
                        # Otherwise, add resonance to the list
                        else:
                            resonances_a.append(delta_f)
                            speed_a.append(abs(Hc_trans[i, j]))
                            if verbose:
                                print("    Resonance from ", ids_j, "to ", ids_i, ", freq", delta_f, ", growth rate=", abs(Hc_trans[i, j]))

        # Append resonances for this qubit to overall list        
        resonances.append(resonances_a)
        speed.append(speed_a)

    # Prepare output for carrier frequencies (om) and growth_rate
    Nfreq = np.zeros(nqubits, dtype=int)
    om = [[0.0] for _ in range(nqubits)]
    growth_rate = [[] for _ in range(nqubits)]
    for q in range(len(resonances)):
        Nfreq[q] = max(1, len(resonances[q]))  # at least one being 0.0
        om[q] = np.zeros(Nfreq[q])
        if len(resonances[q]) > 0:
            om[q] = np.array(resonances[q])
        growth_rate[q] = np.ones(Nfreq[q])
        if len(speed[q]) > 0:
            growth_rate[q] = np.array(speed[q])

    return om, growth_rate


def lowering(n):
    """ Lowering operator of dimension n """
    return np.diag(np.sqrt(np.arange(1, n)), k=1)
def number(n):
    """ Number operator of dimension n """
    return np.diag(np.arange(n))
def map_to_oscillators(id, Ne, Ng):
    """ Return the local energy level of each oscillator for a given global index id """
    # len(Ne) = number of subsystems
    nlevels = [Ne[i]+Ng[i] for i in range(len(Ne))]
    localIDs = []

    index = int(id)
    for iosc in range(len(Ne)):
        postdim = np.prod(nlevels[iosc+1:])
        localIDs.append(int(index / postdim))
        index = index % postdim 

    return localIDs 

def resolve_datadir(datadir: str) -> str:
    """Helper function to resolve the output directory using environment variable

    Parameters:
    ----------
    datadir : Output directory
              If $QUANDARY_BASE_DATADIR is set, relative paths will be resolved against it

    Returns:
    -------
    Resolved absolute or relative path for the data directory

    Raises:
    ------
    ValueError: If QUANDARY_BASE_DATADIR is set but doesn't exist or isn't a directory
    """
    if os.path.isabs(datadir):
        return datadir

    base_dir = os.environ.get("QUANDARY_BASE_DATADIR")
    if base_dir:
        if not os.path.exists(base_dir):
            raise ValueError(f"Environment variable QUANDARY_BASE_DATADIR points to non-existent path: {base_dir}")
        if not os.path.isdir(base_dir):
            raise ValueError(f"Environment variable QUANDARY_BASE_DATADIR is not a directory: {base_dir}")

        datadir = os.path.join(base_dir, datadir)

    return os.path.normpath(datadir)


def hamiltonians(*, N, freq01, selfkerr, crosskerr=[], Jkl = [], rotfreq=[], verbose=True):
    """  
    Create standard Hamiltonian operators to model pulse-driven superconducting qubits.

    Parameters:
    -----------
    N        :  Total levels per oscillator (essential plus guard levels for each qubit)
    freq 01  :  Groundstate frequency for each qubit [GHz]
    selfkerr :  Self-kerr coefficient for each qubit [GHz]

    Optional Parameters:
    ---------------------
    Crosskerr  : ZZ-coupling strength [GHz]
    Jkl        : dipole-dipole coupling strength [GHz]
    rotfreq    : Rotational frequencies for each qubit
    verbose    : Switch to turn on more output. Default: True

    Returns:
    --------
    Hsys        : System Hamiltonian (time-independent), units rad/ns
    Hc_re       : Real parts of control Hamiltonian operators for each qubit (Hc = [ [Hc_qubit1], [Hc_qubit2],... ]). Unit-less.
    Hc_im       : Imag parts of control Hamiltonian operators for each qubit (Hc = [ [Hc_qubit1], [Hc_qubit2],... ]). Unit-less.
    """

    if len(rotfreq) == 0:
        rotfreq=np.zeros(len(N))

    nqubits = len(N)
    assert len(selfkerr) == nqubits
    assert len(freq01) == nqubits
    assert len(rotfreq) == nqubits

    n = np.prod(N)     # System size 

    # Set up lowering operators in full dimension
    Amat = []
    for i in range(len(N)):
        ai = lowering(N[i])
        for j in range(i):
            ai = np.kron(np.identity(N[j]), ai) 
        for j in range(i+1,len(N)):
            ai = np.kron(ai, np.identity(N[j]))
        #print("Amat i =", i)
        #print(ai) 
        Amat.append(ai)

    # Set up system Hamiltonian: Duffing oscillators
    Hsys = np.zeros((n, n))
    for q in range(nqubits):
        domega_radns =  2.0*np.pi * (freq01[q] - rotfreq[q])
        selfkerr_radns = 2.0*np.pi * selfkerr[q]
        Hsys +=  domega_radns * Amat[q].T @ Amat[q]
        Hsys -= selfkerr_radns/2.0 * Amat[q].T @ Amat[q].T @ Amat[q] @ Amat[q]

    # Add cross cerr coupling, if given
    if len(crosskerr)>0:
        idkl = 0 
        for q in range(nqubits):
            for p in range(q + 1, nqubits):
                if abs(crosskerr[idkl]) > 1e-14:
                    crosskerr_radns = 2.0*np.pi * crosskerr[idkl]
                    Hsys -= crosskerr_radns * Amat[q].T @ Amat[q] @ Amat[p].T @ Amat[p]
                idkl += 1
    
    # Add Jkl coupling term. 
    # Note that if the rotating frame frequencies are different amongst oscillators, then this contributes to a *time-dependent* system Hamiltonian. Here, we treat this as time-independent, because this Hamiltonian here is *ONLY* used to compute the time-step size and resonances, and it is NOT passed to the quandary code. Quandary sets up the standard model with a time-dependent system Hamiltonian if the frequencies of rotation differ amongst oscillators.  
    if len(Jkl)>0:
        idkl = 0 
        for q in range(nqubits):
            for p in range(q + 1, nqubits):
                if abs(Jkl[idkl]) > 1e-14:
                    Jkl_radns  = 2.0*np.pi*Jkl[idkl]
                    Hsys += Jkl_radns * (Amat[q].T @ Amat[p] + Amat[q] @ Amat[p].T)
                idkl += 1

    # Set up control Hamiltonians
    Hc_re = [Amat[q] + Amat[q].T for q in range(nqubits)]
    Hc_im = [Amat[q] - Amat[q].T for q in range(nqubits)]

    if verbose:
        print(f"*** {nqubits} coupled quantum systems setup ***")
        print("System Hamiltonian frequencies [GHz]: f01 =", freq01, "rot. freq =", rotfreq)
        print("Selfkerr=", selfkerr)
        print("Coupling: X-Kerr=", crosskerr, ", J-C=", Jkl)

    return Hsys, Hc_re, Hc_im


def plot_pulse(Ne, time, pt, qt):
    """ 
    Plot the control pulse for all qubits
    """
    fig = plt.figure()
    nrows = len(Ne)
    ncols = 1
    for iosc in range(len(Ne)):
        plt.subplot(nrows, ncols, iosc+1)
        plt.plot(time, pt[iosc], "r", label="p(t)")
        plt.plot(time, qt[iosc], "b", label="q(t)")
        plt.xlabel('time (ns)')
        plt.ylabel('Drive strength [MHz]')
        maxp = max(np.abs(pt[iosc]))
        maxq = max(np.abs(qt[iosc]))
        plt.title('Qubit '+str(iosc)+'\n max. drive '+str(round(maxp,1))+", "+str(round(maxq,1))+" MHz")
        plt.legend(loc='lower right')
        plt.xlim([0.0, time[-1]])
    # plt.grid()
    plt.subplots_adjust(hspace=0.6)
    plt.draw()
    plt.show() 

def plot_expectedEnergy(Ne, time, expectedEnergy):
    """ Plot evolution of expected energy levels """

    ninit = len(expectedEnergy[0])
    nplots = ninit                    # one plot for each initial state
    # nplots = np.prod(Ne)                # one plot for each initial state
    ncols = 2 if nplots >= 4 else 1     # 2 rows if more than 3 plots
    nrows = int(np.ceil(nplots/ncols))
    figsizex = 6.4*nrows*0.75 
    figsizey = 4.8*nrows*0.75 
    fig = plt.figure(figsize=(figsizex,figsizey))
    for iplot in range(nplots):
        iinit = iplot 
        plt.subplot(nrows, ncols, iplot+1)
        plt.figsize=(15, 15)
        emax = 1.0
        for iosc in range(len(Ne)):
            label = 'Qubit '+str(iosc) if len(Ne)>1 else ''
            plt.plot(time, expectedEnergy[iosc][iinit], label=label)
            emax_iosc = np.max(expectedEnergy[iosc][iinit])
            emax = max(emax, emax_iosc) # keep track of max energy level for setting ylim
        plt.xlabel('time (ns)')
        plt.ylabel('expected energy')

        plt.ylim([0.0-1e-2, emax + 1e-2])
        plt.xlim([0.0, time[-1]])
        binary_ID = iplot if len(Ne) == 1 else bin(iplot).replace("0b", "").zfill(len(Ne))
        plt.title("from |"+str(binary_ID)+">")
        plt.legend(loc='lower right')
    plt.subplots_adjust(hspace=0.5)
    plt.subplots_adjust(wspace=0.5)
    plt.draw()
    plt.show()

def plot_population(Ne, time, population):
    """ Plot evolution of population """ 

    ninit = len(population[0])
    # nplots = np.prod(Ne)                # one plot for each initial state
    nplots = ninit                      # one plot for each initial state
    ncols = 2 if nplots >= 4 else 1     # 2 rows if more than 3 plots
    nrows = int(np.ceil(nplots/ncols))
    figsizex = 6.4*nrows*0.75 
    figsizey = 4.8*nrows*0.75 
    fig = plt.figure(figsize=(figsizex,figsizey))

    # Iterate over initial conditions (one plot for each)
    for iplot in range(nplots):
        iinit = iplot 
        plt.subplot(nrows, ncols, iplot+1)
        plt.figsize=(15, 15)
        for iosc in range(len(Ne)):
            for istate in range(Ne[iosc]):
                label = 'Qubit '+str(iosc) if len(Ne)>1 else ''
                label = label + " |"+str(istate)+">"
                plt.plot(time, population[iosc][iinit][istate], label=label)
        plt.xlabel('time (ns)')
        plt.ylabel('population')
        plt.ylim([0.0-1e-4, 1.0 + 1e-2])
        plt.xlim([0.0, time[-1]])
        binary_ID = iplot if len(Ne) == 1 else bin(iplot).replace("0b", "").zfill(len(Ne))
        plt.title("from |"+str(binary_ID)+">")
        plt.legend(loc='lower right')
    plt.subplots_adjust(hspace=0.5)
    plt.subplots_adjust(wspace=0.5)
    plt.draw()
    plt.show()

def plot_results_1osc(myconfig, p, q, expectedEnergy, population):
    """ Plot all results of one oscillator """

    fig, ax = plt.subplots(2, 3, figsize=(20,8))
    fig.subplots_adjust(hspace=0.3)

    t = myconfig.time

    # Plot pulses
    ax[0,0].plot(t, p, label='I') # y label: MHz
    ax[0,0].plot(t, q, label='Q') # y label: MHz
    ax[0,0].set_ylabel('Pulse amplitude (MHz)')
    ax[0,0].set_xlabel('Time (ns)')
    ax[0,0].legend()
    ax[0,0].grid()


    # Compute and plot FFT
    zlist = np.array(p)*1e-3 + 1j*np.array(q)*1e-3
    fft = np.fft.fft(zlist)
    dt = myconfig.T / myconfig.nsteps
    fftfr = np.fft.fftfreq(len(zlist), d=dt)

    ax[0,1].scatter(fftfr*1e3, np.abs(fft)**2)
    ax[0,1].set_ylabel('FFT')
    ax[0,1].set_xlabel('Frequency (MHz)')
    ax[0,1].grid()
    ax[0,1].set_title('FFT')
    ax[0,1].set_yscale('log')
    ax[0,1].set_xlim(-500, 500)
    ax[0,1].set_ylim(1e-8, 1e5)

    # Plot Populations for each initial condition 
    for iinit in range(len(population)):  # for each of the initial states
        row = 1
        col = iinit
            
        for istate in range(myconfig.Ne[0]): # for each essential level
            label = "|"+str(istate)+">"
            ax[row, col].plot(t, population[iinit][istate], label=label)
            # ax[row, col+1].plot(np.arange(0, numgate), prob_me_gate[i].real, label=str(i))
            
        ax[row, col].set_xlabel('Time (ns)')
        ax[row, col].set_ylabel('Population')
        ax[row, col].legend()
        ax[row, col].set_title('Populations from |%d>' % iinit)
        ax[row, col].grid()
        
        # ax[row, col+1].set_xlabel('Gate repetition')
        # ax[row, col+1].set_ylabel('Population')
        # ax[row, col+1].legend()
        # ax[row, col+1].set_title('ME solve, starting from %d' % state)

    # Plot expected Energy
    row, col = 0, 2
    for iinit in range(len(expectedEnergy)):
        label = 'from |'+str(iinit)+'>' 
        ax[row, col].plot(t, expectedEnergy[iinit], label=label)
    ax[row, col].set_xlabel('Time (ns)')
    ax[row, col].set_ylabel('Expected Energy Level')
    ax[row, col].legend()
    ax[row, col].set_title('Expected Energy Level')
    ax[row, col].grid()

    plt.draw()
    plt.show()

def timestep_richardson_est(quandary, tol=1e-8, order=2, quandary_exec=""):
    """ Decrease timestep size until Richardson error estimate meets threshold """

    # Factor by which timestep size is decreased
    m = 2 
    
    # Initialize 
    quandary.verbose=False
    t, pt, qt, infidelity, _, _ = quandary.simulate(quandary_exec=quandary_exec, datadir="TS_test")

    Jcurr = infidelity
    uT = quandary.uT.copy()

    # Loop
    errs_J = []
    errs_u = []
    dts = []
    for i in range(10):

        # Update configuration number of timesteps. Note: dt will be set in dump()
        dt_org = quandary.T / quandary.nsteps
        quandary.nsteps= quandary.nsteps* m
        quandary.dT = quandary.T / quandary.nsteps

        # Get u(dt/m) 
        quandary.verbose=False
        t, pt, qt, infidelity, _, _ = quandary.simulate(quandary_exec=quandary_exec, datadir="TS_test")

        # Richardson error estimate 
        err_J = np.abs(Jcurr - infidelity) / (m**order-1.0)
        # err_u = np.abs(uT[1,1]- quandary.uT[1,1]) / (m**order - 1.0)
        err_u = np.linalg.norm(np.subtract(uT, quandary.uT)) / (m**order - 1.0)
        errs_J.append(err_J)
        errs_u.append(err_u)
        dts.append(dt_org)

        # Output 
        print(" -> Error at i=", i, ", dt = ", dt_org, ": err_J = ", err_J, " err_u=", err_u)


        # Stop if tolerance is reached
        if err_J < tol:
            print("\n -> Tolerance reached. N=", quandary.nsteps, ", dt=",dt_org)
            break

        # Update
        Jcurr = infidelity
        uT = np.copy(quandary.uT)

    return errs_J, errs_u, dts


def execute(*, runtype="simulation", ncores=1, config_filename="config.cfg", datadir=".", quandary_exec="", verbose=False, cygwinbash="", mpi_exec="mpirun -np ", batchargs=[]):
    """ 
    Helper function to evoke a subprocess that executes Quandary.

    Optional Parameters:
    -----------
    ncores              (int)       : Number of cores to run on. Default: 1 
    config_filename     (string)    : Name of Quandaries configuration file. Default: "config.cfg"
    datadir             (string)    : Directory for running Quandary and output files. Default: "./"
    quandary_exec       (string)    : Absolute path to quandary's executable. Default: "" (expecting quandary to be in the $PATH)
    verbose             (Bool)      : Flag to print more output. Default: False
    cygwinbash          (string)    : Path to Cygwin bash.exe, if running on Windows machine. Default: None
    mpi_exec            (string)    : MPI launcher prefix, e.g. "mpirun -np" or "srun -n". The string should include the flag for core counts, but not the number of cores itself which will be appended automatically
    batchargs           (List)      : Submit to batch system by setting batchargs= [maxime, accountname, nodes]. Default: []. Compare end of this file. Specify the max. runtime (string), the account name (string) and the number of requested nodes (int). Note, the number of executing *cores* is defined through 'ncores'. 

    Returns:
    ---------
    err     (int)  : Error indicator (didn't run successfully if this this 1. Hence, load the results only if this is 0)
    """

    # Enter the directory where Quandary will be run
    dir_org = os.getcwd() 
    os.chdir(datadir)
 
    # Set up the run command
    exe = "quandary"
    if len(quandary_exec) > 0:
        exe = quandary_exec
    runcommand = f"{exe} {config_filename}"
    if not verbose:
        runcommand += " --quiet"
    # If parallel run, specify runcommand. Default is 'mpirun -np ', unless batch args are given, then currently using 'srun -n', see end of this file for changing that to other batch systems.
    if len(batchargs)>0:
        myrun = batch_run  # currently set to "srun -n"
    else:
        myrun = mpi_exec
    runcommand = f"{myrun} {ncores} " + runcommand
    if verbose:
        print("Running Quandary ... ")

    # If batchargs option is given, submit a batch script to schedule execution of Quandary
    if len(batchargs) > 0:
        maxtime, account, nodes = batchargs
        batch_args = copy.deepcopy(default_batch_args)
        batch_args[batch_args_mapping["NAME"]]            = datadir
        batch_args[batch_args_mapping["ERROR"]]           = datadir+".err"
        batch_args[batch_args_mapping["OUTPUT"]]          = datadir+".out"
        batch_args[batch_args_mapping["NTASKS"]]          = ncores
        batch_args[batch_args_mapping["ACCOUNT"]]         = account
        batch_args[batch_args_mapping["NODES"]]           = nodes
        batch_args[batch_args_mapping["TIME"]]            = maxtime
        assemble_batch_script(datadir+".batch", runcommand, batch_args)
        if True:
            call("sbatch " + datadir+".batch", shell=True)
        err=1
    # If Cygwin option is given, execute on Windows through Cygwin bash
    elif len(cygwinbash)>0: 
        cygwinbash= str(cygwinbash)
        # p = Popen(r"C:/cygwin64/bin/bash.exe", stdin=PIPE, stdout=PIPE, stderr=PIPE)  
        p = Popen(cygwinbash, stdin=PIPE, stdout=PIPE, stderr=PIPE)  
        p.stdin.write(runcommand.encode('ASCII')) 
        p.stdin.close()
        std_out, std_err = p.communicate()
        # Print stdout and stderr
        print(std_out.strip().decode('ascii'))
        print(std_err.strip().decode('ascii'))
        err=0
    # Else: (Default) Execute Quandary locally, for Linux or MacOS 
    else: 
        # Pipe std output to file rather than screen
        # with open(os.path.join(datadir, "out.log"), "w") as stdout_file, \
        #      open(os.path.join(datadir, "err.log"), "w") as stderr_file:
        #         exec = run(runcommand, shell=True, stdout=stdout_file, stderr=stderr_file)
        print("Executing '", runcommand, ". Runtype: ", runtype, "...")
        exec = run(runcommand, shell=True)
        # Check return code
        err = exec.check_returncode()
    if verbose: 
        print("DONE. \n")

    # Return to previous directory
    os.chdir(dir_org)
    return err

####################################
# Batch system settings
####################################

# Define batch argument names, here for SLURM scheduler. You can add others if needed, make sure to 
batch_args_mapping_slurm = {"NAME"  : "--job-name",
                            "ERROR" : "--error",
                            "OUTPUT" : "--output",
                            "TIME"  : "--time",
                            "NTASKS"  : "--ntasks",
                            "NODES"  : "--nodes",
                            "ACCOUNT"  : "--account",
                            "NTASKSPERNODE"  : "--ntasks-per-node",
                            "CONSTRAINT" : "--constraint",
                            "PARTITION" : "--partition"}

# Specify batch system config. Here: Slurm
batch_args_mapping = batch_args_mapping_slurm
batch_run = "srun -n "

# Set default batch arguments
default_batch_args = {batch_args_mapping["NAME"]          : "default",
                      batch_args_mapping["OUTPUT"]        : "default-%j.out",
                      batch_args_mapping["ERROR"]         : "default-%j.err",
                      batch_args_mapping["ACCOUNT"]       : "sqoc",
                      batch_args_mapping["PARTITION"]     : "pbatch",
                      batch_args_mapping["TIME"]          : "00:05:00",
                      batch_args_mapping["NODES"]         : 1,
                      batch_args_mapping["NTASKS"]        : 1}

def assemble_batch_script(name, run_command, batch_args, exclusive=True):
    outfile = open(name, 'w')
    outfile.write("#!/usr/bin/bash\n")
    for arg,value in batch_args.items():
        outfile.write("#SBATCH " + arg + "=" + str(value) + "\n")
    if exclusive:
        outfile.write("#SBATCH --exclusive\n")
    outfile.write(run_command)
    outfile.close()


def infidelity_(A,B):
	dim = int(np.sqrt(A.size))
	return 1.0 - np.abs(np.trace(A.conj().transpose() @ B))**2/dim**2

