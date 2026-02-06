#include <sys/stat.h> 
#include <petscmat.h>
#include <iostream> 
#include "config.hpp"
#include "mastereq.hpp"
#pragma once

/**
 * @brief Output management for quantum control simulations and optimization.
 *
 * This class handles all file output operations for Quandary, including optimization
 * progress logging, control pulse output, state evolution data, and population dynamics.
 * It manages MPI-aware output to avoid duplicate writes in parallel runs.
 */
class Output{
  protected:

  int mpirank_world; ///< Rank of processor in MPI_COMM_WORLD
  int mpirank_petsc; ///< Rank of processor for PETSc parallelization
  int mpisize_petsc; ///< Size of communicator for PETSc parallelization
  int mpirank_init; ///< Rank of processor for initial condition parallelization

  bool quietmode; ///< Flag for reduced screen output
  
  FILE* optimfile; ///< Output file for logging optimization progress
  int output_timestep_stride; ///< Time domain output frequency (write every N time steps)
  std::vector<OutputType> output_observables; ///< List of output types applied to all oscillators

  size_t noscillators; ///< Number of oscillators in the system
  bool writeFullState; ///< Flag to determine if evolution of full state vector should be written to file
  bool writeExpectedEnergy; ///< Flag to determine if evolution of expected energy per oscillator should be written to files
  bool writeExpectedEnergy_comp; ///< Flag to determine if evolution of expected energy of the full composite system should be written to file
  bool writePopulation; ///< Flag to determine if the evolution of the energy level occupations per oscillator should be written to files
  bool writePopulation_comp; ///< Flag to determine if the evolution of the energy level occupations of the full composite system should be written to file
  FILE *ufile; ///< File for writing real part of fullstate evolution
  FILE *vfile; ///< File for writing imaginary part of fullstate evolution
  std::vector<FILE *>expectedfile; ///< Files for expected energy evolution per oscillator
  std::vector<FILE *>populationfile; ///< Files for population evolution per oscillator
  FILE *expectedfile_comp; ///< File for expected energy evolution of the full composite system
  FILE *populationfile_comp; ///< File for population evolution of the full composite system

  // VecScatter scat; ///< PETSc's scatter context for state communication across cores
  // Vec xseq; ///< Sequential vector for I/O operations

  public:
    std::string output_dir; ///< Directory path for output data files
    int output_optimization_stride; ///< Write output files every N optimization iterations

  public:
    Output();

    /**
     * @brief Constructor with configuration and MPI setup.
     *
     * @param config Configuration parameters from input file
     * @param comm_petsc MPI communicator for PETSc parallelization
     * @param comm_init MPI communicator for initial condition parallelization
     * @param quietmode Flag for reduced output (default: false)
     */
    Output(const Config& config, MPI_Comm comm_petsc, MPI_Comm comm_init, bool quietmode=false);

    ~Output();

    /**
     * @brief Writes optimization progress to history file.
     *
     * Called at every optimization iteration to log convergence data. 
     * Optimization history will be written to `<output_dir>/optim_history.dat`.
     *
     * @param optim_iter Current optimization iteration
     * @param objective Total objective function value
     * @param gnorm Gradient norm
     * @param stepsize Optimization step size
     * @param Favg Average fidelity
     * @param cost Final-time cost term
     * @param tikh_regul Tikhonov regularization term
     * @param penalty_leakage Penalty term for leakage
     * @param penalty_weightedcost Weighted cost penalty
     * @param penalty_dpdm Second-order derivative penalty
     * @param penalty_energy Energy penalty term
     * @param penalty_variation Control variation penalty
     */
    void writeOptimFile(int optim_iter, double objective, double gnorm, double stepsize, double Favg, double cost, double tikh_regul,  double penalty_leakage, double penalty_dpdm, double penalty_energy, double penalty_variation, double penalty_weightedcost);

    /**
     * @brief Writes current control pulses per oscillator and control parameters.
     *
     * Called every output_optimization_stride optimization iterations. 
     * Control pulses are written to `<output_dir>/control<ioscillator>.dat`
     * Control parameters are written to `<output_dir>/params.dat`
     *
     * @param params Current parameter vector
     * @param mastereq Pointer to master equation solver
     * @param ntime Total number of time steps
     * @param dt Time step size
     */
    void writeControls(Vec params, MasterEq* mastereq, int ntime, double dt);

    /**
     * @brief Writes gradient vector for debugging adjoint calculations.
     * 
     * Gradient is written to `<output_dir>/grad.dat`
     *
     * @param grad Gradient vector to output
     */
    void writeGradient(Vec grad);

    /**
     * @brief Opens data files for time evolution output for one oscillator.
     *
     * Prepares files for writing full state, expected energy, and population evolution data. Called before timestepping starts. 
     *
     * @param prefix Filename prefix for output files
     * @param initid Initial condition identifier
     */
    void openTrajectoryDataFiles(std::string prefix, int initid);

    /**
     * @brief Writes time evolution data to files.
     *
     * Outputs state vector, expected energies, and populations at current time step. Called at each time step 
     *
     * @param timestep Current time step number
     * @param time Current time value
     * @param state Current state vector
     * @param mastereq Pointer to master equation solver
     */
    void writeTrajectoryDataFiles(int timestep, double time, const Vec state, MasterEq* mastereq);

    /**
     * @brief Closes open time evolution data files.
     *
     * Properly closes and flushes all output files after time-stepping completion.
     */
    void closeTrajectoryDataFiles();

};
