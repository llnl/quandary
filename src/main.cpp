#include "timestepper.hpp"
#include "defs.hpp"
#include <string>
#include "oscillator.hpp" 
#include "mastereq.hpp"
#include "config.hpp"
#include <stdlib.h>
#include <sys/resource.h>
#include <cassert>
#include <algorithm>
#include <cmath>
#include "optimproblem.hpp"
#include "output.hpp"
#include "petsc.h"
#include <random>
#include "util.hpp"
#ifdef WITH_SLEPC
#include <slepceps.h>
#include <slepcsvd.h>
#endif

#define TEST_FD_GRAD 0    // Run Finite Differences gradient test
#define TEST_FD_HESS 0    // Run Finite Differences Hessian test
#define TEST_FD_LINEARIZED_FWD 0 // Run Finite Differences Linearized Forward test
#define TEST_GAUSSNEWTON_LINEARSYSTEM 0
#define TEST_GAUSSNEWTON_LEASTSQUARES 0
#define HESSIAN_DECOMPOSITION 0 // Run eigenvalue analysis for Hessian
#define EPS 1e-5          // Epsilon for Finite Differences

int main(int argc,char **argv)
{
  /* Parse command line arguments */
  ParsedArgs args = parseArguments(argc, argv);

  char filename[255];
  PetscErrorCode ierr;

  /* Initialize MPI */
  int mpisize_world, mpirank_world;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank_world);
  MPI_Comm_size(MPI_COMM_WORLD, &mpisize_world);

  bool quietmode = args.quietmode;
  if (mpirank_world == 0 && !quietmode) printf("Running on %d cores.\n", mpisize_world);

  MPILogger logger(mpirank_world, quietmode);
  std::string config_file = args.config_filename;
  Config config = Config::fromFile(config_file, logger);
  std::stringstream config_log;
  config.printConfig(config_log);

  /* Initialize random number generator: Check if rand_seed is provided from config file, otherwise set random. */
  int rand_seed = config.getRandSeed();
  MPI_Bcast(&rand_seed, 1, MPI_INT, 0, MPI_COMM_WORLD); // Broadcast from rank 0 to all.
  std::mt19937 rand_engine{}; // Use Mersenne Twister for cross-platform reproducibility
  rand_engine.seed(rand_seed);

  /* Get type and the total number of initial conditions */
  int ninit = config.getNInitialConditions();

  /* --- Split communicators for distributed initial conditions, distributed linear algebra, parallel optimization --- */
  int mpirank_init, mpisize_init;
  int mpirank_optim, mpisize_optim;
  int mpirank_petsc, mpisize_petsc;
  MPI_Comm comm_optim, comm_init, comm_petsc;

  /* Get the size of communicators  */
  // Number of cores for optimization. Under development, set to 1 for now. 
  // int np_optim= config.GetIntParam("np_optim", 1);
  // np_optim= min(np_optim, mpisize_world); 
  // int np_optim= 1;
  // Number of cores for initial condition distribution. Since this gives perfect speedup, choose maximum.
  int np_init = std::min(ninit, mpisize_world); 
  // Number of cores for Petsc: All the remaining ones. 
  // int np_petsc = mpisize_world / (np_init * np_optim);
  int np_petsc = 1;
  int np_optim = mpisize_world / (np_init * np_petsc);

  /* Sanity check for communicator sizes */ 
  if (mpisize_world % ninit != 0 && ninit % mpisize_world != 0) {
    if (mpirank_world == 0) printf("ERROR: Number of threads (%d) must be integer multiplier or divisor of the number of initial conditions (%d)!\n", mpisize_world, ninit);
    exit(1);
  }

  /* Split communicators */
  // Distributed initial conditions 
  int color_init = mpirank_world % (np_petsc * np_optim);
  MPI_Comm_split(MPI_COMM_WORLD, color_init, mpirank_world, &comm_init);
  MPI_Comm_rank(comm_init, &mpirank_init);
  MPI_Comm_size(comm_init, &mpisize_init);

  // Time-parallel Optimization
  int color_optim = mpirank_world % np_petsc + mpirank_init * np_petsc;
  MPI_Comm_split(MPI_COMM_WORLD, color_optim, mpirank_world, &comm_optim);
  MPI_Comm_rank(comm_optim, &mpirank_optim);
  MPI_Comm_size(comm_optim, &mpisize_optim);

  // Distributed Linear algebra: Petsc
  int color_petsc = mpirank_world / np_petsc;
  MPI_Comm_split(MPI_COMM_WORLD, color_petsc, mpirank_world, &comm_petsc);
  MPI_Comm_rank(comm_petsc, &mpirank_petsc);
  MPI_Comm_size(comm_petsc, &mpisize_petsc);

  /* Set Petsc using petsc's communicator */
  PETSC_COMM_WORLD = comm_petsc;

  if (mpirank_world == 0 && !quietmode)  std::cout<< "Parallel distribution: " << mpisize_init << " np_init  X  " << mpisize_petsc<< " np_petsc  X " << mpisize_optim << " np_optim" << std::endl;

  char** petsc_argv = args.petsc_argv.data();
#ifdef WITH_SLEPC
  ierr = SlepcInitialize(&args.petsc_argc, &petsc_argv, (char*)0, NULL);if (ierr) return ierr;
#else
  ierr = PetscInitialize(&args.petsc_argc, &petsc_argv, (char*)0, NULL);if (ierr) return ierr;
#endif
  PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD, 	PETSC_VIEWER_ASCII_MATLAB );

  size_t num_osc = config.getNumOsc();

  /* --- Initialize the Oscillators --- */
  Oscillator** oscil_vec = new Oscillator*[num_osc];
  int param_offset = 0;
  for (size_t i = 0; i < num_osc; i++){
    oscil_vec[i] = new Oscillator(config, i, rand_engine, param_offset, quietmode);
    param_offset += oscil_vec[i]->getNParams();
  }


  /* --- Initialize the Master Equation  --- */
  // Sanity check for matrix free solver
  if (config.getUseMatFree() && mpisize_petsc > 1) {
    if (mpirank_world == 0) printf("ERROR: No Petsc-parallel version for the matrix free solver available!");
    exit(1);
  }

  MasterEq* mastereq = new MasterEq(config, oscil_vec, quietmode);

  /* Output */
  Output* output = new Output(config, comm_petsc, comm_init, quietmode);

  /* --- Initialize the time-stepper --- */
  TimeStepperType timesteppertype = config.getTimestepperType();
  TimeStepper* timestepper = nullptr;
  int ninit_local = ninit / mpisize_init; 
  switch (timesteppertype) {
    case TimeStepperType::IMR:
      timestepper = new ImplMidpoint(config, mastereq, output, ninit_local);
      break;
    case TimeStepperType::IMR4:
      timestepper = new CompositionalImplMidpoint(config, mastereq, output, ninit_local, 4);
      break;
    case TimeStepperType::IMR8:
      timestepper = new CompositionalImplMidpoint(config, mastereq, output, ninit_local, 8);
      break;
    case TimeStepperType::EE:
      timestepper = new ExplEuler(config, mastereq, output, ninit_local);
      break;
    case TimeStepperType::PETSCTS:
      timestepper = new PetscTS(config, mastereq, output, ninit_local);
      break;
    default:
      logger.exitWithError("Unknown timestepper type\n");
  }

  // Some screen output 
  if (mpirank_world == 0 && !quietmode) {
    std::cout<< "System: ";
    for (size_t i=0; i<num_osc; i++) {
      std::cout<< config.getNLevels(i);
      if (i < num_osc-1) std::cout<< "x";
    }
    std::cout<<"  (essential levels: ";
    for (size_t i=0; i<num_osc; i++) {
      std::cout<< config.getNEssential(i);
      if (i < num_osc-1) std::cout<< "x";
    }
    std::cout << ") " << std::endl;

    std::cout<<"State dimension (complex): " << mastereq->getDim() << std::endl;
    std::cout << "Time domain: [0:" << config.getTotalTime() << "]" << std::endl;
    std::cout << "Timestepping type: " << enumToString(config.getTimestepperType(), TIME_STEPPER_TYPE_MAP);
    if (config.getTimestepperType() != TimeStepperType::PETSCTS)
      std::cout << ", N="<< config.getNTime()<< ", dt=" << config.getDt();
    std::cout << std::endl;
  }

  /* --- Initialize optimization --- */
  // Create optimization target
  OptimTarget* optim_target = new OptimTarget(config, mastereq, quietmode);
  timestepper->setOptimTarget(optim_target); // Pass pointer to optimization target to timestepper for objective function evaluation.

  // Create optimization problem context 
  OptimProblem* optimctx = new OptimProblem(config, optim_target, timestepper, mastereq, comm_init, comm_optim, output, quietmode);

  /* Set upt solution and gradient vector */
  Vec xinit;
  VecCreateSeq(PETSC_COMM_SELF, optimctx->getNdesign(), &xinit);
  VecSetFromOptions(xinit);
  Vec grad;
  VecCreateSeq(PETSC_COMM_SELF, optimctx->getNdesign(), &grad);
  VecSetUp(grad);
  VecZeroEntries(grad);
  Vec opt;

  /* Some output */
  if (mpirank_world == 0)
  {
    /* Print parameters to file */
    snprintf(filename, 254, "%s/config_log.toml", output->output_dir.c_str());
    std::ofstream logfile(filename);
    if (logfile.is_open()){
      logfile << config_log.str();
      logfile.close();
      if (!quietmode) printf("File written: %s\n", filename);
    }
    else std::cerr << "Unable to open " << filename;
  }

  /* Start timer */
  double StartTime = MPI_Wtime();
  double objective;
  double gnorm = 0.0;
  /* --- Solve primal --- */
  if (config.getRuntype() == RunType::SIMULATION) {
    optimctx->getStartingPoint(xinit);
    output->writeControlParams(xinit); // Write params to file

    if (mpirank_world == 0 && !quietmode) printf("\nStarting primal solver... \n");
    bool writeTrajectoryDataFiles = true;
    objective = optimctx->evalF(xinit, writeTrajectoryDataFiles);
    if (mpirank_world == 0 && !quietmode) printf("\nTotal objective = %1.14e, \n", objective);
    optimctx->getSolution(&opt);

    // Write control pulses to file
    output->writeControls(xinit, mastereq, config.getTotalTime(), config.getDt(), timestepper->getMinTimestepSize()); // Write the control pulses 
  } 

  /* Test Gauss-Newton eigenvalue computation */
  if (config.getRuntype() == RunType::GAUSSNEWTON_EVALS) {
    if (mpirank_world == 0 && !quietmode) printf("\nStarting Gauss-Newton eigenvalue computation...\n");
    optimctx->getStartingPoint(xinit);

    Mat evecs;
    std::vector<double> evals = optimctx->computeGaussNewtonNormalEqEvals(xinit, &evecs);
    // Print the eigenvalues to file:
    snprintf(filename, 254, "%s/eigenvalues.dat", output->output_dir.c_str());
    std::ofstream evalfile(filename);
    if (evalfile.is_open()){
      for (int i=0; i<evals.size(); i++) {
        evalfile << evals[i] << "\n";
      }
      evalfile.close();
      if (mpirank_world == 0 && !quietmode) printf("Eigenvalues written to file: %s\n", filename);
    }
    else std::cerr << "Unable to open " << filename;

    // Print all eigenvectors to file
    snprintf(filename, 254, "%s/eigenvectors.dat", output->output_dir.c_str());
    std::ofstream evecfile(filename);
    if (evecfile.is_open()){
      PetscInt nrows, ncols;
      MatGetSize(evecs, &nrows, &ncols);
      for (PetscInt i = 0; i < nrows; i++) {
        for (PetscInt j = 0; j < ncols; j++) {
          PetscScalar val;
          MatGetValues(evecs, 1, &i, 1, &j, &val);
          evecfile << PetscRealPart(val) << " ";
        }
        evecfile << "\n";
      }
      evecfile.close();
      if (mpirank_world == 0 && !quietmode) printf("Eigenvectors written to file: %s\n", filename);
    }
    else std::cerr << "Unable to open " << filename;

    // // Test if evecs are orthonormal
    // Mat evecsT;
    // MatTranspose(evecs, MAT_INITIAL_MATRIX, &evecsT);
    // Mat product;
    // MatMatMult(evecsT, evecs, MAT_INITIAL_MATRIX, 1.0, &product);
    // // Check if product is approximately the identity matrix
    // PetscInt nrows, ncols;
    // MatGetSize(product, &nrows, &ncols);
    // for (PetscInt i = 0; i < nrows; i++) {
    //   for (PetscInt j = 0; j < ncols; j++) {
    //     PetscScalar val;
    //     MatGetValues(product, 1, &i, 1, &j, &val);
    //     if (mpirank_world == 0) printf("product(%d,%d) = %1.14e\n", i, j, PetscRealPart(val));
    //   }
    // }
    // MatDestroy(&evecsT);
    // MatDestroy(&product);
    MatDestroy(&evecs);

  }

  /* Test Gauss-Newton linear system solve */
  if (config.getRuntype() == RunType::GAUSSNEWTON_LS) {
    if (mpirank_world == 0 && !quietmode) {
      printf("\nStarting Gauss-Newton solves ...\n");
    }
    optimctx->getStartingPoint(xinit);
    // Example: Override max iterations for all Gauss-Newton solvers at once
    // (All three solvers use gn_maxiter from config by default)
    // optimctx->setGaussNewtonMaxiter(200);

    // Do one gradient evaluation first to store the forward states and get the right hand side
    bool writeTrajectoryDataFiles = true;
    optimctx->evalGradF(xinit, grad, writeTrajectoryDataFiles);

    Vec v_LeastSquares, v_result, v_zero;
    VecDuplicate(xinit, &v_LeastSquares);
    VecDuplicate(xinit, &v_zero);
    VecDuplicate(xinit, &v_result); // for testing
    VecZeroEntries(v_zero); // Zero initial guess

    if (mpirank_world == 0 && !quietmode) printf("\nCalling solveGaussNewtonLeastSquares\n");
    optimctx->solveGaussNewtonLeastSquares(xinit, v_zero, v_LeastSquares);
    double v_LeastSquares_norm;
    VecNorm(v_LeastSquares, NORM_2, &v_LeastSquares_norm);
    if (mpirank_world == 0 && !quietmode) {
      printf("Norm of LeastSquares solution: %1.14e\n", v_LeastSquares_norm);
    }

    // call the Gauss-Newton solver again to evaluate the residual
    // Save current max iterations, set to 1 for residual check, then restore
    int saved_maxiter = optimctx->getGaussNewtonMaxiter();
    optimctx->setGaussNewtonMaxiter(1); // just for checking the residual
    if (mpirank_world == 0 && !quietmode) printf("\nTesting: call solveGaussNewtonLeastSquares again with converged search direction\n");
    optimctx->solveGaussNewtonLeastSquares(xinit, v_LeastSquares, v_result);
    // check if the solution norm has changed
    VecNorm(v_result, NORM_2, &v_LeastSquares_norm);
    if (mpirank_world == 0 && !quietmode) {
      printf("Norm of re-evaluated LeastSquares solution: %1.14e\n", v_LeastSquares_norm);
      printf("End test\n\n");
    }
    optimctx->setGaussNewtonMaxiter(saved_maxiter); // reset max iterations
    // exit(1);

    // Linear systems: Set right hand side
    Vec gnrhs; 
    VecDuplicate(grad, &gnrhs); 
    VecCopy(grad, gnrhs);
    VecScale(gnrhs, -1.0);
    
    // Solve Gauss-Newton Normal Equation with KSP
    Vec v_KSP;
    VecDuplicate(grad, &v_KSP);
    if (mpirank_world == 0 && !quietmode) printf("\nCalling solveGaussNewtonNormalEqKSP\n");
    optimctx->solveGaussNewtonNormalEqKSP(xinit, v_zero, gnrhs, v_KSP);

    // call the Gauss-Newton solver again to re-evaluate the residual
    optimctx->setGaussNewtonMaxiter(1); // just for checking the residual
    if (mpirank_world == 0 && !quietmode) printf("\nTesting: call solveGaussNewtonNormalEqKSP again with converged search direction\n");
    optimctx->solveGaussNewtonNormalEqKSP(xinit, v_KSP, gnrhs, v_result);
    // check if the solution norm has changed
    VecNorm(v_result, NORM_2, &v_LeastSquares_norm);
    if (mpirank_world == 0 && !quietmode) {
      printf("Norm of re-evaluated NormalEq solution: %1.14e\n", v_LeastSquares_norm);
    }

    // check if the solution from the least squares routine gives a small residual in KSP and vv?
    if (mpirank_world == 0 && !quietmode){
      printf("\n");
      printf("Testing: call solveGaussNewtonLeastSquares with converged search direction from solveGaussNewtonNormalEqKSP\n");
    }
    optimctx->solveGaussNewtonLeastSquares(xinit, v_KSP, v_result);
  
    // check if the solution from the least squares routine gives a small residual?
    if (mpirank_world == 0 && !quietmode){
        printf("\n");
        printf("Testing: call solveGaussNewtonNormalEqKSP with converged search direction from solveGaussNewtonLeastSquares\n");
    }
    optimctx->solveGaussNewtonNormalEqKSP(xinit, v_LeastSquares, gnrhs, v_result);

    // Evaluate residuals from Least Squares solution
    Vec v_test_1;
    VecDuplicate(grad, &v_test_1);
    MatMult(optimctx->getGN_NormalEq_MatShell(), v_LeastSquares, v_test_1);
    VecAXPY(v_test_1, -1.0, gnrhs);
    VecNorm(v_test_1, NORM_2, &v_LeastSquares_norm);
    if (mpirank_world == 0 && !quietmode) {
      printf("Norm of Normal equation residual on the LeastSquares solution: %1.14e\n", v_LeastSquares_norm);
      printf("End test\n\n");
    }

    optimctx->setGaussNewtonMaxiter(saved_maxiter); // reset max iterations


    // Solve Gauss-Newton via SVD
    Vec v_EPS;
    VecDuplicate(grad, &v_EPS);
    optimctx->solveGaussNewtonNormalEqEPS(xinit, gnrhs, v_EPS);

    // Compare the solutions from KSP and EPS and v_LeastSquares
    if (mpirank_world == 0 && !quietmode) {
      Vec diff;
      VecDuplicate(grad, &diff);
      VecCopy(v_KSP, diff);
      VecAXPY(diff, -1.0, v_EPS);
      double diff_norm;
      VecNorm(diff, NORM_2, &diff_norm);
      double vnorm;
      VecNorm(v_KSP, NORM_2, &vnorm);
      printf("\n");
      printf("GN Least-squares solver type: %s\n", config.getGnLeastSquaresSolver().c_str());
      printf("GN NormalEq solver type: %s", config.getGnNormaleqSolver().c_str());
      if (config.getGnNormaleqSolver() == "minres" && config.getGnNormaleqMinresQlp()) {
        printf(" (QLP variant)");
      }
      printf("\n");
      printf("GN NormalEq damping: %1.4e\n", config.getGnNormaleqDamping());
      printf("\n");
      printf("Relative difference norm between KSP and EPS solutions: %1.14e (absolute: %1.14e)\n", diff_norm/vnorm, diff_norm);
      VecCopy(v_KSP, diff);
      VecAXPY(diff, -1.0, v_LeastSquares);
      VecNorm(diff, NORM_2, &diff_norm);
      VecNorm(v_KSP, NORM_2, &vnorm);
      printf("Relative difference norm between KSP and LeastSquares solutions: %1.14e (absolute: %1.14e)\n", diff_norm/vnorm, diff_norm);
      VecCopy(v_EPS, diff);
      VecAXPY(diff, -1.0, v_LeastSquares);
      VecNorm(diff, NORM_2, &diff_norm);
      VecNorm(v_EPS, NORM_2, &vnorm);
      printf("Relative difference norm between EPS and LeastSquares solutions: %1.14e (absolute: %1.14e)\n", diff_norm/vnorm, diff_norm);
      VecDestroy(&diff);
    }
    
    // Check if v_KSP is a descent direction
    double dot_ksp, dot_eps, dot_ls;
    VecDot(grad, v_KSP, &dot_ksp);
    VecDot(grad, v_EPS, &dot_eps);
    VecDot(grad, v_LeastSquares, &dot_ls);
    if (mpirank_world == 0 && !quietmode) {
      printf("\n Dot product of gradient and KSP solution (should be negative for descent): %1.14e\n", dot_ksp);
      printf(" Dot product of gradient and EPS solution (should be negative for descent): %1.14e\n", dot_eps);
      printf(" Dot product of gradient and LeastSquares solution (should be negative for descent): %1.14e\n", dot_ls);
    }

    VecDestroy(&v_KSP);
    VecDestroy(&v_EPS);
    VecDestroy(&v_LeastSquares);
    VecDestroy(&v_result);
    VecDestroy(&v_zero);
    VecDestroy(&gnrhs);
  }

  /* --- Solve adjoint --- */
  if (config.getRuntype() == RunType::GRADIENT) {
    optimctx->getStartingPoint(xinit);
    output->writeControlParams(xinit); // Write params to file

    if (mpirank_world == 0 && !quietmode) printf("\nStarting adjoint solver...\n");
    bool writeTrajectoryDataFiles=true;
    optimctx->evalGradF(xinit, grad, writeTrajectoryDataFiles);
    VecNorm(grad, NORM_2, &gnorm);
    // VecView(grad, PETSC_VIEWER_STDOUT_WORLD);
    if (mpirank_world == 0 && !quietmode) {
      printf("\nGradient norm: %1.14e\n", gnorm);
    }
    output->writeGradient(grad);

    // Write control pulses to file
    output->writeControls(xinit, mastereq, config.getTotalTime(), config.getDt(), timestepper->getMinTimestepSize()); // Write the control pulses 
  }

  /* --- Solve the optimization  --- */
  if (config.getRuntype() == RunType::OPTIMIZATION) {
    /* Set initial starting point */
    optimctx->getStartingPoint(xinit);
    output->writeControlParams(xinit); // Write params to file

    if (mpirank_world == 0 && !quietmode) printf("\nStarting Optimization solver ... \n");
    optimctx->solve(xinit);
    optimctx->getSolution(&opt);

    // Write control and parameters to file. 
    output->writeControlParams(opt);
    output->writeControls(opt, mastereq, config.getTotalTime(), config.getDt(), timestepper->getMinTimestepSize());

    // Do one last forward evaluation while writing trajectory files
    optimctx->evalF(opt, true); 
  }
  
  /* Only evaluate and write control pulses (no propagation) */
  if (config.getRuntype() == RunType::EVALCONTROLS) {
    std::vector<double> pt, qt;
    if (mpirank_world == 0 && !quietmode) printf("\nEvaluating current controls ... \n");
    optimctx->getStartingPoint(xinit);
    output->writeControlParams(xinit); // Write params to file
    output->writeControls(xinit, mastereq, config.getTotalTime(), config.getDt(), timestepper->getMinTimestepSize()); // Write the control pulses 
  }

  /* Output */
  if (config.getRuntype() != RunType::OPTIMIZATION) {
    output->writeOptimFile(0, optimctx->getObjective(), gnorm, 0.0, optimctx->getFidelity(), optimctx->getCostT(), optimctx->getRegul(), optimctx->getPenaltyLeakage(), optimctx->getPenaltyDpDm(), optimctx->getPenaltyEnergy(), optimctx->getPenaltyVariation(), optimctx->getPenaltyWeightedCost(), 0);
  }

  /* --- Finalize --- */

  /* Get timings */
  // #ifdef WITH_MPI
  double UsedTime = MPI_Wtime() - StartTime;
  // #else
  // double UsedTime = 0.0; // TODO
  // #endif
  /* Get memory usage */
  struct rusage r_usage;
  getrusage(RUSAGE_SELF, &r_usage);
  double myMB;
  #ifdef __APPLE__
      // On macOS, ru_maxrss is in bytes
      myMB = (double)r_usage.ru_maxrss / (1024.0 * 1024.0);
  #else
      // On Linux, ru_maxrss is in kilobytes
      myMB = (double)r_usage.ru_maxrss / 1024.0;
  #endif
  double globalMB = myMB;
  MPI_Allreduce(&myMB, &globalMB, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  /* Print statistics */
  if (mpirank_world == 0 && !quietmode) {
    printf("\n");
    printf(" Used Time:        %.2f seconds\n", UsedTime);
    printf(" Processors used:  %d\n", mpisize_world);
    printf(" Global Memory:    %.2f MB    [~ %.2f MB per proc]\n", globalMB, globalMB / mpisize_world);
    printf("\n");
  }
  // printf("Rank %d: %.2fMB\n", mpirank_world, myMB );

  /* Print timing to file */
  if (mpirank_world == 0) {
    snprintf(filename, 254, "%s/timing.dat", output->output_dir.c_str());
    FILE* timefile = fopen(filename, "w");
    fprintf(timefile, "%d  %1.8e\n", mpisize_world, UsedTime);
    fclose(timefile);
  }


#if TEST_FD_GRAD
  if (mpirank_world == 0)  {
    printf("\n\n#########################\n");
    printf(" FD Testing for Gradient ... \n");
    printf("#########################\n\n");
  }

  if (config.getTimestepperType() == TimeStepperType::PETSCTS) {
    if (mpirank_world == 0) printf("WARNING: Finite Difference test with PETSc's adaptive timestepper gives weird results when EPS gets small! Better to switch to TSAdapt=NONE for finite differences testing.\n");
  }

  double obj_org;
  double obj_pert1, obj_pert2;

  optimctx->getStartingPoint(xinit);
  output->writeControlParams(xinit); // Write params to file

  /* --- Solve primal --- */
  if (mpirank_world == 0) printf("\nRunning optimizer eval_f... ");
  obj_org = optimctx->evalF(xinit);
  if (mpirank_world == 0) printf(" Obj_orig %1.14e\n", obj_org);

  /* --- Solve adjoint --- */
  if (mpirank_world == 0) printf("\nRunning optimizer eval_grad_f...\n");
  optimctx->evalGradF(xinit, grad);
  // VecView(grad, PETSC_VIEWER_STDOUT_WORLD);
  

  /* --- Finite Differences --- */
  if (mpirank_world == 0) printf("\nFinite Difference testing...\n");
  double max_err = 0.0;
  double max_abs_err = 0.0;
  for (PetscInt i=0; i<optimctx->getNdesign(); i++){
  // {int i=0;

    double xi = 0.0;
    VecGetValues(xinit, 1, &i, &xi);
    const double eps_i = EPS * std::max(1.0, std::abs(xi));

    /* Evaluate f(p+eps)*/
    VecSetValue(xinit, i, eps_i, ADD_VALUES);
    VecAssemblyBegin(xinit); VecAssemblyEnd(xinit);
    obj_pert1 = optimctx->evalF(xinit);

    /* Evaluate f(p-eps)*/
    VecSetValue(xinit, i, -2*eps_i, ADD_VALUES);
    VecAssemblyBegin(xinit); VecAssemblyEnd(xinit);
    obj_pert2 = optimctx->evalF(xinit);

    /* Eval FD and error */
    double fd = (obj_pert1 - obj_pert2) / (2.*eps_i);
    double gradi; 
    VecGetValues(grad, 1, &i, &gradi);
    const double abs_err = std::abs(gradi - fd);
    const double rel_denom = std::max({1.0, std::abs(fd), std::abs(gradi)});
    const double err = abs_err / rel_denom;
    if (mpirank_world == 0) printf(" %d: eps_i %1.14e, obj %1.14e, obj_pert1 %1.14e, obj_pert2 %1.14e, fd %1.14e, grad %1.14e, abs_err %1.14e, rel_err %1.14e\n", i, eps_i, obj_org, obj_pert1, obj_pert2, fd, gradi, abs_err, err);
    if (abs(err) > max_err) max_err = err;
    if (abs_err > max_abs_err) max_abs_err = abs_err;

    /* Restore parameter */
    VecSetValue(xinit, i, eps_i, ADD_VALUES);
    VecAssemblyBegin(xinit); VecAssemblyEnd(xinit);
  }

  printf("\nMax. Finite Difference relative error: %1.14e\n", max_err);
  printf("Max. Finite Difference absolute error: %1.14e\n\n", max_abs_err);
  
#endif

#if TEST_FD_LINEARIZED_FWD
  if (mpirank_world == 0)  {
    printf("\n\n#########################\n");
    printf(" FD Testing for linearized forward solve... \n");
    printf("#########################\n\n");
  }

  // Point of evaluation
  optimctx->getStartingPoint(xinit);
  output->writeControlParams(xinit); // Write params to file

  // one forward just to get a state of correct dimension
  optimctx->evalF(xinit);
  Vec state;
  VecDuplicate(timestepper->getFinalState(0), &state);
  VecAssemblyBegin(state); VecAssemblyEnd(state);
 
  // Create storate
  Vec FD_approx, FD_err;
  VecDuplicate(state, &FD_approx);
  VecDuplicate(state, &FD_err);
  std::vector<Vec> states_plus (ninit_local);
  std::vector<Vec> states_minus (ninit_local);
  std::vector<Vec> linearized_state (ninit_local);
  for (int i =0; i<ninit_local; i++){
    VecDuplicate(state, &states_plus[i]);
    VecDuplicate(state, &states_minus[i]);
    VecDuplicate(state, &linearized_state[i]);
    VecAssemblyBegin(states_plus[i]); VecAssemblyEnd(states_plus[i]);
    VecAssemblyBegin(states_minus[i]); VecAssemblyEnd(states_minus[i]);
    VecAssemblyBegin(linearized_state[i]); VecAssemblyEnd(linearized_state[i]);
    VecZeroEntries(states_plus[i]);
    VecZeroEntries(states_minus[i]);
    VecZeroEntries(linearized_state[i]);
  }
  Vec v;
  VecDuplicate(xinit, &v);
  VecAssemblyBegin(v); VecAssemblyEnd(v);

  /* --- Finite Differences --- */
  double max_abs_err = 0.0;
  double abs_err = 0.0;
  double rel_err = 0.0;

  for (PetscInt ix=0; ix<optimctx->getNdesign(); ix++){
  // PetscInt i=5; {
    double xi = 0.0;
    VecGetValues(xinit, 1, &ix, &xi);
    const double eps_i = EPS * std::max(1.0, std::abs(xi));
    printf("Testing finite difference for parameter index %d with eps_i = %e\n", ix, eps_i);

    // Set linearization direction to i-th unit vector
    VecZeroEntries(v);
    VecSetValue(v, ix, 1.0, ADD_VALUES);
    VecAssemblyBegin(v); VecAssemblyEnd(v);

    // Get linearized forward results
    optimctx->evalLinearizedForward(xinit, v);
    for (int i =0; i<ninit_local; i++){
      VecCopy(timestepper->getLinearizedState(i, config.getNTime()), linearized_state[i]);
    }

    /* Evaluate perturbed state U(p+eps)*/
    VecSetValue(xinit, ix, eps_i, ADD_VALUES);
    VecAssemblyBegin(xinit); VecAssemblyEnd(xinit);
    optimctx->evalF(xinit);
    for (int i =0; i<ninit_local; i++){
      VecCopy(timestepper->getFinalState(i), states_plus[i]);
    }

    /* Evaluate U(p-eps)*/
    VecSetValue(xinit, ix, -2*eps_i, ADD_VALUES);
    VecAssemblyBegin(xinit); VecAssemblyEnd(xinit);
    optimctx->evalF(xinit);
    for (int i =0; i<ninit_local; i++){
      VecCopy(timestepper->getFinalState(i), states_minus[i]);
    }

    /* Restore original parameters xinit */
    VecSetValue(xinit, ix, eps_i, ADD_VALUES);
    VecAssemblyBegin(xinit); VecAssemblyEnd(xinit);

    /* Evaluate finite difference and error */
    // FD : dU/dalpha_k = 1/(2EPS)* (Uplus - Uminus)
    // error = norm(DU_FD - DU_exact)
    for (int iinit=0; iinit<ninit_local; iinit++){

      // FD_approx = (U(p+eps) - U(p-eps)) / 2eps
      VecCopy(states_plus[iinit], FD_approx);
      VecAXPY(FD_approx, -1.0, states_minus[iinit]);
      VecScale(FD_approx, 1./(2.*eps_i)); 

      // FD_error = states_minus = Exact - FDapprox = exact - states_plus
      VecCopy(linearized_state[iinit], FD_err);
      VecAXPY(FD_err, -1.0, FD_approx);

      // error: norm(DU_exact - FD_approx) 
      VecNorm(FD_err, NORM_2, &abs_err);

      if (mpirank_world == 0)
      printf(" %d: %d/%d iinit %d abs_err %1.14e \n", mpirank_world, ix, optimctx->getNdesign(), iinit, abs_err);

      max_abs_err = std::max(abs_err, max_abs_err);
    }
  }

  printf("\n Max. absolute error = %1.14e\n", max_abs_err);

  // Cleanup
  VecDestroy(&v);
  VecDestroy(&FD_approx);
  VecDestroy(&FD_err);
  for (int i =0; i<ninit_local; i++){
    VecDestroy(&states_plus[i]);
    VecDestroy(&states_minus[i]);
    VecDestroy(&linearized_state[i]);
  }

#endif

#if TEST_GAUSSNEWTON_LEASTSQUARES
  optimctx->getStartingPoint(xinit);
  output->writeControlParams(xinit); // Write params to file

  // One forward and backward first to populate final_states.
  bool writeTrajectoryDataFiles = true;
  optimctx->evalGradF(xinit, grad, writeTrajectoryDataFiles);

  // Set point of evaluation for Gauss-Newton matrix
  optimctx->setXevalGN(xinit);

  Vec v, Av;
  int ndesign = optimctx->getNdesign();

  // Create vectors with correct dimensions using MatCreateVecs
  MatCreateVecs(optimctx->getGN_LeastSquares_MatShell(), &v, &Av);

  // Get the dimension of the output vector (nested vector with ninit subvectors)
  PetscInt nrows;
  VecGetSize(Av, &nrows);

  // Create matrix to store all columns: nrows x ndesign
  Mat A_LeastSquares;
  MatCreate(PETSC_COMM_SELF, &A_LeastSquares);
  MatSetSizes(A_LeastSquares, PETSC_DECIDE, PETSC_DECIDE, nrows, ndesign);
  MatSetType(A_LeastSquares, MATDENSE);
  MatSetUp(A_LeastSquares);
  MatZeroEntries(A_LeastSquares);

  // Apply GN_LeastSquares_MatShell operator to each unit vector, storing the resulting matrix
  for (int ix = 0; ix < ndesign; ix++){
    VecZeroEntries(v);
    VecSetValue(v, ix, 1.0, INSERT_VALUES);
    VecAssemblyBegin(v); VecAssemblyEnd(v);

    MatMult(optimctx->getGN_LeastSquares_MatShell(), v, Av);

    // Store Av in ix-th column of A_LeastSquares
    const PetscScalar *Av_ptr;
    VecGetArrayRead(Av, &Av_ptr);
    for (PetscInt row = 0; row < nrows; row++){
      MatSetValue(A_LeastSquares, row, ix, Av_ptr[row], INSERT_VALUES);
    }
    VecRestoreArrayRead(Av, &Av_ptr);
  }
  MatAssemblyBegin(A_LeastSquares, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A_LeastSquares, MAT_FINAL_ASSEMBLY);

  // Write the full matrix to file in Python-friendly format
  if (mpirank_world == 0) {
    snprintf(filename, 254, "%s/gaussnewton_leastsquares_matrix.dat", output->output_dir.c_str());
    FILE* matfile = fopen(filename, "w");
    if (matfile) {
      // Write header comment with dimensions
      fprintf(matfile, "# Gauss-Newton Least Squares matrix (Jacobian L)\n");
      fprintf(matfile, "# Dimensions: %d x %d (rows x cols)\n", (int)nrows, ndesign);

      // Write matrix row by row
      for (PetscInt i = 0; i < nrows; i++) {
        for (int j = 0; j < ndesign; j++) {
          double Aij;
          MatGetValue(A_LeastSquares, i, j, &Aij);
          fprintf(matfile, "%.16e", Aij);
          if (j < ndesign - 1) {
            fprintf(matfile, " ");
          }
        }
        fprintf(matfile, "\n");
      }
      fclose(matfile);
      printf("File written: %s\n", filename);
    } else {
      printf("ERROR: Could not open file %s for writing\n", filename);
    }
  }

  // Clean up
  VecDestroy(&v);
  VecDestroy(&Av);
  MatDestroy(&A_LeastSquares);

#endif

#if TEST_GAUSSNEWTON_LINEARSYSTEM
  /*  ---- TEST: Evaluate GaussNewton matrix columns ---- */
  optimctx->getStartingPoint(xinit);
  output->writeControlParams(xinit); // Write params to file

  // One forward and backward first to populate final_states.
  bool writeTrajectoryDataFiles_GN = true;
  optimctx->evalGradF(xinit, grad, writeTrajectoryDataFiles_GN);

  // Set point of evaluation for Gauss-Newton matrix
  optimctx->setXevalGN(xinit);

  Vec v_GN, Av_GN;
  VecDuplicate(xinit, &v_GN);
  VecDuplicate(xinit, &Av_GN);
  VecZeroEntries(v_GN);
  VecZeroEntries(Av_GN);
  int ndesign_GN = optimctx->getNdesign();

  // Storage for Gauss-Newton matrix A = L^T * L
  Mat A_GaussNewton;
  MatCreate(PETSC_COMM_SELF, &A_GaussNewton);
  MatSetSizes(A_GaussNewton, PETSC_DECIDE, PETSC_DECIDE, ndesign_GN, ndesign_GN);
  MatSetType(A_GaussNewton, MATDENSE);
  MatSetUp(A_GaussNewton);
  MatZeroEntries(A_GaussNewton);

  // Number of local matrix columns for this processor
  int ncols_local = ndesign_GN / mpisize_optim;
  // If not integer divisible, let the last processor handle the remainder
  if (mpirank_optim == mpisize_optim - 1) {
    ncols_local = ndesign_GN - mpirank_optim * ncols_local;
  }

  // Iterate over local columns of the Gauss-Newton matrix
  for (int ix_local = 0; ix_local < ncols_local; ix_local++) {
    int ix = mpirank_optim * ncols_local + ix_local;

    // Set v_GN to the ix-th unit vector
    VecZeroEntries(v_GN);
    VecSetValue(v_GN, ix, 1.0, INSERT_VALUES);
    VecAssemblyBegin(v_GN); VecAssemblyEnd(v_GN);

    // Evaluate Av_GN = A * e_ix
    MatMult(optimctx->getGN_NormalEq_MatShell(), v_GN, Av_GN);

    // Store Av_GN in ix-th column of A_GaussNewton
    const PetscScalar *Av_ptr;
    VecGetArrayRead(Av_GN, &Av_ptr);
    for (int row = 0; row < ndesign_GN; row++){
      MatSetValue(A_GaussNewton, row, ix, Av_ptr[row], INSERT_VALUES);
    }
    VecRestoreArrayRead(Av_GN, &Av_ptr);
    MatAssemblyBegin(A_GaussNewton, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A_GaussNewton, MAT_FINAL_ASSEMBLY);
  }

  // Need to sum up the columns of A_GaussNewton from all optim_comm processors
  PetscScalar *A_data;
  MatDenseGetArray(A_GaussNewton, &A_data);
  int size = ndesign_GN * ndesign_GN;
  MPI_Allreduce(MPI_IN_PLACE, A_data, size, MPIU_SCALAR, MPI_SUM, comm_optim);
  MatDenseRestoreArray(A_GaussNewton, &A_data);

  // Write the full matrix to file in Python-friendly format
  if (mpirank_world == 0) {
    snprintf(filename, 254, "%s/gaussnewton_linearsystem_matrix.dat", output->output_dir.c_str());
    FILE* matfile = fopen(filename, "w");
    if (matfile) {
      // Write header comment with dimensions
      fprintf(matfile, "# Gauss-Newton matrix A = L^T * L\n");
      fprintf(matfile, "# Dimensions: %d x %d (rows x cols)\n", ndesign_GN, ndesign_GN);

      // Write matrix row by row
      for (int i = 0; i < ndesign_GN; i++) {
        for (int j = 0; j < ndesign_GN; j++) {
          double Aij;
          MatGetValue(A_GaussNewton, i, j, &Aij);
          fprintf(matfile, "%.16e", Aij);
          if (j < ndesign_GN - 1) {
            fprintf(matfile, " ");
          }
        }
        fprintf(matfile, "\n");
      }
      fclose(matfile);
      printf("File written: %s\n", filename);
    } else {
      printf("ERROR: Could not open file %s for writing\n", filename);
    }
  }

  // Clean up
  VecDestroy(&v_GN);
  VecDestroy(&Av_GN);
  MatDestroy(&A_GaussNewton);

#endif


#if TEST_FD_HESS
  if (mpirank_world == 0)  {
    printf("\n\n#########################\n");
    printf(" FD Testing for Hessian... \n");
    printf("#########################\n\n");
  }
  optimctx->getStartingPoint(xinit);
  output->writeControlParams(xinit); // Write params to file

  /* Figure out which parameters are hitting bounds */
  double bound_tol = 1e-3;
  std::vector<int> Ihess; // Index set for all elements that do NOT hit a bound
  for (PetscInt i=0; i<optimctx->getNdesign(); i++){
    // get x_i and bounds for x_i
    double xi, blower, bupper;
    VecGetValues(xinit, 1, &i, &xi);
    VecGetValues(optimctx->xlower, 1, &i, &blower);
    VecGetValues(optimctx->xupper, 1, &i, &bupper);
    // compare 
    if (fabs(xi - blower) < bound_tol || 
        fabs(xi - bupper) < bound_tol  ) {
          printf("Parameter %d hits bound: x=%f\n", i, xi);
    } else {
      Ihess.push_back(i);
    }
  }

  double grad_org;
  double grad_pert1, grad_pert2;
  Mat Hess;
  int nhess = Ihess.size();
  MatCreateSeqDense(PETSC_COMM_SELF, nhess, nhess, NULL, &Hess);
  MatSetUp(Hess);

  Vec grad1, grad2;
  VecDuplicate(grad, &grad1);
  VecDuplicate(grad, &grad2);


  /* Iterate over all params that do not hit a bound */
  for (PetscInt k=0; k< Ihess.size(); k++){
    PetscInt j = Ihess[k];
    printf("Computing column %d\n", j);

    /* Evaluate \nabla_x J(x + eps * e_j) */
    VecSetValue(xinit, j, EPS, ADD_VALUES); 
    optimctx->evalGradF(xinit, grad);        
    VecCopy(grad, grad1);

    /* Evaluate \nabla_x J(x - eps * e_j) */
    VecSetValue(xinit, j, -2.*EPS, ADD_VALUES); 
    optimctx->evalGradF(xinit, grad);
    VecCopy(grad, grad2);

    for (PetscInt l=0; l<Ihess.size(); l++){
      PetscInt i = Ihess[l];

      /* Get the derivative wrt parameter i */
      VecGetValues(grad1, 1, &i, &grad_pert1);   // \nabla_x_i J(x+eps*e_j)
      VecGetValues(grad2, 1, &i, &grad_pert2);    // \nabla_x_i J(x-eps*e_j)

      /* Finite difference for element Hess(l,k) */
      double fd = (grad_pert1 - grad_pert2) / (2.*EPS);
      MatSetValue(Hess, l, k, fd, INSERT_VALUES);
    }

    /* Restore parameters xinit */
    VecSetValue(xinit, j, EPS, ADD_VALUES);
  }
  /* Assemble the Hessian */
  MatAssemblyBegin(Hess, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(Hess, MAT_FINAL_ASSEMBLY);
  
  /* Clean up */
  VecDestroy(&grad1);
  VecDestroy(&grad2);


  /* Epsilon test: compute ||1/2(H-H^T)||_F  */
  MatScale(Hess, 0.5);
  Mat HessT, Htest;
  MatDuplicate(Hess, MAT_COPY_VALUES, &Htest);
  MatTranspose(Hess, MAT_INITIAL_MATRIX, &HessT);
  MatAXPY(Htest, -1.0, HessT, SAME_NONZERO_PATTERN);
  double fnorm;
  MatNorm(Htest, NORM_FROBENIUS, &fnorm);
  printf("EPS-test: ||1/2(H-H^T)||= %1.14e\n", fnorm);

  /* symmetrize H_symm = 1/2(H+H^T) */
  MatAXPY(Hess, 1.0, HessT, SAME_NONZERO_PATTERN);

  /* --- Print Hessian to file */
  
  snprintf(filename, 254, "%s/hessian.dat", output->output_dir.c_str());
  printf("File written: %s.\n", filename);
  PetscViewer viewer;
  PetscViewerCreate(MPI_COMM_WORLD, &viewer);
  PetscViewerSetType(viewer, PETSCVIEWERASCII);
  PetscViewerFileSetMode(viewer, FILE_MODE_WRITE);
  PetscViewerFileSetName(viewer, filename);
  // PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_DENSE);
  MatView(Hess, viewer);
  PetscViewerPopFormat(viewer);
  PetscViewerDestroy(&viewer);

  // write again in binary
  snprintf(filename, 254, "%s/hessian_bin.dat", output->output_dir.c_str());
  printf("File written: %s.\n", filename);
  PetscViewerBinaryOpen(MPI_COMM_WORLD, filename, FILE_MODE_WRITE, &viewer);
  MatView(Hess, viewer);
  PetscViewerDestroy(&viewer);

  MatDestroy(&Hess);

#endif

#if HESSIAN_DECOMPOSITION 
  /* --- Compute eigenvalues of Hessian --- */
  printf("\n\n#########################\n");
  printf(" Eigenvalue analysis... \n");
  printf("#########################\n\n");

  /* Load Hessian from file */
  Mat Hess;
  MatCreate(PETSC_COMM_SELF, &Hess);
  snprintf(filename, 254, "%s/hessian_bin.dat", output->output_dir.c_str());
  printf("Reading file: %s\n", filename);
  PetscViewer viewer;
  PetscViewerCreate(MPI_COMM_WORLD, &viewer);
  PetscViewerSetType(viewer, PETSCVIEWERBINARY);
  PetscViewerFileSetMode(viewer, FILE_MODE_READ);
  PetscViewerFileSetName(viewer, filename);
  PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_DENSE);
  MatLoad(Hess, viewer);
  PetscViewerPopFormat(viewer);
  PetscViewerDestroy(&viewer);
  int nrows, ncols;
  MatGetSize(Hess, &nrows, &ncols);


  /* Set the percentage of eigenpairs that should be computed */
  double frac = 1.0;  // 1.0 = 100%
  int neigvals = nrows * frac;     // hopefully rounds to closest int 
  printf("\nComputing %d eigenpairs now...\n", neigvals);
  
  /* Compute eigenpair */
  std::vector<double> eigvals;
  std::vector<Vec> eigvecs;
  getEigvals(Hess, neigvals, eigvals, eigvecs);

  /* Print eigenvalues to file. */
  FILE *file;
  snprintf(filename, 254, "%s/eigvals.dat", output->output_dir.c_str());
  file =fopen(filename,"w");
  for (int i=0; i<eigvals.size(); i++){
      fprintf(file, "% 1.8e\n", eigvals[i]);  
  }
  fclose(file);
  printf("File written: %s.\n", filename);

  /* Print eigenvectors to file. Columns wise */
  snprintf(filename, 254, "%s/eigvecs.dat", output->output_dir.c_str());
  file =fopen(filename,"w");
  for (PetscInt j=0; j<nrows; j++){  // rows
    for (PetscInt i=0; i<eigvals.size(); i++){
      double val;
      VecGetValues(eigvecs[i], 1, &j, &val); // j-th row of eigenvalue i
      fprintf(file, "% 1.8e  ", val);  
    }
    fprintf(file, "\n");
  }
  fclose(file);
  printf("File written: %s.\n", filename);


#endif

#ifdef SANITY_CHECK
  printf("\n\n Sanity checks have been performed. Check output for warnings and errors!\n\n");
#endif

  /* Clean up */
  for (size_t i=0; i<num_osc; i++){
    delete oscil_vec[i];
  }
  delete [] oscil_vec;
  delete mastereq;
  delete timestepper;
  delete optimctx;
  delete optim_target;
  delete output;

  VecDestroy(&xinit);
  VecDestroy(&grad);


  /* Finallize Petsc */
#ifdef WITH_SLEPC
  ierr = SlepcFinalize();
#else
  PetscOptionsSetValue(NULL, "-options_left", "no"); // Remove warning about unused options.
  ierr = PetscFinalize();
#endif


  MPI_Finalize();
  return ierr;
}
