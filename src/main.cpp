#include "timestepper.hpp"
#include "defs.hpp"
#include <string>
#include "oscillator.hpp" 
#include "mastereq.hpp"
#include "config.hpp"
#include <stdlib.h>
#include <sys/resource.h>
#include <cassert>
#include "optimproblem.hpp"
#include <filesystem>
#include "output.hpp"
#include "petsc.h"
#include <random>
#include "util.hpp"
#ifdef WITH_SLEPC
#include <slepceps.h>
#endif

#define TEST_FD_GRAD 0    // Run Finite Differences gradient test
#define TEST_FD_HESS 0    // Run Finite Differences Hessian test
#define HESSIAN_DECOMPOSITION 0 // Run eigenvalue analysis for Hessian
#define HESSIAN_COMPUTATION 0
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
  // int np_optim= 1;
  int np_petsc= 1;
  // Number of cores for initial condition distribution. Since this gives perfect speedup, choose maximum.
  int np_init = std::min(ninit, mpisize_world); 
  // Number of cores for Optim: All the remaining ones. 
  int np_optim= mpisize_world / (np_init * np_petsc);

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

  if (mpirank_world == 0 && !quietmode)  std::cout<< "Parallel distribution: " << mpisize_init << " np_init  X  " << mpisize_petsc<< " np_petsc  X " << mpisize_optim << " np_optim "  << std::endl;

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

  int ntime = config.getNTime();
  double dt = config.getDt();
  double total_time = config.getTotalTime();

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
    std::cout << "Time: [0:" << total_time << "], ";
    std::cout << "N="<< ntime << ", dt=" << dt << std::endl;
  }

  /* --- Initialize the time-stepper --- */
  LinearSolverType linsolvetype = config.getLinearSolverType();
  int linsolve_maxiter = config.getLinearSolverMaxiter();

  /* My time stepper */
  bool storeFWD = true; // TODO: dont store always.
  RunType runtype = config.getRuntype();
  if (mastereq->decoherence_type != DecoherenceType::NONE &&   
     (runtype == RunType::GRADIENT || runtype == RunType::OPTIMIZATION) ) storeFWD = true;  // if NOT Schroedinger solver and running gradient optim: store forward states. Otherwise, they will be recomputed during gradient. 

  TimeStepperType timesteppertype = config.getTimestepperType();
  TimeStepper* mytimestepper = nullptr;
  switch (timesteppertype) {
    case TimeStepperType::IMR:
      mytimestepper = new ImplMidpoint(mastereq, ntime, total_time, linsolvetype, linsolve_maxiter, output, storeFWD, ninit/mpisize_init, comm_init);
      break;
    case TimeStepperType::IMR4:
      mytimestepper = new CompositionalImplMidpoint(4, mastereq, ntime, total_time, linsolvetype, linsolve_maxiter, output, storeFWD, ninit/mpisize_init, comm_init);
      break;
    case TimeStepperType::IMR8:
      mytimestepper = new CompositionalImplMidpoint(8, mastereq, ntime, total_time, linsolvetype, linsolve_maxiter, output, storeFWD, ninit/mpisize_init, comm_init);
      break;
    case TimeStepperType::EE:
      mytimestepper = new ExplEuler(mastereq, ntime, total_time, output, storeFWD, ninit/mpisize_init, comm_init);
      break;
    default:
      logger.exitWithError("Unknown timestepper type\n");
  }

  /* --- Initialize optimization --- */
  OptimProblem* optimctx = new OptimProblem(config, mytimestepper, comm_init, comm_optim, output, quietmode);

  /* Set up initial control vector and gradient for TAO */
  Vec xinit;
  VecCreateSeq(PETSC_COMM_SELF, optimctx->getNdesign(), &xinit);
  VecSetFromOptions(xinit);
  optimctx->getStartingPoint(xinit);
  VecCopy(xinit, optimctx->xinit); // Store the initial guess
  Vec grad;
  VecCreateSeq(PETSC_COMM_SELF, optimctx->getNdesign(), &grad);
  VecSetUp(grad);
  VecZeroEntries(grad);
  // Vec opt;n

  /* Some screen output */
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
  if (runtype == RunType::SIMULATION) {
    if (mpirank_world == 0 && !quietmode) printf("\nStarting primal solver... \n");
    optimctx->timestepper->writeTrajectoryDataFiles = true;
    // if (optimsolvertype==OptimSolverType::TAO) {
    objective = optimctx->evalF(xinit);
    
    if (mpirank_world == 0 && !quietmode) printf("\nTotal objective = %1.14e, \n", objective);
  } 
  
  /* --- Solve adjoint --- */
  if (runtype == RunType::GRADIENT) {
    if (mpirank_world == 0 && !quietmode) printf("\nStarting adjoint solver...\n");
    optimctx->timestepper->writeTrajectoryDataFiles = true;

    optimctx->evalGradF(xinit, grad);
    VecNorm(grad, NORM_2, &gnorm);
    optimctx->output->writeGradient(grad);

    if (mpirank_world == 0 && !quietmode) {
      printf("\nGradient norm: %1.14e\n", gnorm);
    }

    if (config.getOptimSolverType() == OptimSolverType::TAO_HESSIAN) {
      // TEST HESSIAN FUNCTION
      optimctx->evalHessian(xinit, optimctx->Hessian, NULL);
    }
  }

  /* --- Solve the optimization  --- */
  if (runtype == RunType::OPTIMIZATION) {
    /* Set initial starting point */
    if (mpirank_world == 0 && !quietmode) printf("\nStarting Optimization solver ... \n");
    optimctx->timestepper->writeTrajectoryDataFiles = false;

    StartTime = MPI_Wtime();
    optimctx->solve(xinit);

    // // write optimal controls 
    // const myVec& ex = dynamic_cast<const myVec&>(*x); 
    // output->writeControls(ex.getVector(), mytimestepper->mastereq, mytimestepper->ntime, mytimestepper->dt);
    // // one last forward evaluation while writing trajectory data
    // mytimestepper->writeTrajectoryDataFiles = true;
    // optimctx->evalF(ex.getVector()); 
  }

  /* Only evaluate and write control pulses (no propagation) */
  if (runtype == RunType::EVALCONTROLS) {
    std::vector<double> pt, qt;
    optimctx->getStartingPoint(xinit);
    if (mpirank_world == 0 && !quietmode) printf("\nEvaluating current controls ... \n");
    output->writeControls(xinit, mastereq, ntime, dt);
  }

  /* Output */
  if (runtype != RunType::OPTIMIZATION) {
    optimctx->output->writeOptimFile(0, optimctx->getObjective(), gnorm, 0.0, optimctx->getFidelity(), optimctx->getCostT(), optimctx->getRegul(), optimctx->getPenaltyLeakage(), optimctx->getPenaltyDpDm(), optimctx->getPenaltyEnergy(), optimctx->getPenaltyVariation(), optimctx->getPenaltyWeightedCost());
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

  double obj_org;
  double obj_pert1, obj_pert2;

  optimctx->getStartingPoint(xinit);

  /* --- Solve primal --- */
  if (mpirank_world == 0) printf("\nRunning optimizer eval_f... ");
  obj_org = optimctx->evalF(xinit);
  if (mpirank_world == 0) printf(" Obj_orig %1.14e\n", obj_org);


  // Run once without testing to store all adjoint states, then set testnow to true so that the next evalGradF will compare to the stored adjoint states.
  optimctx->evalGradF(xinit, grad);
  mytimestepper->testnow = true;

  /* --- Solve adjoint --- */
  if (mpirank_world == 0) printf("\nRunning optimizer eval_grad_f...\n");
  optimctx->evalGradF(xinit, grad);
  VecView(grad, PETSC_VIEWER_STDOUT_WORLD);
  

  /* --- Finite Differences --- */
  if (mpirank_world == 0) printf("\nFinite Difference testing...\n");
  double max_err = 0.0;
  for (PetscInt i=0; i<optimctx->getNdesign(); i++){
  // {int i=0;

    /* Evaluate f(p+eps)*/
    VecSetValue(xinit, i, EPS, ADD_VALUES);
    obj_pert1 = optimctx->evalF(xinit);

    /* Evaluate f(p-eps)*/
    VecSetValue(xinit, i, -2*EPS, ADD_VALUES);
    obj_pert2 = optimctx->evalF(xinit);

    /* Eval FD and error */
    double fd = (obj_pert1 - obj_pert2) / (2.*EPS);
    double err = 0.0;
    double gradi; 
    VecGetValues(grad, 1, &i, &gradi);
    if (fd != 0.0) err = (gradi - fd) / fd;
    if (mpirank_world == 0) printf(" %d: obj %1.14e, obj_pert1 %1.14e, obj_pert2 %1.14e, fd %1.14e, grad %1.14e, err %1.14e\n", i, obj_org, obj_pert1, obj_pert2, fd, gradi, err);
    if (abs(err) > max_err) max_err = err;

    /* Restore parameter */
    VecSetValue(xinit, i, EPS, ADD_VALUES);
  }

  printf("\nMax. Finite Difference error: %1.14e\n\n", max_err);
  
#endif


#if HESSIAN_COMPUTATION

    /* Hessian computation */
    int sizex = 0;
    VecGetSize(xinit, &sizex);


    for (int i=0; i<sizex; i++) {
    // for (int i=1085; i<sizex; i++) {
      printf("Hessian for index %d\n", i);

      Vec ei, Hei;
      // Create unit vector ei
      VecCreateSeq(PETSC_COMM_SELF, sizex, &ei);
      VecSetFromOptions(ei);
      VecZeroEntries(ei);
      VecSetValue(ei, i, 1.0, INSERT_VALUES);
      VecAssemblyBegin(ei);
      VecAssemblyEnd(ei);

      VecCreateSeq(PETSC_COMM_SELF, sizex, &Hei);
      VecSetUp(Hei);
      VecZeroEntries(Hei);

      optimctx->evalHessVec(xinit, ei, Hei);

      // Output Hessian column to file 
      if (mpirank_world == 0) {
        snprintf(filename, 254, "%s/hessian_col_%04d.dat", output->datadir.c_str(), i);
        FILE* hessfile = fopen(filename, "w");
        double hessval;
        for (int j=0; j<sizex; j++) {
          VecGetValues(Hei, 1, &j, &hessval);
          fprintf(hessfile, "%1.14e\n", hessval);
        }
        fclose(hessfile);
      }
      VecDestroy(&ei);
      VecDestroy(&Hei);
    }
#endif


#if TEST_FD_HESS
  if (mpirank_world == 0)  {
    printf("\n\n#########################\n");
    printf(" FD Testing for Hessian... \n");
    printf("#########################\n\n");
  }
  optimctx->getStartingPoint(xinit);

  Vec v, hessv;
  VecDuplicate(xinit, &v);
  VecDuplicate(xinit, &hessv);
  Vec hessv_fd;
  VecDuplicate(xinit, &hessv_fd);
  Vec gplus, gminus;
  VecDuplicate(xinit, &gplus);
  VecDuplicate(xinit, &gminus);
  Vec xplus, xminus;
  VecDuplicate(xinit, &xplus);
  VecDuplicate(xinit, &xminus);

  double epsilon = 1e-7;
  // for (int ieps = 0; ieps < 8; ieps++) {
  for (int ieps = 0; ieps < 1; ieps++) {
    epsilon *= 0.1; // Decrease epsilon by factor of 10
    printf("\nTesting Hessian vector product with epsilon = %1.14e\n", epsilon);

    // Iterate over design variables testing H*e_i against finite differences
    // for (int itest=0; itest < optimctx->getNdesign(); itest++){
    { int itest = 3; // fixed column

      // Choose the direction v = e_i
      VecZeroEntries(v);
      const int i = itest;
      VecSetValue(v, i, 1.0, INSERT_ALL_VALUES);

      // Evaluate HessianVector product
      VecZeroEntries(hessv);
      optimctx->evalHessVec(xinit, v, hessv);

      // if (mpirank_init==0) {
      //   printf("Hess vec hv = \n");
      //   VecView(hessv, NULL);
      // }
      // exit(1);

      // Perturb xinit by epsilon in direction v
      VecCopy(xinit, xplus);
      VecAXPY(xplus, epsilon, v); // xplus = xinit + EPS * v
      VecCopy(xinit, xminus);
      VecAXPY(xminus, -epsilon, v); // xminus = xinit - EPS * v

      // Evaluate gradient at xplus and xminus
      optimctx->evalGradF(xplus, gplus);
      optimctx->evalGradF(xminus, gminus);

      // Compute finite differences Hessian vector product
      VecCopy(gplus, hessv_fd);
      VecAXPY(hessv_fd, -1.0, gminus); // hessv_fd = gplus - gminus
      VecScale(hessv_fd, 1.0 / (2.0 * epsilon)); // hessv_fd = (gplus - gminus) / (2*EPS)

      // Compute error 
      const double* hessv_ptr, *hessv_fd_ptr;
      VecGetArrayRead(hessv, &hessv_ptr);
      VecGetArrayRead(hessv_fd, &hessv_fd_ptr);
      for (int j=0; j<optimctx->getNdesign(); j++) {
        printf("Diff at i=%d, j=%d: Hv= %1.14e, Hv_FD=%1.14e, rel. err=%1.14e\n", i, j, hessv_ptr[j], hessv_fd_ptr[j], (hessv_ptr[j] - hessv_fd_ptr[j])/hessv_fd_ptr[j]);
      }
      VecRestoreArrayRead(hessv, &hessv_ptr);
      VecRestoreArrayRead(hessv_fd, &hessv_fd_ptr);
    }
  }

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
  delete mytimestepper;
  delete optimctx;
  delete output;

  // ROL might be doing this??
  // VecDestroy(&xinit);
  // VecDestroy(&grad);


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
