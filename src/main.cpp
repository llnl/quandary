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
  int np_optim= 1;
  // Number of cores for initial condition distribution. Since this gives perfect speedup, choose maximum.
  int np_init = std::min(ninit, mpisize_world); 
  // Number of cores for Petsc: All the remaining ones. 
  int np_petsc = mpisize_world / (np_init * np_optim);

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

  if (mpirank_world == 0 && !quietmode)  std::cout<< "Parallel distribution: " << mpisize_init << " np_init  X  " << mpisize_petsc<< " np_petsc  " << std::endl;

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
  bool storeFWD = false;
  RunType runtype = config.getRuntype();
  if (mastereq->decoherence_type != DecoherenceType::NONE &&   
     (runtype == RunType::GRADIENT || runtype == RunType::OPTIMIZATION) ) storeFWD = true;  // if NOT Schroedinger solver and running gradient optim: store forward states. Otherwise, they will be recomputed during gradient. 

  TimeStepperType timesteppertype = config.getTimestepperType();
  TimeStepper* mytimestepper = nullptr;
  switch (timesteppertype) {
    case TimeStepperType::IMR:
      mytimestepper = new ImplMidpoint(mastereq, ntime, total_time, linsolvetype, linsolve_maxiter, output, storeFWD);
      break;
    case TimeStepperType::IMR4:
      mytimestepper = new CompositionalImplMidpoint(4, mastereq, ntime, total_time, linsolvetype, linsolve_maxiter, output, storeFWD);
      break;
    case TimeStepperType::IMR8:
      mytimestepper = new CompositionalImplMidpoint(8, mastereq, ntime, total_time, linsolvetype, linsolve_maxiter, output, storeFWD);
      break;
    case TimeStepperType::EE:
      mytimestepper = new ExplEuler(mastereq, ntime, total_time, output, storeFWD);
      break;
    default:
      logger.exitWithError("Unknown timestepper type\n");
  }

  /* --- Initialize optimization --- */
  OptimProblem* optimctx = new OptimProblem(config, mytimestepper, comm_init, comm_optim, output, quietmode);

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
  if (runtype == RunType::SIMULATION) {
    optimctx->getStartingPoint(xinit);
    VecCopy(xinit, optimctx->xinit); // Store the initial guess
    if (mpirank_world == 0 && !quietmode) printf("\nStarting primal solver... \n");
    optimctx->timestepper->writeTrajectoryDataFiles = true;
    objective = optimctx->evalF(xinit);
    if (mpirank_world == 0 && !quietmode) printf("\nTotal objective = %1.14e, \n", objective);
    optimctx->getSolution(&opt);
  } 
  
  /* --- Solve adjoint --- */
  if (runtype == RunType::GRADIENT) {
    optimctx->getStartingPoint(xinit);
    VecCopy(xinit, optimctx->xinit); // Store the initial guess
    if (mpirank_world == 0 && !quietmode) printf("\nStarting adjoint solver...\n");
    optimctx->timestepper->writeTrajectoryDataFiles = true;
    optimctx->evalGradF(xinit, grad);
    VecNorm(grad, NORM_2, &gnorm);
    // VecView(grad, PETSC_VIEWER_STDOUT_WORLD);
    if (mpirank_world == 0 && !quietmode) {
      printf("\nGradient norm: %1.14e\n", gnorm);
    }
    optimctx->output->writeGradient(grad);
  }

  /* --- Solve the optimization  --- */
  if (runtype == RunType::OPTIMIZATION) {
    /* Set initial starting point */
    optimctx->getStartingPoint(xinit);
    VecCopy(xinit, optimctx->xinit); // Store the initial guess
    if (mpirank_world == 0 && !quietmode) printf("\nStarting Optimization solver ... \n");
    optimctx->timestepper->writeTrajectoryDataFiles = false;
    optimctx->solve(xinit);
    optimctx->getSolution(&opt);
  }

  /* Only evaluate and write control pulses (no propagation) */
  if (runtype == RunType::EVALCONTROLS) {
    std::vector<double> pt, qt;
    optimctx->getStartingPoint(xinit);
    if (mpirank_world == 0 && !quietmode) printf("\nEvaluating current controls ... \n");
    output->writeControls(xinit, mastereq);
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


#if TEST_FD_HESS
  if (mpirank_world == 0)  {
    printf("\n\n#########################\n");
    printf(" FD Testing for Hessian... \n");
    printf("#########################\n\n");
  }
  optimctx->getStartingPoint(xinit);

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
  delete mytimestepper;
  delete optimctx;
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
