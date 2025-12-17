#include "optimproblem.hpp"

OptimProblem::OptimProblem(Config config, TimeStepper* timestepper_, MPI_Comm comm_init_, MPI_Comm comm_optim_, int ninit_, Output* output_, bool quietmode_){

  timestepper = timestepper_;
  ninit = ninit_;
  output = output_;
  quietmode = quietmode_;
  /* Reset */
  objective = 0.0;

  /* Store communicators */
  comm_init = comm_init_;
  comm_optim = comm_optim_;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank_world);
  MPI_Comm_size(MPI_COMM_WORLD, &mpisize_world);
  MPI_Comm_rank(PETSC_COMM_WORLD, &mpirank_petsc);
  MPI_Comm_size(PETSC_COMM_WORLD, &mpisize_petsc);
  MPI_Comm_rank(comm_init, &mpirank_init);
  MPI_Comm_size(comm_init, &mpisize_init);
  MPI_Comm_rank(comm_optim, &mpirank_optim);
  MPI_Comm_size(comm_optim, &mpisize_optim);

  /* Store number of initial conditions per init-processor group */
  ninit_local = ninit / mpisize_init; 

  /*  If Schroedingers solver, allocate storage for the final states at time T for each initial condition. Schroedinger's solver does not store the time-trajectories during forward ODE solve, but instead recomputes the primal states during the adjoint solve. Therefore we need to store the terminal condition for the backwards primal solve. Be aware that the final states stored here will be overwritten during backwards computation!! */
  if (timestepper->mastereq->lindbladtype == LindbladType::NONE) {
    for (int i = 0; i < ninit_local; i++) {

      PetscInt globalsize = 2 * timestepper->mastereq->getDim();  // Global state vector: 2 for real and imaginary part
      PetscInt localsize = globalsize / mpisize_petsc;  // Local vector per processor
      Vec state;
      VecCreate(PETSC_COMM_WORLD, &state);
      VecSetSizes(state, localsize, globalsize);
      VecSetFromOptions(state);
      store_finalstates.push_back(state);
    }
  }

  // Check for new Riemannian objective function:
  use_new_objective = config.GetBoolParam("use_new_objective", false, true);
  phase_invariant = config.GetBoolParam("phase_invariant", false, true);
  if (timestepper->mastereq->lindbladtype != LindbladType::NONE) {
    use_new_objective = false;
  }
  if (use_new_objective) {
    if (mpisize_petsc>1){
      printf("New objective function only works with one petsc core right now.\n");
      exit(1);
    }
    if (timestepper->mastereq->getDimEss() != timestepper->mastereq->getDim()){
      printf("New objective function does not work with guard levels right now.\n");
      exit(1);
    }
  }

  // Allocate storage for final-time unitary if new objective function is used
  if (use_new_objective) {
    PetscInt globalsize_rows = timestepper->mastereq->getDim();
    PetscInt globalsize_cols = ninit;;
    PetscInt localsize_rows = globalsize_rows / mpisize_petsc;
    PetscInt localsize_cols = globalsize_cols / mpisize_petsc;
    MatCreateDense(PETSC_COMM_WORLD, localsize_rows, localsize_cols, globalsize_rows, globalsize_cols, NULL, &U_final_re);
    MatCreateDense(PETSC_COMM_WORLD, localsize_rows, localsize_cols, globalsize_rows, globalsize_cols, NULL, &U_final_im);
    MatCreateDense(PETSC_COMM_WORLD, localsize_rows, localsize_cols, globalsize_rows, globalsize_cols, NULL, &U_final_re_bar);
    MatCreateDense(PETSC_COMM_WORLD, localsize_rows, localsize_cols, globalsize_rows, globalsize_cols, NULL, &U_final_im_bar);
    MatSetUp(U_final_re);
    MatSetUp(U_final_im);
    MatSetUp(U_final_re_bar);
    MatSetUp(U_final_im_bar);
    MatZeroEntries(U_final_re);
    MatZeroEntries(U_final_im);
    MatZeroEntries(U_final_re_bar);
    MatZeroEntries(U_final_im_bar);
    MatAssemblyBegin(U_final_re, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(U_final_im, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(U_final_re_bar, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(U_final_im_bar, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(U_final_re, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(U_final_im, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(U_final_re_bar, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(U_final_im_bar, MAT_FINAL_ASSEMBLY);
  }


  /* Store number of design parameters */
  int n = 0;
  for (size_t ioscil = 0; ioscil < timestepper->mastereq->getNOscillators(); ioscil++) {
      n += timestepper->mastereq->getOscillator(ioscil)->getNParams(); 
  }
  ndesign = n;
  if (mpirank_world == 0 && !quietmode) std::cout<< "Number of control parameters: " << ndesign << std::endl;

  /* Allocate the initial condition vector and adjoint terminal state */
  VecCreate(PETSC_COMM_WORLD, &rho_t0); 
  PetscInt globalsize = 2 * timestepper->mastereq->getDim();  // Global state vector: 2 for real and imaginary part
  PetscInt localsize = globalsize / mpisize_petsc;  // Local vector per processor
  VecSetSizes(rho_t0, localsize, globalsize);
  VecSetFromOptions(rho_t0);
  VecZeroEntries(rho_t0);
  VecAssemblyBegin(rho_t0); VecAssemblyEnd(rho_t0);
  VecDuplicate(rho_t0, &rho_t0_bar);
  VecZeroEntries(rho_t0_bar);
  VecAssemblyBegin(rho_t0_bar); VecAssemblyEnd(rho_t0_bar);

  /* Initialize the optimization target, including setting of initial state rho_t0 if read from file or pure state or ensemble */
  std::vector<std::string> target_str;
  config.GetVecStrParam("optim_target", target_str, "pure");
  std::string objective_str = config.GetStrParam("optim_objective", "Jfrobenius");
  std::vector<double> read_gate_rot;
  config.GetVecDoubleParam("gate_rot_freq", read_gate_rot, 1e20, true, false); 
  std::vector<std::string> initcond_str;
  config.GetVecStrParam("initialcondition", initcond_str, "none", false);
  optim_target = new OptimTarget(target_str, objective_str, initcond_str, timestepper->mastereq, timestepper->total_time, read_gate_rot, rho_t0, quietmode);

  /* Get weights for the objective function (weighting the different initial conditions */
  config.GetVecDoubleParam("optim_weights", obj_weights, 1.0);
  int nfill = 0;
  if (obj_weights.size() < ninit) nfill = ninit - obj_weights.size();
  double val = obj_weights[obj_weights.size()-1];
  if (obj_weights.size() < ninit){
    for (int i = 0; i < nfill; i++) 
      obj_weights.push_back(val);
  }
  assert(obj_weights.size() >= ninit);
  // Scale the weights such that they sum up to one: beta_i <- beta_i / (\sum_i beta_i)
  double scaleweights = 0.0;
  for (size_t i=0; i<ninit; i++) scaleweights += obj_weights[i];
  for (size_t i=0; i<ninit; i++) obj_weights[i] = obj_weights[i] / scaleweights;
  // Distribute over mpi_init processes 
  std::vector<double> sendbuf = obj_weights;
  std::vector<double> recvbuf = obj_weights;
  int nscatter = ninit_local;
  MPI_Scatter(sendbuf.data(), nscatter, MPI_DOUBLE, recvbuf.data(), nscatter,  MPI_DOUBLE, 0, comm_init);
  for (int i = 0; i < nscatter; i++) obj_weights[i] = recvbuf[i];
  for (size_t i=nscatter; i < obj_weights.size(); i++) obj_weights[i] = 0.0;


  /* Store other optimization parameters */
  gamma_tik = config.GetDoubleParam("optim_regul", 1e-4);
  gatol = config.GetDoubleParam("optim_atol", 1e-8);
  fatol = config.GetDoubleParam("optim_ftol", 1e-8);
  inftol = config.GetDoubleParam("optim_inftol", 1e-5);
  grtol = config.GetDoubleParam("optim_rtol", 1e-4);
  maxiter = config.GetIntParam("optim_maxiter", 200);
  gamma_penalty = config.GetDoubleParam("optim_penalty", 0.0);
  penalty_param = config.GetDoubleParam("optim_penalty_param", 0.5);
  gamma_penalty_energy = config.GetDoubleParam("optim_penalty_energy", 0.0);
  // Check for new parameter name first, then old name for backward compatibility
  if (config.count("optim_regul_tik0") > 0) {
    gamma_tik_interpolate = config.GetBoolParam("optim_regul_tik0", false, false);
  } else if (config.count("optim_regul_interpolate") > 0) {
    gamma_tik_interpolate = config.GetBoolParam("optim_regul_interpolate", false, false);
    if (mpirank_world == 0 && !quietmode) {
      printf("Warning: Parameter 'optim_regul_interpolate' is deprecated. Please use 'optim_regul_tik0' instead.\n");
    }
  } else {
    gamma_tik_interpolate = false;
  }
  gamma_penalty_dpdm = config.GetDoubleParam("optim_penalty_dpdm", 0.0);
  gamma_penalty_variation = config.GetDoubleParam("optim_penalty_variation", 0.01); 
  

  if (gamma_penalty_dpdm > 1e-13 && timestepper->mastereq->lindbladtype != LindbladType::NONE){
    if (mpirank_world == 0 && !quietmode) {
      printf("Warning: Disabling DpDm penalty term because it is not implemented for the Lindblad solver.\n");
    }
    gamma_penalty_dpdm = 0.0;
  }

  /* Pass information on objective function to the time stepper needed for penalty objective function */
  timestepper->penalty_param = penalty_param;
  timestepper->gamma_penalty = gamma_penalty;
  timestepper->gamma_penalty_dpdm = gamma_penalty_dpdm;
  timestepper->gamma_penalty_energy = gamma_penalty_energy;
  timestepper->optim_target = optim_target;

  /* Store optimization bounds */
  VecCreateSeq(PETSC_COMM_SELF, ndesign, &xlower);
  VecSetFromOptions(xlower);
  VecDuplicate(xlower, &xupper);
  int col = 0;
  for (size_t iosc = 0; iosc < timestepper->mastereq->getNOscillators(); iosc++){
    std::vector<std::string> bound_str;
    config.GetVecStrParam("control_bounds" + std::to_string(iosc), bound_str, "10000.0");
    for (size_t iseg = 0; iseg < timestepper->mastereq->getOscillator(iosc)->getNSegments(); iseg++){
      double boundval = 0.0;
      if (bound_str.size() <= iseg) boundval =  atof(bound_str[bound_str.size()-1].c_str());
      else boundval = atof(bound_str[iseg].c_str());
      // Scale bounds by the number of carrier waves, and convert to radians */
      boundval = boundval / (sqrt(2) * timestepper->mastereq->getOscillator(iosc)->getNCarrierfrequencies());
      boundval = boundval * 2.0*M_PI;
      for (int i=0; i<timestepper->mastereq->getOscillator(iosc)->getNSegParams(iseg); i++){
        VecSetValue(xupper, col + i, boundval, INSERT_VALUES);
        VecSetValue(xlower, col + i, -1. * boundval, INSERT_VALUES);
      }
      // Disable bound for phase if this is spline_amplitude control
      if (timestepper->mastereq->getOscillator(iosc)->getControlType() == ControlType::BSPLINEAMP) {
        for (size_t f = 0; f < timestepper->mastereq->getOscillator(iosc)->getNCarrierfrequencies(); f++){
          int nsplines = timestepper->mastereq->getOscillator(iosc)->getNSplines();
          boundval = 1e+10;
          VecSetValue(xupper, col + f*(nsplines+1) + nsplines, boundval, INSERT_VALUES);
          VecSetValue(xlower, col + f*(nsplines+1) + nsplines, -1.*boundval, INSERT_VALUES);
        }
      }
      col = col + timestepper->mastereq->getOscillator(iosc)->getNSegParams(iseg);
    }
  }
  VecAssemblyBegin(xlower); VecAssemblyEnd(xlower);
  VecAssemblyBegin(xupper); VecAssemblyEnd(xupper);

  /* Store the initial guess if read from file */
  std::vector<std::string> controlinit_str;
  config.GetVecStrParam("control_initialization0", controlinit_str, "constant, 0.0");
  if ( controlinit_str.size() > 0 && controlinit_str[0].compare("file") == 0 ) {
    assert(controlinit_str.size() >=2);
    for (int i=0; i<ndesign; i++) initguess_fromfile.push_back(0.0);
    if (mpirank_world == 0) read_vector(controlinit_str[1].c_str(), initguess_fromfile.data(), ndesign, quietmode);
    MPI_Bcast(initguess_fromfile.data(), ndesign, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
 
  /* Create Petsc's optimization solver */
  TaoCreate(PETSC_COMM_WORLD, &tao);
  /* Set optimization type and parameters */
  TaoSetType(tao,TAOBQNLS);         // Optim type: taoblmvm vs BQNLS ??
  TaoSetMaximumIterations(tao, maxiter);
  TaoSetTolerances(tao, gatol, PETSC_DEFAULT, grtol);
  TaoMonitorSet(tao, TaoMonitor, (void*)this, NULL);
  TaoSetVariableBounds(tao, xlower, xupper);
  TaoSetFromOptions(tao);
  /* Set user-defined objective and gradient evaluation routines */
  TaoSetObjective(tao, TaoEvalObjective, (void *)this);
  TaoSetGradient(tao, NULL, TaoEvalGradient,(void *)this);
  TaoSetObjectiveAndGradient(tao, NULL, TaoEvalObjectiveAndGradient, (void*) this);

  /* Allocate auxiliary vector */
  mygrad = new double[ndesign];

  /* Allocat xinit, xtmp */
  VecCreateSeq(PETSC_COMM_SELF, ndesign, &xinit);
  VecSetFromOptions(xinit);
  VecZeroEntries(xinit);
  VecCreateSeq(PETSC_COMM_SELF, ndesign, &xtmp);
  VecSetFromOptions(xtmp);
  VecZeroEntries(xtmp);
}


OptimProblem::~OptimProblem() {
  delete [] mygrad;
  delete optim_target;
  VecDestroy(&rho_t0);
  VecDestroy(&rho_t0_bar);

  VecDestroy(&xlower);
  VecDestroy(&xupper);
  VecDestroy(&xinit);
  VecDestroy(&xtmp);

  for (size_t i = 0; i < store_finalstates.size(); i++) {
    VecDestroy(&(store_finalstates[i]));
  }
  if (use_new_objective) {
    MatDestroy(&U_final_re);
    MatDestroy(&U_final_im);
    MatDestroy(&U_final_re_bar);
    MatDestroy(&U_final_im_bar);
  }

  TaoDestroy(&tao);
}



double OptimProblem::evalF(const Vec x) {

  MasterEq* mastereq = timestepper->mastereq;

  if (mpirank_world == 0 && !quietmode) printf("EVAL F... \n");

  /* Pass design vector x to oscillators */
  mastereq->setControlAmplitudes(x); 

  /*  Iterate over initial condition */
  obj_cost  = 0.0;
  obj_regul = 0.0;
  obj_penal = 0.0;
  obj_penal_dpdm = 0.0;
  obj_penal_energy = 0.0;
  obj_penal_variation = 0.0;
  fidelity = 0.0;
  double obj_cost_re = 0.0;
  double obj_cost_im = 0.0;
  double fidelity_re = 0.0;
  double fidelity_im = 0.0;
  for (int iinit = 0; iinit < ninit_local; iinit++) {
      
    /* Prepare the initial condition in [rank * ninit_local, ... , (rank+1) * ninit_local - 1] */
    int iinit_global = mpirank_init * ninit_local + iinit;
    int initid = optim_target->prepareInitialState(iinit_global, ninit, timestepper->mastereq->nlevels, timestepper->mastereq->nessential, rho_t0);
    if (mpirank_optim == 0 && !quietmode) printf("%d: Initial condition id=%d ...\n", mpirank_init, initid);

    /* If gate optimiztion, compute the target state rho^target = Vrho(0)V^dagger */
    optim_target->prepareTargetState(rho_t0);

    /* Run forward with initial condition initid */
    Vec finalstate = timestepper->solveODE(initid, rho_t0);

    /* Store the final state for Riemannian objective function */
    if (use_new_objective) {
      const PetscScalar *finalstate_array;
      VecGetArrayRead(finalstate, &finalstate_array);
      for (size_t row = 0; row < timestepper->mastereq->getDim(); row++) {
        int id_re = row;
        int id_im = row + timestepper->mastereq->getDim();
        MatSetValue(U_final_re, row, iinit_global, finalstate_array[id_re], INSERT_VALUES);
        MatSetValue(U_final_im, row, iinit_global, finalstate_array[id_im], INSERT_VALUES);
      }
      VecRestoreArrayRead(finalstate, &finalstate_array);
    }

    /* Add to integral penalty term */
    obj_penal += obj_weights[iinit] * gamma_penalty * timestepper->penalty_integral;

    /* Add to second derivative penalty term */
    obj_penal_dpdm += obj_weights[iinit] * gamma_penalty_dpdm * timestepper->penalty_dpdm;
    
    /* Add to energy integral penalty term */
    obj_penal_energy += obj_weights[iinit] * gamma_penalty_energy* timestepper->energy_penalty_integral;

    /* Evaluate J(finalstate) and add to final-time cost */
    double obj_iinit_re = 0.0;
    double obj_iinit_im = 0.0;
    optim_target->evalJ(finalstate,  &obj_iinit_re, &obj_iinit_im);
    obj_cost_re += obj_weights[iinit] * obj_iinit_re;
    obj_cost_im += obj_weights[iinit] * obj_iinit_im;

    /* Add to final-time fidelity */
    double fidelity_iinit_re = 0.0;
    double fidelity_iinit_im = 0.0;
    optim_target->HilbertSchmidtOverlap(finalstate, false, &fidelity_iinit_re, &fidelity_iinit_im);
    fidelity_re += 1./ ninit * fidelity_iinit_re;
    fidelity_im += 1./ ninit * fidelity_iinit_im;

    // printf("%d, %d: iinit obj_iinit: %f * (%1.14e + i %1.14e, Overlap=%1.14e + i %1.14e\n", mpirank_world, mpirank_init, obj_weights[iinit], obj_iinit_re, obj_iinit_im, fidelity_iinit_re, fidelity_iinit_im);
  }
  if (use_new_objective) {
    MatAssemblyBegin(U_final_re, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(U_final_im, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(U_final_re, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(U_final_im, MAT_FINAL_ASSEMBLY);
  }

  /* Sum up from initial conditions processors */
  double mypen = obj_penal;
  double mypen_dpdm = obj_penal_dpdm;
  double mypenen = obj_penal_energy;
  double mycost_re = obj_cost_re;
  double mycost_im = obj_cost_im;
  double myfidelity_re = fidelity_re;
  double myfidelity_im = fidelity_im;
  MPI_Allreduce(&mypen, &obj_penal, 1, MPI_DOUBLE, MPI_SUM, comm_init);
  MPI_Allreduce(&mypen_dpdm, &obj_penal_dpdm, 1, MPI_DOUBLE, MPI_SUM, comm_init);
  MPI_Allreduce(&mypenen, &obj_penal_energy, 1, MPI_DOUBLE, MPI_SUM, comm_init);
  MPI_Allreduce(&mycost_re, &obj_cost_re, 1, MPI_DOUBLE, MPI_SUM, comm_init);
  MPI_Allreduce(&mycost_im, &obj_cost_im, 1, MPI_DOUBLE, MPI_SUM, comm_init);
  MPI_Allreduce(&myfidelity_re, &fidelity_re, 1, MPI_DOUBLE, MPI_SUM, comm_init);
  MPI_Allreduce(&myfidelity_im, &fidelity_im, 1, MPI_DOUBLE, MPI_SUM, comm_init);

  /* Set the fidelity: If Schroedinger, need to compute the absolute value: Fid= |\sum_i \phi^\dagger \phi_target|^2 */
  if (timestepper->mastereq->lindbladtype == LindbladType::NONE) {
    fidelity = pow(fidelity_re, 2.0) + pow(fidelity_im, 2.0);
  } else {
    fidelity = fidelity_re; 
  }
 
  /* Finalize the objective function */
  obj_cost = optim_target->finalizeJ(obj_cost_re, obj_cost_im);

  /* New objective function: Riemannian distance */
  if (use_new_objective) {

    MPI_Barrier(MPI_COMM_WORLD);

    /* allreduce the U_final matrix */
    PetscScalar *data;
    MatDenseGetArray(U_final_re, &data);
    int size = timestepper->mastereq->getDim();
    MPI_Allreduce(MPI_IN_PLACE, data, size * size, MPIU_SCALAR, MPI_SUM, comm_init);
    MatDenseRestoreArray(U_final_re, &data);

    MatDenseGetArray(U_final_im, &data);
    MPI_Allreduce(MPI_IN_PLACE, data, size * size, MPIU_SCALAR, MPI_SUM, comm_init);
    MatDenseRestoreArray(U_final_im, &data);

    double obj_riemannian = optim_target->RiemannianDistance(U_final_re, U_final_im, phase_invariant);

    // if (mpirank_world == 0) printf("\nRiemannian distance objective: %1.14e\n\n", obj_riemannian);
    obj_cost = obj_riemannian;
  }

  /* Evaluate Tikhonov regularization term: gamma/2 * ||x-x0||^2*/
  double xnorm;
  if (!gamma_tik_interpolate){  // ||x||^2
    VecNorm(x, NORM_2, &xnorm);
  } else {
    VecCopy(x, xtmp);
    VecAXPY(xtmp, -1.0, xinit);    // xtmp =  x - x_0
    VecNorm(xtmp, NORM_2, &xnorm);
  }
  obj_regul = gamma_tik / 2. * pow(xnorm,2.0);

  /* Evaluate penality term for control variation */
  double var_reg = 0.0;
  for (size_t iosc = 0; iosc < timestepper->mastereq->getNOscillators(); iosc++){
    var_reg += timestepper->mastereq->getOscillator(iosc)->evalControlVariation(); // uses Oscillator::params instead of 'x'
  }
  obj_penal_variation = 0.5*gamma_penalty_variation*var_reg; 

  /* Sum, store and return objective value */
  objective = obj_cost + obj_regul + obj_penal + obj_penal_dpdm + obj_penal_energy + obj_penal_variation;

  /* Output */
  if (mpirank_world == 0 && !quietmode) {
    std::cout<< "Objective = " << std::scientific<<std::setprecision(14) << obj_cost << " + " << obj_regul << " + " << obj_penal << " + " << obj_penal_dpdm << " + " << obj_penal_energy << " + " << obj_penal_variation << std::endl;
    std::cout<< "Fidelity = " << fidelity  << std::endl;
  }

  return objective;
}



void OptimProblem::evalGradF(const Vec x, Vec G){

  MasterEq* mastereq = timestepper->mastereq;

  if (mpirank_world == 0 && !quietmode) std::cout<< "EVAL GRAD F... " << std::endl;

  /* Pass design vector x to oscillators */
  mastereq->setControlAmplitudes(x); 

  /* Reset Gradient */
  VecZeroEntries(G);

  /* Derivative of regulatization terms (ADD ON ONE PROC ONLY!) */
  // if (mpirank_init == 0 && mpirank_optim == 0) { // TODO: Which one?? 
  if (mpirank_init == 0 ) {

    // Derivative of Tikhonov 0.5*gamma * ||x||^2 
    VecAXPY(G, gamma_tik, x);   // + gamma_tik * x
    if (gamma_tik_interpolate){
      VecAXPY(G, -1.0*gamma_tik, xinit); // -gamma_tik * xinit
    }

    // Derivative of penalization of control variation 
    double var_reg_bar = 0.5*gamma_penalty_variation;
    int skip_to_oscillator = 0;
    for (size_t iosc = 0; iosc < timestepper->mastereq->getNOscillators(); iosc++){
      Oscillator* osc = timestepper->mastereq->getOscillator(iosc);
      osc->evalControlVariationDiff(G, var_reg_bar, skip_to_oscillator);
      skip_to_oscillator += osc->getNParams();
    }
  }

  /*  Iterate over initial condition */
  obj_cost = 0.0;
  obj_regul = 0.0;
  obj_penal = 0.0;
  obj_penal_dpdm = 0.0;
  obj_penal_energy = 0.0;
  obj_penal_variation = 0.0;
  fidelity = 0.0;
  double obj_cost_re = 0.0;
  double obj_cost_im = 0.0;
  double fidelity_re = 0.0;
  double fidelity_im = 0.0;
  for (int iinit = 0; iinit < ninit_local; iinit++) {

    /* Prepare the initial condition */
    int iinit_global = mpirank_init * ninit_local + iinit;
    int initid = optim_target->prepareInitialState(iinit_global, ninit, timestepper->mastereq->nlevels, timestepper->mastereq->nessential, rho_t0);

    /* If gate optimiztion, compute the target state rho^target = Vrho(0)V^dagger */
    optim_target->prepareTargetState(rho_t0);

    /* --- Solve primal --- */
    // if (mpirank_optim == 0) printf("%d: %d FWD. ", mpirank_init, initid);

    /* Run forward with initial condition rho_t0 */
    Vec finalstate = timestepper->solveODE(initid, rho_t0);

    /* Store the final state for Riemannian objective function */
    if (use_new_objective) {
      const PetscScalar *finalstate_array;
      VecGetArrayRead(finalstate, &finalstate_array);
      for (size_t row = 0; row < timestepper->mastereq->getDim(); row++) {
        int id_re = row;
        int id_im = row + timestepper->mastereq->getDim();
        MatSetValue(U_final_re, row, iinit_global, finalstate_array[id_re], INSERT_VALUES);
        MatSetValue(U_final_im, row, iinit_global, finalstate_array[id_im], INSERT_VALUES);
      }
      VecRestoreArrayRead(finalstate, &finalstate_array);
    }

    /* Store the final state for the Schroedinger solver */
    if (timestepper->mastereq->lindbladtype == LindbladType::NONE) VecCopy(finalstate, store_finalstates[iinit]);

    /* Add to integral penalty term */
    obj_penal += obj_weights[iinit] * gamma_penalty * timestepper->penalty_integral;

    /* Add to second derivative dpdm integral penalty term */
    obj_penal_dpdm += obj_weights[iinit] * gamma_penalty_dpdm * timestepper->penalty_dpdm;
    /* Add to energy integral penalty term */
    obj_penal_energy += obj_weights[iinit] * gamma_penalty_energy * timestepper->energy_penalty_integral;

    /* Evaluate J(finalstate) and add to final-time cost */
    double obj_iinit_re = 0.0;
    double obj_iinit_im = 0.0;
    optim_target->evalJ(finalstate,  &obj_iinit_re, &obj_iinit_im);
    obj_cost_re += obj_weights[iinit] * obj_iinit_re;
    obj_cost_im += obj_weights[iinit] * obj_iinit_im;

    /* Add to final-time fidelity */
    double fidelity_iinit_re = 0.0;
    double fidelity_iinit_im = 0.0;
    optim_target->HilbertSchmidtOverlap(finalstate, false, &fidelity_iinit_re, &fidelity_iinit_im);
    fidelity_re += 1./ ninit * fidelity_iinit_re;
    fidelity_im += 1./ ninit * fidelity_iinit_im;

    /* If Lindblas solver, compute adjoint for this initial condition. Otherwise (Schroedinger solver), compute adjoint only after all initial conditions have been propagated through (separate loop below) */
    if (timestepper->mastereq->lindbladtype != LindbladType::NONE) {
      // if (mpirank_optim == 0) printf("%d: %d BWD.", mpirank_init, initid);

      /* Reset adjoint */
      VecZeroEntries(rho_t0_bar);

      /* Terminal condition for adjoint variable: Derivative of final time objective J */
      double obj_cost_re_bar, obj_cost_im_bar;
      optim_target->finalizeJ_diff(obj_cost_re, obj_cost_im, &obj_cost_re_bar, &obj_cost_im_bar);
      optim_target->evalJ_diff(finalstate, rho_t0_bar, obj_weights[iinit]*obj_cost_re_bar, obj_weights[iinit]*obj_cost_im_bar);

      /* Derivative of time-stepping */
      timestepper->solveAdjointODE(rho_t0_bar, finalstate, obj_weights[iinit] * gamma_penalty, obj_weights[iinit]*gamma_penalty_dpdm, obj_weights[iinit]*gamma_penalty_energy);

      /* Add to optimizers's gradient */
      VecAXPY(G, 1.0, timestepper->redgrad);
    }
  }
  if (use_new_objective) {
    MatAssemblyBegin(U_final_re, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(U_final_im, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(U_final_re, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(U_final_im, MAT_FINAL_ASSEMBLY);
  }

  /* Sum up from initial conditions processors */
  double mypen = obj_penal;
  double mypen_dpdm = obj_penal_dpdm;
  double mypenen = obj_penal_energy;
  double mycost_re = obj_cost_re;
  double mycost_im = obj_cost_im;
  double myfidelity_re = fidelity_re;
  double myfidelity_im = fidelity_im;
  MPI_Allreduce(&mypen, &obj_penal, 1, MPI_DOUBLE, MPI_SUM, comm_init);
  MPI_Allreduce(&mypen_dpdm, &obj_penal_dpdm, 1, MPI_DOUBLE, MPI_SUM, comm_init);
  MPI_Allreduce(&mypenen, &obj_penal_energy, 1, MPI_DOUBLE, MPI_SUM, comm_init);
  MPI_Allreduce(&mycost_re, &obj_cost_re, 1, MPI_DOUBLE, MPI_SUM, comm_init);
  MPI_Allreduce(&mycost_im, &obj_cost_im, 1, MPI_DOUBLE, MPI_SUM, comm_init);
  MPI_Allreduce(&myfidelity_re, &fidelity_re, 1, MPI_DOUBLE, MPI_SUM, comm_init);
  MPI_Allreduce(&myfidelity_im, &fidelity_im, 1, MPI_DOUBLE, MPI_SUM, comm_init);

  /* Set the fidelity: If Schroedinger, need to compute the absolute value: Fid= |\sum_i \phi^\dagger \phi_target|^2 */
  if (timestepper->mastereq->lindbladtype == LindbladType::NONE) {
    fidelity = pow(fidelity_re, 2.0) + pow(fidelity_im, 2.0);
  } else {
    fidelity = fidelity_re; 
  }
 
  /* Finalize the objective function Jtrace to get the infidelity. 
     If Schroedingers solver, need to take the absolute value */
  obj_cost = optim_target->finalizeJ(obj_cost_re, obj_cost_im);

  /* New objective function: Riemannian distance */
  if (use_new_objective) {

    MPI_Barrier(MPI_COMM_WORLD);

    /* allreduce the U_final matrix */
    PetscScalar *data;
    MatDenseGetArray(U_final_re, &data);
    int size = timestepper->mastereq->getDim();
    MPI_Allreduce(MPI_IN_PLACE, data, size * size, MPIU_SCALAR, MPI_SUM, comm_init);
    MatDenseRestoreArray(U_final_re, &data);

    MatDenseGetArray(U_final_im, &data);
    MPI_Allreduce(MPI_IN_PLACE, data, size * size, MPIU_SCALAR, MPI_SUM, comm_init);
    MatDenseRestoreArray(U_final_im, &data);

    double obj_riemannian = optim_target->RiemannianDistance(U_final_re, U_final_im, phase_invariant);

    // if (mpirank_world == 0) printf("\nRiemannian distance objective: %1.14e\n\n", obj_riemannian);
    obj_cost = obj_riemannian;
  }

  /* Evaluate Tikhonov regularization term += gamma/2 * ||x||^2*/
  double xnorm;
  if (!gamma_tik_interpolate){  // ||x||^2
    VecNorm(x, NORM_2, &xnorm);
  } else {
    VecCopy(x, xtmp);
    VecAXPY(xtmp, -1.0, xinit);    // xtmp =  x_k - x_0
    VecNorm(xtmp, NORM_2, &xnorm);
  }
  obj_regul = gamma_tik / 2. * pow(xnorm,2.0);

  /* Evaluate penalty term for control parameter variation */
  double var_reg = 0.0;
  for (size_t iosc = 0; iosc < timestepper->mastereq->getNOscillators(); iosc++){
    var_reg += timestepper->mastereq->getOscillator(iosc)->evalControlVariation(); // uses Oscillator::params instead of 'x'
  }
  obj_penal_variation = 0.5*gamma_penalty_variation*var_reg; 

  /* Sum, store and return objective value */
  objective = obj_cost + obj_regul + obj_penal + obj_penal_dpdm + obj_penal_energy + obj_penal_variation;

  /* Derivative of new objective function */
  if (use_new_objective) {
    optim_target->RiemannianDistance_diff(U_final_re, U_final_im, U_final_re_bar, U_final_im_bar, phase_invariant);
  }


  /* For Schroedinger solver: Solve adjoint equations for all initial conditions here. */
  if (timestepper->mastereq->lindbladtype == LindbladType::NONE) {

    // Iterate over all initial conditions 
    for (int iinit = 0; iinit < ninit_local; iinit++) {
      int iinit_global = mpirank_init * ninit_local + iinit;

      /* Recompute the initial state and target */
      optim_target->prepareInitialState(iinit_global, ninit, timestepper->mastereq->nlevels, timestepper->mastereq->nessential, rho_t0);
      optim_target->prepareTargetState(rho_t0);
     
      /* Reset adjoint */
      VecZeroEntries(rho_t0_bar);

      /* Terminal condition for adjoint variable: Derivative of final time objective J */
      double obj_cost_re_bar, obj_cost_im_bar;
      optim_target->finalizeJ_diff(obj_cost_re, obj_cost_im, &obj_cost_re_bar, &obj_cost_im_bar);
      optim_target->evalJ_diff(store_finalstates[iinit], rho_t0_bar, obj_weights[iinit]*obj_cost_re_bar, obj_weights[iinit]*obj_cost_im_bar);

      /* Derivative of time-stepping */
      timestepper->solveAdjointODE(rho_t0_bar, store_finalstates[iinit], obj_weights[iinit] * gamma_penalty, obj_weights[iinit]*gamma_penalty_dpdm, obj_weights[iinit]*gamma_penalty_energy);

      /* Add to optimizers's gradient */
      VecAXPY(G, 1.0, timestepper->redgrad);
    } // end of initial condition loop 
  } // end of adjoint for Schroedinger

  /* Sum up the gradient from all initial condition processors */
  PetscScalar* grad; 
  VecGetArray(G, &grad);
  for (int i=0; i<ndesign; i++) {
    mygrad[i] = grad[i];
  }
  MPI_Allreduce(mygrad, grad, ndesign, MPI_DOUBLE, MPI_SUM, comm_init);
  VecRestoreArray(G, &grad);

  /* Compute and store gradient norm */
  VecNorm(G, NORM_2, &(gnorm));

  /* Output */
  // if (mpirank_world == 0 && !quietmode) {
  //   std::cout<< "Objective = " << std::scientific<<std::setprecision(14) << obj_cost << " + " << obj_regul << " + " << obj_penal << " + " << obj_penal_dpdm << " + " << obj_penal_energy << " + " << obj_penal_variation << std::endl;
  //   std::cout<< "Fidelity = " << fidelity << std::endl;
  // }
}

void OptimProblem::solve(Vec xinit) {
  TaoSetSolution(tao, xinit);
  TaoSolve(tao);
}

void OptimProblem::getStartingPoint(Vec xinit){
  MasterEq* mastereq = timestepper->mastereq;

  if (initguess_fromfile.size() > 0) {
    /* Set the initial guess from file */
    for (size_t i=0; i<initguess_fromfile.size(); i++) {
      VecSetValue(xinit, i, initguess_fromfile[i], INSERT_VALUES);
    }

  } else { // copy from initialization in oscillators contructor
    PetscScalar* xptr;
    VecGetArray(xinit, &xptr);
    int shift = 0;
    for (size_t ioscil = 0; ioscil<mastereq->getNOscillators(); ioscil++){
      mastereq->getOscillator(ioscil)->getParams(xptr + shift);
      shift += mastereq->getOscillator(ioscil)->getNParams();
    }
    VecRestoreArray(xinit, &xptr);
  }

  /* Assemble initial guess */
  VecAssemblyBegin(xinit);
  VecAssemblyEnd(xinit);

  /* Pass to oscillator */
  timestepper->mastereq->setControlAmplitudes(xinit);
  
  /* Write initial control functions to file */
  output->writeControls(xinit, timestepper->mastereq, timestepper->ntime, timestepper->dt);

}


void OptimProblem::getSolution(Vec* param_ptr){
  
  /* Get ref to optimized parameters */
  Vec params;
  TaoGetSolution(tao, &params);
  *param_ptr = params;
}

PetscErrorCode TaoMonitor(Tao tao,void*ptr){
  OptimProblem* ctx = (OptimProblem*) ptr;

  /* Get information from Tao optimization */
  PetscInt iter;
  PetscScalar deltax;
  Vec params;
  TaoConvergedReason reason;
  PetscScalar f, gnorm;
  TaoGetSolutionStatus(tao, &iter, &f, &gnorm, NULL, &deltax, &reason);
  TaoGetSolution(tao, &params);

  /* Grab some output stuff */
  double obj_cost = ctx->getCostT();
  double obj_regul = ctx->getRegul();
  double obj_penal = ctx->getPenalty();
  double obj_penal_dpdm = ctx->getPenaltyDpDm();
  double obj_penal_energy = ctx->getPenaltyEnergy();
  double obj_penal_variation= ctx->getPenaltyVariation();
  double F_avg = ctx->getFidelity();

  /* Additional Stopping criteria */
  bool lastIter = false;
  std::string finalReason_str = "";
  if (1.0 - F_avg <= ctx->getInfTol()) {
    finalReason_str = "Optimization converged with small infidelity.";
    TaoSetConvergedReason(tao, TAO_CONVERGED_USER);
    lastIter = true;
  } else if (obj_cost <= ctx->getFaTol()) {
    finalReason_str = "Optimization converged with small final time cost.";
    TaoSetConvergedReason(tao, TAO_CONVERGED_USER);
    lastIter = true;
  } else if (iter == ctx->getMaxIter()) {
    finalReason_str = "Optimization stopped at maximum number of iterations.";
    lastIter = true;
  } else if (gnorm < ctx->getGaTol()) {
    finalReason_str = "OPtimization converged with small gradient norm.";
    lastIter=true;
  }

  /* First iteration: Header for screen output of optimization history */
  if (iter == 0 && ctx->getMPIrank_world() == 0) {
    std::cout<<  "    Objective             Tikhonov                Penalty-Leakage        Penalty-StateVar       Penalty-TotalEnergy    Penalty-CtrlVar" << std::endl;
  }

  /* Every <optim_monitor_freq> iterations: Output of optimization history */
  if (iter % ctx->output->optim_monitor_freq == 0 ||lastIter) {
    // Add to optimization history file 
    ctx->output->writeOptimFile(iter, f, gnorm, deltax, F_avg, obj_cost, obj_regul, obj_penal, obj_penal_dpdm, obj_penal_energy, obj_penal_variation);
    // Screen output 
    if (ctx->getMPIrank_world() == 0) {
      std::cout<< iter <<  "  " << std::scientific<<std::setprecision(14) << obj_cost << " + " << obj_regul << " + " << obj_penal << " + " << obj_penal_dpdm << " + " << obj_penal_energy << " + " << obj_penal_variation;
      std::cout<< "  Fidelity = " << F_avg;
      std::cout<< "  ||Grad|| = " << gnorm;
      std::cout<< std::endl;
    }
  }

  /* Last iteration: Print solution, controls and trajectory data to files */
  if (lastIter) {
    ctx->output->writeControls(params, ctx->timestepper->mastereq, ctx->timestepper->ntime, ctx->timestepper->dt);

    // do one last forward evaluation while writing trajectory files
    ctx->timestepper->writeTrajectoryDataFiles = true;
    ctx->evalF(params); 

    // Print stopping reason to screen
    if (ctx->getMPIrank_world() == 0){
      std::cout<< finalReason_str << std::endl;
    }
  }


  return 0;
}


PetscErrorCode TaoEvalObjectiveAndGradient(Tao tao, Vec x, PetscReal *f, Vec G, void*ptr){

  TaoEvalGradient(tao, x, G, ptr);
  OptimProblem* ctx = (OptimProblem*) ptr;
  *f = ctx->getObjective();

  return 0;
}

PetscErrorCode TaoEvalObjective(Tao /*tao*/, Vec x, PetscReal *f, void*ptr){

  OptimProblem* ctx = (OptimProblem*) ptr;
  *f = ctx->evalF(x);
  
  return 0;
}


PetscErrorCode TaoEvalGradient(Tao /*tao*/, Vec x, Vec G, void*ptr){

  OptimProblem* ctx = (OptimProblem*) ptr;
  ctx->evalGradF(x, G);
  
  return 0;
}
