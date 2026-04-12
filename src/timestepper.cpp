#include "timestepper.hpp"
#include "petscvec.h"

TimeStepper::TimeStepper() {
  mastereq = NULL;
  ntime = 0;
  total_time = 0.0;
  dt = 0.0;
  storeFWD = false;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank_world);
  MPI_Comm_rank(PETSC_COMM_WORLD, &mpirank_petsc);
  MPI_Comm_size(PETSC_COMM_WORLD, &mpisize_petsc);
  writeTrajectoryDataFiles = false;
}

TimeStepper::TimeStepper(size_t ninit_local, MasterEq* mastereq_, int ntime_, double total_time_, Output* output_, bool storeFWD_) : TimeStepper() {
  mastereq = mastereq_;
  ntime = ntime_;
  total_time = total_time_;
  output = output_;
  storeFWD = storeFWD_;
  eval_leakage = false;
  eval_energy = false;
  eval_dpdm = false;
  eval_weightedcost = false;
  weightedcost_width = 0.0;

  // Set local sizes of subvectors u,v in state x=[u,v]
  localsize_u = mastereq->getDim() / mpisize_petsc; 
  ilow = mpirank_petsc * localsize_u;
  iupp = ilow + localsize_u;         

  /* Set the time-step size */
  dt = total_time / ntime;

  /* Allocate storage of primal state */
  if (storeFWD) { 
    for (int n = 0; n <=ntime; n++) {
      Vec state;
      VecCreate(PETSC_COMM_WORLD, &state);
      PetscInt globalsize = 2 * mastereq->getDim();  // 2 for real and imaginary part
      PetscInt localsize = globalsize / mpisize_petsc;  // Local vector per processor
      VecSetSizes(state,localsize,globalsize);
      VecSetFromOptions(state);
      store_states.push_back(state);
    }
  }

  /* Allocate auxiliary state vector */
  VecCreate(PETSC_COMM_WORLD, &x);

  PetscInt globalsize = 2 * mastereq->getDim(); 
  PetscInt localsize = globalsize / mpisize_petsc;  // Local vector per processor
  VecSetSizes(x,localsize,globalsize);
  VecSetFromOptions(x);
  VecZeroEntries(x);
  VecDuplicate(x, &xadj);
  VecDuplicate(x, &xprimal);

  /* Allocate the reduced gradient */
  int ndesign = 0;
  for (size_t ioscil = 0; ioscil < mastereq->getNOscillators(); ioscil++) {
      ndesign += mastereq->getOscillator(ioscil)->getNParams(); 
  }
  VecCreateSeq(PETSC_COMM_SELF, ndesign, &redgrad);
  VecSetFromOptions(redgrad);
  VecAssemblyBegin(redgrad);
  VecAssemblyEnd(redgrad);

}


TimeStepper::~TimeStepper() {
  for (size_t n = 0; n < store_states.size(); n++) {
    VecDestroy(&(store_states[n]));
  }
  VecDestroy(&x);
  VecDestroy(&xadj);
  VecDestroy(&xprimal);
  VecDestroy(&redgrad);
}



Vec TimeStepper::getState(size_t tindex){
  
  if (tindex >= store_states.size()) {
    printf("ERROR: Time-stepper requested state at time index %zu, but didn't store it.\n", tindex);
    exit(1);
  }

  return store_states[tindex];
}

Vec TimeStepper::solveODE(int initid, int iinit_local, Vec rho_t0){

  /* Prepare storage for trajectory output data */
  if (writeTrajectoryDataFiles) {
    output->resetTrajectoryData(initid);
  }

  /* Set initial condition  */
  VecCopy(rho_t0, x);


  /* Store initial state for dpdm integral term */
  if (eval_dpdm){
    for (int i = 0; i < 2; i++) {
      Vec state;
      VecCreate(PETSC_COMM_WORLD, &state);
      PetscInt globalsize = 2 * mastereq->getDim(); 
      PetscInt localsize = globalsize / mpisize_petsc;  // Local vector per processor
      VecSetSizes(state,localsize,globalsize);
      VecSetFromOptions(state);
      dpdm_states.push_back(state);
    }
    VecCopy(x, dpdm_states[0]);
  }

  /* --- Loop over time interval --- */
  leakage_integral = 0.0;
  weightedcost_integral = 0.0;
  energy_integral = 0.0;
  dpdm_integral = 0.0;
  for (int n = 0; n < ntime; n++){

    /* current time */
    double tstart = n * dt;
    double tstop  = (n+1) * dt;

    /* store and write current state. */
    if (storeFWD) VecCopy(x, store_states[n]);
    if (writeTrajectoryDataFiles) {
      output->evalTrajectoryData(n, tstart, x, mastereq);
    }

    /* Take one time step */
    evolveFWD(tstart, tstop, x);

    /* Add to integral leakage term */
    if (eval_leakage) leakage_integral += evalLeakage(x);

    /* Add to integral weighted cost function term */
    if (eval_weightedcost) weightedcost_integral += evalWeightedCost(tstop, x);

    /* Add to second derivative variation term */
    if (eval_dpdm) {
      if (n > 0) dpdm_integral += evalDpDm(x, dpdm_states[n%2], dpdm_states[(n+1)%2]);  // uses x, x_n, x_n-1
      // Update storage of primal states. Should build a history of 3 states.
      VecCopy(x, dpdm_states[(n+1)%2]);
    }

    /* Add to energy integral term */
    if (eval_energy) energy_integral += evalEnergy(tstop);

#ifdef SANITY_CHECK
    SanityTests(x, tstart);
#endif
  }

  /* Average integral terms over time steps */
  dpdm_integral = dpdm_integral/ntime;
  leakage_integral = leakage_integral / ntime;
  energy_integral = energy_integral / ntime;
  weightedcost_integral = weightedcost_integral * dt;

  /* Store last time step */
  if (storeFWD) VecCopy(x, store_states[ntime]);

  /* Clear out dpdm storage */
  if (eval_dpdm) {
    for (size_t i=0; i<dpdm_states.size(); i++) {
      VecDestroy(&(dpdm_states[i]));
    }
    dpdm_states.clear();
  }

  /* Write last time step and close files */
  if (writeTrajectoryDataFiles) {
    output->evalTrajectoryData(ntime, ntime*dt, x, mastereq);
    output->writeTrajectoryData();
  }
  

  return x;
}


void TimeStepper::solveAdjointODE(int iinit_local, Vec rho_t0_bar, Vec finalstate, double Jbar_leakage, double Jbar_weightedcost, double Jbar_dpdm, double Jbar_energy) {

  /* Reset gradient */
  VecZeroEntries(redgrad);

  /* Set terminal adjoint condition */
  VecCopy(rho_t0_bar, xadj);

  /* Set terminal primal state */
  VecCopy(finalstate, xprimal);

  /* Store states at N, N-1, N-2 for dpdm integral term */
  if (eval_dpdm){
    for (int i = 0; i < 5; i++) {
      Vec state;
      VecCreate(PETSC_COMM_WORLD, &state);
      PetscInt globalsize = 2 * mastereq->getDim(); 
      PetscInt localsize = globalsize / mpisize_petsc;  // Local vector per processor
      VecSetSizes(state,localsize,globalsize);
      VecSetFromOptions(state);
      dpdm_states.push_back(state);
    }
    VecCopy(xprimal, dpdm_states[2]);  // x_N
    VecCopy(dpdm_states[2], dpdm_states[1]);
    evolveFWD(ntime*dt, (ntime-1)*dt, dpdm_states[1]);  // x_N-1
    VecCopy(dpdm_states[1], dpdm_states[0]);
    evolveFWD((ntime-1)*dt, (ntime-2)*dt, dpdm_states[0]);  // x_N-2
  }
 

  /* Loop over time interval */
  for (int n = ntime; n > 0; n--){
    double tstop  = n * dt;
    double tstart = (n-1) * dt;
    // printf("Backwards %d -> %d ... \n", n, n-1);

    /* Derivative of energy integral term */
    if (eval_energy) evalEnergy_diff(tstop, Jbar_energy/ntime, redgrad);

    /* Derivative of dpdm term */
    if (eval_dpdm) evalDpDm_diff(n, xadj, Jbar_dpdm/ntime);

    /* Derivative of weighted cost term */
    if (eval_weightedcost) evalWeightedCost_diff(tstop, xprimal, xadj, Jbar_weightedcost*dt);

    /* Derivative of leakage term */
    if (eval_leakage) evalLeakage_diff(xprimal, xadj, Jbar_leakage/ntime);

    /* Get the state at n-1. If Schroedinger solver, recompute it by taking a step backwards with the forward solver, otherwise get it from storage. */
    if (storeFWD) VecCopy(getState(n-1), xprimal);
    else evolveFWD(tstop, tstart, xprimal);

    /* Take one time step backwards for the adjoint */
    evolveBWD(tstop, tstart, xprimal, xadj, redgrad, true);

    /* Update dpdm storage */
    if (eval_dpdm) {
      int k = ntime - n;
      int idx = 4-((k+4)%5);
      int idx1 = 4-(k%5);
      VecCopy(dpdm_states[idx], dpdm_states[idx1]);
      if (n>2) evolveFWD((n-2)*dt, (n-3)*dt, dpdm_states[idx1]);
    }
  }

  /* Clear out dpdm storage */
  if (eval_dpdm) {
    for (size_t i=0; i<dpdm_states.size(); i++) {
      VecDestroy(&(dpdm_states[i]));
    }
    dpdm_states.clear();
  }
}

double TimeStepper::evalLeakage(const Vec x){
  PetscInt dim_rho = mastereq->getDimRho(); // N
  double x_re, x_im;

  /* Add guard-level occupation to prevent leakage. A guard level is the LAST NON-ESSENTIAL energy level of an oscillator */
  double leakage = 0.0;
  /* Sum over all diagonal elements that correspond to a non-essential guard level. */
  for (PetscInt i=0; i<dim_rho; i++) {
    if ( isGuardLevel(i, mastereq->nlevels, mastereq->nessential) ) {
      // vectorize if lindblad
      PetscInt vecID = i;
      if (mastereq->decoherence_type != DecoherenceType::NONE) vecID = getVecID(i,i,dim_rho);
      x_re = 0.0; 
      x_im = 0.0;
      if (ilow <= vecID && vecID < iupp) {
        PetscInt id_global_x = vecID + mpirank_petsc*localsize_u; 
        VecGetValues(x, 1, &id_global_x, &x_re);
        id_global_x += localsize_u; 
        VecGetValues(x, 1, &id_global_x, &x_im); 
      }
      leakage += (x_re * x_re + x_im * x_im);
    }
  }
  double mine = leakage;
  MPI_Allreduce(&mine, &leakage, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

  return leakage;
}

void TimeStepper::evalLeakage_diff(const Vec x, Vec xbar, double Jbar){
  PetscInt dim_rho = mastereq->getDimRho();  // N

  /* If gate optimization: Derivative of adding guard-level occupation */
  double x_re, x_im;
  for (PetscInt i=0; i<dim_rho; i++) {
    if ( isGuardLevel(i, mastereq->nlevels, mastereq->nessential) ) {
      PetscInt  vecID = i;
      if (mastereq->decoherence_type != DecoherenceType::NONE) vecID = getVecID(i,i,dim_rho);
      x_re = 0.0; 
      x_im = 0.0;
      if (ilow <= vecID && vecID < iupp) {
        PetscInt id_global_x = vecID + mpirank_petsc*localsize_u; 
        VecGetValues(x, 1, &id_global_x, &x_re);
        VecSetValue(xbar, id_global_x, 2.*x_re*Jbar, ADD_VALUES);
        id_global_x += localsize_u; 
        VecGetValues(x, 1, &id_global_x, &x_im);
        VecSetValue(xbar, id_global_x, 2.*x_im*Jbar, ADD_VALUES);
      }
    }
  }
  VecAssemblyBegin(xbar);
  VecAssemblyEnd(xbar);
}

double TimeStepper::evalWeightedCost(double time, const Vec x){

  if (weightedcost_width <= 0.0) {
    return 0.0;
  }
  /* weighted integral of the objective function */
  double weight = 1./weightedcost_width* exp(- pow((time - total_time)/weightedcost_width, 2));
  double obj_re = 0.0;
  double obj_im = 0.0;
  optim_target->evalJ(x, &obj_re, &obj_im);
  double obj_cost = optim_target->finalizeJ(obj_re, obj_im);

  return weight * obj_cost;
}

void TimeStepper::evalWeightedCost_diff(double time, const Vec x, Vec xbar, double Jbar){
  if (weightedcost_width <= 0.0) {
    return;
  }
  /* Derivative of weighted integral of the objective function */
  double weight = 1./weightedcost_width* exp(- pow((time - total_time)/weightedcost_width, 2));
  
  double obj_cost_re = 0.0;
  double obj_cost_im = 0.0;
  optim_target->evalJ(x, &obj_cost_re, &obj_cost_im);

  double obj_cost_re_bar = 0.0; 
  double obj_cost_im_bar = 0.0;
  optim_target->finalizeJ_diff(obj_cost_re, obj_cost_im, &obj_cost_re_bar, &obj_cost_im_bar);
  optim_target->evalJ_diff(x, xbar, weight*obj_cost_re_bar*Jbar, weight*obj_cost_im_bar*Jbar);
  
}


double TimeStepper::evalDpDm(Vec x, Vec xm1, Vec xm2){

    // Get local data pointers
    const PetscScalar *xptr, *xm1ptr, *xm2ptr;
    VecGetArrayRead(x, &xptr);
    VecGetArrayRead(xm1, &xm1ptr);
    VecGetArrayRead(xm2, &xm2ptr);
     
    /* precompute 1/dt^4 */
    double dtinv = 1.0 / (dt*dt*dt*dt);

    double dpdm = 0.0;
    for (PetscInt i=0; i<localsize_u; i++) {
        PetscInt vecID_re = i;
        PetscInt vecID_im = i + localsize_u;

        double tmp1 = xptr[vecID_re]*xptr[vecID_re] - 2.0*xm1ptr[vecID_re]*xm1ptr[vecID_re] + xm2ptr[vecID_re]*xm2ptr[vecID_re];
        double tmp2 = xptr[vecID_im]*xptr[vecID_im] - 2.0*xm1ptr[vecID_im]*xm1ptr[vecID_im] + xm2ptr[vecID_im]*xm2ptr[vecID_im];

        dpdm += dtinv * (tmp1 + tmp2)*(tmp1 + tmp2);
    }

    VecRestoreArrayRead(x, &xptr);
    VecRestoreArrayRead(xm1, &xm1ptr);
    VecRestoreArrayRead(xm2, &xm2ptr);

    return dpdm;
}


void TimeStepper::evalDpDm_diff(int n, Vec xbar, double Jbar){

    const PetscScalar *xptr, *xm1ptr, *xm2ptr, *xp1ptr, *xp2ptr;
    Vec x = nullptr, xm1 = nullptr, xm2 = nullptr, xp1 = nullptr, xp2 = nullptr;

    int k = ntime - n;
    int idx = 4-(k+4)%5;
    if (n > 1)         xm2 = dpdm_states[idx];
    if (n > 0 )        xm1 = dpdm_states[(idx+1)%5];
    x = dpdm_states[(idx+2)%5];
    if (n < ntime)     xp1 = dpdm_states[(idx+3)%5];
    if (n < ntime-1)   xp2 = dpdm_states[(idx+4)%5];

    VecGetArrayRead(x, &xptr);
    if (n > 0)        VecGetArrayRead(xm1, &xm1ptr);
    if (n > 1)        VecGetArrayRead(xm2, &xm2ptr);
    if (n < ntime)    VecGetArrayRead(xp1, &xp1ptr);
    if (n < ntime-1)  VecGetArrayRead(xp2, &xp2ptr);

    /* precompute 1/dt^4 */
    double dtinv = 1.0 / (dt*dt*dt*dt);

    for (PetscInt i=0; i<localsize_u; i++) {
        PetscInt vecID_re = i;
        PetscInt vecID_im = i + localsize_u;

        // first term
        if (n > 1) {
          // if (i==0) printf("DPDM BWD, update %d from on first f(%d %d %d) \n", n, n-2, n-1, n);

            double tmp1 = xm2ptr[vecID_re]*xm2ptr[vecID_re] - 2.0*xm1ptr[vecID_re]*xm1ptr[vecID_re] + xptr[vecID_re]*xptr[vecID_re];
            double tmp2 = xm2ptr[vecID_im]*xm2ptr[vecID_im] - 2.0*xm1ptr[vecID_im]*xm1ptr[vecID_im] + xptr[vecID_im]*xptr[vecID_im];
            double pop = tmp1+tmp2;
            double dpdphi_re = 2.0*xptr[vecID_re];
            double dpdphi_im = 2.0*xptr[vecID_im];
            VecSetValue(xbar,vecID_re,  2.0 * pop * dpdphi_re * dtinv * Jbar, ADD_VALUES);
            VecSetValue(xbar,vecID_im,  2.0 * pop * dpdphi_im * dtinv * Jbar, ADD_VALUES);
        }

        // center term
        if (n > 0 && n < ntime) {
            // if (i==0) printf("DPDM BWD, update %d from on second f(%d %d %d) \n", n, n-1, n, n+1);
              double tmp1 = xm1ptr[vecID_re]*xm1ptr[vecID_re] - 2.0*xptr[vecID_re]*xptr[vecID_re] + xp1ptr[vecID_re]*xp1ptr[vecID_re];
              double tmp2 = xm1ptr[vecID_im]*xm1ptr[vecID_im] - 2.0*xptr[vecID_im]*xptr[vecID_im] + xp1ptr[vecID_im]*xp1ptr[vecID_im];
              double pop = tmp1 + tmp2;
              double dpdphi_re = 2.0*xptr[vecID_re];
              double dpdphi_im = 2.0*xptr[vecID_im];
              VecSetValue(xbar, vecID_re, - 4.0 * pop * dpdphi_re * dtinv * Jbar, ADD_VALUES);
              VecSetValue(xbar, vecID_im, - 4.0 * pop * dpdphi_im * dtinv * Jbar, ADD_VALUES);
        }
        
        // last term 
        if (n < ntime-1) {
            // if (i==0) printf("DPDM BWD, update %d from on third f(%d %d %d) \n", n, n, n+1, n+2);
              double tmp1 = xptr[vecID_re]*xptr[vecID_re] - 2.0*xp1ptr[vecID_re]*xp1ptr[vecID_re] + xp2ptr[vecID_re]*xp2ptr[vecID_re];
              double tmp2 = xptr[vecID_im]*xptr[vecID_im] - 2.0*xp1ptr[vecID_im]*xp1ptr[vecID_im] + xp2ptr[vecID_im]*xp2ptr[vecID_im];
              double pop = tmp1 + tmp2;
              double dpdphi_re = 2.0*xptr[vecID_re];
              double dpdphi_im = 2.0*xptr[vecID_im];
              VecSetValue(xbar, vecID_re, 2.0* pop * dpdphi_re * dtinv * Jbar, ADD_VALUES);
              VecSetValue(xbar, vecID_im, 2.0* pop * dpdphi_im * dtinv * Jbar, ADD_VALUES);
        }
    }

    VecRestoreArrayRead(x, &xptr);
    if (n > 0 )       VecRestoreArrayRead(xm1, &xm1ptr);
    if (n > 1)        VecRestoreArrayRead(xm2, &xm2ptr);
    if (n < ntime)    VecRestoreArrayRead(xp1, &xp1ptr);
    if (n < ntime-1)  VecRestoreArrayRead(xp2, &xp2ptr);
 
}

double TimeStepper::evalEnergy(double time){
  double pen = 0.0;

  /* Loop over oscillators */
  for (size_t iosc = 0; iosc < mastereq->getNOscillators(); iosc++) {
    double p,q;
    mastereq->getOscillator(iosc)->evalControl(time, &p, &q); 
    pen += p*p + q*q;
  }

  return pen;
}


void TimeStepper::evalEnergy_diff(double time, double Jbar, Vec redgrad){

  PetscInt col_shift = 0;
  double* grad_ptr;
  VecGetArray(redgrad, &grad_ptr);

  for (size_t iosc = 0; iosc < mastereq->getNOscillators(); iosc++){

    /* Reevaluate the controls to set pbar, qbar */
    double p,q;
    mastereq->getOscillator(iosc)->evalControl(time, &p, &q); 
    double pbar = Jbar * 2.0 * p;
    double qbar = Jbar * 2.0 * q;
    /* Derivative of evalControls */
    double* grad_for_this_oscillator = grad_ptr + col_shift;
    mastereq->getOscillator(iosc)->evalControl_diff(time, grad_for_this_oscillator, pbar, qbar);

    // Skip in gradient for next oscillator
    col_shift += mastereq->getOscillator(iosc)->getNParams();
  } 
  VecRestoreArray(redgrad, &grad_ptr);
}

void TimeStepper::evolveBWD(const double /*tstart*/, const double /*tstop*/, const Vec /*x_stop*/, Vec /*x_adj*/, Vec /*grad*/, bool /*compute_gradient*/){}

ExplEuler::ExplEuler(size_t ninit_local, MasterEq* mastereq_, int ntime_, double total_time_, Output* output_, bool storeFWD_) : TimeStepper(ninit_local, mastereq_, ntime_, total_time_, output_, storeFWD_) {
  MatCreateVecs(mastereq->getRHS(), &stage, NULL);
  VecZeroEntries(stage);
}

ExplEuler::~ExplEuler() {
  VecDestroy(&stage);
}

void ExplEuler::evolveFWD(const double tstart,const  double tstop, Vec x) {

  double dt = tstop - tstart;

   /* Compute A(tstart) */
  mastereq->assemble_RHS(tstart);
  Mat A = mastereq->getRHS(); 

  /* update x = x + hAx */
  MatMult(A, x, stage);
  VecAXPY(x, dt, stage);
}

void ExplEuler::evolveBWD(const double tstop,const  double tstart,const  Vec x, Vec x_adj, Vec grad, bool compute_gradient){
  double dt = tstop - tstart;

  /* Add to reduced gradient */
  if (compute_gradient) {
    mastereq->compute_dRHS_dParams(tstop, x, x_adj, dt, grad);
  }

  /* update x_adj = x_adj + hA^Tx_adj */
  mastereq->assemble_RHS(tstop);
  Mat A = mastereq->getRHS(); 
  MatMultTranspose(A, x_adj, stage);
  VecAXPY(x_adj, dt, stage);

}

ImplMidpoint::ImplMidpoint(size_t ninit_local, MasterEq* mastereq_, int ntime_, double total_time_, LinearSolverType linsolve_type_, int linsolve_maxiter_, Output* output_, bool storeFWD_) : TimeStepper(ninit_local, mastereq_, ntime_, total_time_, output_, storeFWD_) {

  /* Create and reset the intermediate vectors */
  MatCreateVecs(mastereq->getRHS(), &stage, NULL);
  VecDuplicate(stage, &stage_adj);
  VecDuplicate(stage, &rhs);
  VecDuplicate(stage, &rhs_adj);
  VecZeroEntries(stage);
  VecZeroEntries(stage_adj);
  VecZeroEntries(rhs);
  VecZeroEntries(rhs_adj);
  linsolve_type = linsolve_type_;
  linsolve_maxiter = linsolve_maxiter_;
  linsolve_reltol = 1.e-20;
  linsolve_abstol = 1.e-10;
  linsolve_iterstaken_avg = 0;
  linsolve_counter = 0;
  linsolve_error_avg = 0.0;

  if (linsolve_type == LinearSolverType::GMRES) {
    /* Create Petsc's linear solver */
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPGetPC(ksp, &preconditioner);
    PCSetType(preconditioner, PCNONE);
    KSPSetTolerances(ksp, linsolve_reltol, linsolve_abstol, PETSC_DEFAULT, linsolve_maxiter);
    KSPSetType(ksp, KSPGMRES);
    KSPSetOperators(ksp, mastereq->getRHS(), mastereq->getRHS());
    KSPSetFromOptions(ksp);
  }
  else {
    /* For Neumann iterations, allocate a temporary vector */
    MatCreateVecs(mastereq->getRHS(), &tmp, NULL);
    MatCreateVecs(mastereq->getRHS(), &err, NULL);
  }
}


ImplMidpoint::~ImplMidpoint(){

  /* Print linear solver statistics */
  if (linsolve_counter <= 0) linsolve_counter = 1;
  linsolve_iterstaken_avg = (int) linsolve_iterstaken_avg / linsolve_counter;
  linsolve_error_avg = linsolve_error_avg / linsolve_counter;
  // MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  // if (myrank == 0) printf("Linear solver type %d: Average iterations = %d, average error = %1.2e\n", linsolve_type, linsolve_iterstaken_avg, linsolve_error_avg);

  /* Free up Petsc's linear solver */
  if (linsolve_type == LinearSolverType::GMRES) {
    KSPDestroy(&ksp);
  } else {
    VecDestroy(&tmp);
    VecDestroy(&err);
  }

  /* Free up intermediate vectors */
  VecDestroy(&stage_adj);
  VecDestroy(&stage);
  VecDestroy(&rhs_adj);
  VecDestroy(&rhs);

}

void ImplMidpoint::evolveFWD(const double tstart,const  double tstop, Vec x) {

  /* Compute time step size */
  double dt = tstop - tstart;  

  /* Compute A(t_n+h/2) */
  mastereq->assemble_RHS( (tstart + tstop) / 2.0);
  Mat A = mastereq->getRHS(); 

  /* Compute rhs = A x */
  MatMult(A, x, rhs);

  /* Solve for the stage variable (I-dt/2 A) k1 = Ax */
  switch (linsolve_type) {
    case LinearSolverType::GMRES:
      /* Set up I-dt/2 A, then solve */
      MatScale(A, - dt/2.0);
      MatShift(A, 1.0);  
      KSPSolve(ksp, rhs, stage);

      /* Monitor error */
      double rnorm;
      PetscInt iters_taken;
      KSPGetResidualNorm(ksp, &rnorm);
      KSPGetIterationNumber(ksp, &iters_taken);
      // printf("Residual norm %d: %1.5e\n", iters_taken, rnorm);
      linsolve_iterstaken_avg += iters_taken;
      linsolve_error_avg += rnorm;
      if (rnorm > 1e-3)  {
        printf("WARNING: Linear solver residual norm: %1.5e\n", rnorm);
      }
 
      /* Revert the scaling and shifting if gmres solver */
      MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
      break;

    case LinearSolverType::NEUMANN:
      linsolve_iterstaken_avg += NeumannSolve(A, rhs, stage, dt/2.0, false);
      break;
  }
  linsolve_counter++;

  /* --- Update state x += dt * stage --- */
  VecAXPY(x, dt, stage);
}

void ImplMidpoint::evolveBWD(const double tstop, const double tstart, const Vec x, Vec x_adj, Vec grad, bool compute_gradient){
  Mat A;

  /* Compute time step size */
  double dt = tstop - tstart;
  double thalf = (tstart + tstop) / 2.0;

  /* Assemble RHS(t_1/2) */
  mastereq->assemble_RHS( (tstart + tstop) / 2.0);
  A = mastereq->getRHS();

  /* Get Ax_n for use in gradient */
  if (compute_gradient) {
    MatMult(A, x, rhs);
  }

  /* Solve for adjoint stage variable */
  switch (linsolve_type) {
    case LinearSolverType::GMRES:
      MatScale(A, - dt/2.0);
      MatShift(A, 1.0);  // WARNING: this can be very slow if some diagonal elements are missing.
      KSPSolveTranspose(ksp, x_adj, stage_adj);
      double rnorm;
      KSPGetResidualNorm(ksp, &rnorm);
      PetscInt iters_taken;
      KSPGetIterationNumber(ksp, &iters_taken);
      if (rnorm > 1e-3)  {
        printf("WARNING: Linear solver residual norm: %1.5e\n", rnorm);
      }
      break;

    case LinearSolverType::NEUMANN: 
      NeumannSolve(A, x_adj, stage_adj, dt/2.0, true);
      break;
  }

  // k_bar = h*k_bar 
  VecScale(stage_adj, dt);

  /* Add to reduced gradient */
  if (compute_gradient) {
    switch (linsolve_type) {
      case LinearSolverType::GMRES: 
        KSPSolve(ksp, rhs, stage);
        break;
      case LinearSolverType::NEUMANN:
        NeumannSolve(A, rhs, stage, dt/2.0, false);
        break;
    }
    VecAYPX(stage, dt / 2.0, x);
    mastereq->compute_dRHS_dParams(thalf, stage, stage_adj, 1.0, grad);
  }

  /* Revert changes to RHS from above, if gmres solver */
  A = mastereq->getRHS();
  if (linsolve_type == LinearSolverType::GMRES) {
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  }

  /* Update adjoint state x_adj += dt * A^Tstage_adj --- */
  MatMultTransposeAdd(A, stage_adj, x_adj, x_adj);

}


int ImplMidpoint::NeumannSolve(Mat A, Vec b, Vec y, double alpha, bool transpose){

  double errnorm, errnorm0;

  // Initialize y = b
  VecCopy(b, y);

  int iter;
  for (iter = 0; iter < linsolve_maxiter; iter++) {
    VecCopy(y, err);

    // y = b + alpha * A *  y
    if (!transpose) MatMult(A, y, tmp);
    else            MatMultTranspose(A, y, tmp);
    VecAXPBYPCZ(y, 1.0, alpha, 0.0, b, tmp);

    /* Error approximation  */
    VecAXPY(err, -1.0, y); // err = yprev - y 
    VecNorm(err, NORM_2, &errnorm);

    /* Stopping criteria */
    if (iter == 0) errnorm0 = errnorm;
    if (errnorm < linsolve_abstol) break;
    if (errnorm / errnorm0 < linsolve_reltol) break;
  }

  // printf("Neumann error: %1.14e\n", errnorm);
  linsolve_error_avg += errnorm;

  return iter;
}



CompositionalImplMidpoint::CompositionalImplMidpoint(size_t ninit_local, int order_, MasterEq* mastereq_, int ntime_, double total_time_, LinearSolverType linsolve_type_, int linsolve_maxiter_, Output* output_, bool storeFWD_): ImplMidpoint(ninit_local, mastereq_, ntime_, total_time_, linsolve_type_, linsolve_maxiter_, output_, storeFWD_) {

  order = order_;

  // coefficients for order 8, stages s=15
  gamma.clear();
  if (order == 8){
    gamma.push_back(0.74167036435061295344822780);
    gamma.push_back(-0.40910082580003159399730010);
    gamma.push_back(0.19075471029623837995387626);
    gamma.push_back(-0.57386247111608226665638773);
    gamma.push_back(0.29906418130365592384446354);
    gamma.push_back(0.33462491824529818378495798);
    gamma.push_back(0.31529309239676659663205666);
    gamma.push_back(-0.79688793935291635401978884); // 8
    gamma.push_back(0.31529309239676659663205666);
    gamma.push_back(0.33462491824529818378495798);
    gamma.push_back(0.29906418130365592384446354);
    gamma.push_back(-0.57386247111608226665638773);
    gamma.push_back(0.19075471029623837995387626);
    gamma.push_back(-0.40910082580003159399730010);
    gamma.push_back(0.74167036435061295344822780);
  } else if (order == 4) {
    gamma.push_back(1./(2. - pow(2., 1./3.)));
    gamma.push_back(- pow(2., 1./3.)*gamma[0] );
    gamma.push_back(1./(2. - pow(2., 1./3.)));
  }

  if (mpirank_world == 0) printf("Timestepper: Compositional Impl. Midpoint, order %d, %lu stages\n", order, gamma.size());

  // Allocate storage of stages for backward process 
  PetscInt globalsize = 2 * mastereq->getDim(); 
  PetscInt localsize = globalsize / mpisize_petsc;  // Local vector per processor
  for (size_t i = 0; i <gamma.size(); i++) {
    Vec state;
    VecCreate(PETSC_COMM_WORLD, &state);
    VecSetSizes(state,localsize,globalsize);
    VecSetFromOptions(state);
    x_stage.push_back(state);
  }
  VecCreate(PETSC_COMM_WORLD, &aux);
  VecSetSizes(aux,localsize,globalsize);
  VecSetFromOptions(aux);
}

CompositionalImplMidpoint::~CompositionalImplMidpoint(){
  for (size_t i = 0; i <gamma.size(); i++) {
    VecDestroy(&(x_stage[i]));
  }
  VecDestroy(&aux);
}


void CompositionalImplMidpoint::evolveFWD(const double tstart,const  double tstop, Vec x) {

  double dt = tstop - tstart;
  double tcurr = tstart;

  // Loop over stages
  for (size_t istage = 0; istage < gamma.size(); istage++) {
    // time-step size and tstart,tstop for compositional step
    double dt_stage = gamma[istage] * dt;

    // Evolve 'tcurr -> tcurr + gamma*dt' using ImpliMidpointrule
    ImplMidpoint::evolveFWD(tcurr, tcurr + dt_stage, x);

    // Update current time
    tcurr = tcurr + dt_stage;
  }
  assert(fabs(tcurr - tstop) < 1e-12);

}

void CompositionalImplMidpoint::evolveBWD(const double tstop, const double tstart, const Vec x, Vec x_adj, Vec grad, bool compute_gradient){
  
  double dt = tstop - tstart;

  // Run forward again to store the (primal) stages
  double tcurr = tstart;
  VecCopy(x, aux);
  for (size_t istage = 0; istage < gamma.size(); istage++) {
    VecCopy(aux, x_stage[istage]);
    double dt_stage = gamma[istage] * dt;
    ImplMidpoint::evolveFWD(tcurr, tcurr + dt_stage, aux);
    tcurr = tcurr + dt_stage;
  }
  assert(fabs(tcurr - tstop) < 1e-12);

  // Run backwards while updating adjoint and gradient
  for (int istage = gamma.size()-1; istage >=0; istage--){
    double dt_stage = gamma[istage] * dt;
    ImplMidpoint::evolveBWD(tcurr, tcurr-dt_stage, x_stage[istage], x_adj, grad, compute_gradient);
    tcurr = tcurr - gamma[istage]*dt;
  }
  assert(fabs(tcurr - tstart) < 1e-12);
}

PetscTS::PetscTS(size_t ninit_local, MasterEq* mastereq_, int ntime_, double total_time_, Output* output_, bool storeFWD_) : TimeStepper(ninit_local, mastereq_, ntime_, total_time_, output_, storeFWD_) {

  ninit_pool = ninit_local; 
  ts_pool.resize(ninit_pool, nullptr);
  q_pool.resize(ninit_pool, nullptr);

  // Prepare a shared dRHSdp MatShell used by all TS objects.
  PetscInt nstate_global, nparam_global;
  VecGetSize(x, &nstate_global);
  VecGetSize(redgrad, &nparam_global);
  MatCreateShell(PETSC_COMM_WORLD, PETSC_DECIDE, nparam_global, nstate_global, nparam_global, this, &dRHSdp);
  MatShellSetOperation(dRHSdp, MATOP_MULT_TRANSPOSE, (void(*)(void)) computedRHSdp);

  // TSSetCostGradients is configured with numcost=1, so quadrature Jacobians
  // must provide a single combined running-cost channel for adjoint.
  const PetscInt ncost_terms = 1;
  MatCreateDense(PETSC_COMM_SELF, PETSC_DECIDE, PETSC_DECIDE, ncost_terms, nstate_global, NULL, &dIntegralCostdY);
  MatSetUp(dIntegralCostdY);
  MatCreateDense(PETSC_COMM_SELF, PETSC_DECIDE, PETSC_DECIDE, ncost_terms, nparam_global, NULL, &dIntegralCostdP);
  MatSetUp(dIntegralCostdP);

  dRHSdp_time = 0.0;
  adj_scale_leakage = 0.0;
  adj_scale_weightedcost = 0.0;
  adj_scale_energy = 0.0;

  // Helper function to create TS objects
  auto configureTS = [&](TS tsi) {
    TSSetProblemType(tsi, TS_LINEAR);

    /* Explicit RK */
    TSSetType(tsi, TSRK);
    TSRKSetType(tsi, TSRK5DP);

    // Pass the RHS function and Jacobian to Petsc.
    TSSetRHSFunction(tsi, NULL, TSComputeRHSFunctionLinear, mastereq->getRHSctx());
    TSSetRHSJacobian(tsi, mastereq->getRHS(), mastereq->getRHS(), RHSMatrixUpdate, mastereq->getRHSctx());
    TSSetRHSJacobianP(tsi, dRHSdp, dRHSdpMatrixUpdate, this);

    // Set time domain and adaptivity.
    TSSetTime(tsi, 0.0);
    TSSetMaxTime(tsi, total_time);
    TSSetExactFinalTime(tsi, TS_EXACTFINALTIME_MATCHSTEP);
    TSAdapt adapt;
    TSGetAdapt(tsi, &adapt);
    TSAdaptSetType(adapt, TSADAPTBASIC);
    // TSAdaptSetType(adapt, TSADAPTNONE);
    // TSSetTimeStep(tsi, total_time / ntime);
    double atol = 1e-7;
    double rtol = 1e-7;
    TSSetTolerances(tsi, atol, NULL, rtol, NULL);

    TSMonitorSet(tsi, PetscTS::monitorTrajectory, this, NULL);
    TSSetFromOptions(tsi);

    // Enable in-memory trajectory storage for adjoint solves.
    TSSetSaveTrajectory(tsi);
    TSTrajectory tj;
    TSGetTrajectory(tsi, &tj);
    TSTrajectorySetType(tj, tsi, TSTRAJECTORYMEMORY);
  };

  // Create a TS object for each initial condition. 
  for (int i = 0; i < ninit_pool; i++) {
    TSCreate(PETSC_COMM_SELF, &ts_pool[i]);
    configureTS(ts_pool[i]);

    // Attach quadrature integrator and a integral state per TS instance.
    TS ts_quad_i;
    TSCreateQuadratureTS(ts_pool[i], PETSC_TRUE, &ts_quad_i);
    TSSetRHSFunction(ts_quad_i, NULL, IntegralCosts, this);
    TSSetRHSJacobian(ts_quad_i, dIntegralCostdY, dIntegralCostdY, dIntegralCostdYUpdate, this);
    TSSetRHSJacobianP(ts_quad_i, dIntegralCostdP, dIntegralCostdPUpdate, this);
    VecCreate(PETSC_COMM_SELF, &q_pool[i]);
    VecSetSizes(q_pool[i], PETSC_DECIDE, 3);
    VecSetFromOptions(q_pool[i]);
    VecSet(q_pool[i], 0.0);
    TSSetSolution(ts_quad_i, q_pool[i]);
  }

  // Keep existing member as alias to first TS instance.
  ts = ts_pool[0];

  // Create an internal TS gradient 
  VecCreate(PETSC_COMM_SELF, &redgrad_ts);
  VecSetSizes(redgrad_ts, PETSC_DECIDE, nparam_global);
  VecSetFromOptions(redgrad_ts);
  VecZeroEntries(redgrad_ts);

  // // Gradient of quadrature register callbacks for r, dr/dy, dr/dp
  // TSSetCostIntegrand(ts, 1, q, IntegralCosts, DRDYFunction,DRDPFunction,ctx);

}

PetscTS::~PetscTS() {
  for (auto &qi : q_pool) {
    if (qi) VecDestroy(&qi);
  }
  for (auto &tsi : ts_pool) {
    if (tsi) TSDestroy(&tsi);
  }
  MatDestroy(&dRHSdp);
  MatDestroy(&dIntegralCostdY);
  MatDestroy(&dIntegralCostdP);
  VecDestroy(&redgrad_ts);
}

Vec PetscTS::solveODE(int initid, int iinit_local, Vec rho_t0){
  // Grab the timestepper for this initial condit
  const int iinit = iinit_local;;
  TS ts_run = getTSForInit(iinit_local);
  Vec q_run = q_pool[iinit];

  /* Prepare storage for trajectory output data */
  if (writeTrajectoryDataFiles) {
    output->resetTrajectoryData(initid);
  }

  /* Reset the time stepper */
  TSSetTime(ts_run, 0.0);
  TSSetStepNumber(ts_run, 0);
  TSSetTimeStep(ts_run, 0.1);  // This is Petsc's default initial guess. Will be adpted during TSSolve().

  /* Reset integral terms for this initial condition. */
  VecSet(q_run, 0.0);
  
  /* Set initial condition for timestepping */
  VecCopy(rho_t0, x);
  TSSetSolution(ts_run, x);

  // Reset/setup trajectory for this forward solve after solution is known.
  TSResetTrajectory(ts_run);

  /* Solve the ODE */
  TSSolve(ts_run, x);

  /* Store integral cost terms */
  const PetscScalar *terms;
  VecGetArrayRead(q_run, &terms);
  leakage_integral = terms[0] / total_time;
  weightedcost_integral = terms[1] / total_time;
  energy_integral = terms[2] / total_time;
  VecRestoreArrayRead(q_run, &terms);

  /* Write trajectory data files. Data was collected in monitorTrajectory() during TSSolve. */
  PetscInt nsteps;
  TSGetStepNumber(ts_run, &nsteps);
  if (writeTrajectoryDataFiles) {
    output->writeTrajectoryData();
  }

  return x;
}

void PetscTS::solveAdjointODE(int iinit_local, Vec rho_t0_bar, Vec finalstate, double Jbar_leakage, double Jbar_weightedcost, double Jbar_dpdm, double Jbar_energy) {
  TS ts_run = getTSForInit(iinit_local);
  (void)Jbar_dpdm;

  // Scaling.
  adj_scale_leakage = Jbar_leakage / total_time;
  adj_scale_weightedcost = Jbar_weightedcost / total_time;
  adj_scale_energy = Jbar_energy / total_time;

  /* Build terminal adjoint condition lambda(T) and terminal parameter gradient mu(T). */
  VecCopy(finalstate, xprimal);
  TSSetSolution(ts_run, xprimal);

  VecCopy(rho_t0_bar, xadj);
  VecZeroEntries(redgrad_ts);
  TSSetCostGradients(ts_run, 1, &xadj, &redgrad_ts);

  // backward solve
  TSAdjointSolve(ts_run);

  // Copy gradient to timestepper:
  VecCopy(redgrad_ts, redgrad);

}


PetscErrorCode PetscTS::RHSMatrixUpdate(TS, PetscReal t, Vec, Mat, Mat, void *ptr){
  MatShellCtx *ctx = (MatShellCtx *)ptr;
  if (!ctx->assembled || fabs(t - ctx->time) > 1e-8) {
    ctx->mastereq->assemble_RHS(t);
  }
  return 0;
};


PetscErrorCode PetscTS::dRHSdpMatrixUpdate(TS, PetscReal t, Vec xstate, Mat, void *ptr){

  PetscTS *self = static_cast<PetscTS *>(ptr);
  self->dRHSdp_time = t;
  VecCopy(xstate, self->xprimal);

  return 0;
}


PetscErrorCode PetscTS::computedRHSdp(Mat A, Vec xbar, Vec grad){
  PetscTS *self;
  MatShellGetContext(A, (void**)&self);

  VecZeroEntries(grad); // Need to reset here because compute_dRHS_dParams adds to grad. 
  self->mastereq->compute_dRHS_dParams(self->dRHSdp_time, self->xprimal, xbar, 1.0, grad);

  return 0;
}

PetscErrorCode PetscTS::dIntegralCostdYUpdate(TS, PetscReal t, Vec xstate, Mat, Mat, void *ptr){
  PetscTS *self = static_cast<PetscTS *>(ptr);
  self->dRHSdp_time = t;
  VecCopy(xstate, self->xprimal);

  Mat dRdy = self->dIntegralCostdY;
  MatZeroEntries(dRdy);

  Vec row;
  VecDuplicate(self->xprimal, &row);

  PetscInt ilow, iupp;
  VecGetOwnershipRange(row, &ilow, &iupp);
  const PetscScalar *vals = NULL;

  // Build a single combined derivative row for all active running costs.
  VecZeroEntries(row);
  if (self->eval_leakage) {
    self->evalLeakage_diff(self->xprimal, row, self->adj_scale_leakage);
  }
  if (self->eval_weightedcost) {
    self->evalWeightedCost_diff(self->dRHSdp_time, self->xprimal, row, self->adj_scale_weightedcost);
  }
  VecGetArrayRead(row, &vals);
  for (PetscInt j = ilow; j < iupp; ++j) {
    MatSetValue(dRdy, 0, j, vals[j - ilow], INSERT_VALUES);
  }
  VecRestoreArrayRead(row, &vals);

  VecDestroy(&row);

  MatAssemblyBegin(dRdy, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(dRdy, MAT_FINAL_ASSEMBLY);
  return 0;
}

PetscErrorCode PetscTS::dIntegralCostdPUpdate(TS, PetscReal t, Vec xstate, Mat, void *ptr){
  PetscTS *self = static_cast<PetscTS *>(ptr);
  self->dRHSdp_time = t;
  VecCopy(xstate, self->xprimal);

  MatZeroEntries(self->dIntegralCostdP);

  Vec prow;
  VecDuplicate(self->redgrad_ts, &prow);
  VecZeroEntries(prow);
  if (self->eval_energy) {
    self->evalEnergy_diff(self->dRHSdp_time, self->adj_scale_energy, prow);
  }

  PetscInt ilow, iupp;
  VecGetOwnershipRange(prow, &ilow, &iupp);
  const PetscScalar *vals = NULL;
  VecGetArrayRead(prow, &vals);
  for (PetscInt j = ilow; j < iupp; ++j) {
    MatSetValue(self->dIntegralCostdP, 0, j, vals[j - ilow], INSERT_VALUES);
  }
  VecRestoreArrayRead(prow, &vals);
  VecDestroy(&prow);

  MatAssemblyBegin(self->dIntegralCostdP, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(self->dIntegralCostdP, MAT_FINAL_ASSEMBLY);
  return 0;
}



PetscErrorCode PetscTS::IntegralCosts(TS, PetscReal t, Vec x, Vec F, void *ctx){
  PetscTS *self = static_cast<PetscTS *>(ctx); // The base Timestepper class.
  PetscScalar *integral_terms;
  VecGetArray(F, &integral_terms);

  integral_terms[0] = 0.0;
  integral_terms[1] = 0.0;
  integral_terms[2] = 0.0;

  /* Evaluate the cost integrand at time t, given state x. Store result in F. */
  if (self->eval_leakage)       integral_terms[0] = self->evalLeakage(x);
  if (self->eval_weightedcost)  integral_terms[1] = self->evalWeightedCost(t, x);
  if (self->eval_energy)        integral_terms[2] = self->evalEnergy(t);
  // TODO: DPDM integral term

  VecRestoreArray(F, &integral_terms);
  return 0;
}

PetscErrorCode PetscTS::monitorTrajectory(TS ts, PetscInt step, PetscReal time, Vec state, void *ctx){
  (void)ts;
  PetscTS *self = static_cast<PetscTS *>(ctx); // The base Timestepper class.

  // evaluate trajectory output 
  if (self->writeTrajectoryDataFiles) {
    self->output->evalTrajectoryData(step, time, state, self->mastereq);
  }

  return 0;
}