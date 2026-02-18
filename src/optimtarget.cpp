#include "optimtarget.hpp"
#include "config.hpp"
#include "defs.hpp"
#include <string>

OptimTarget::OptimTarget(){
  dim = 0;
  dim_rho = 0;
  dim_ess = 0;
  noscillators = 0;
  target_type = TargetType::NONE;
  objective_type = ObjectiveType::JTRACE;
  decoherence_type = DecoherenceType::NONE;
  targetgate = NULL;
  purity_rho0 = 1.0;
  purestateID = -1;
  targetstate = NULL;
  mpisize_petsc=0;
  mpirank_petsc=0;
}


OptimTarget::OptimTarget(const Config& config, MasterEq* mastereq, double total_time, Vec rho_t0, bool quietmode_) : OptimTarget() {

  // initialize
  dim = mastereq->getDim();
  dim_rho = mastereq->getDimRho();
  dim_ess = mastereq->getDimEss();
  quietmode = quietmode_;
  decoherence_type = mastereq->decoherence_type;
  MPI_Comm_size(PETSC_COMM_WORLD, &mpisize_petsc);
  MPI_Comm_rank(PETSC_COMM_WORLD, &mpirank_petsc);
  int mpirank_world;
  int mpisize_petsc;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank_world);
  MPI_Comm_size(PETSC_COMM_WORLD, &mpisize_petsc);
  // Set local sizes of subvectors u,v in state x=[u,v]
  localsize_u = dim / mpisize_petsc; 
  ilow = mpirank_petsc * localsize_u;
  iupp = ilow + localsize_u;         

  /* Get initial condition type and IDs */
  initcond = config.getInitialCondition();

  /* Prepare initial state rho_t0 if PRODUCT_STATE or FROMFILE or ENSEMBLE initialization. Otherwise they are set within prepareInitialState during evalF. */
  if (initcond.type == InitialConditionType::PRODUCT_STATE) {
    const auto& level_indices = initcond.levels.value();
    // Find the id within the global composite system 
    PetscInt diag_id = 0;
    for (size_t k=0; k < level_indices.size(); k++) {
      PetscInt dim_postkron = 1;
      for (size_t m=k+1; m < level_indices.size(); m++) {
        dim_postkron *= mastereq->getOscillator(m)->getNLevels();
      }
      diag_id += level_indices[k] * dim_postkron;
    }
    // Vectorize if lindblad solver
    PetscInt vec_id = diag_id;
    if (decoherence_type != DecoherenceType::NONE) vec_id = getVecID( diag_id, diag_id, dim_rho); 
    // Set 1.0 on the processor who owns this index
    if (ilow <= vec_id && vec_id < iupp) {
      PetscInt id_global_x =  vec_id + mpirank_petsc*localsize_u; // Global index of u_i in x=[u,v]
      VecSetValue(rho_t0, id_global_x, 1.0, INSERT_VALUES);
    }
  }
  else if (initcond.type == InitialConditionType::FROMFILE) {
    /* Read initial condition from file */
    int nelems = 0;
    // Always read the state vector; density matrix will be formed here for Lindblad.
    nelems = 2 * dim_ess;
    double * vec = new double[nelems];
    if (mpirank_world == 0) {
      std::string filename = initcond.filename.value();
      read_vector(filename.c_str(), vec, nelems, quietmode);
    }
    MPI_Bcast(vec, nelems, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (decoherence_type != DecoherenceType::NONE) { // Lindblad solver, build density matrix from state vector
      for (PetscInt i = 0; i < dim_ess; i++) {
        double psi_i_re = vec[i];
        double psi_i_im = vec[i + dim_ess];
        for (PetscInt j = 0; j < dim_ess; j++) {
          double psi_j_re = vec[j];
          double psi_j_im = vec[j + dim_ess];
          // rho_ij = psi_i * conj(psi_j)
          double rho_re = psi_i_re * psi_j_re + psi_i_im * psi_j_im;
          double rho_im = psi_i_im * psi_j_re - psi_i_re * psi_j_im;
          PetscInt k = i;
          PetscInt l = j;
          if (dim_ess*dim_ess < mastereq->getDim()) {
            k = mapEssToFull(k, mastereq->nlevels, mastereq->nessential);
            l = mapEssToFull(l, mastereq->nlevels, mastereq->nessential);
          }
          PetscInt elemid = getVecID(k, l, dim_rho);
          if (ilow <= elemid && elemid < iupp) {
            PetscInt id_global_x =  elemid + mpirank_petsc*localsize_u; // Global index of u_i in x=[u,v]
            VecSetValue(rho_t0, id_global_x, rho_re, INSERT_VALUES);  // RealPart
            VecSetValue(rho_t0, id_global_x + localsize_u, rho_im, INSERT_VALUES); // Imaginary Part
          }
        }
      }
    } else { // Schroedinger solver, fill vector
      for (PetscInt i = 0; i < dim_ess; i++) {
        PetscInt k = i;
        if (dim_ess < mastereq->getDim()) 
          k = mapEssToFull(i, mastereq->nlevels, mastereq->nessential);
        PetscInt elemid = k;
        if (ilow <= elemid && elemid < iupp) {
          PetscInt id_global_x =  elemid + mpirank_petsc*localsize_u; // Global index of u_i in x=[u,v]
          VecSetValue(rho_t0, id_global_x, vec[i], INSERT_VALUES);  // RealPart
          VecSetValue(rho_t0, id_global_x + localsize_u, vec[i + dim_ess], INSERT_VALUES); // Imaginary Part
        }
      }
    }
    delete [] vec;
  } else if (initcond.type == InitialConditionType::ENSEMBLE) {
    const auto& osc_IDs = initcond.subsystem.value();

    // get dimension of subsystems defined by ensemble_init.level_indices, as well as the one before and after. Span in essential levels only.
    PetscInt dimpost = 1;
    PetscInt dimsub = 1;
    for (size_t i=0; i<mastereq->getNOscillators(); i++){
      if (osc_IDs[0] <= i && i <= osc_IDs[osc_IDs.size()-1]) dimsub *= mastereq->nessential[i];
      else dimpost *= mastereq->nessential[i];
    }
    PetscInt dimrho = mastereq->getDimRho();
    PetscInt dimrhoess = mastereq->getDimEss();
    // Loop over ensemble state elements in essential level dimensions of the subsystem defined by the initcond_ids:
    for (PetscInt i=0; i < dimsub; i++){
      for (PetscInt j=i; j < dimsub; j++){
        PetscInt ifull = i * dimpost; // account for the system behind
        PetscInt jfull = j * dimpost;
        if (dimrhoess < dimrho) ifull = mapEssToFull(ifull, mastereq->nlevels, mastereq->nessential);
        if (dimrhoess < dimrho) jfull = mapEssToFull(jfull, mastereq->nlevels, mastereq->nessential);
        // printf(" i=%d j=%d ifull %d, jfull %d\n", i, j, ifull, jfull);
        if (i == j) { 
          // diagonal element: 1/N_sub
          PetscInt elemid = getVecID(ifull, jfull, dimrho);
          if (ilow <= elemid && elemid < iupp) {
            PetscInt id_global_x =  elemid + mpirank_petsc*localsize_u; // Global index of u_i in x=[u,v]
            VecSetValue(rho_t0, id_global_x, 1./dimsub, INSERT_VALUES);
          }
        } else {
          // upper diagonal (0.5 + 0.5*i) / (N_sub^2)
          PetscInt elemid = getVecID(ifull, jfull, dimrho);
          if (ilow <= elemid && elemid < iupp) {
            PetscInt id_global_x =  elemid + mpirank_petsc*localsize_u; // Global index of u_i in x=[u,v]
            VecSetValue(rho_t0, id_global_x, 0.5/(dimsub*dimsub), INSERT_VALUES);
            VecSetValue(rho_t0, id_global_x + localsize_u, 0.5/(dimsub*dimsub), INSERT_VALUES);
          }
          // lower diagonal (0.5 - 0.5*i) / (N_sub^2)
          elemid = getVecID(jfull, ifull, dimrho);
          if (ilow <= elemid && elemid < iupp) {
            PetscInt id_global_x =  elemid + mpirank_petsc*localsize_u; // Global index of u_i in x=[u,v]
            VecSetValue(rho_t0, id_global_x,  0.5/(dimsub*dimsub), INSERT_VALUES);
            VecSetValue(rho_t0, id_global_x + localsize_u, -0.5/(dimsub*dimsub), INSERT_VALUES);
          }
        } 
      }
    }
  }
  VecAssemblyBegin(rho_t0); VecAssemblyEnd(rho_t0);

  /* Allocate storage for the target state */
  VecCreate(PETSC_COMM_WORLD, &targetstate); 
  PetscInt globalsize = 2 * mastereq->getDim();  // Global state vector: 2 for real and imaginary part
  PetscInt localsize = globalsize / mpisize_petsc;  // Local vector per processor
  VecSetSizes(targetstate,localsize,globalsize);
  VecSetFromOptions(targetstate);
  VecZeroEntries(targetstate);

  /* Store objective function type */
  objective_type = config.getOptimObjective();

  /* Get target type */  
  purestateID = -1;
  const auto& target = config.getOptimTarget();
  target_type = target.type;

  if (target_type == TargetType::GATE) {
    // Get optional gate rotation frequencies
    const std::vector<double>& gate_rot_freq = target.gate_rot_freq.value_or(std::vector<double>(mastereq->getNOscillators(), 0.0));

    /* Initialize the targetgate, either from file or using default set of gates */
    targetgate = initTargetGate(target.gate_type.value(), target.filename.value_or(""), mastereq->nlevels, mastereq->nessential, total_time, decoherence_type, gate_rot_freq, quietmode);

  } else if (target_type == TargetType::STATE) {
    // Initialize a the target state, either pure state prep (store only the purestateID) or read state from file
    if (target.levels.has_value()){
      purestateID = 0;
      const std::vector<size_t>& purestate_levels = target.levels.value();
      for (size_t i=0; i < mastereq->getNOscillators(); i++) {
        purestateID += purestate_levels[i] * mastereq->getOscillator(i)->dim_postOsc;
      }
    } else if (target.filename.has_value()) {
      /* Read the target state from file into vec */
      auto target_filename = target.filename.value();
      PetscInt nelems = 0;
      if (mastereq->decoherence_type != DecoherenceType::NONE) nelems = 2*dim_ess*dim_ess;
      else nelems = 2 * dim_ess;
      double* vec = new double[nelems];
      if (mpirank_world == 0) 
        read_vector(target_filename.c_str(), vec, nelems, quietmode);
      MPI_Bcast(vec, nelems, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      if (decoherence_type != DecoherenceType::NONE) { // Lindblad solver, fill density matrix
        for (PetscInt i = 0; i < dim_ess*dim_ess; i++) {
          PetscInt k = i % dim_ess;
          PetscInt j = i / dim_ess;
          if (dim_ess*dim_ess < mastereq->getDim()) {
            k = mapEssToFull(k, mastereq->nlevels, mastereq->nessential);
            j = mapEssToFull(j, mastereq->nlevels, mastereq->nessential);
          }
          PetscInt elemid = getVecID(k,j,dim_rho);
          if (ilow <= elemid && elemid < iupp) {
            PetscInt id_global_x =  elemid + mpirank_petsc*localsize_u; // Global index of u_i in x=[u,v]
            VecSetValue(targetstate, id_global_x, vec[i],       INSERT_VALUES); // RealPart
            VecSetValue(targetstate, id_global_x + localsize_u, vec[i + dim_ess*dim_ess], INSERT_VALUES); // Imaginary Part
          }
        }
      } else {  // Schroedinger solver, fill vector
        for (int i = 0; i < dim_ess; i++) {
          int k = i;
          if (dim_ess < mastereq->getDim()) 
            k = mapEssToFull(i, mastereq->nlevels, mastereq->nessential);
          PetscInt elemid = k;
          if (ilow <= elemid && elemid < iupp) {
            PetscInt id_global_x =  elemid + mpirank_petsc*localsize_u; // Global index of u_i in x=[u,v]
            VecSetValue(targetstate, id_global_x, vec[i], INSERT_VALUES);        // RealPart
            VecSetValue(targetstate, id_global_x + localsize_u, vec[i + dim_ess], INSERT_VALUES); // Imaginary Part
          }
        }
      }
      VecAssemblyBegin(targetstate); VecAssemblyEnd(targetstate);
      delete [] vec;
    }
  }

  /* Allocate an auxiliary vec needed for evaluating the frobenius norm */
  if (objective_type == ObjectiveType::JFROBENIUS) {
    VecCreate(PETSC_COMM_WORLD, &aux); 
    PetscInt globalsize = 2 * mastereq->getDim();  // 2 for real and imaginary part
    PetscInt localsize = globalsize / mpisize_petsc;  // Local vector per processor
    VecSetSizes(aux,localsize, globalsize);
    VecSetFromOptions(aux);
  }
}

OptimTarget::~OptimTarget(){
  if (objective_type == ObjectiveType::JFROBENIUS) VecDestroy(&aux);
  VecDestroy(&targetstate);
  delete targetgate;
}

double OptimTarget::FrobeniusDistance(const Vec state){
  // Frobenius distance F = 1/2 || targetstate - state ||^2_F  = 1/2 || vec(targetstate-state)||^2_2
  double norm;
  VecAYPX(aux, 0.0, targetstate);    // aux = targetstate
  VecAXPY(aux, -1.0, state);   // aux = targetstate - state
  VecNorm(aux, NORM_2, &norm);
  double J = norm * norm;

  return J;
}

void OptimTarget::FrobeniusDistance_diff(const Vec state, Vec statebar, const double Jbar){
  
  // Derivative of frobenius distance : statebar += 2 * (targetstate - state) * (-1) * Jbar 
  VecAXPY(statebar,  2.0*Jbar, state);
  VecAXPY(statebar, -2.0*Jbar, targetstate);  
}

void OptimTarget::HilbertSchmidtOverlap(const Vec state, const bool scalebypurity, double* HS_re_ptr, double* HS_im_ptr ){
  /* Lindblas solver: Tr(state * target^\dagger) = vec(target)^dagger * vec(state), will be real!
   * Schroedinger:    Tr(state * target^\dagger) = target^\dag * state, will be complex!*/

  // Reset
  double HS_re = 0.0;
  double HS_im = 0.0;

  /* Simplify computation if the target is PRODUCT_STATE, i.e. target = e_m or e_m * e_m^\dag */
  /* Tr(...) = phi_m if Schroedinger, or \rho_mm if Lindblad */
  if (target_type == TargetType::STATE && purestateID >= 0){

    // Vectorize pure state ID if Lindblad
    PetscInt idm = purestateID;
    if (decoherence_type != DecoherenceType::NONE) idm = getVecID(purestateID, purestateID, (PetscInt)sqrt(dim));

    // Get real and imag values from the processor who owns the subvector index.
    if (ilow <= idm && idm < iupp) {
      PetscInt id_global_x = idm + mpirank_petsc*localsize_u; // Global index of u_i in x=[u,v]
      VecGetValues(state, 1, &id_global_x, &HS_re);
      id_global_x += localsize_u; // Imaginary part
      VecGetValues(state, 1, &id_global_x, &HS_im); // Should be 0.0 if Lindblad!
    }
    if (decoherence_type != DecoherenceType::NONE) assert(fabs(HS_im) <= 1e-14);

    // Communicate over all petsc processors.
    double myre = HS_re;
    double myim = HS_im;
    MPI_Allreduce(&myre, &HS_re, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
    MPI_Allreduce(&myim, &HS_im, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

  } else { // Target is not of the form e_m (schroedinger) or e_m e_m^\dagger (lindblad).

    if (decoherence_type != DecoherenceType::NONE) // Lindblad solver. HS overlap is real!
      VecTDot(targetstate, state, &HS_re);  
    else {  // Schroedinger solver. target^\dagger * state
      // Get local data pointers
      const PetscScalar* target_ptr;
      const PetscScalar* state_ptr;
      VecGetArrayRead(targetstate, &target_ptr); 
      VecGetArrayRead(state, &state_ptr);
      for (PetscInt i=0; i<localsize_u; i++){
        PetscInt idre = i;
        PetscInt idim = i + localsize_u;
        HS_re +=  target_ptr[idre]*state_ptr[idre] + target_ptr[idim]*state_ptr[idim];
        HS_im += -target_ptr[idim]*state_ptr[idre] + target_ptr[idre]*state_ptr[idim];
      } 
      VecRestoreArrayRead(targetstate, &target_ptr);
      VecRestoreArrayRead(state, &state_ptr);
      // The above computation was local, so have to sum up here.
      double re=HS_re;
      double im=HS_im;
      MPI_Allreduce(&re, &HS_re, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
      MPI_Allreduce(&im, &HS_im, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
    }
  }

  // scale by purity Tr(rho(0)^2). 
  if (scalebypurity){ 
    HS_re = HS_re / purity_rho0;
  }

  // return
  *HS_re_ptr = HS_re;
  *HS_im_ptr = HS_im;
}

void OptimTarget::HilbertSchmidtOverlap_diff(Vec statebar, bool scalebypurity, const double HS_re_bar, const double HS_im_bar){

  double scale = 1.0;
  if (scalebypurity){ 
    scale = 1./purity_rho0;
  }

  // Simplified computation if target is product state
  if (target_type == TargetType::STATE && purestateID >= 0){
    PetscInt idm = purestateID;
    if (decoherence_type != DecoherenceType::NONE) idm = getVecID(purestateID, purestateID, (PetscInt)sqrt(dim));

    if (ilow <= idm && idm < iupp) {
      PetscInt id_global_x = idm + mpirank_petsc*localsize_u; // Global index of u_i in x=[u,v]
      VecSetValue(statebar, id_global_x, HS_re_bar*scale, ADD_VALUES);
      VecSetValue(statebar, id_global_x + localsize_u, HS_im_bar, ADD_VALUES);
    }

  } else { // Target is not of the form e_m or e_m*e_m^\dagger 

    if (decoherence_type != DecoherenceType::NONE)
      VecAXPY(statebar, HS_re_bar*scale, targetstate);
    else {
      const PetscScalar* target_ptr;
      PetscScalar* statebar_ptr;
      VecGetArrayRead(targetstate, &target_ptr); 
      VecGetArray(statebar, &statebar_ptr);
      for (PetscInt i=0; i<localsize_u; i++){
        PetscInt idre = i;
        PetscInt idim = i + localsize_u;
        statebar_ptr[idre] += target_ptr[idre] * HS_re_bar*scale  - target_ptr[idim] * HS_im_bar;
        statebar_ptr[idim] += target_ptr[idim] * HS_re_bar*scale  + target_ptr[idre] * HS_im_bar;
      }
      VecRestoreArrayRead(targetstate, &target_ptr);
      VecRestoreArray(statebar, &statebar_ptr);
    }
  }
}


int OptimTarget::prepareInitialState(const int iinit, const int ninit, const std::vector<size_t>& nlevels, const std::vector<size_t>& nessential, Vec rho0){

  PetscInt elemID;
  double val;
  PetscInt dim_post;
  int initID = 0;    // Output: ID for this initial condition */

  /* Conditionals over type of initial condition */
  if (initcond.type == InitialConditionType::PERFORMANCE) {
    /* Set up Input state psi = 1/sqrt(2N)*(Ones(N) + im*Ones(N)) or rho = psi*psi^\dag */
    VecZeroEntries(rho0);

    for (PetscInt i=0; i<dim_rho; i++){
      if (decoherence_type == DecoherenceType::NONE) {
        double val = 1./ sqrt(2.*dim_rho);
        if (ilow <= i && i < iupp) {
          PetscInt id_global_x =  i + mpirank_petsc*localsize_u; // Global index of u_i in x=[u,v]
          VecSetValue(rho0, id_global_x, val, INSERT_VALUES);
          VecSetValue(rho0, id_global_x + localsize_u, val, INSERT_VALUES);
        }
      } else {
        PetscInt elem_re = getVecID(i, i, dim_rho);
        double val = 1./ dim_rho;
        if (ilow <= elem_re && elem_re < iupp) {
          PetscInt id_global_x =  i + mpirank_petsc*localsize_u; // Global index of u_i in x=[u,v]
          VecSetValue(rho0, id_global_x, val, INSERT_VALUES);
        }
      }
    }
  } else if(initcond.type == InitialConditionType::FROMFILE) {
    /* Do nothing. Init cond is already stored */
  } else if(initcond.type == InitialConditionType::PRODUCT_STATE) {
    /* Do nothing. Init cond is already stored */
  } else if(initcond.type == InitialConditionType::ENSEMBLE) {
    /* Do nothing. Init cond is already stored */
  } else if (initcond.type == InitialConditionType::THREESTATES) {
    assert(decoherence_type != DecoherenceType::NONE);
    VecZeroEntries(rho0);

    /* Set the <iinit>'th initial state */
    if (iinit == 0) {
      // 1st initial state: rho(0)_IJ = 2(N-i)/(N(N+1)) Delta_IJ
      initID = 1;
      for (PetscInt i_full = 0; i_full<dim_rho; i_full++) {
        PetscInt diagID = getVecID(i_full,i_full,dim_rho);
        double val = 2.*(dim_rho - i_full) / (dim_rho * (dim_rho + 1));
        if (ilow <= diagID && diagID < iupp) {
          PetscInt id_global_x =  diagID + mpirank_petsc*localsize_u; // Global index of u_i in x=[u,v]
          VecSetValue(rho0, id_global_x, val, INSERT_VALUES);
        }
      }
    } else if (iinit == 1) {
      // 2nd initial state: rho(0)_IJ = 1/N
      initID = 2;
      for (PetscInt i_full = 0; i_full<dim_rho; i_full++) {
        for (PetscInt j_full = 0; j_full<dim_rho; j_full++) {
          double val = 1./dim_rho;
          PetscInt index = getVecID(i_full,j_full,dim_rho);
          if (ilow <= index && index < iupp) {
            PetscInt id_global_x =  index + mpirank_petsc*localsize_u;
            VecSetValue(rho0, id_global_x, val, INSERT_VALUES); 
          }
        }
      }
    } else if (iinit == 2) {
      // 3rd initial state: rho(0)_IJ = 1/N Delta_IJ
      initID = 3;
      for (PetscInt i_full = 0; i_full<dim_rho; i_full++) {
        PetscInt diagID = getVecID(i_full,i_full,dim_rho);
        double val = 1./ dim_rho;
        if (ilow <= diagID && diagID < iupp) {
          PetscInt id_global_x =  diagID + mpirank_petsc*localsize_u;
          VecSetValue(rho0, id_global_x, val, INSERT_VALUES);
        }
      }
    } else {
      printf("ERROR: Wrong initial condition setting! Should never happen.\n");
      exit(1);
    }
    VecAssemblyBegin(rho0); VecAssemblyEnd(rho0);
  } else if (initcond.type == InitialConditionType::NPLUSONE) {
    assert(decoherence_type != DecoherenceType::NONE);

    if (iinit < dim_rho) {// Diagonal e_j e_j^\dag
      VecZeroEntries(rho0);
      elemID = getVecID(iinit, iinit, dim_rho);
      val = 1.0;
      if (ilow <= elemID && elemID < iupp) {
        PetscInt id_global_x = elemID+ mpirank_petsc*localsize_u;
        VecSetValues(rho0, 1, &id_global_x, &val, INSERT_VALUES);
      }
    }
    else if (iinit == dim_rho) { // fully rotated 1/d*Ones(d)
      for (PetscInt i=0; i<dim_rho; i++){
        for (PetscInt j=0; j<dim_rho; j++){
          elemID = getVecID(i,j,dim_rho);
          val = 1.0 / dim_rho;
          if (ilow <= elemID && elemID < iupp) {
            PetscInt id_global_x = elemID + mpirank_petsc*localsize_u;
            VecSetValues(rho0, 1, &id_global_x, &val, INSERT_VALUES);
          }
        }
      }
    }
    else {
      printf("Wrong initial condition index. Should never happen!\n");
      exit(1);
    }
    initID = iinit;
    VecAssemblyBegin(rho0); VecAssemblyEnd(rho0);
  } else if (initcond.type == InitialConditionType::DIAGONAL) {
    const auto& initcond_IDs = initcond.subsystem.value();
    PetscInt diagelem;
    VecZeroEntries(rho0);

    /* Get dimension of partial system behind last oscillator ID (essential levels only) */
    dim_post = 1;
    for (size_t k = initcond_IDs[initcond_IDs.size()-1] + 1; k < nessential.size(); k++) {
      // dim_post *= getOscillator(k)->getNLevels();
      dim_post *= nessential[k];
    }

    /* Compute index of the nonzero element in rho_m(0) = E_pre \otimes |m><m| \otimes E_post */
    diagelem = iinit * dim_post;
    if (dim_ess < dim_rho)  diagelem = mapEssToFull(diagelem, nlevels, nessential);

    // Vectorize if Lindblad
    elemID = diagelem;
    if (decoherence_type != DecoherenceType::NONE) elemID = getVecID(diagelem, diagelem, dim_rho); 
    val = 1.0;
    if (ilow <= elemID && elemID < iupp) {
      PetscInt id_global_x =  elemID + mpirank_petsc*localsize_u; 
      VecSetValues(rho0, 1, &id_global_x, &val, INSERT_VALUES);
    }
    VecAssemblyBegin(rho0); VecAssemblyEnd(rho0);

    /* Set initial conditon ID */
    if (decoherence_type != DecoherenceType::NONE) initID = iinit * ninit + iinit;
    else initID = iinit;

  } else if (initcond.type == InitialConditionType::BASIS) {
    const auto& initcond_IDs = initcond.subsystem.value();

    assert(decoherence_type != DecoherenceType::NONE); // should never happen. For Schroedinger: BASIS equals DIAGONAL, and should go into the above switch case. 

    /* Reset the initial conditions */
    VecZeroEntries(rho0);

    /* Get dimension of partial system behind last oscillator ID (essential levels only) */
    dim_post = 1;
    for (size_t k = initcond_IDs[initcond_IDs.size()-1] + 1; k < nessential.size(); k++) {
      dim_post *= nessential[k];
    }

    /* Get index (k,j) of basis element B_{k,j} for this initial condition index iinit */
    PetscInt k, j;
    k = iinit % ( (int) sqrt(ninit) );
    j = iinit / ( (int) sqrt(ninit) );

    /* Set initial condition ID */
    initID = j * ( (int) sqrt(ninit)) + k;

    /* Set position in rho */
    k = k*dim_post;
    j = j*dim_post;
    if (dim_ess < dim_rho) { 
      k = mapEssToFull(k, nlevels, nessential);
      j = mapEssToFull(j, nlevels, nessential);
    }

    if (k == j) {
      /* B_{kk} = E_{kk} -> set only one element at (k,k) */
      elemID = getVecID(k, k, dim_rho); 
      double val = 1.0;
      if (ilow <= elemID && elemID < iupp) {
        PetscInt id_global_x =  elemID + mpirank_petsc*localsize_u; 
        VecSetValues(rho0, 1, &id_global_x, &val, INSERT_VALUES);
      }
    } else {
    //   /* B_{kj} contains four non-zeros, two per row */
      PetscInt* rows = new PetscInt[4];
      PetscScalar* vals = new PetscScalar[4];

      /* Get storage index of Re(x) */
      rows[0] = getVecID(k, k, dim_rho); // (k,k)
      rows[1] = getVecID(j, j, dim_rho); // (j,j)
      rows[2] = getVecID(k, j, dim_rho); // (k,j)
      rows[3] = getVecID(j, k, dim_rho); // (j,k)

      if (k < j) { // B_{kj} = 1/2(E_kk + E_jj) + 1/2(E_kj + E_jk)
        vals[0] = 0.5;
        vals[1] = 0.5;
        vals[2] = 0.5;
        vals[3] = 0.5;
        for (int i=0; i<4; i++) {
          if (ilow <= rows[i] && rows[i] < iupp) {
            PetscInt id_global_x =  rows[i]+ mpirank_petsc*localsize_u; 
            VecSetValues(rho0, 1, &id_global_x, &(vals[i]), INSERT_VALUES);
          }
        }
      } else {  // B_{kj} = 1/2(E_kk + E_jj) + i/2(E_jk - E_kj)
        vals[0] = 0.5;
        vals[1] = 0.5;
        for (int i=0; i<2; i++) {
          if (ilow <= rows[i] && rows[i] < iupp) {
            PetscInt id_global_x =  rows[i]+ mpirank_petsc*localsize_u; 
            VecSetValues(rho0, 1, &id_global_x, &(vals[i]), INSERT_VALUES);
          }
        }
        vals[2] = -0.5;
        vals[3] = 0.5;
        rows[2] = getVecID(k, j, dim_rho); // (k,j)
        rows[3] = getVecID(j, k, dim_rho); // (j,k)
        for (int i=2; i<4; i++) {
          if (ilow <= rows[i] && rows[i] < iupp) {
            PetscInt id_global_x =  rows[i]+ mpirank_petsc*localsize_u + localsize_u; 
            VecSetValues(rho0, 1, &id_global_x, &(vals[i]), INSERT_VALUES);
          }
        }
      }
      delete [] rows;
      delete [] vals;
    }

    /* Assemble rho0 */
    VecAssemblyBegin(rho0); VecAssemblyEnd(rho0);
  } else {
    printf("ERROR! Wrong initial condition type.\n This should never happen!\n");
    exit(1);
}

  return initID;
}


void OptimTarget::prepareTargetState(const Vec rho_t0){

  // If no target specified, set target state to zero.
  if (target_type == TargetType::NONE) VecZeroEntries(targetstate);

  // If gate optimization, apply the gate and store targetstate for later use. Else, do nothing.
  if (target_type == TargetType::GATE) targetgate->applyGate(rho_t0, targetstate);

  /* Compute and store the purity of rho(0), Tr(rho(0)^2), so that it can be used by JTrace (HS overlap) */
  VecNorm(rho_t0, NORM_2, &purity_rho0);
  purity_rho0 = purity_rho0 * purity_rho0;
}



void OptimTarget::evalJ(const Vec state, double* J_re_ptr, double* J_im_ptr){
  // Don't evaluate any objective function if the target type is NONE
  if (target_type == TargetType::NONE) {
    *J_re_ptr = 0.0;
    *J_im_ptr = 0.0;
    return;
  }

  double J_re = 0.0;
  double J_im = 0.0;
  double sum, rhoii, rhoii_re, rhoii_im, lambdai, norm;
  PetscInt dimsq;

  switch(objective_type) {

    /* J_Frob = 1/2 * || rho_target - rho(T)||^2_F  */
    case ObjectiveType::JFROBENIUS:

      // Simplify computation if target is a product state, otherwise, just evaluate the frobenius distance.
      if (target_type == TargetType::STATE && purestateID >= 0){ 
        // substract 1.0 from m-th diagonal element then take the vector norm 
        PetscInt diagID = purestateID;
        if (decoherence_type != DecoherenceType::NONE) diagID = getVecID(purestateID,purestateID,(PetscInt)sqrt(dim));
        if (ilow <= diagID && diagID < iupp) {
          PetscInt id_global_x = diagID + mpirank_petsc*localsize_u; 
          VecSetValue(state, id_global_x, -1.0, ADD_VALUES);
        }
        VecAssemblyBegin(state); VecAssemblyEnd(state);
        norm = 0.0;
        VecNorm(state, NORM_2, &norm);
        J_re = pow(norm, 2.0) / 2.0;
        if (ilow <= diagID && diagID < iupp) {
          PetscInt id_global_x = diagID + mpirank_petsc*localsize_u; 
          VecSetValue(state, id_global_x, +1.0, ADD_VALUES); // restore original state!
        }
        VecAssemblyBegin(state); VecAssemblyEnd(state);

      } else {
        J_re = FrobeniusDistance(state) / 2.0;
      }
      break;  // case Frobenius

    /* J_Trace:  1 / purity * Tr(state * target^\dagger)  =  HilbertSchmidtOverlap(target, state) is real if Lindblad, and complex if Schroedinger! */
    case ObjectiveType::JTRACE:

      HilbertSchmidtOverlap(state, true, &J_re, &J_im); // is real if Lindblad solver. 
      break; // case J_Trace

    /* J_Measure = Tr(O_m rho(T)) = \sum_i |i-m| rho_ii(T) if Lindblad and \sum_i |i-m| |phi_i(T)|^2  if Schroedinger */
    case ObjectiveType::JMEASURE:
      // Sanity check
      if (target_type != TargetType::STATE || purestateID < 0) {
        printf("Error: Wrong setting for objective function. Jmeasure can only be used for pure product state targets.\n");
        exit(1);
      }

      dimsq = dim;   // Schroedinger solver: dim = N
      if (decoherence_type != DecoherenceType::NONE) dimsq = (PetscInt)sqrt(dim); // Lindblad solver: dim = N^2

      // iterate over diagonal elements 
      sum = 0.0;
      for (PetscInt i=0; i<dimsq; i++){
        if (decoherence_type != DecoherenceType::NONE) {
          PetscInt diagID = getVecID(i,i,dimsq);
          rhoii = 0.0;
          if (ilow <= diagID && diagID < iupp) {
            PetscInt id_global_x =  diagID + mpirank_petsc*localsize_u;
            VecGetValues(state, 1, &id_global_x, &rhoii);
          }
        } else  {
          PetscInt diagID = i;
          rhoii_re = 0.0;
          rhoii_im = 0.0;
          if (ilow <= diagID && diagID < iupp) {
            PetscInt id_global_x =  diagID + mpirank_petsc*localsize_u;
            VecGetValues(state, 1, &id_global_x, &rhoii_re);
            id_global_x += localsize_u;
            VecGetValues(state, 1, &id_global_x, &rhoii_im);
          }
          rhoii = pow(rhoii_re, 2.0) + pow(rhoii_im, 2.0);
        }
        lambdai = fabs(i - purestateID);
        sum += lambdai * rhoii;
      }
      J_re = sum;
      MPI_Allreduce(&sum, &J_re, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
      break; // case J_MEASURE
  }

  // return
  *J_re_ptr = J_re;
  *J_im_ptr = J_im;
}


void OptimTarget::evalJ_diff(const Vec state, Vec statebar, const double J_re_bar, const double J_im_bar){
  // Do nothing if no target is specified.
  if (target_type == TargetType::NONE) {
    return;
  }

  double lambdai, val, rhoii_re, rhoii_im;

  switch (objective_type) {

    case ObjectiveType::JFROBENIUS:

      if (target_type == TargetType::STATE && purestateID >= 0) {
        // Derivative of J = 1/2||x||^2 is xbar += x * Jbar, where x = rho(t) - E_mm
        VecAXPY(statebar, J_re_bar, state);
        // now substract 1.0*Jbar from m-th diagonal element
        PetscInt diagID = purestateID;
        if (decoherence_type != DecoherenceType::NONE) diagID = getVecID(purestateID,purestateID,(PetscInt)sqrt(dim));
        if (ilow <= diagID && diagID < iupp) {
          PetscInt id_global_x = diagID + mpirank_petsc*localsize_u;
          VecSetValue(statebar, id_global_x, -1.0*J_re_bar, ADD_VALUES);
        }
      } else {
        FrobeniusDistance_diff(state, statebar, J_re_bar/ 2.0);
      }
      break; // case JFROBENIUS

    case ObjectiveType::JTRACE:
      HilbertSchmidtOverlap_diff(statebar, true, J_re_bar, J_im_bar);
    break;

    case ObjectiveType::JMEASURE:

      PetscInt dimsq = dim;   // Schroedinger solver: dim = N
      if (decoherence_type != DecoherenceType::NONE) dimsq = (PetscInt)sqrt(dim); // Lindblad solver: dim = N^2

      // iterate over diagonal elements 
      for (PetscInt i=0; i<dimsq; i++){
        lambdai = fabs(i - purestateID);
        if (decoherence_type != DecoherenceType::NONE) {
          PetscInt diagID = getVecID(i,i,dimsq);
          val = lambdai * J_re_bar;
          if (ilow <= diagID && diagID < iupp) {
            PetscInt id_global_x =  diagID + mpirank_petsc*localsize_u;
            VecSetValue(statebar, id_global_x, val, ADD_VALUES);
          }
        } else {
          PetscInt diagID = i;
          rhoii_re = 0.0;
          rhoii_im = 0.0;
          if (ilow <= diagID && diagID < iupp) {
            PetscInt id_global_x =  diagID + mpirank_petsc*localsize_u;
            VecGetValues(state, 1, &id_global_x, &rhoii_re);
            VecSetValue(statebar, id_global_x, 2.*J_re_bar*lambdai*rhoii_re, ADD_VALUES);
            id_global_x += localsize_u;
            VecGetValues(state, 1, &id_global_x, &rhoii_im);
            VecSetValue(statebar, id_global_x, 2.*J_re_bar*lambdai*rhoii_im, ADD_VALUES);
          }
        }
      }
    break;
  }
  VecAssemblyBegin(statebar); VecAssemblyEnd(statebar);
}

double OptimTarget::finalizeJ(const double obj_cost_re, const double obj_cost_im) {
  if (target_type == TargetType::NONE) {
    return 0.0;
  }

  double obj_cost = 0.0;
  if (objective_type == ObjectiveType::JTRACE) {
    if (decoherence_type == DecoherenceType::NONE) {
      obj_cost = 1.0 - (pow(obj_cost_re,2.0) + pow(obj_cost_im, 2.0));
    } else {
      obj_cost = 1.0 - obj_cost_re;
    }
  } else {
    obj_cost = obj_cost_re;
    assert(obj_cost_im <= 1e-14);
  }

  return obj_cost;
}


void OptimTarget::finalizeJ_diff(const double obj_cost_re, const double obj_cost_im, double* obj_cost_re_bar, double* obj_cost_im_bar){
  if (target_type == TargetType::NONE) {
    *obj_cost_re_bar = 0.0;
    *obj_cost_im_bar = 0.0;
    return;
  }

  if (objective_type == ObjectiveType::JTRACE) {
    if (decoherence_type == DecoherenceType::NONE) {
      // obj_cost = 1.0 - (pow(obj_cost_re,2.0) + pow(obj_cost_im, 2.0));
      *obj_cost_re_bar = -2.*obj_cost_re;
      *obj_cost_im_bar = -2.*obj_cost_im;
    } else {
      *obj_cost_re_bar = -1.0;
      *obj_cost_im_bar = 0.0;
    }
  } else {
    *obj_cost_re_bar = 1.0;
    *obj_cost_im_bar = 0.0;
  }
}
