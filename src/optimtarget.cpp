#include "optimtarget.hpp"

OptimTarget::OptimTarget(){
  dim = 0;
  dim_rho = 0;
  dim_ess = 0;
  noscillators = 0;
  target_type = TargetType::GATE;
  objective_type = ObjectiveType::JTRACE;
  initcond_type = InitialConditionType::BASIS;
  lindbladtype = LindbladType::NONE;
  targetgate = NULL;
  purity_rho0 = 1.0;
  purestateID = -1;
  target_filename = "";
  targetstate = NULL;
  mpisize_petsc=0;
  mpirank_petsc=0;
}


OptimTarget::OptimTarget(std::vector<std::string> target_str, const std::string& objective_str, std::vector<std::string> initcond_str, MasterEq* mastereq, double total_time, std::vector<double> read_gate_rot, Vec rho_t0, bool quietmode_) : OptimTarget() {

  // initialize
  dim = mastereq->getDim();
  dim_rho = mastereq->getDimRho();
  dim_ess = mastereq->getDimEss();
  quietmode = quietmode_;
  lindbladtype = mastereq->lindbladtype;
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

  /* Get initial condition type */
  if (initcond_str[0].compare("file") == 0)              initcond_type = InitialConditionType::FROMFILE;
  else if  (initcond_str[0].compare("pure") == 0)        initcond_type = InitialConditionType::PURE;
  else if  (initcond_str[0].compare("ensemble") == 0)    initcond_type = InitialConditionType::ENSEMBLE;
  else if  (initcond_str[0].compare("performance") == 0) initcond_type = InitialConditionType::PERFORMANCE;
  else if  (initcond_str[0].compare("3states") == 0)     initcond_type = InitialConditionType::THREESTATES;
  else if  (initcond_str[0].compare("Nplus1") == 0)      initcond_type = InitialConditionType::NPLUSONE;
  else if  (initcond_str[0].compare("diagonal") == 0)    initcond_type = InitialConditionType::DIAGONAL;
  else if  (initcond_str[0].compare("basis") == 0)       initcond_type = InitialConditionType::BASIS;
  else {
    printf("ERROR: Unknown initial condition type %s!\n", initcond_str[0].c_str());
    exit(1);
  }
  /* Sanity check for Schrodinger solver initial conditions */
  if (lindbladtype == LindbladType::NONE){
    if (initcond_type == InitialConditionType::ENSEMBLE ||
        initcond_type == InitialConditionType::THREESTATES ||
        initcond_type == InitialConditionType::NPLUSONE ){
          printf("\n\n ERROR for initial condition setting: \n When running Schroedingers solver (collapse_type == NONE), the initial condition needs to be either 'pure' or 'from file' or 'diagonal' or 'basis'. Note that 'diagonal' and 'basis' in the Schroedinger case are the same (all unit vectors).\n\n");
          exit(1);
    } else if (initcond_type == InitialConditionType::BASIS) {
      // DIAGONAL and BASIS initial conditions in the Schroedinger case are the same. Overwrite it to DIAGONAL
      initcond_type = InitialConditionType::DIAGONAL;  
    }
  }
  /* Get list of involved oscillators */
  if (initcond_str.size() < 2) 
    for (size_t j=0; j<mastereq->getNOscillators(); j++)   
      initcond_str.push_back(std::to_string(j)); // Default: all oscillators
  for (size_t i=1; i<initcond_str.size(); i++) 
    initcond_IDs.push_back(atoi(initcond_str[i].c_str())); // Overwrite with config option, if given.

  /* Prepare initial state rho_t0 if PURE or FROMFILE or ENSEMBLE initialization. Otherwise they are set within prepareInitialState during evalF. */
  if (initcond_type == InitialConditionType::PURE) { 
    /* Initialize with tensor product of unit vectors. */
    if (initcond_IDs.size() != mastereq->getNOscillators()) {
      printf("ERROR during pure-state initialization: List of IDs must contain %zu elements!\n", mastereq->getNOscillators());
      exit(1);
    }
    // Find the id within the global composite system 
    PetscInt diag_id = 0;
    for (size_t k=0; k < initcond_IDs.size(); k++) {
      if (initcond_IDs[k] > mastereq->getOscillator(k)->getNLevels()-1){
        printf("ERROR in config setting. The requested pure state initialization |%zu> exceeds the number of allowed levels for that oscillator (%zu).\n", initcond_IDs[k], mastereq->getOscillator(k)->getNLevels());
        exit(1);
      }
      assert (initcond_IDs[k] < mastereq->getOscillator(k)->getNLevels());
      PetscInt dim_postkron = 1;
      for (size_t m=k+1; m < initcond_IDs.size(); m++) {
        dim_postkron *= mastereq->getOscillator(m)->getNLevels();
      }
      diag_id += initcond_IDs[k] * dim_postkron;
    }
    // Vectorize if lindblad solver
    PetscInt vec_id = diag_id;
    if (lindbladtype != LindbladType::NONE) vec_id = getVecID( diag_id, diag_id, dim_rho); 
    // Set 1.0 on the processor who owns this index
    if (ilow <= vec_id && vec_id < iupp) {
      PetscInt id_global_x =  vec_id + mpirank_petsc*localsize_u; // Global index of u_i in x=[u,v]
      VecSetValue(rho_t0, id_global_x, 1.0, INSERT_VALUES);
    }
  }
  else if (initcond_type == InitialConditionType::FROMFILE) { 
    /* Read initial condition from file */
    int nelems = 0;
    if (mastereq->lindbladtype != LindbladType::NONE) nelems = 2*dim_ess*dim_ess;
    else nelems = 2 * dim_ess;
    double * vec = new double[nelems];
    if (mpirank_world == 0) {
      assert (initcond_str.size()==2);
      std::string filename = initcond_str[1];
      read_vector(filename.c_str(), vec, nelems, quietmode);
    }
    MPI_Bcast(vec, nelems, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (lindbladtype != LindbladType::NONE) { // Lindblad solver, fill density matrix
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
          VecSetValue(rho_t0, id_global_x, vec[i], INSERT_VALUES);  // RealPart
          VecSetValue(rho_t0, id_global_x + localsize_u, vec[i + dim_ess*dim_ess], INSERT_VALUES); // Imaginary Part
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
  } else if (initcond_type == InitialConditionType::ENSEMBLE) {
    // Sanity check for the list in initcond_IDs!
    assert(initcond_IDs.size() >= 1); // at least one element 
    assert(initcond_IDs[initcond_IDs.size()-1] < mastereq->getNOscillators()); // last element can't exceed total number of oscillators
    for (size_t i=0; i < initcond_IDs.size()-1; i++){ // list should be consecutive!
      if (initcond_IDs[i]+1 != initcond_IDs[i+1]) {
        printf("ERROR: List of oscillators for ensemble initialization should be consecutive!\n");
        exit(1);
      }
    }
    // get dimension of subsystems defined by initcond_IDs, as well as the one before and after. Span in essential levels only.
    PetscInt dimpost = 1;
    PetscInt dimsub = 1;
    for (size_t i=0; i<mastereq->getNOscillators(); i++){
      if (initcond_IDs[0] <= i && i <= initcond_IDs[initcond_IDs.size()-1]) dimsub *= mastereq->nessential[i];
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

  /* Get target type */  
  purestateID = -1;
  target_filename = "";
  if ( target_str[0].compare("gate") ==0 ) {
    target_type = TargetType::GATE;

    /* Get gate rotation frequencies. */
    int noscillators = mastereq->nlevels.size();
    copyLast(read_gate_rot, noscillators);
    std::vector<double> gate_rot_freq(noscillators); 
    for (int iosc=0; iosc<noscillators; iosc++) {
      if (read_gate_rot[0] < 1e20) // the config option exists, use it, else use default 0.0 
        gate_rot_freq[iosc] = read_gate_rot[iosc];
      else
        gate_rot_freq[iosc] = 0.0;
    }
    /* Initialize the targetgate */
    targetgate = initTargetGate(target_str, mastereq->nlevels, mastereq->nessential, total_time, lindbladtype, gate_rot_freq, quietmode);
  }  
  else if (target_str[0].compare("pure")==0) {
    target_type = TargetType::PURE;
    purestateID = 0;
    if (target_str.size() < 2) {
      printf("# Warning: You want to prepare a pure state, but didn't specify which one. Taking default: ground-state |0...0> \n");
    } else {
      /* Compute the index m for preparing e_m e_m^\dagger. Note that the input is given for pure states PER OSCILLATOR such as |m_1 m_2 ... m_Q> and hence m = m_1 * dimPost(oscil 1) + m_2 * dimPost(oscil 2) + ... + m_Q */
      if (target_str.size() - 1 < mastereq->getNOscillators()) {
        copyLast(target_str, mastereq->getNOscillators()+1);
      }
      for (size_t i=0; i < mastereq->getNOscillators(); i++) {
        size_t Qi_state = atoi(target_str[i+1].c_str());
        if (Qi_state >= mastereq->getOscillator(i)->getNLevels()) {
          printf("ERROR in config setting. The requested pure state target |%zu> exceeds the number of modeled levels for that oscillator (%zu).\n", Qi_state, mastereq->getOscillator(i)->getNLevels());
          exit(1);
        }
        purestateID += Qi_state * mastereq->getOscillator(i)->dim_postOsc;
      }
    }
  } 
  else if (target_str[0].compare("file")==0) { 
    // Get the name of the file and pass it to the OptimTarget type later.
    target_type = TargetType::FROMFILE;
    assert(target_str.size() >= 2);
    target_filename = target_str[1];
  }
  else {
      printf("\n\n ERROR: Unknown optimization target: %s\n", target_str[0].c_str());
      exit(1);
  }

  /* Get the objective function */
  if (objective_str.compare("Jfrobenius")==0)     objective_type = ObjectiveType::JFROBENIUS;
  else if (objective_str.compare("Jtrace")==0)    objective_type = ObjectiveType::JTRACE;
  else if (objective_str.compare("Jmeasure")==0)  objective_type = ObjectiveType::JMEASURE;
  else  {
    printf("\n\n ERROR: Unknown objective function: %s\n", objective_str.c_str());
    exit(1);
  }

  /* Allocate target state, if it is read from file, or if target is a gate transformation VrhoV. If pure target, only store the ID. */
  if (target_type == TargetType::GATE || target_type == TargetType::FROMFILE) {
    VecCreate(PETSC_COMM_WORLD, &targetstate); 
    PetscInt globalsize = 2 * mastereq->getDim();  // Global state vector: 2 for real and imaginary part
    PetscInt localsize = globalsize / mpisize_petsc;  // Local vector per processor
    VecSetSizes(targetstate,localsize,globalsize);
    VecSetFromOptions(targetstate);
  }

  /* Read the target state from file into vec */
  if (target_type == TargetType::FROMFILE) {
    PetscInt nelems = 0;
    if (mastereq->lindbladtype != LindbladType::NONE) nelems = 2*dim_ess*dim_ess;
    else nelems = 2 * dim_ess;
    double* vec = new double[nelems];
    if (mpirank_world == 0) 
      read_vector(target_filename.c_str(), vec, nelems, quietmode);
    MPI_Bcast(vec, nelems, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (lindbladtype != LindbladType::NONE) { // Lindblad solver, fill density matrix
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
  if (target_type == TargetType::GATE || target_type == TargetType::FROMFILE)  VecDestroy(&targetstate);

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

  /* Simplify computation if the target is PURE, i.e. target = e_m or e_m * e_m^\dag */
  /* Tr(...) = phi_m if Schroedinger, or \rho_mm if Lindblad */
  if (target_type == TargetType::PURE){

    // Vectorize pure state ID if Lindblad
    PetscInt idm = purestateID;
    if (lindbladtype != LindbladType::NONE) idm = getVecID(purestateID, purestateID, (PetscInt)sqrt(dim));

    // Get real and imag values from the processor who owns the subvector index.
    if (ilow <= idm && idm < iupp) {
      PetscInt id_global_x = idm + mpirank_petsc*localsize_u; // Global index of u_i in x=[u,v]
      VecGetValues(state, 1, &id_global_x, &HS_re);
      id_global_x += localsize_u; // Imaginary part
      VecGetValues(state, 1, &id_global_x, &HS_im); // Should be 0.0 if Lindblad!
    }
    if (lindbladtype != LindbladType::NONE) assert(fabs(HS_im) <= 1e-14);

    // Communicate over all petsc processors.
    double myre = HS_re;
    double myim = HS_im;
    MPI_Allreduce(&myre, &HS_re, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
    MPI_Allreduce(&myim, &HS_im, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

  } else { // Target is not of the form e_m (schroedinger) or e_m e_m^\dagger (lindblad).

    if (lindbladtype != LindbladType::NONE) // Lindblad solver. HS overlap is real!
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

  // Simplified computation if target is pure 
  if (target_type == TargetType::PURE){
    PetscInt idm = purestateID;
    if (lindbladtype != LindbladType::NONE) idm = getVecID(purestateID, purestateID, (PetscInt)sqrt(dim));

    if (ilow <= idm && idm < iupp) {
      PetscInt id_global_x = idm + mpirank_petsc*localsize_u; // Global index of u_i in x=[u,v]
      VecSetValue(statebar, id_global_x, HS_re_bar*scale, ADD_VALUES);
      VecSetValue(statebar, id_global_x + localsize_u, HS_im_bar, ADD_VALUES);
    }

  } else { // Target is not of the form e_m or e_m*e_m^\dagger 

    if (lindbladtype != LindbladType::NONE)
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

  /* Switch over type of initial condition */
  switch (initcond_type) {

    case InitialConditionType::PERFORMANCE:
      /* Set up Input state psi = 1/sqrt(2N)*(Ones(N) + im*Ones(N)) or rho = psi*psi^\dag */
      VecZeroEntries(rho0);

      for (PetscInt i=0; i<dim_rho; i++){
        if (lindbladtype == LindbladType::NONE) {
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
      break;

    case InitialConditionType::FROMFILE:
      /* Do nothing. Init cond is already stored */
      break;

    case InitialConditionType::PURE:
      /* Do nothing. Init cond is already stored */
      break;

    case InitialConditionType::ENSEMBLE:
      /* Do nothing. Init cond is already stored */
      break;

    case InitialConditionType::THREESTATES:
      assert(lindbladtype != LindbladType::NONE);
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
      break;

    case InitialConditionType::NPLUSONE:
      assert(lindbladtype != LindbladType::NONE);

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
      break;

    case InitialConditionType::DIAGONAL:
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
      if (lindbladtype != LindbladType::NONE) elemID = getVecID(diagelem, diagelem, dim_rho); 
      val = 1.0;
      if (ilow <= elemID && elemID < iupp) {
        PetscInt id_global_x =  elemID + mpirank_petsc*localsize_u; 
        VecSetValues(rho0, 1, &id_global_x, &val, INSERT_VALUES);
      }
      VecAssemblyBegin(rho0); VecAssemblyEnd(rho0);

      /* Set initial conditon ID */
      if (lindbladtype != LindbladType::NONE) initID = iinit * ninit + iinit;
      else initID = iinit;

      break;

    case InitialConditionType::BASIS:
      assert(lindbladtype != LindbladType::NONE); // should never happen. For Schroedinger: BASIS equals DIAGONAL, and should go into the above switch case. 

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

      break;

    default:
      printf("ERROR! Wrong initial condition type: %d\n This should never happen!\n", static_cast<int>(initcond_type));
      exit(1);
  }

  return initID;
}


void OptimTarget::prepareTargetState(const Vec rho_t0){
  // If gate optimization, apply the gate and store targetstate for later use. Else, do nothing.
  if (target_type == TargetType::GATE) targetgate->applyGate(rho_t0, targetstate);

  /* Compute and store the purity of rho(0), Tr(rho(0)^2), so that it can be used by JTrace (HS overlap) */
  VecNorm(rho_t0, NORM_2, &purity_rho0);
  purity_rho0 = purity_rho0 * purity_rho0;
}



void OptimTarget::evalJ(const Vec state, double* J_re_ptr, double* J_im_ptr){
  double J_re = 0.0;
  double J_im = 0.0;
  double sum, rhoii, rhoii_re, rhoii_im, lambdai, norm;
  PetscInt dimsq;

  switch(objective_type) {

    /* J_Frob = 1/2 * || rho_target - rho(T)||^2_F  */
    case ObjectiveType::JFROBENIUS:

      if (target_type == TargetType::GATE || target_type == TargetType::FROMFILE ) {
        // target state is already set. Either \rho_target = Vrho(0)V^\dagger or read from file. Just eval norm.
        J_re = FrobeniusDistance(state) / 2.0;
      } 
      else {  // target = e_me_m^\dagger ( or target = e_m for Schroedinger)
        assert(target_type == TargetType::PURE);

        // substract 1.0 from m-th diagonal element then take the vector norm 
        PetscInt diagID = purestateID;
        if (lindbladtype != LindbladType::NONE) diagID = getVecID(purestateID,purestateID,(PetscInt)sqrt(dim));
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
      }
      break;  // case Frobenius

    /* J_Trace:  1 / purity * Tr(state * target^\dagger)  =  HilbertSchmidtOverlap(target, state) is real if Lindblad, and complex if Schroedinger! */
    case ObjectiveType::JTRACE:

      HilbertSchmidtOverlap(state, true, &J_re, &J_im); // is real if Lindblad solver. 
      break; // case J_Trace

    /* J_Measure = Tr(O_m rho(T)) = \sum_i |i-m| rho_ii(T) if Lindblad and \sum_i |i-m| |phi_i(T)|^2  if Schroedinger */
    case ObjectiveType::JMEASURE:
      // Sanity check
      if (target_type != TargetType::PURE) {
        printf("Error: Wrong setting for objective function. Jmeasure can only be used for 'pure' targets.\n");
        exit(1);
      }

      dimsq = dim;   // Schroedinger solver: dim = N
      if (lindbladtype != LindbladType::NONE) dimsq = (PetscInt)sqrt(dim); // Lindblad solver: dim = N^2

      // iterate over diagonal elements 
      sum = 0.0;
      for (PetscInt i=0; i<dimsq; i++){
        if (lindbladtype != LindbladType::NONE) {
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
  double lambdai, val, rhoii_re, rhoii_im;

  switch (objective_type) {

    case ObjectiveType::JFROBENIUS:

      if (target_type == TargetType::GATE || target_type == TargetType::FROMFILE ) {
        FrobeniusDistance_diff(state, statebar, J_re_bar/ 2.0);
      } else {
        assert(target_type == TargetType::PURE);         
        // Derivative of J = 1/2||x||^2 is xbar += x * Jbar, where x = rho(t) - E_mm
        VecAXPY(statebar, J_re_bar, state);
        // now substract 1.0*Jbar from m-th diagonal element
        PetscInt diagID = purestateID;
        if (lindbladtype != LindbladType::NONE) diagID = getVecID(purestateID,purestateID,(PetscInt)sqrt(dim));
        if (ilow <= diagID && diagID < iupp) {
          PetscInt id_global_x = diagID + mpirank_petsc*localsize_u;
          VecSetValue(statebar, id_global_x, -1.0*J_re_bar, ADD_VALUES);
        }
      }
      break; // case JFROBENIUS

    case ObjectiveType::JTRACE:
      HilbertSchmidtOverlap_diff(statebar, true, J_re_bar, J_im_bar);
    break;

    case ObjectiveType::JMEASURE:
      assert(target_type == TargetType::PURE);         

      PetscInt dimsq = dim;   // Schroedinger solver: dim = N
      if (lindbladtype != LindbladType::NONE) dimsq = (PetscInt)sqrt(dim); // Lindblad solver: dim = N^2

      // iterate over diagonal elements 
      for (PetscInt i=0; i<dimsq; i++){
        lambdai = fabs(i - purestateID);
        if (lindbladtype != LindbladType::NONE) {
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
  double obj_cost = 0.0;

  if (objective_type == ObjectiveType::JTRACE) {
    if (lindbladtype == LindbladType::NONE) {
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

  if (objective_type == ObjectiveType::JTRACE) {
    if (lindbladtype == LindbladType::NONE) {
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


double OptimTarget::RiemannianDistance(const Mat U_final_re, const Mat U_final_im){

  /* First, allreduce the U_final matrix */
  if (lindbladtype == LindbladType::NONE) {
    PetscScalar *data;
    MatDenseGetArray(U_final_re, &data);
    MPI_Allreduce(MPI_IN_PLACE, data, 2*dim * 2*dim, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD);
    MatDenseRestoreArray(U_final_re, &data);

    MatDenseGetArray(U_final_im, &data);
    MPI_Allreduce(MPI_IN_PLACE, data, 2*dim * 2*dim, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD);
    MatDenseRestoreArray(U_final_im, &data);
  }

  /* Compute the Riemannian distance */
  // Use the util getEigenvalues function. 

  // Here, for testing, reimplement the final fidelity using U and VxV from the targetgate, just to make sure that U_final is correct. 

  // Compute fidelity = 1/N | Tr(U_final^\dagger V) |^2, where U_final = U_final_re + i U_final_im and V = VxV_re + i VxV_im
  double trace_re = 0.0;
  double trace_im = 0.0;
  for (PetscInt i=0; i<dim; i++){
    for (PetscInt j=0; j<dim; j++){
      double U_re_ji, U_im_ji, V_re_ij, V_im_ij;
      MatGetValue(U_final_re, j, i, &U_re_ji);
      MatGetValue(U_final_im, j, i, &U_im_ji);
      MatGetValue(targetgate->VxV_re, i, j, &V_re_ij);
      MatGetValue(targetgate->VxV_im, i, j, &V_im_ij);
      // Conjugate U_final
      trace_re += U_re_ji * V_re_ij + U_im_ji * V_im_ij;
      trace_im += -U_im_ji * V_re_ij + U_re_ji * V_im_ij;
    } 
  }
  double obj = (trace_re * trace_re + trace_im * trace_im) / (dim * dim);

  return obj; 
}