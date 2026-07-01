#include "oscillator.hpp"
#include "config.hpp"
#include "defs.hpp"
#include "mpi_logger.hpp"
#include <stdexcept>

const double TOLERANCE = 1e-10; // Tolerance for avoiding numerical precision issues when comparing floating point numbers in evalControl.

Oscillator::Oscillator(){
  nlevels = 0;
  total_time = 0;
  ground_freq = 0.0;
  control_zero_boundary_condition = true;
  control_flux_zero_boundary_condition = true;
}

Oscillator::Oscillator(const Config& config, size_t id, std::mt19937 rand_engine, int param_offset, bool quietmode){

  myid = id;

  // Extract parameters from config
  const std::vector<size_t>& nlevels_all_ = config.getNLevels();
  nlevels = nlevels_all_[id];

  total_time = config.getTotalTime();

  const std::vector<double>& trans_freq = config.getTransitionFrequency();
  const std::vector<double>& rot_freq = config.getRotationFrequency();
  const std::vector<double>& selfkerr_config = config.getSelfKerr();

  ground_freq = trans_freq[id] * 2.0 * M_PI;
  selfkerr = selfkerr_config[id] * 2.0 * M_PI;
  detuning_freq = 2.0 * M_PI * (trans_freq[id] - rot_freq[id]);

  const std::vector<double>& carrier_freq_config = config.getCarrierFrequencies(id);
  carrier_freq = carrier_freq_config;
  for (size_t i=0; i<carrier_freq.size(); i++) {
    carrier_freq[i] *= 2.0*M_PI;
  }

  decoherence_type = config.getDecoherenceType();
  const std::vector<double>& decay_time_config = config.getDecayTime();
  const std::vector<double>& dephase_time_config = config.getDephaseTime();
  decay_time = decay_time_config[id];
  dephase_time = dephase_time_config[id];

  MPI_Comm_rank(PETSC_COMM_WORLD, &mpirank_petsc);
  MPI_Comm_size(PETSC_COMM_WORLD, &mpisize_petsc);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank_world);

  MPILogger logger(mpirank_world);

  // Get system dimension N (Schroedinger) or N^2 (Lindblad)
  PetscInt dim = 1;
  for (size_t ioscil = 0; ioscil < nlevels_all_.size(); ioscil++) {
    dim *= nlevels_all_[ioscil];
  }
  if (decoherence_type != DecoherenceType::NONE) dim *= dim; // Lindblad: N^2

  // Set local sizes of subvectors u,v in state x=[u,v]
  localsize_u = dim / mpisize_petsc; 
  ilow = mpirank_petsc * localsize_u;
  iupp = ilow + localsize_u;         

  /* Check if boundary conditions for controls should be enfored (default: yes). */
  control_zero_boundary_condition = config.getControlZeroBoundaryCondition();
  control_flux_zero_boundary_condition = config.getControlFluxZeroBoundaryCondition();

  // Initialize the control parameterization basis functions. Note: Currently only one control parameterization is supported. nsegments <= 1!
  int nparams_per_seg = 0;
  // for (auto drive_parameterization : config.getControlParameterizations(id)) {
  const auto& drive_parameterization = config.getControlParameterizations(id);
  auto tstart = drive_parameterization.tstart.value_or(0.0);
  auto tstop  = drive_parameterization.tstop.value_or(total_time);
 
  switch (drive_parameterization.type) {
    case ControlType::BSPLINE: {
      ControlBasis* mysplinebasis = new BSpline2nd(*drive_parameterization.nspline, tstart, tstop, control_zero_boundary_condition);
      mysplinebasis->setSkip(nparams_per_seg);
      nparams_per_seg += mysplinebasis->getNparams() * carrier_freq.size();
      drive_basisfunctions.push_back(mysplinebasis);
      break;
    }
    case ControlType::BSPLINE0: {
      ControlBasis* mysplinebasis = new BSpline0(*drive_parameterization.nspline, tstart, tstop, control_zero_boundary_condition);
      mysplinebasis->setSkip(nparams_per_seg);
      nparams_per_seg += mysplinebasis->getNparams() * carrier_freq.size();
      drive_basisfunctions.push_back(mysplinebasis);
      break;
    }
    case ControlType::BSPLINEAMP: {
      ControlBasis* mysplinebasis = new BSpline2ndAmplitude(*drive_parameterization.nspline, *drive_parameterization.scaling, tstart, tstop, control_zero_boundary_condition);
      mysplinebasis->setSkip(nparams_per_seg);
      nparams_per_seg += mysplinebasis->getNparams() * carrier_freq.size();
      drive_basisfunctions.push_back(mysplinebasis);
      break;
    }
    case ControlType::NONE: {
      // logger.exitWithError("Control type 'none' not supported.");
      // Do nothing.
    }
    // } 
  } // end switch

  // Initialize optional flux control parameterization. Independent from drive controls.
  int nparams_flux = 0;
  const auto& flux_parameterization = config.getControlFluxParameterizations(id);
  auto flux_tstart = flux_parameterization.tstart.value_or(0.0);
  auto flux_tstop  = flux_parameterization.tstop.value_or(total_time);
  switch (flux_parameterization.type) {
    case ControlType::BSPLINE: {
      ControlBasis* fluxbasis = new BSpline2nd(*flux_parameterization.nspline, flux_tstart, flux_tstop, control_flux_zero_boundary_condition);
      fluxbasis->setSkip(nparams_flux);
      nparams_flux += fluxbasis->getNparams(); // flux has no carrier expansion
      flux_basisfunctions.push_back(fluxbasis);
      break;
    }
    case ControlType::BSPLINE0: {
      ControlBasis* fluxbasis = new BSpline0(*flux_parameterization.nspline, flux_tstart, flux_tstop, control_flux_zero_boundary_condition);
      fluxbasis->setSkip(nparams_flux);
      nparams_flux += fluxbasis->getNparams(); // flux has no carrier expansion
      flux_basisfunctions.push_back(fluxbasis);
      break;
    }
    case ControlType::BSPLINEAMP: {
      logger.exitWithError("Flux control does not support 'spline_amplitude' parameterization.");
      break;
    }
    case ControlType::NONE: {
      break;
    }
  }

  /* Initialization of the control parameters.  */
  for (size_t iseg = 0; iseg < drive_basisfunctions.size(); iseg++) { // NOTE: Currently only one control parameterization supported: iseg = 0!

    // const auto& drive_initialization = config.getControlInitializations(id)[seg];
    const auto& drive_initialization = config.getControlInitializations(id);
    if (drive_initialization.type == ControlInitializationType::FILE) { // read from file 

      size_t nparams = drive_basisfunctions[iseg]->getNparams() * carrier_freq.size();
      drive_params.resize(nparams);
      if (mpirank_world == 0) {
        read_vector(drive_initialization.filename.value().c_str(), drive_params.data(), nparams, quietmode, param_offset);
      }
      MPI_Bcast(drive_params.data(), nparams, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    } else if (drive_initialization.type == ControlInitializationType::CONSTANT || drive_initialization.type == ControlInitializationType::RANDOM) {

      // Note, the config amplitude is multiplied by 2pi here!!
      double initval = drive_initialization.amplitude.value()*2.0*M_PI;
      
      for (size_t f = 0; f<carrier_freq.size(); f++) {
        for (int i=0; i<drive_basisfunctions[iseg]->getNparams(); i++){

          double val; 
          if (drive_initialization.type == ControlInitializationType::CONSTANT) {
            val = initval;
          } else if (drive_initialization.type == ControlInitializationType::RANDOM) {
            // Uniform distribution [-a,a)
            std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);
            double randval = uniform_dist(rand_engine);  // random in [0,1)
            // scale to chosen amplitude [-a,a]
            val = initval*randval;
            val = 2*val - initval;
          } else {
            logger.exitWithError("Unknown control initialization type.");
          }

          // Push the value to the parameter storage
          drive_params.push_back(val);
        }

        // if BSPLINEAMP: Two values can be provided: First one for the amplitude (set above), second one for the phase which otherwise is set to 0.0 (overwrite here)
        if (drive_basisfunctions[iseg]->getType() == ControlType::BSPLINEAMP) {
          drive_params[drive_params.size()-1] = drive_initialization.phase.value();
        }
      }
    }
  }

  // Flux parameter initialization (single channel, no carrier expansion)
  for (size_t iseg = 0; iseg < flux_basisfunctions.size(); iseg++) {
    const auto& flux_initialization = config.getControlFluxInitializations(id);
    if (flux_initialization.type == ControlInitializationType::FILE) {
      size_t nparams = flux_basisfunctions[iseg]->getNparams();
      flux_params.resize(nparams);
      if (mpirank_world == 0) {
        read_vector(flux_initialization.filename.value().c_str(), flux_params.data(), nparams, quietmode, 0);
      }
      MPI_Bcast(flux_params.data(), nparams, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    } else if (flux_initialization.type == ControlInitializationType::CONSTANT || flux_initialization.type == ControlInitializationType::RANDOM) {
      double initval = flux_initialization.amplitude.value() * 2.0 * M_PI;
      for (int i = 0; i < flux_basisfunctions[iseg]->getNparams(); i++) {
        double val;
        if (flux_initialization.type == ControlInitializationType::CONSTANT) {
          val = initval;
        } else {
          std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);
          double randval = uniform_dist(rand_engine);
          val = initval * randval;
          val = 2 * val - initval;
        }
        flux_params.push_back(val);
      }
    }
  }

  /* Make sure the initial guess satisfies the boundary conditions, if needed */
  if (drive_params.size() > 0 && control_zero_boundary_condition){
    for (size_t bs = 0; bs < drive_basisfunctions.size(); bs++){
      for (size_t f=0; f < carrier_freq.size(); f++) {
        drive_basisfunctions[bs]->enforceBoundary(drive_params.data(), f);
      }
    }
  }

  if (flux_params.size() > 0 && control_flux_zero_boundary_condition) {
    for (size_t bs = 0; bs < flux_basisfunctions.size(); bs++) {
      // Flux is a single channel with one logical carrier id (0)
      flux_basisfunctions[bs]->enforceBoundary(flux_params.data(), 0);
    }
  }

  /* Compute and store dimension of preceding and following oscillators */
  dim_preOsc = 1;
  dim_postOsc = 1;
  for (size_t j=0; j<nlevels_all_.size(); j++) {
    if (j < id) dim_preOsc  *= nlevels_all_[j];
    if (j > id) dim_postOsc *= nlevels_all_[j];
  }
}


Oscillator::~Oscillator(){
  if (drive_params.size() > 0) {
    for (size_t i=0; i<drive_basisfunctions.size(); i++) 
      delete drive_basisfunctions[i];
  }
  if (flux_params.size() > 0) {
    for (size_t i = 0; i < flux_basisfunctions.size(); i++)
      delete flux_basisfunctions[i];
  }
}

void Oscillator::setParams(const double* x){

  // copy x into the oscillators parameter storage
  for (size_t i=0; i<drive_params.size(); i++) {
    drive_params[i] = x[i]; 
  }
  for (size_t i = 0; i < flux_params.size(); i++) {
    flux_params[i] = x[drive_params.size() + i];
  }
}

void Oscillator::getParams(double* x){
  for (size_t i=0; i<drive_params.size(); i++) {
    x[i] = drive_params[i]; 
  }
  for (size_t i = 0; i < flux_params.size(); i++) {
    x[drive_params.size() + i] = flux_params[i];
  }
}

int Oscillator::getNSegParams(int parameterizationID){
  int n = 0;
  if (drive_params.size()>0) {
    assert(drive_basisfunctions.size() > static_cast<size_t>(parameterizationID));
    n = drive_basisfunctions[parameterizationID]->getNparams()*carrier_freq.size();
  }
  return n; 
}

int Oscillator::getNFluxSegParams(int parameterizationID) {
  int n = 0;
  if (flux_params.size() > 0) {
    assert(flux_basisfunctions.size() > static_cast<size_t>(parameterizationID));
    n = flux_basisfunctions[parameterizationID]->getNparams();
  }
  return n;
}

double Oscillator::evalControlVariation(){
  // NOTE: drive_params holds the relevant copy of the optimizers 'x' vector 
  double var_reg = 0.0;
  if (drive_params.size()>0) {
    // Iterate over control parameterizations. NOTE: Currently only one parameterization segment is supported. iseg = 0!
    for (size_t iseg= 0; iseg< drive_basisfunctions.size(); iseg++){
      /* Iterate over carrier frequencies */
      for (size_t f=0; f < carrier_freq.size(); f++) {
        var_reg += drive_basisfunctions[iseg]->computeVariation(drive_params, f);
      }
    }
  } 
  return var_reg;
}

void Oscillator::evalControlVariationDiff(Vec G, double var_reg_bar, int skip_to_oscillator){
  // NOTE: drive_params holds the relevant copy of the 'x' array

  if (drive_params.size()>0) {
    PetscScalar* grad; 
    VecGetArray(G, &grad);

    // Iterate over basis parameterizations??? 
    for (size_t iseg = 0; iseg< drive_basisfunctions.size(); iseg++){
      /* Iterate over carrier frequencies */
      for (size_t f=0; f < carrier_freq.size(); f++) {
        // pass the portion of the gradient that corresponds to this oscillator
        drive_basisfunctions[iseg]->computeVariation_diff(grad+skip_to_oscillator, drive_params, var_reg_bar, f);
     }
    }
    VecRestoreArray(G, &grad);
  } 
}

int Oscillator::evalDriveControl(const double t, double* p_ptr, double* q_ptr){

  // Default: Non controllable oscillator. Will typically be overwritten below. 
  *p_ptr = 0.0;
  *q_ptr = 0.0;

  /* Evaluate p(t) and q(t) using the parameters */
  if (drive_params.size()>0) {
    // Iterate over control parameterizations. Only one will be used, see the break-statement. 
    for (size_t bs = 0; bs < drive_basisfunctions.size(); bs++){ 
     if (drive_basisfunctions[bs]->getTstart() - TOLERANCE <= t && 
          drive_basisfunctions[bs]->getTstop() + TOLERANCE >= t ) {
        /* Iterate over carrier frequencies */
        double sum_p = 0.0;
        double sum_q = 0.0;
        for (size_t f=0; f < carrier_freq.size(); f++) {
          /* Evaluate the Bspline for this carrier wave */
          double Blt1 = 0.0; // Sums over alpha^1 * basisfunction(t) (real)
          double Blt2 = 0.0; // Sums over alpha^2 * basisfunction(t) (imag)
          drive_basisfunctions[bs]->evaluate(t, drive_params, f, &Blt1, &Blt2);
          if (drive_basisfunctions[bs]->getType() == ControlType::BSPLINEAMP) {
            double cos_omt = cos(carrier_freq[f]*t + Blt2);
            double sin_omt = sin(carrier_freq[f]*t + Blt2);
            sum_p += cos_omt * Blt1; 
            sum_q += sin_omt * Blt1;
          } else {
            double cos_omt = cos(carrier_freq[f]*t);
            double sin_omt = sin(carrier_freq[f]*t);
            sum_p += cos_omt * Blt1 - sin_omt * Blt2; 
            sum_q += sin_omt * Blt1 + cos_omt * Blt2;
          }
        }
        *p_ptr = sum_p;
        *q_ptr = sum_q;
        break;
      }
    }
  } 

  return 0;
}

int Oscillator::evalFluxControl(const double t, double* flux_ptr) {
  *flux_ptr = 0.0;

  if (flux_params.size() > 0) {
    for (size_t bs = 0; bs < flux_basisfunctions.size(); bs++) {
      if (flux_basisfunctions[bs]->getTstart() - TOLERANCE <= t &&
          flux_basisfunctions[bs]->getTstop() + TOLERANCE >= t) {
        double Blt1 = 0.0;
        double Blt2 = 0.0;
        // Flux uses a single scalar channel. We take the first basis channel.
        flux_basisfunctions[bs]->evaluate(t, flux_params, 0, &Blt1, &Blt2);
        *flux_ptr = Blt1;
        break;
      }
    }
  }

  return 0;
}

int Oscillator::evalControl(const double t, double* p_ptr, double* q_ptr, double* flux_ptr) {
  evalDriveControl(t, p_ptr, q_ptr);
  evalFluxControl(t, flux_ptr);
  return 0;
}

int Oscillator::evalDriveControl_diff(const double t, double* grad, const double pbar, const double qbar) {

  if (drive_params.size()>0) {

    // Iterate over control parameterizations. Only one is active, see break statement.
    for (size_t bs = 0; bs < drive_basisfunctions.size(); bs++){
      if (drive_basisfunctions[bs]->getTstart() - TOLERANCE <= t && 
          drive_basisfunctions[bs]->getTstop() + TOLERANCE >= t ) {
        /* Iterate over carrier frequencies */
        for (size_t f=0; f < carrier_freq.size(); f++) {

          if (drive_basisfunctions[bs]->getType() == ControlType::BSPLINEAMP) {
            // drive_basisfunctions[bs]->derivative(t, drive_params, dpdalpha, carrier_freq[f], 1.0, f);  // +/-1.0 is used as a flag inside Bsline2ndAmplitude->evaluate() to determine whether this is for p (1.0) or for q (-1.0)
            // drive_basisfunctions[bs]->derivative(t, drive_params, dqdalpha, carrier_freq[f], -1.0, f);
            // drive_basisfunctions[bs]->derivative(t, drive_params, dpdalpha, carrier_freq[f], 1.0, f);  // +/-1.0 is used as a flag inside Bsline2ndAmplitude->evaluate() to determine whether this is for p (1.0) or for q (-1.0)
            printf("Gradient for BsplineAmp parameterization is currently not implemented. Reach out to guenther5@llnl.gov if you need this. TODO.\n");
            exit(1);
          } else {

            double cos_omt = cos(carrier_freq[f]*t);
            double sin_omt = sin(carrier_freq[f]*t);
            double Blt1bar = sin_omt*qbar + cos_omt*pbar;
            double Blt2bar = cos_omt*qbar - sin_omt*pbar;

            /* Derivative wrt control alpha */
            drive_basisfunctions[bs]->derivative(t, drive_params, grad, Blt1bar, Blt2bar, f); // dp(t) / dalpha
          }
        }
        break;
      }
    }
  } 

  return 0;
}

int Oscillator::evalControl_diff(const double t, double* grad, const double pbar, const double qbar, const double fbar) {
  // First, accumulate drive-channel sensitivity
  evalDriveControl_diff(t, grad, pbar, qbar);

  // Then, accumulate flux-channel sensitivity in the flux parameter block
  if (flux_params.size() > 0) {
    double* flux_grad = grad + drive_params.size();
    for (size_t bs = 0; bs < flux_basisfunctions.size(); bs++) {
      if (flux_basisfunctions[bs]->getTstart() - TOLERANCE <= t &&
          flux_basisfunctions[bs]->getTstop() + TOLERANCE >= t) {
        flux_basisfunctions[bs]->derivative(t, flux_params, flux_grad, fbar, 0.0, 0);
        break;
      }
    }
  }

  return 0;
}

double Oscillator::expectedEnergy(const Vec x) {
 
  PetscInt dim;
  VecGetSize(x, &dim);
  PetscInt dimmat;
  if (decoherence_type != DecoherenceType::NONE)  dimmat = (PetscInt) sqrt(dim/2);
  else dimmat = (PetscInt) dim/2;

  /* Iterate over diagonal elements to add up expected energy level */
  double expected = 0.0;
  for (PetscInt i=0; i<dimmat; i++) {
    /* Get diagonal element in number operator */
    PetscInt num_diag = i % (nlevels*dim_postOsc);
    num_diag = num_diag / dim_postOsc;
    // Vectorize if Lindblad 
    PetscInt idx_diag = i;
    if (decoherence_type != DecoherenceType::NONE) idx_diag = getVecID(i,i,dimmat);
    
    double xdiag = 0.0;
    if (decoherence_type != DecoherenceType::NONE){ // Lindblad solver: += i * rho_ii
      if (ilow <= idx_diag && idx_diag < iupp) {
        PetscInt id_global_x = idx_diag + mpirank_petsc*localsize_u; 
        VecGetValues(x, 1, &id_global_x, &xdiag);
      }
      expected += num_diag * xdiag;
    }
    else { // Schoedinger solver: += i * | psi_i |^2
      if (ilow <= idx_diag && idx_diag < iupp) {
        PetscInt id_global_x = idx_diag + mpirank_petsc*localsize_u; 
        VecGetValues(x, 1, &id_global_x, &xdiag);
        expected += num_diag * xdiag * xdiag;
        id_global_x += localsize_u; 
        VecGetValues(x, 1, &id_global_x, &xdiag);
        expected += num_diag * xdiag * xdiag;
      }
    }
  }
  
  /* Sum up from all Petsc processors */
  double myexp = expected;
  MPI_Allreduce(&myexp, &expected, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

  return expected;
}


void Oscillator::expectedEnergy_diff(const Vec x, Vec x_bar, const double obj_bar) {
  PetscInt dim;
  VecGetSize(x, &dim);
  PetscInt dimmat;
  if (decoherence_type != DecoherenceType::NONE) dimmat = (PetscInt) sqrt(dim/2);
  else dimmat = dim/2;
  double xdiag, val;

  /* Derivative of projective measure */
  for (PetscInt i=0; i<dimmat; i++) {
    PetscInt num_diag = i % (nlevels*dim_postOsc);
    num_diag = num_diag / dim_postOsc;
    if (decoherence_type != DecoherenceType::NONE) { // Lindblas solver
      val = num_diag * obj_bar;
      PetscInt idx_diag = getVecID(i, i, dimmat);
      if (ilow <= idx_diag && idx_diag < iupp) {
        PetscInt id_global_x = idx_diag + mpirank_petsc*localsize_u; 
        VecSetValues(x_bar, 1, &id_global_x, &val, ADD_VALUES);
      }
    }
    else {
      // Real part
      PetscInt idx_diag = i;
      xdiag = 0.0;
      if (ilow <= idx_diag && idx_diag < iupp) {
        PetscInt id_global_x = idx_diag + mpirank_petsc*localsize_u; 
        VecGetValues(x, 1, &id_global_x, &xdiag);
        val = num_diag * xdiag * obj_bar;
        VecSetValues(x_bar, 1, &id_global_x, &val, ADD_VALUES);
        // Imaginary part
        id_global_x += localsize_u; 
        VecGetValues(x, 1, &id_global_x, &xdiag);
        val = - num_diag * xdiag * obj_bar; // TODO: Is this a minus or a plus?? 
        VecSetValues(x_bar, 1, &id_global_x, &val, ADD_VALUES);
      }
    }
  }
  VecAssemblyBegin(x_bar); VecAssemblyEnd(x_bar);

}


void Oscillator::population(const Vec x, std::vector<double> &pop) {

  PetscInt dimN = dim_preOsc * nlevels * dim_postOsc;
  double val;

  assert (pop.size() == nlevels);

  std::vector<double> mypop(nlevels, 0.0);

  /* Iterate over diagonal elements of the reduced density matrix for this oscillator */
  for (size_t i=0; i < nlevels; i++) {
    PetscInt identitystartID = i * dim_postOsc;
    /* Sum up elements from all dim_preOsc blocks of size (n_k * dim_postOsc) */
    double sum = 0.0;
    for (PetscInt j=0; j < dim_preOsc; j++) {
      PetscInt blockstartID = j * nlevels * dim_postOsc; // Go to the block
      /* Iterate over identity */
      for (PetscInt l=0; l < dim_postOsc; l++) {
        /* Get diagonal element */
        PetscInt rhoID = blockstartID + identitystartID + l; // Diagonal element of rho
        if (decoherence_type != DecoherenceType::NONE) { // Lindblad solver
          PetscInt diagID = getVecID(rhoID, rhoID, dimN);  // Position in vectorized rho
          double val = 0.0;
          if (ilow <= diagID && diagID < iupp)  {
            PetscInt id_global_x = diagID+ mpirank_petsc*localsize_u; 
            VecGetValues(x, 1, &id_global_x, &val);
          }
          sum += val;
        } else {
          PetscInt diagID = rhoID;
          val = 0.0;
          if (ilow <= diagID && diagID < iupp) {
            PetscInt id_global_x = diagID+ mpirank_petsc*localsize_u; 
            VecGetValues(x, 1, &id_global_x, &val);
            sum += val * val;
            id_global_x += localsize_u; 
            VecGetValues(x, 1, &id_global_x, &val);
            sum += val * val;
          }
        }
      }
    }
    mypop[i] = sum;
  } 

  /* Gather poppulation from all Petsc processors */
  for (size_t i=0; i<mypop.size(); i++) {pop[i] = mypop[i];}
  MPI_Allreduce(mypop.data(), pop.data(), static_cast<int>(nlevels), MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
}
