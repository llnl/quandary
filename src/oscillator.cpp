#include "oscillator.hpp"
#include "config.hpp"
#include "defs.hpp"
#include "mpi_logger.hpp"
#include <stdexcept>


Oscillator::Oscillator(){
  myid = 0;
  nlevels = 0;
  total_time = 0;
  ground_freq = 0.0;
  selfkerr = 0.0;
  detuning_freq = 0.0;
  decay_time = 0.0;
  dephase_time = 0.0;
  drive_basisfunctions_re = nullptr;
  drive_basisfunctions_im = nullptr;
  flux_basisfunctions = nullptr;
}

Oscillator::Oscillator(const Config& config, size_t id, std::mt19937& rand_engine, int param_offset, bool quietmode){
  
  MPI_Comm_rank(PETSC_COMM_WORLD, &mpirank_petsc);
  MPI_Comm_size(PETSC_COMM_WORLD, &mpisize_petsc);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank_world);
  MPILogger logger(mpirank_world);

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

  decoherence_type = config.getDecoherenceType();
  const std::vector<double>& decay_time_config = config.getDecayTime();
  const std::vector<double>& dephase_time_config = config.getDephaseTime();
  decay_time = decay_time_config[id];
  dephase_time = dephase_time_config[id];

  carrier_freq = config.getCarrierFrequencies(id);
  for (size_t i=0; i<carrier_freq.size(); i++) {
    carrier_freq[i] *= 2.0*M_PI;
  }

  // Get system dimension N (Schroedinger) or N^2 (Lindblad)
  PetscInt dim = 1;
  for (size_t ioscil = 0; ioscil < nlevels_all_.size(); ioscil++) {
    dim *= nlevels_all_[ioscil];
  }
  if (decoherence_type != DecoherenceType::NONE) dim *= dim; // Lindblad: N^2

  // Set local sizes of subvectors u,v in real-valued state x=[u,v]
  localsize_u = dim / mpisize_petsc; 
  ilow = mpirank_petsc * localsize_u;
  iupp = ilow + localsize_u;         

  /* Compute and store dimension of preceding and following oscillators */
  dim_preOsc = 1;
  dim_postOsc = 1;
  for (size_t j=0; j<nlevels_all_.size(); j++) {
    if (j < id) dim_preOsc  *= nlevels_all_[j];
    if (j > id) dim_postOsc *= nlevels_all_[j];
  }

  // Create p/q drive control parameterizations with carrier waves.
  const auto& pq_drive_settings = config.getControlParameterizations(id);
  auto nspline = pq_drive_settings.nspline.value_or(0);
  auto tstart = pq_drive_settings.tstart.value_or(0.0);
  auto tstop  = pq_drive_settings.tstop.value_or(total_time);
  const bool drive_zero_bc = config.getControlZeroBoundaryCondition();
  switch (pq_drive_settings.type) {
    case ControlType::BSPLINE: {
      drive_basisfunctions_re = new BSpline2nd(static_cast<int>(nspline),static_cast<int>(carrier_freq.size()),tstart, tstop, drive_zero_bc);
      drive_basisfunctions_im = new BSpline2nd(static_cast<int>(nspline), static_cast<int>(carrier_freq.size()), tstart, tstop, drive_zero_bc);
      break;
    }
    case ControlType::BSPLINE0: {
      drive_basisfunctions_re = new BSpline0(static_cast<int>(nspline), static_cast<int>(carrier_freq.size()), tstart, tstop, drive_zero_bc);
      drive_basisfunctions_im = new BSpline0(static_cast<int>(nspline), static_cast<int>(carrier_freq.size()), tstart, tstop, drive_zero_bc);
      break;
    }
    case ControlType::NONE: {
      // Keep drive basis function pointers null.
      break;
    }
  } 

  // Initialize parameters of the p/q drive 
  const auto& drive_init_settings = config.getControlInitializations(id);
  std::vector<double> params_drive(getNDriveParams(), 0.0);
  if (drive_init_settings.type == ControlInitializationType::FILE) {
    if (mpirank_world == 0) {
      read_vector(drive_init_settings.filename.value().c_str(), params_drive.data(), params_drive.size(), quietmode, param_offset);
    }
    MPI_Bcast(params_drive.data(), getNDriveParams(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  } else if (drive_init_settings.type == ControlInitializationType::CONSTANT) {
    const double initval = drive_init_settings.amplitude.value_or(0.0) * 2.0 * M_PI;
    std::fill(params_drive.begin(), params_drive.end(), initval);
  } else if (drive_init_settings.type == ControlInitializationType::RANDOM) {
    const double initval = drive_init_settings.amplitude.value_or(0.0) * 2.0 * M_PI;
    std::uniform_real_distribution<double> uniform_dist(-initval, initval);
    for (double& param : params_drive) {
      param = uniform_dist(rand_engine);
    }
  } else {
    logger.exitWithError("Unknown control initialization type for drive controls.");
  }
  setControlParams(params_drive.data());
  // Make sure initial parameters satisfy boundary conditions.
  if (drive_basisfunctions_re) drive_basisfunctions_re->enforceBoundary();
  if (drive_basisfunctions_im) drive_basisfunctions_im->enforceBoundary();

  // Create flux control parameterization
  const auto& flux_settings= config.getControlFluxParameterizations(id);
  const bool flux_zero_bc = config.getControlFluxZeroBoundaryCondition();
  auto flux_tstart = flux_settings.tstart.value_or(0.0);
  auto flux_tstop  = flux_settings.tstop.value_or(total_time);
  auto flux_nspline = flux_settings.nspline.value_or(0);
  switch (flux_settings.type) {
    case ControlType::BSPLINE: {
      flux_basisfunctions = new BSpline2nd(static_cast<int>(flux_nspline), 1, flux_tstart, flux_tstop, flux_zero_bc);
      break;
    }
    case ControlType::BSPLINE0: {
      flux_basisfunctions = new BSpline0(static_cast<int>(flux_nspline), 1, flux_tstart, flux_tstop, flux_zero_bc);
      break;
    }
    case ControlType::NONE: {
      break;
    }
  }

  // Initialize flux control parameter
  const auto& flux_init_settings = config.getControlFluxInitializations(id);
  const int param_offset_flux = param_offset + getNDriveParams();
  std::vector<double> params_flux(getNFluxParams(), 0.0);
  if (flux_init_settings.type == ControlInitializationType::FILE) {
    if (mpirank_world == 0) {
      read_vector(flux_init_settings.filename.value().c_str(), params_flux.data(), params_flux.size(), quietmode, param_offset_flux);
    }
    MPI_Bcast(params_flux.data(), getNFluxParams(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  } else if (flux_init_settings.type == ControlInitializationType::CONSTANT) {
    const double initval = flux_init_settings.amplitude.value_or(0.0) * 2.0 * M_PI;
    std::fill(params_flux.begin(), params_flux.end(), initval);
  } else if (flux_init_settings.type == ControlInitializationType::RANDOM) {
    const double initval = flux_init_settings.amplitude.value_or(0.0) * 2.0 * M_PI;
    std::uniform_real_distribution<double> uniform_dist(-initval, initval);
    for (double& param : params_flux) {
      param = uniform_dist(rand_engine);
    }
  } else {
    logger.exitWithError("Unknown control initialization type for flux controls.");
  }
  if (flux_basisfunctions) {
    flux_basisfunctions->setParams(params_flux.data(), 0);
    flux_basisfunctions->enforceBoundary();
  }
}


Oscillator::~Oscillator(){
  if (drive_basisfunctions_re) delete drive_basisfunctions_re;
  if (drive_basisfunctions_im) delete drive_basisfunctions_im;
  if (flux_basisfunctions) delete flux_basisfunctions;
}

void Oscillator::setControlParams(const double* x) {

  // Copy global parameter block x into p/q and flux parameterizations.
  int skip = 0;
  for (size_t f=0; f<carrier_freq.size(); f++) {
    if (drive_basisfunctions_re) drive_basisfunctions_re->setParams(x + skip, f);
    skip += drive_basisfunctions_re ? drive_basisfunctions_re->getNparams(f) : 0;
    if (drive_basisfunctions_im) drive_basisfunctions_im->setParams(x + skip, f);
    skip += drive_basisfunctions_im ? drive_basisfunctions_im->getNparams(f) : 0;
  }
  if (flux_basisfunctions) flux_basisfunctions->setParams(x + getNDriveParams(), 0);

}

void Oscillator::getControlParams(double* x){

  // Copy p/q and flux parameterizations into the global parameter block x.
  int skip = 0;
  for (size_t f=0; f<carrier_freq.size(); f++) {
    if (drive_basisfunctions_re) drive_basisfunctions_re->getParams(x + skip, f);
    skip += drive_basisfunctions_re ? drive_basisfunctions_re->getNparams(f) : 0;
    if (drive_basisfunctions_im) drive_basisfunctions_im->getParams(x + skip, f);
    skip += drive_basisfunctions_im ? drive_basisfunctions_im->getNparams(f) : 0;
  }
  if (flux_basisfunctions) flux_basisfunctions->getParams(x + getNDriveParams(), 0);
}

int Oscillator::evalDriveControl(const double t, double* p_ptr, double* q_ptr) {
  double flux_dummy = 0.0;
  return evalControl(t, p_ptr, q_ptr, &flux_dummy);
}


double Oscillator::evalControlVariation(){

  double var_reg = 0.0;
  if (drive_basisfunctions_re) var_reg += drive_basisfunctions_re->computeVariation();
  if (drive_basisfunctions_im) var_reg += drive_basisfunctions_im->computeVariation();
  if (flux_basisfunctions) var_reg += flux_basisfunctions->computeVariation();

  return var_reg;
}

void Oscillator::evalControlVariationDiff(Vec G, double var_reg_bar, int skip_to_oscillator){

    PetscScalar* grad; 
    VecGetArray(G, &grad);

    // pass the portion of the gradient that corresponds to this oscillator
    int skip = skip_to_oscillator;
    if (drive_basisfunctions_re) {
      drive_basisfunctions_re->computeVariation_diff(grad + skip, var_reg_bar);
      skip += drive_basisfunctions_re->getNparams();
    }
    if (drive_basisfunctions_im) {
      drive_basisfunctions_im->computeVariation_diff(grad + skip, var_reg_bar);
      skip += drive_basisfunctions_im->getNparams();
    }
    if (flux_basisfunctions) {
      flux_basisfunctions->computeVariation_diff(grad + skip, var_reg_bar);
    }
    VecRestoreArray(G, &grad);
}

int Oscillator::evalControl(const double t, double* p_ptr, double* q_ptr, double* flux_ptr) {
  
  // Evaluate p/q drives.
  double p_drive = 0.0;
  double q_drive = 0.0;
  // Iterate over carrier frequencies.
  for (size_t f=0; f < carrier_freq.size(); f++) {
    // Evaluate the spline functions. Blt1= sum_i alpha_i^1 B_i(t), Blt2 = sum_i alpha_i^2 B_i(t)
    double Blt1 = drive_basisfunctions_re ? drive_basisfunctions_re->evaluate(t, static_cast<int>(f)) : 0.0;
    double Blt2 = drive_basisfunctions_im ? drive_basisfunctions_im->evaluate(t, static_cast<int>(f)) : 0.0;
    /* Mix in the carrier waves */
    double cos_omt = cos(carrier_freq[f]*t);
    double sin_omt = sin(carrier_freq[f]*t);
    p_drive += cos_omt * Blt1 - sin_omt * Blt2; 
    q_drive += sin_omt * Blt1 + cos_omt * Blt2;
  }
  *p_ptr = p_drive;
  *q_ptr = q_drive;


  // Evaluate flux control
  *flux_ptr = flux_basisfunctions ? flux_basisfunctions->evaluate(t, 0) : 0.0;

  return 0;
}

int Oscillator::evalControl_diff(const double t, double* grad, const double pbar, const double qbar, const double fbar) {

  // First, accumulate drive-channel sensitivity
  int skip = 0;
  for (size_t f = 0; f < carrier_freq.size(); f++) {

    double cos_omt = cos(carrier_freq[f]*t);
    double sin_omt = sin(carrier_freq[f]*t);
    double Blt1bar = sin_omt*qbar + cos_omt*pbar;
    double Blt2bar = cos_omt*qbar - sin_omt*pbar;

    /* Derivative with respect to control coefficients. */
    if (drive_basisfunctions_re) {
      drive_basisfunctions_re->derivative(t, static_cast<int>(f), grad + skip, Blt1bar);
      skip += drive_basisfunctions_re->getNparams(f);
    }
    if (drive_basisfunctions_im) {
      drive_basisfunctions_im->derivative(t, static_cast<int>(f), grad + skip, Blt2bar);
      skip += drive_basisfunctions_im->getNparams(f);
    }
  }

  // Then, accumulate flux-channel sensitivity in the flux parameter block
  if (flux_basisfunctions) {
    int skip = 0;
    if (drive_basisfunctions_re) skip += drive_basisfunctions_re->getNparams();
    if (drive_basisfunctions_im) skip += drive_basisfunctions_im->getNparams();
    flux_basisfunctions->derivative(t, 0, grad + skip, fbar);
  }

  return 0;
}

int Oscillator::evalDriveControl_diff(const double t, double* grad_for_this_oscillator, const double pbar, const double qbar) {
  return evalControl_diff(t, grad_for_this_oscillator, pbar, qbar, 0.0);
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
    else { // Schroedinger solver: += i * |psi_i|^2
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
    if (decoherence_type != DecoherenceType::NONE) { // Lindblad solver
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

  /* Gather population from all PETSc processors. */
  for (size_t i=0; i<mypop.size(); i++) {pop[i] = mypop[i];}
  MPI_Allreduce(mypop.data(), pop.data(), static_cast<int>(nlevels), MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
}
