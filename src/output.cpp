#include "output.hpp"
#include "defs.hpp"
#include <vector>

Output::Output(){
  mpirank_world = -1;
  mpirank_petsc = -1;
  mpirank_init  = -1;
  output_timestep_stride = 0;
  quietmode = false;
  trajectory_initid = -1;
}

Output::Output(const Config& config, MPI_Comm comm_petsc, MPI_Comm comm_init, bool quietmode_) : Output() {

  /* Get communicator ranks */
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank_world);
  MPI_Comm_rank(comm_petsc, &mpirank_petsc);
  MPI_Comm_size(comm_petsc, &mpisize_petsc);
  MPI_Comm_rank(comm_init, &mpirank_init);

  /* Reduced output */
  quietmode = quietmode_;

  /* Store number of oscillators */
  noscillators = config.getNumOsc();

  /* Create Data directory */
  output_dir = config.getOutputDirectory();
  if (mpirank_world == 0) {
    mkdir(output_dir.c_str(), 0777);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  /* Prepare output for optimizer */
  output_optimization_stride = config.getOutputOptimizationStride();
  output_timestep_stride = config.getOutputTimestepStride();
  if (mpirank_world == 0) {
    char filename[255];
    snprintf(filename, 254, "%s/optim_history.dat", output_dir.c_str());
    optimfile = fopen(filename, "w");
    if (optimfile == nullptr) {
      printf("ERROR: Could not open file %s\n", filename);
      exit(1);
    }
    fprintf(optimfile, "#\"iter\"    \"Objective\"           \"||Pr(grad)||\"           \"LS step\"           \"F_avg\"           \"Terminal cost\"         \"Tikhonov-regul\"        \"Penalty-term\"          \"State variation\"        \"Energy-term\"           \"Control variation\"\n");
  } 

  /* Check which output should be written to files (applies them to all oscillators) */
  writeFullState = false;
  writeExpectedEnergy_comp = false;
  writePopulation_comp = false;
  writeExpectedEnergy = false;
  writePopulation = false;
  output_observables = config.getOutputObservables();
  for (auto type : output_observables) { // iterates over output types
    switch (type) {
      case OutputType::EXPECTED_ENERGY:
        writeExpectedEnergy = true;
        break;
      case OutputType::EXPECTED_ENERGY_COMPOSITE:
        writeExpectedEnergy_comp = true;
        break;
      case OutputType::POPULATION:
        writePopulation = true;
        break;
      case OutputType::POPULATION_COMPOSITE:
        writePopulation_comp = true;
        break;
      case OutputType::FULLSTATE:
        writeFullState = true;
        break;
    }
  }

  // Flag for writing control trajectory. Currently always true.
  writeControls = true;
}


Output::~Output(){
  if (mpirank_world == 0 && !quietmode) printf("Output directory: %s\n", output_dir.c_str());
  if (mpirank_world == 0) fclose(optimfile);
}


void Output::writeOptimFile(int optim_iter, double objective, double gnorm, double stepsize, double Favg, double costT, double tikh_regul, double penalty_leakage, double penalty_dpdm, double penalty_energy, double penalty_variation, double penalty_weightedcost){

  if (mpirank_world == 0){
    fprintf(optimfile, "%05d  %1.14e  %1.14e  %.8f  %1.14e  %1.14e  %1.14e  %1.14e  %1.14e  %1.14e  %1.14e\n", optim_iter, objective, gnorm, stepsize, Favg, costT, tikh_regul, penalty_leakage + penalty_weightedcost, penalty_dpdm, penalty_energy, penalty_variation);
    fflush(optimfile);
  } 
}

void Output::writeGradient(Vec grad){
  char filename[255];  
  PetscInt ngrad;
  VecGetSize(grad, &ngrad);

  if (mpirank_world == 0) {
    /* Print current gradients to file */
    FILE *file;
    // sprintf(filename, "%s/grad_iter%04d.dat", output_dir.c_str(), optim_iter);
    snprintf(filename, 254, "%s/grad.dat", output_dir.c_str());
    file = fopen(filename, "w");
    if (file == nullptr) {
      printf("ERROR: Could not open file %s\n", filename);
      exit(1);
    }
    const PetscScalar* grad_ptr;
    VecGetArrayRead(grad, &grad_ptr);
    for (int i=0; i<ngrad; i++){
      fprintf(file, "%1.14e\n", grad_ptr[i]);
    }
    fclose(file);
    VecRestoreArrayRead(grad, &grad_ptr);
    // if (!quietmode) printf("File written: %s\n", filename);
  }
}

void Output::writeControlParams(Vec params){

  /* Write control parameters */
  if ( mpirank_world == 0 ) { 

    char filename[255];
    PetscInt ndesign;
    VecGetSize(params, &ndesign);

    /* Print current parameters to file */
    FILE *file;
    snprintf(filename, 254, "%s/params.dat", output_dir.c_str());
    file = fopen(filename, "w");
    if (file == nullptr) {
      printf("ERROR: Could not open file %s\n", filename);
      exit(1);
    }
    const PetscScalar* params_ptr;
    VecGetArrayRead(params, &params_ptr);
    for (int i=0; i<ndesign; i++){
      fprintf(file, "%1.14e\n", params_ptr[i]);
    }
    fclose(file);
    VecRestoreArrayRead(params, &params_ptr);
    // if (!quietmode) printf("File written: %s\n", filename);
  }
}


void Output::resetTrajectoryData(int initid){
  trajectory_initid = initid;

  // Reset all trajectory buffers
  trajectory_times.clear();
  controls_Re_buffer.clear();
  controls_Im_buffer.clear();
  expected_energy_buffer.clear();
  expected_energy_comp_buffer.clear();
  population_buffer.clear();
  population_comp_buffer.clear();
  fullstate_re_buffer.clear();
  fullstate_im_buffer.clear();
}

void Output::evalTrajectoryData(int timestep, double time, const Vec state, MasterEq* mastereq){

  if (output_timestep_stride <= 0 ) {
    return;
  }

  /* Write output only every <num> time-steps */
  if (timestep % output_timestep_stride == 0) {

    if (mpirank_petsc == 0) {
      trajectory_times.push_back(time);
    }

    /* Grab controls from the master equation  */
    if (writeControls) {
      controls_Re_buffer.push_back(mastereq->getRHSctx()->control_Re);
      controls_Im_buffer.push_back(mastereq->getRHSctx()->control_Im);
    }

    /* Compute expected energy level */
    if (writeExpectedEnergy) {
      std::vector<double> vals(noscillators, 0.0);
      for (size_t iosc = 0; iosc < noscillators; iosc++) {
        vals[iosc] = mastereq->getOscillator(iosc)->expectedEnergy(state);
      }
      if (mpirank_petsc == 0) expected_energy_buffer.push_back(vals);
    }
    if (writeExpectedEnergy_comp) {
      double expected_comp = mastereq->expectedEnergy(state);
      if (mpirank_petsc == 0) expected_energy_comp_buffer.push_back(expected_comp);
    }

    /* Write population to file */
    if (writePopulation) {
      std::vector<std::vector<double>> vals(noscillators);
      for (size_t iosc = 0; iosc < noscillators; iosc++) {
        // std::vector<double> pop (mastereq->getOscillator(iosc)->getNLevels(), 0.0);
        vals[iosc].resize(mastereq->getOscillator(iosc)->getNLevels(), 0.0);
        mastereq->getOscillator(iosc)->population(state, vals[iosc]);
      }
      if (mpirank_petsc == 0) population_buffer.push_back(vals);
    }
    if (writePopulation_comp) {
      std::vector<double> population_comp; 
      mastereq->population(state, population_comp);
      if (mpirank_petsc == 0) population_comp_buffer.push_back(population_comp);
    }

    /* Write full state to file. Currently not available if Petsc-parallel */
    if (writeFullState && mpisize_petsc == 1) {
      /* TODO: Make this work in parallel! */
      /* Gather the vector from all petsc processors onto the first one */
      // VecScatterCreateToZero(x, &scat, &xseq);
      // VecScatterBegin(scat, u->x, xseq, INSERT_VALUES, SCATTER_FORWARD);
      // VecScatterEnd(scat, u->x, xseq, INSERT_VALUES, SCATTER_FORWARD);

      /* On first petsc rank, write full state vector to file */
      if (mpirank_petsc == 0) {
        const PetscScalar *x;
        VecGetArrayRead(state, &x);
        std::vector<double> re(mastereq->getDim());
        std::vector<double> im(mastereq->getDim());
        for (int i=0; i<mastereq->getDim(); i++) {
          re[i] = x[i];
          im[i] = x[i + mastereq->getDim()];
        }
        fullstate_re_buffer.push_back(std::move(re));
        fullstate_im_buffer.push_back(std::move(im));
        VecRestoreArrayRead(state, &x);
      }
      /* Destroy scatter context and vector */
      // VecScatterDestroy(&scat);
      // VecDestroy(&xseq); // TODO create and destroy scatter and xseq in contructor/destructor
    }
  }
}

void Output::writeTrajectoryData(){
  char filename[255];

  if (mpirank_petsc == 0) {
    size_t ntimepoints = trajectory_times.size();

    if (writeControls) {
      if (mpirank_world == 0) {// Only one proc write controls
        for (size_t iosc = 0; iosc < noscillators; iosc++) {
          snprintf(filename, 254, "%s/control%zu.dat", output_dir.c_str(), iosc);
          FILE* controlfile = fopen(filename, "w");
          if (controlfile == nullptr) {
            printf("ERROR: Could not open file %s\n", filename);
            exit(1);
          }
          fprintf(controlfile, "#\"time\"         \"p(t) (rotating)\"          \"q(t) (rotating)\"\n");
          for (size_t s = 0; s < ntimepoints && s < controls_Re_buffer.size() && s < controls_Im_buffer.size(); s++) {
            fprintf(controlfile, "% 1.8f   % 1.14e   % 1.14e \n", trajectory_times[s], controls_Re_buffer[s][iosc]/(2.0*M_PI), controls_Im_buffer[s][iosc]/(2.0*M_PI));
          }
          fclose(controlfile);
        }
      }
    }

    if (writeExpectedEnergy) {
      for (size_t i=0; i<noscillators; i++) {
        snprintf(filename, 254, "%s/expected%zu.iinit%04d.dat", output_dir.c_str(), i, trajectory_initid);
        FILE* expectedfile = fopen(filename, "w");
        if (expectedfile == nullptr) {
          printf("ERROR: Could not open file %s\n", filename);
          exit(1);
        }
        fprintf(expectedfile, "#\"time\"      \"expected energy level\"\n");
        for (size_t s = 0; s < ntimepoints; s++) {
          fprintf(expectedfile, "%.8f %1.14e\n", trajectory_times[s], expected_energy_buffer[s][i]);
        }
        fclose(expectedfile);
      }
    }

    if (writeExpectedEnergy_comp) {
      snprintf(filename, 254, "%s/expected_composite.iinit%04d.dat", output_dir.c_str(), trajectory_initid);
      FILE* expectedfile_comp = fopen(filename, "w");
      if (expectedfile_comp == nullptr) {
        printf("ERROR: Could not open file %s\n", filename);
        exit(1);
      }
      fprintf(expectedfile_comp, "#\"time\"      \"expected energy level\"\n");
      for (size_t s = 0; s < ntimepoints; s++) {
        fprintf(expectedfile_comp, "%.8f %1.14e\n", trajectory_times[s], expected_energy_comp_buffer[s]);
      }
      fclose(expectedfile_comp);
    }

    if (writePopulation) {
      for (size_t i=0; i<noscillators; i++) {
        snprintf(filename, 254, "%s/population%zu.iinit%04d.dat", output_dir.c_str(), i, trajectory_initid);
        FILE* populationfile = fopen(filename, "w");
        if (populationfile == nullptr) {
          printf("ERROR: Could not open file %s\n", filename);
          exit(1);
        }
        fprintf(populationfile, "#\"time\"      \"diagonal of the density matrix\"\n");
        for (size_t s = 0; s < ntimepoints; s++) {
          fprintf(populationfile, "%.8f ", trajectory_times[s]);
          for (size_t j = 0; j < population_buffer[s][i].size(); j++) {
            fprintf(populationfile, " %1.14e", population_buffer[s][i][j]);
          }
          fprintf(populationfile, "\n");
        }
        fclose(populationfile);
      }
    }

    if (writePopulation_comp) {
      snprintf(filename, 254, "%s/population_composite.iinit%04d.dat", output_dir.c_str(), trajectory_initid);
      FILE* populationfile_comp = fopen(filename, "w");
      if (populationfile_comp == nullptr) {
        printf("ERROR: Could not open file %s\n", filename);
        exit(1);
      }
      fprintf(populationfile_comp, "#\"time\"      \"population\"\n");
      for (size_t s = 0; s < ntimepoints && s < population_comp_buffer.size(); s++) {
        fprintf(populationfile_comp, "%.8f  ", trajectory_times[s]);
        for (size_t i=0; i<population_comp_buffer[s].size(); i++){
          fprintf(populationfile_comp, "%1.14e  ", population_comp_buffer[s][i]);
        }
        fprintf(populationfile_comp, "\n");
      }
      fclose(populationfile_comp);
    }

    if (writeFullState && mpisize_petsc == 1) {
      snprintf(filename, 254, "%s/rho_Re.iinit%04d.dat", output_dir.c_str(), trajectory_initid);
      FILE* ufile = fopen(filename, "w");
      if (ufile == nullptr) {
        printf("ERROR: Could not open file %s\n", filename);
        exit(1);
      }
      snprintf(filename, 254, "%s/rho_Im.iinit%04d.dat", output_dir.c_str(), trajectory_initid);
      FILE* vfile = fopen(filename, "w");
      if (vfile == nullptr) {
        printf("ERROR: Could not open file %s\n", filename);
        exit(1);
      }
      for (size_t s = 0; s < ntimepoints && s < fullstate_re_buffer.size() && s < fullstate_im_buffer.size(); s++) {
        fprintf(ufile,  "%.8f  ", trajectory_times[s]);
        fprintf(vfile,  "%.8f  ", trajectory_times[s]);
        for (size_t i=0; i<fullstate_re_buffer[s].size(); i++) {
          fprintf(ufile, "%1.10e  ", fullstate_re_buffer[s][i]);  
          fprintf(vfile, "%1.10e  ", fullstate_im_buffer[s][i]);  
        }
        fprintf(ufile, "\n");
        fprintf(vfile, "\n");
      }
      fclose(ufile);
      fclose(vfile);
    }
  }
}
