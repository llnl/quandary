#include "output.hpp"
#include "defs.hpp"
#include <vector>

Output::Output(){
  mpirank_world = -1;
  mpirank_petsc = -1;
  mpirank_init  = -1;
  output_timestep_stride = 0;
  quietmode = false;
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

  /* Reset flags and data file pointers */
  ufile = NULL;
  vfile = NULL;
  for (size_t i=0; i< noscillators; i++) expectedfile.push_back (NULL);
  for (size_t i=0; i< noscillators; i++) populationfile.push_back (NULL);

  /* Check which output should be written to files (applies them to all oscillators) */
  expectedfile_comp=NULL;
  populationfile_comp=NULL;
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
}


Output::~Output(){
  if (mpirank_world == 0 && !quietmode) printf("Output directory: %s\n", output_dir.c_str());
  if (mpirank_world == 0) fclose(optimfile);
}


void Output::writeOptimFile(int optim_iter, double objective, double gnorm, double stepsize, double Favg, double costT, double tikh_regul, double penalty_leakage, double penalty_dpdm, double penalty_energy, double penalty_variation, double penalty_weightedcost, double obj_robust){

  if (mpirank_world == 0){
    if (obj_robust > 0.0) fprintf(optimfile, "%05d  %1.14e  %1.14e  %.8f  %1.14e  %1.14e  %1.14e  %1.14e  %1.14e  %1.14e  %1.14e  %1.14e\n", optim_iter, objective, gnorm, stepsize, Favg, costT, tikh_regul, penalty_leakage + penalty_weightedcost, penalty_dpdm, penalty_energy, penalty_variation, obj_robust);
    else fprintf(optimfile, "%05d  %1.14e  %1.14e  %.8f  %1.14e  %1.14e  %1.14e  %1.14e  %1.14e  %1.14e  %1.14e\n", optim_iter, objective, gnorm, stepsize, Favg, costT, tikh_regul, penalty_leakage + penalty_weightedcost, penalty_dpdm, penalty_energy, penalty_variation);
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

void Output::writeControls(Vec params, MasterEq* mastereq, int ntime, double dt){

  /* Write controls every <outfreq> iterations */
  if ( mpirank_world == 0 ) { 

    char filename[255];
    PetscInt ndesign;
    VecGetSize(params, &ndesign);

    /* Print current parameters to file */
    FILE *file, *file_c;
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
    if (!quietmode) printf("File written: %s\n", filename);

    /* Print control to file for each oscillator */
    mastereq->setControlAmplitudes(params);
    for (size_t ioscil = 0; ioscil < mastereq->getNOscillators(); ioscil++) {
      snprintf(filename, 254, "%s/control%zu.dat", output_dir.c_str(), ioscil);
      file_c = fopen(filename, "w");
      if (file_c == nullptr) {
        printf("ERROR: Could not open file %s\n", filename);
        exit(1);
      }
      fprintf(file_c, "#\"time\"         \"p(t) (rotating)\"          \"q(t) (rotating)\"         \"f(t) (labframe)\"\n");

      /* Write every <num> timestep to file */
      for (int i=0; i<=ntime; i+=output_timestep_stride) {
        double time = i*dt; 

        double ReI, ImI, LabI;
        mastereq->getOscillator(ioscil)->evalControl(time, &ReI, &ImI);
        mastereq->getOscillator(ioscil)->evalControl_Labframe(time, &LabI);
        // Write control drives
        fprintf(file_c, "% 1.8f   % 1.14e   % 1.14e   % 1.14e \n", time, ReI/(2.0*M_PI), ImI/(2.0*M_PI), LabI/(2.0*M_PI));
     } // end of time loop 

      fclose(file_c);
      if (!quietmode) printf("File written: %s\n", filename);
    } // end of oscillator loop
  }
}


void Output::openTrajectoryDataFiles(std::string prefix, int initid){
  char filename[255];

  // On the first petsc rank, open required files and print header information
  if (mpirank_petsc == 0) {

    // Open files for expected energy per oscillator  
    if (writeExpectedEnergy) {
      for (size_t i=0; i<noscillators; i++) { 
        snprintf(filename, 254, "%s/expected%zu.iinit%04d.dat", output_dir.c_str(), i, initid);
        expectedfile[i] = fopen(filename, "w");
        if (expectedfile[i] == nullptr) {
          printf("ERROR: Could not open file %s\n", filename);
          exit(1);
        }
        fprintf(expectedfile[i], "#\"time\"      \"expected energy level\"\n");
      }
    }
    // Open file for expected energy of the full composite system
    if (writeExpectedEnergy_comp) {
      snprintf(filename, 254, "%s/expected_composite.iinit%04d.dat", output_dir.c_str(), initid);
      expectedfile_comp = fopen(filename, "w");
      if (expectedfile_comp == nullptr) {
        printf("ERROR: Could not open file %s\n", filename);
        exit(1);
      }
      fprintf(expectedfile_comp, "#\"time\"      \"expected energy level\"\n");
    }
    // Open files for populations per oscillator
    if (writePopulation) {
      for (size_t i=0; i<noscillators; i++) { 
        snprintf(filename, 254, "%s/population%zu.iinit%04d.dat", output_dir.c_str(), i, initid);
        populationfile[i] = fopen(filename, "w");
        if (populationfile[i] == nullptr) {
          printf("ERROR: Could not open file %s\n", filename);
          exit(1);
        }
        fprintf(populationfile[i], "#\"time\"      \"diagonal of the density matrix\"\n");
      }
    }
    // Open file for population for full composite system 
    if (writePopulation_comp) {
      snprintf(filename, 254, "%s/population_composite.iinit%04d.dat", output_dir.c_str(), initid);
      populationfile_comp = fopen(filename, "w");
      if (populationfile_comp == nullptr) {
        printf("ERROR: Could not open file %s\n", filename);
        exit(1);
      }
      fprintf(populationfile_comp, "#\"time\"      \"population\"\n");
    }
    // Open file for full vectorized state 
    if (writeFullState) {
      snprintf(filename, 254, "%s/%s_Re.iinit%04d.dat", output_dir.c_str(), prefix.c_str(), initid);
      ufile = fopen(filename, "w");
      if (ufile == nullptr) {
        printf("ERROR: Could not open file %s\n", filename);
        exit(1);
      }
      snprintf(filename, 254, "%s/%s_Im.iinit%04d.dat", output_dir.c_str(), prefix.c_str(), initid);
      vfile = fopen(filename, "w"); 
      if (vfile == nullptr) {
        printf("ERROR: Could not open file %s\n", filename);
        exit(1);
      }
    }
  }
}

void Output::writeTrajectoryDataFiles(int timestep, double time, const Vec state, MasterEq* mastereq){

  /* Write output only every <num> time-steps */
  if (timestep % output_timestep_stride == 0) {

    /* Write expected energy levels to file */
    if (writeExpectedEnergy) {
      for (size_t iosc = 0; iosc < expectedfile.size(); iosc++) {
        double expected = mastereq->getOscillator(iosc)->expectedEnergy(state);
        if (mpirank_petsc==0) fprintf(expectedfile[iosc], "%.8f %1.14e\n", time, expected);
      }
    }
    if (writeExpectedEnergy_comp) {
      double expected_comp = mastereq->expectedEnergy(state);
      if (mpirank_petsc==0) fprintf(expectedfile_comp, "%.8f %1.14e\n", time, expected_comp);
    }

    /* Write population to file */
    if (writePopulation) {
      for (size_t iosc = 0; iosc < populationfile.size(); iosc++) {
        std::vector<double> pop (mastereq->getOscillator(iosc)->getNLevels(), 0.0);
        mastereq->getOscillator(iosc)->population(state, pop);
        if (mpirank_petsc == 0) {
          fprintf(populationfile[iosc], "%.8f ", time);
          for (size_t i = 0; i<pop.size(); i++) {
            fprintf(populationfile[iosc], " %1.14e", pop[i]);
          }
          fprintf(populationfile[iosc], "\n");
        }
      }
    }
    if (writePopulation_comp) {
      std::vector<double> population_comp; 
      mastereq->population(state, population_comp);
      if (mpirank_petsc == 0) {
        fprintf(populationfile_comp, "%.8f  ", time);
        for (size_t i=0; i<population_comp.size(); i++){
          fprintf(populationfile_comp, "%1.14e  ", population_comp[i]);
        }
        fprintf(populationfile_comp, "\n");
      }
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
        fprintf(ufile,  "%.8f  ", time);
        fprintf(vfile,  "%.8f  ", time);
        const PetscScalar *x;
        VecGetArrayRead(state, &x);
        for (int i=0; i<mastereq->getDim(); i++) {
          fprintf(ufile, "%1.10e  ", x[i]);  
          fprintf(vfile, "%1.10e  ", x[i + mastereq->getDim()]);  
        }
        fprintf(ufile, "\n");
        fprintf(vfile, "\n");
        VecRestoreArrayRead(state, &x);
      }
      /* Destroy scatter context and vector */
      // VecScatterDestroy(&scat);
      // VecDestroy(&xseq); // TODO create and destroy scatter and xseq in contructor/destructor
    }
  }
}

void Output::closeTrajectoryDataFiles(){

  /* Close output data files */
  if (ufile != NULL) {
    fclose(ufile);
    ufile = NULL;
  }
  if (vfile != NULL) {
    fclose(vfile);
    vfile = NULL;
  }
  for (size_t i=0; i< expectedfile.size(); i++) {
    if (expectedfile[i] != NULL) {
      fclose(expectedfile[i]);
      expectedfile[i] = NULL;
    }
  }
  if (expectedfile_comp != NULL) fclose(expectedfile_comp);
  expectedfile_comp = NULL;
  for (size_t i=0; i< populationfile.size(); i++) {
    if (populationfile[i] != NULL) {
      fclose(populationfile[i]);
      populationfile[i] = NULL;
    }
  }
  if (populationfile_comp != NULL) fclose(populationfile_comp);
  populationfile_comp = NULL;
}
