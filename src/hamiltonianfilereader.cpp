#include "hamiltonianfilereader.hpp"
#include <optional>

HamiltonianFileReader::HamiltonianFileReader(){
}


HamiltonianFileReader::HamiltonianFileReader(std::optional<std::string> hamiltonian_file_Hsys_, std::optional<std::string> hamiltonian_file_Hc_, DecoherenceType decoherence_type_, PetscInt dim_rho_, bool quietmode_) {

  decoherence_type = decoherence_type_;
  dim_rho = dim_rho_;
  hamiltonian_file_Hsys = hamiltonian_file_Hsys_;
  hamiltonian_file_Hc = hamiltonian_file_Hc_;
  quietmode=quietmode_;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank_world);
  MPI_Comm_size(MPI_COMM_WORLD, &mpisize_world);
}

HamiltonianFileReader::~HamiltonianFileReader(){
}

void HamiltonianFileReader::receiveHsys(Mat& Ad, Mat& Bd){

  if (!hamiltonian_file_Hsys.has_value()) return;

  MatSetOption(Bd, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);


  // Read and broadcast matrix data 
  std::vector<PetscInt> rows, cols;
  std::vector<PetscScalar> real_vals, imag_vals;

  // Only rank 0 reads the file 
  if (mpirank_world == 0) {
    std::ifstream infile(hamiltonian_file_Hsys.value());
    if (!infile.is_open()) {
        std::cerr << "Could not open " << hamiltonian_file_Hsys.value() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        PetscInt row, col;
        double real, imag;
        if (!(iss >> row >> col >> real >> imag)) continue;
        rows.push_back(row);
        cols.push_back(col);
        real_vals.push_back(real);
        imag_vals.push_back(imag);
    }
    infile.close();
  }

  // Broadcast the size and data to all ranks
  int num_entries = rows.size();
  MPI_Bcast(&num_entries, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (mpirank_world != 0) {
    rows.resize(num_entries);
    cols.resize(num_entries);
    real_vals.resize(num_entries);
    imag_vals.resize(num_entries);
  }

  MPI_Bcast(rows.data(), num_entries, MPIU_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(cols.data(), num_entries, MPIU_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(real_vals.data(), num_entries, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(imag_vals.data(), num_entries, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // Now all ranks set their local values
  PetscInt ilow, iupp;
  MatGetOwnershipRange(Ad, &ilow, &iupp); 
  for (int i = 0; i < num_entries; i++) {
    PetscInt row = rows[i]; 
    PetscInt col = cols[i];
    double real = real_vals[i]; 
    double imag = imag_vals[i];

    // Assemble: Ad = Real(-i*Hsys) = Imag(Hsys)
    // Assemble: Bd = Imag(-i*Hsys) = -Real(Hsys)
    if (decoherence_type == DecoherenceType::NONE) {
      // Schroedinger 
      if (fabs(imag) > 1e-15 && ilow <= row && row < iupp) MatSetValue(Ad, row, col, imag, INSERT_VALUES);
      if (fabs(real) > 1e-15 && ilow <= row && row < iupp) MatSetValue(Bd, row, col, -real, INSERT_VALUES);
    } else {
      // Lindblad: 
      // Vectorized Ad = Real(I_N \kron (-iH) - (-iH)^T \kron I_N) = I_N \kron Im(Hsys) - Im(Hsys)^T \kron I_N
      // Vectorized Bd = Imag(I_N \kron (-iH) - (-iH)^T \kron I_N) = - I_n \kron Real(Hsys) + Real(Hsys)^T \kron I_N
      for (PetscInt k = 0; k < dim_rho; k++) {
        PetscInt rowk = row + dim_rho * k;
        PetscInt colk = col + dim_rho * k;
        if (fabs(imag) > 1e-15 && ilow <= rowk && rowk < iupp) MatSetValue(Ad, rowk, colk, imag, ADD_VALUES);
        if (fabs(real) > 1e-15 && ilow <= rowk && rowk < iupp) MatSetValue(Bd, rowk, colk, -real, ADD_VALUES);
        rowk = col * dim_rho + k;
        colk = row * dim_rho + k;
        if (fabs(imag) > 1e-15 && ilow <= rowk && rowk < iupp) MatSetValue(Ad, rowk, colk, -imag, ADD_VALUES);
        if (fabs(real) > 1e-15 && ilow <= rowk && rowk < iupp) MatSetValue(Bd, rowk, colk, real, ADD_VALUES);
      }
    }
  }
}

void HamiltonianFileReader::receiveHc(std::vector<Mat>& Ac_vec, std::vector<Mat>& Bc_vec){

  if (!hamiltonian_file_Hc.has_value()) return;

  for (size_t i=0; i<Ac_vec.size(); ++i) {
    MatSetOption(Ac_vec[i], MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    MatSetOption(Bc_vec[i], MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  }

  // Read and broadcast matrix data 
  std::vector<PetscInt> rows, cols, oscs;
  std::vector<PetscScalar> real_vals, imag_vals;

  // Only rank 0 reads the file
  if (mpirank_world == 0) {
    std::ifstream infile(hamiltonian_file_Hc.value());
    if (!infile.is_open()) {
        std::cerr << "Could not open " << hamiltonian_file_Hc.value() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        int osc;
        PetscInt row, col;
        double real, imag;
        if (!(iss >> osc >> row >> col >> real >> imag)) continue;
        oscs.push_back(osc);
        rows.push_back(row);
        cols.push_back(col);
        real_vals.push_back(real);
        imag_vals.push_back(imag);
    }
    infile.close();
  }

  // Broadcast the size and data to all ranks
  int num_entries = rows.size();
  MPI_Bcast(&num_entries, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (mpirank_world != 0) {
    rows.resize(num_entries);
    cols.resize(num_entries);
    oscs.resize(num_entries);
    real_vals.resize(num_entries);
    imag_vals.resize(num_entries);
  }

  MPI_Bcast(rows.data(), num_entries, MPIU_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(cols.data(), num_entries, MPIU_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(oscs.data(), num_entries, MPIU_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(real_vals.data(), num_entries, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(imag_vals.data(), num_entries, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // Now all ranks set their local values
  PetscInt ilow, iupp;
  MatGetOwnershipRange(Ac_vec[0], &ilow, &iupp); 
  for (int i = 0; i < num_entries; i++) {
    PetscInt row = rows[i]; 
    PetscInt col = cols[i];
    PetscInt osc = oscs[i];
    double real = real_vals[i]; 
    double imag = imag_vals[i];
    // Assemble: Ac = Real(-i*Hc) = Imag(Hc)
    // Assemble: Bc = Imag(-i*Hc) = -Real(Hc)
    if (decoherence_type == DecoherenceType::NONE) {
      // Schroedinger
      if (fabs(imag) > 1e-15 && ilow <= row && row < iupp) MatSetValue(Ac_vec[osc], row, col, imag, ADD_VALUES);
      if (fabs(real) > 1e-15 && ilow <= row && row < iupp) MatSetValue(Bc_vec[osc], row, col, -real, ADD_VALUES);
    } else {
      // Lindblad: 
      // Vectorized Ac = Real(I_N \kron (-iH) - (-iH)^T \kron I_N) = I_N \kron Im(Hsys) - Im(Hsys)^T \kron I_N
      // Vectorized Bc = Imag(I_N \kron (-iH) - (-iH)^T \kron I_N) = - I_n \kron Real(Hsys) + Real(Hsys)^T \kron I_N
      for (PetscInt k = 0; k < dim_rho; k++) {
        PetscInt rowk = row + dim_rho * k;
        PetscInt colk = col + dim_rho * k;
        if (fabs(imag) > 1e-15 && ilow <= rowk && rowk < iupp) MatSetValue(Ac_vec[osc], rowk, colk, imag, ADD_VALUES);
        if (fabs(real) > 1e-15 && ilow <= rowk && rowk < iupp) MatSetValue(Bc_vec[osc], rowk, colk, -real, ADD_VALUES);
        rowk = col * dim_rho + k;
        colk = row * dim_rho + k;
        if (fabs(imag) > 1e-15 && ilow <= rowk && rowk < iupp) MatSetValue(Ac_vec[osc], rowk, colk, -imag, ADD_VALUES);
        if (fabs(real) > 1e-15 && ilow <= rowk && rowk < iupp) MatSetValue(Bc_vec[osc], rowk, colk, real, ADD_VALUES);
      }
    }
  }
}