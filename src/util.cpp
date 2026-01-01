#include "util.hpp"

// Suppress compiler warnings about unused parameters in code with #ifdef
#define UNUSED(expr) (void)(expr)

double sigmoid(double width, double x){
  return 1.0 / ( 1.0 + exp(-width*x) );
}

double sigmoid_diff(double width, double x){
  return sigmoid(width, x) * (1.0 - sigmoid(width, x)) * width;
}

double getRampFactor(const double time, const double tstart, const double tstop, const double tramp){

    // double eps = 1e-4; // Cutoff for sigmoid ramp 
    // double steep = log(1./eps - 1.) * 2. / tramp; // steepness of sigmoid such that ramp(x) small for x < - tramp/2
    // printf("steep eval %f\n", steep);

    double rampfactor = 0.0;
    if (time <= tstart + tramp) { // ramp up
      // double center = tstart + tramp/2.0;
      // rampfactor = sigmoid(steep, time - center);
      rampfactor =  1.0/tramp * time - tstart/ tramp;
    }
    else if (tstart + tramp <= time && 
            time <= tstop - tramp) { // center
      rampfactor = 1.0;
    }
    else if (time >= tstop - tramp && time <= tstop) { // down
      // double center = tstop - tramp/2.0;
      // rampfactor = sigmoid(steep, -(time - center));
      // steep = 1842.048073;
      // steep = 1000.0;
      rampfactor =  -1.0/tramp * time + tstop / tramp;
    }

    // If ramp time is larger than total amount of time, turn off control:
    if (tstop < tstart + 2*tramp) rampfactor = 0.0;

    return rampfactor;
}

double getRampFactor_diff(const double time, const double tstart, const double tstop, const double tramp){

    // double eps = 1e-4; // Cutoff for sigmoid ramp 
    // double steep = log(1./eps - 1.) * 2. / tramp; // steepness of sigmoid such that ramp(x) small for x < - tramp/2
    // printf("steep der %f\n", steep);

    double dramp_dtstop= 0.0;
    if (time <= tstart + tramp) { // ramp up
      dramp_dtstop = 0.0;
    }
    else if (tstart + tramp <= time && 
            time <= tstop - tramp) { // center
      dramp_dtstop = 0.0;
    }
    else if (time >= tstop - tramp && time <= tstop) { // down
      // double center = tstop - tramp/2.0;
      // dramp_dtstop = sigmoid_diff(steep, -(time - center));
      // steep = 1842.048073;
      dramp_dtstop = 1.0/tramp;
    }
    
    // If ramp time is larger than total amount of time, turn off control:
    if (tstop < tstart + 2*tramp) dramp_dtstop= 0.0;

    return dramp_dtstop;
}


PetscInt getVecID(const PetscInt row, const PetscInt col, const PetscInt dim){
  return row + col * dim;  
} 


PetscInt mapEssToFull(const PetscInt i, const std::vector<size_t> &nlevels, const std::vector<size_t> &nessential){

  PetscInt id = 0;
  PetscInt index = i;
  for (size_t iosc = 0; iosc<nlevels.size()-1; iosc++){
    PetscInt postdim = 1;
    PetscInt postdim_ess = 1;
    for (size_t j = iosc+1; j<nlevels.size(); j++){
      postdim *= nlevels[j];
      postdim_ess *= nessential[j];
    }
    PetscInt iblock = index / postdim_ess;
    index = index % postdim_ess;
    // move id to that block
    id += iblock * postdim;  
  }
  // move to index in last block
  id += index;

  return id;
}

PetscInt mapFullToEss(const PetscInt i, const std::vector<size_t> &nlevels, const std::vector<size_t> &nessential){

  PetscInt id = 0;
  PetscInt index = i;
  for (size_t iosc = 0; iosc<nlevels.size(); iosc++){
    PetscInt postdim = 1;
    PetscInt postdim_ess = 1;
    for (size_t j = iosc+1; j<nlevels.size(); j++){
      postdim *= nlevels[j];
      postdim_ess *= nessential[j];
    }
    PetscInt iblock = index / postdim;
    index = index % postdim;
    if (iblock >= static_cast<PetscInt>(nessential[iosc])) return -1; // this row/col belongs to a guard level, no mapping defined. 
    // move id to that block
    id += iblock * postdim_ess;  
  }

  return id;
}



// void projectToEss(Vec state,const std::vector<int> &nlevels, const std::vector<int> &nessential){

//   /* Get dimensions */
//   int dim_rho = 1;
//   for (int i=0; i<nlevels.size(); i++){
//     dim_rho *=nlevels[i];
//   }

//   /* Get local ownership of the state */
//   PetscInt ilow, iupp;
//   VecGetOwnershipRange(state, &ilow, &iupp);

//   /* Iterate over rows of system matrix, check if it corresponds to an essential level, and if not, set this row and colum to zero */
//   int reID, imID;
//   for (int i=0; i<dim_rho; i++) {
//     // zero out row and column if this does not belong to an essential level
//     if (!isEssential(i, nlevels, nessential)) { 
//       for (int j=0; j<dim_rho; j++) {
//         // zero out row
//         reID = getIndexReal(getVecID(i,j,dim_rho));
//         imID = getIndexImag(getVecID(i,j,dim_rho));
//         if (ilow <= reID && reID < iupp) VecSetValue(state, reID, 0.0, INSERT_VALUES);
//         if (ilow <= imID && imID < iupp) VecSetValue(state, imID, 0.0, INSERT_VALUES);
//         // zero out colum
//         reID = getIndexReal(getVecID(j,i,dim_rho));
//         imID = getIndexImag(getVecID(j,i,dim_rho));
//         if (ilow <= reID && reID < iupp) VecSetValue(state, reID, 0.0, INSERT_VALUES);
//         if (ilow <= imID && imID < iupp) VecSetValue(state, imID, 0.0, INSERT_VALUES);
//       }
//     } 
//   }
//   VecAssemblyBegin(state);
//   VecAssemblyEnd(state);


// }

int isEssential(const int i, const std::vector<size_t> &nlevels, const std::vector<size_t> &nessential) {

  int isEss = 1;
  int index = i;
  for (size_t iosc = 0; iosc < nlevels.size(); iosc++){

    PetscInt postdim = 1;
    for (size_t j = iosc+1; j<nlevels.size(); j++){
      postdim *= nlevels[j];
    }
    int itest = (int) index / postdim;
    // test if essential for this oscillator
    if (itest >= static_cast<int>(nessential[iosc])) {
      isEss = 0;
      break;
    }
    index = index % postdim;
  }

  return isEss; 
}

int isGuardLevel(const int i, const std::vector<size_t> &nlevels, const std::vector<size_t> &nessential){
  int isGuard =  0;
  int index = i;
  for (size_t iosc = 0; iosc < nlevels.size(); iosc++){

    PetscInt postdim = 1;
    for (size_t j = iosc+1; j<nlevels.size(); j++){
      postdim *= nlevels[j];
    }
    int itest = (int) index / postdim;   // floor(i/n_post)
    // test if this is a guard level for this oscillator
    if (itest == static_cast<int>(nlevels[iosc]) - 1 && itest >= static_cast<int>(nessential[iosc])) {  // last energy level for system 'iosc'
      isGuard = 1;
      break;
    }
    index = index % postdim;
  }

  return isGuard;
}

PetscErrorCode Ikron(const Mat A,const  int dimI, const double alpha, Mat *Out, InsertMode insert_mode){

    PetscInt ierr;
    PetscInt ncols;
    const PetscInt* cols; 
    const PetscScalar* Avals;
    PetscInt* shiftcols;
    PetscScalar* vals;
    PetscInt dimA;
    // PetscInt dimOut;
    // PetscInt nonzeroOut;
    PetscInt rowID;

    MatGetSize(A, &dimA, NULL);

    ierr = PetscMalloc1(dimA, &shiftcols); CHKERRQ(ierr);
    ierr = PetscMalloc1(dimA, &vals); CHKERRQ(ierr);

    /* Loop over dimension of I */
    for (PetscInt i = 0; i < dimI; i++){

        /* Set the diagonal block (i*dimA)::(i+1)*dimA */
        for (PetscInt j=0; j<dimA; j++){
            MatGetRow(A, j, &ncols, &cols, &Avals);
            rowID = i*dimA + j;
            for (int k=0; k<ncols; k++){
                shiftcols[k] = cols[k] + i*dimA;
                vals[k] = Avals[k] * alpha;
            }
            MatSetValues(*Out, 1, &rowID, ncols, shiftcols, vals, insert_mode);
            MatRestoreRow(A, j, &ncols, &cols, &Avals);
        }

    }
    // MatAssemblyBegin(*Out, MAT_FINAL_ASSEMBLY);
    // MatAssemblyEnd(*Out, MAT_FINAL_ASSEMBLY);

    PetscFree(shiftcols);
    PetscFree(vals);
    return 0;
}

PetscErrorCode kronI(const Mat A, const int dimI, const double alpha, Mat *Out, InsertMode insert_mode){
    
    PetscInt ierr;
    PetscInt dimA;
    const PetscInt* cols; 
    const PetscScalar* Avals;
    PetscInt rowid;
    PetscInt colid;
    PetscScalar insertval;
    // PetscInt dimOut;
    // PetscInt nonzeroOut;
    PetscInt ncols;
    // MatInfo Ainfo;
    MatGetSize(A, &dimA, NULL);

    ierr = PetscMalloc1(dimA, &cols); CHKERRQ(ierr);
    ierr = PetscMalloc1(dimA, &Avals);

    /* Loop over rows in A */
    for (PetscInt i = 0; i < dimA; i++){
        MatGetRow(A, i, &ncols, &cols, &Avals);

        /* Loop over non negative columns in row i */
        for (PetscInt j = 0; j < ncols; j++){
            //printf("A: row = %d, col = %d, val = %f\n", i, cols[j], Avals[j]);
            
            // dimI rows. global row indices: i, i+dimI
            for (PetscInt k=0; k<dimI; k++) {
               rowid = i*dimI + k;
               colid = cols[j]*dimI + k;
               insertval = Avals[j] * alpha;
               MatSetValues(*Out, 1, &rowid, 1, &colid, &insertval, insert_mode);
              //  printf("Setting %d,%d %f\n", rowid, colid, insertval);
            }
        }
        MatRestoreRow(A, i, &ncols, &cols, &Avals);
    }

    // MatAssemblyBegin(*Out, MAT_FINAL_ASSEMBLY);
    // MatAssemblyEnd(*Out, MAT_FINAL_ASSEMBLY);

    PetscFree(cols);
    PetscFree(Avals);

    return 0;
}



PetscErrorCode AkronB(const Mat A, const Mat B, const double alpha, Mat *Out, InsertMode insert_mode){
    PetscInt Adim1, Adim2, Bdim1, Bdim2;
    MatGetSize(A, &Adim1, &Adim2);
    MatGetSize(B, &Bdim1, &Bdim2);

    PetscInt ncolsA, ncolsB;
    const PetscInt *colsA, *colsB;
    const double *valsA, *valsB;
    // Iterate over rows of A 
    for (PetscInt irowA = 0; irowA < Adim1; irowA++){
        // Iterate over non-zero columns in this row of A
        MatGetRow(A, irowA, &ncolsA, &colsA, &valsA);
        for (PetscInt j=0; j<ncolsA; j++) {
            PetscInt icolA = colsA[j];
            PetscScalar valA = valsA[j];
            /* put a B-block at position (irowA*Bdim1, icolA*Bdim2): */
            // Iterate over rows of B 
            for (PetscInt irowB = 0; irowB < Bdim1; irowB++){
                // Iterate over non-zero columns in this B-row
                MatGetRow(B, irowB, &ncolsB, &colsB, &valsB);
                for (PetscInt k=0; k< ncolsB; k++) {
                    PetscInt icolB = colsB[k];
                    PetscScalar valB = valsB[k];
                    /* Insert values in Out */
                    PetscInt rowOut = irowA*Bdim1 + irowB;
                    PetscInt colOut = icolA*Bdim2 + icolB;
                    PetscScalar valOut = valA * valB * alpha; 
                    MatSetValue(*Out, rowOut, colOut, valOut, insert_mode);
                }
                MatRestoreRow(B, irowB, &ncolsB, &colsB, &valsB);
            }
        }   
        MatRestoreRow(A, irowA, &ncolsA, &colsA, &valsA);
    }  

  return 0;
}


PetscErrorCode MatIsAntiSymmetric(Mat A, PetscReal tol, PetscBool *flag) {
  
  int ierr; 

  /* Create B = -A */
  Mat B;
  ierr = MatConvert(A, MATSAME, MAT_INITIAL_MATRIX, &B); CHKERRQ(ierr);
  ierr = MatScale(B, -1.0); CHKERRQ(ierr);

  /* Test if B^T = A */
  ierr = MatIsTranspose(B, A, tol, flag); CHKERRQ(ierr);

  /* Cleanup */
  ierr = MatDestroy(&B); CHKERRQ(ierr);

  return ierr;
}



PetscErrorCode StateIsHermitian(Vec x, PetscReal tol, PetscBool *flag) {
  int ierr;
  PetscInt i, j;

  /* TODO: Either make this work in Petsc-parallel, or add error exit if this runs in parallel. */
  
  /* Get u and v from x */
  PetscInt dim;
  ierr = VecGetSize(x, &dim); CHKERRQ(ierr);
  dim = dim/2;
  Vec u, v;
  IS isu, isv;

  PetscInt dimis = dim;
  ierr = ISCreateStride(PETSC_COMM_WORLD, dimis, 0, 1, &isu); CHKERRQ(ierr);
  ierr = ISCreateStride(PETSC_COMM_WORLD, dimis, dimis, 1, &isv); CHKERRQ(ierr);
  ierr = VecGetSubVector(x, isu, &u); CHKERRQ(ierr);
  ierr = VecGetSubVector(x, isv, &v); CHKERRQ(ierr);

  /* Init flags*/
  *flag = PETSC_TRUE;

  /* Check for symmetric u and antisymmetric v */
  const double *u_array;
  const double *v_array;
  double u_diff, v_diff;
  ierr = VecGetArrayRead(u, &u_array); CHKERRQ(ierr);
  ierr = VecGetArrayRead(v, &v_array); CHKERRQ(ierr);
  PetscInt N = sqrt(dim);
  for (i=0; i<N; i++) {
    for (j=i; j<N; j++) {
      u_diff = u_array[i*N+j] - u_array[j*N+i];
      v_diff = v_array[i*N+j] + v_array[j*N+i];
      if (fabs(u_diff) > tol || fabs(v_diff) > tol ) {
        *flag = PETSC_FALSE;
        break;
      }
    }
  }

  ierr = VecRestoreArrayRead(u, &u_array);
  ierr = VecRestoreArrayRead(v, &v_array);
  ierr = VecRestoreSubVector(x, isu, &u);
  ierr = VecRestoreSubVector(x, isv, &v);
  ISDestroy(&isu);
  ISDestroy(&isv);

  return ierr;
}



PetscErrorCode StateHasTrace1(Vec x, PetscReal tol, PetscBool *flag) {

  int ierr;
  PetscInt i;

  /* Get u and v from x */
  PetscInt dim;
  ierr = VecGetSize(x, &dim); CHKERRQ(ierr);
  PetscInt dimis = dim/2;
  Vec u, v;
  IS isu, isv;
  ierr = ISCreateStride(PETSC_COMM_WORLD, dimis, 0, 1, &isu); CHKERRQ(ierr);
  ierr = ISCreateStride(PETSC_COMM_WORLD, dimis, dimis, 1, &isv); CHKERRQ(ierr);
  ierr = VecGetSubVector(x, isu, &u); CHKERRQ(ierr);
  ierr = VecGetSubVector(x, isv, &v); CHKERRQ(ierr);

  /* Init flags*/
  *flag = PETSC_FALSE;
  PetscBool u_hastrace1 = PETSC_FALSE;
  PetscBool v_hasdiag0  = PETSC_FALSE;

  /* Check if diagonal of u sums to 1, and diagonal elements of v are 0 */ 
  const double *u_array;
  const double *v_array;
  double u_sum = 0.0;
  double v_sum = 0.0;
  ierr = VecGetArrayRead(u, &u_array); CHKERRQ(ierr);
  ierr = VecGetArrayRead(v, &v_array); CHKERRQ(ierr);
  PetscInt N = sqrt(dimis);
  for (i=0; i<N; i++) {
    u_sum += u_array[i*N+i];
    v_sum += fabs(v_array[i*N+i]);
  }
  if ( fabs(u_sum - 1.0) < tol ) u_hastrace1 = PETSC_TRUE;
  if ( fabs(v_sum      ) < tol ) v_hasdiag0  = PETSC_TRUE;

  /* Restore vecs */
  ierr = VecRestoreArrayRead(u, &u_array);
  ierr = VecRestoreArrayRead(v, &v_array);
  ierr = VecRestoreSubVector(x, isu, &u);
  ierr = VecRestoreSubVector(x, isv, &v);

  /* Answer*/
  if (u_hastrace1 && v_hasdiag0) {
    *flag = PETSC_TRUE;
  }
  
  /* Destroy vector strides */
  ISDestroy(&isu);
  ISDestroy(&isv);


  return ierr;
}



PetscErrorCode SanityTests(Vec x, double time){

  /* Sanity check. Be careful: This is costly! */
  printf("Trace check %f ...\n", time);
  PetscBool check;
  double tol = 1e-10;
  StateIsHermitian(x, tol, &check);
  if (!check) {
    printf("WARNING at t=%f: rho is not hermitian!\n", time);
    printf("\n rho :\n");
    VecView(x, PETSC_VIEWER_STDOUT_WORLD);
    exit(1);
  }
  else printf("IsHermitian check passed.\n");
  StateHasTrace1(x, tol, &check);
  if (!check) {
    printf("WARNING at t=%f: Tr(rho) is NOT one!\n", time);
    printf("\n rho :\n");
    VecView(x, PETSC_VIEWER_STDOUT_WORLD);
    exit(1);
  }
  else printf("Trace1 check passed.\n");

  return 0;
}


int read_vector(const char *filename, double *var, int dim, bool quietmode, int skiplines, const std::string testheader) {

  FILE *file;
  int success = 0;

  file = fopen(filename, "r");

  if (file != NULL) {
    if (!quietmode) printf("Reading file %s, starting from line %d.\n", filename, skiplines+1);

    /* Scip first <skiplines> lines */
    char buffer[51]; // need one extra element because fscanf adds a '\0' at the end
    for (int ix = 0; ix < skiplines; ix++) {
      int ret = fscanf(file, "%50[^\n]%*c", buffer); // // NOTE: &buffer[50] is a pointer to buffer[50] (i.e. its last element)
      if (ret == EOF) {
        printf("ERROR: EOF reached while skipping lines in file %s.\n", filename);
        fclose(file);
        return success;
      }
      // printf("Skipping %d lines: %s \n:", skiplines, buffer);
    }

    // Test the header, if given, and set vals to zero if header doesn't match.
    if (testheader.size()>0) {
      if (!quietmode) printf("Compare to Header '%s': ", testheader.c_str());
      // read one (first) line
      int ret = fscanf(file, "%50[^\n]%*c", buffer); // NOTE: &buffer[50] is a pointer to buffer[50] (i.e. its last element)
      std::string header = buffer;
      // Compare to testheader, return if it doesn't match 
      if (ret==EOF || header.compare(0,testheader.size(),testheader) != 0) {
        // printf("Header not found: %s != %s\n", header.c_str(), testheader.c_str());
        printf("Header not found.\n");
        for (int ix = 0; ix < dim; ix++) var[ix] = 0.0;
        fclose(file);
        return success;
      } else {
        if (!quietmode) printf(" Header correct! Reading now.\n");
      }
    }
 
    // printf("Either matching header, or no header given. Now reading lines \n");

    /* Read <dim> lines from file */
    for (int ix = 0; ix < dim; ix++) {
      double myval = 0.0;
      // read the next line
      int ret = fscanf(file, "%lf", &myval); 
      // if end of file, set remaining vars to zero
      if (ret == EOF){ 
        for (int j = ix; j<dim; j++) var[j] = 0.0;
        break;
      } else { // otherwise, set the value
        var[ix] = myval;
      }
      success = 1;
    }
  } else {
    printf("ERROR: Can't open file %s\n", filename);
    exit(1);
  }

  fclose(file);
  return success;
}


/* Compute eigenvalues */
int getEigvals(const Mat A, const int neigvals, std::vector<double>& eigvals, std::vector<Vec>& eigvecs){

UNUSED(A);
UNUSED(neigvals);
UNUSED(eigvals);
UNUSED(eigvecs);

int nconv = 0;
#ifdef WITH_SLEPC

  /* Create Slepc's eigensolver */
  EPS eigensolver;       
  EPSCreate(PETSC_COMM_WORLD, &eigensolver);
  EPSSetOperators(eigensolver, A, NULL);
  EPSSetProblemType(eigensolver, EPS_NHEP);
  EPSSetFromOptions(eigensolver);

  /* Number of requested eigenvalues */
  EPSSetDimensions(eigensolver,neigvals,PETSC_DEFAULT,PETSC_DEFAULT);

  // Solve eigenvalue problem
  int ierr = EPSSolve(eigensolver); CHKERRQ(ierr);

  /* Get information about convergence */
  int its, nev, maxit;
  EPSType type;
  double tol;
  EPSGetIterationNumber(eigensolver,&its);
  EPSGetType(eigensolver,&type);
  EPSGetDimensions(eigensolver,&nev,NULL,NULL);
  EPSGetTolerances(eigensolver,&tol,&maxit);
  EPSGetConverged(eigensolver, &nconv );

  PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n",type);
  PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenpairs: %D\n",nev);
  PetscPrintf(PETSC_COMM_WORLD," Number of iterations taken: %D / %D\n",its, maxit);
  PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g\n",(double)tol);
  PetscPrintf(PETSC_COMM_WORLD," Number of converged eigenpairs: %D\n\n",nconv);

  /* Allocate eigenvectors */
  Vec eigvec;
  MatCreateVecs(A, &eigvec, NULL);

  // Get the result
  double kr, ki, error;
  // printf("Eigenvalues: \n");
  for (int j=0; j<nconv; j++) {
      EPSGetEigenpair( eigensolver, j, &kr, &ki, eigvec, NULL);
      EPSComputeError( eigensolver, j, EPS_ERROR_RELATIVE, &error );
      // printf("%f + i%f (err %f)\n", kr, ki, error);

      /* Store the eigenpair */
      eigvals.push_back(kr);
      eigvecs.push_back(eigvec);
      if (ki != 0.0) printf("Warning: eigenvalue imaginary! : %f", ki);
  }
  // printf("\n");
  // EPSView(eigensolver, PETSC_VIEWER_STDOUT_WORLD);

  /* Clean up*/
  EPSDestroy(&eigensolver);
#endif
  return nconv;
}

// test if A+iB is a unitary matrix: (A+iB)^\dag (A+iB) = I!
bool isUnitary(const Mat V_re, const Mat V_im){
  Mat C, D;
  double norm;
  bool isunitary = true;

  // test: C=V_re^T V_re + Vim^TVim should be the identity!
  MatTransposeMatMult(V_re, V_re, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &C);
  MatTransposeMatMult(V_im, V_im, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &D);
  MatAXPY(C, 1.0, D, DIFFERENT_NONZERO_PATTERN); 
  MatShift(C, -1.0);
  MatNorm(C, NORM_FROBENIUS, &norm);
  if (norm > 1e-12) {
    printf("Unitary Test: V_re^TVre+Vim^TVim is not the identity! %1.14e\n", norm);
    // MatView(C, NULL);
    isunitary = false;
  } 
  MatDestroy(&C);
  MatDestroy(&D);

  // test: C=V_re^T V_im - Vre^TVim should be zero!
  MatTransposeMatMult(V_re, V_im, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &C);
  MatTransposeMatMult(V_im, V_re, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &D);
  MatAXPY(C, -1.0, D, DIFFERENT_NONZERO_PATTERN); 
  MatNorm(C, NORM_FROBENIUS, &norm);
  if (norm > 1e-12) {
    printf("Unitary Test: Vre^TVim - Vim^TVre is not zero! %1.14e\n", norm);
    // MatView(C,NULL);
    isunitary = false;
  }
  MatDestroy(&C);
  MatDestroy(&D);

  return isunitary;
}


void getEigenComplex(const Mat A_re, const Mat A_im, std::unique_ptr<std::vector<double>>& eigvals_re, std::unique_ptr<std::vector<double>>& eigvals_im, Mat& Evecs_re, Mat& Evecs_im, bool print){

  PetscInt dim;
  MatGetSize(A_re, &dim, NULL);

  /* Set up the larger matrix A = [A_re -A_im; A_im A_re] */
  Mat A;
  MatCreate(PETSC_COMM_WORLD, &A);
  MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, 2*dim, 2*dim);
  MatSetFromOptions(A);
  MatSetUp(A);

  // Number of eigenvalues to compute of the larger matrix A=[A_re -A_im; A_im A_re]
  int neigvals = 2*dim; 

  // Create storage of Evecs_re and Evecs_im
  if (Evecs_re) MatDestroy(&Evecs_re);
  if (Evecs_im) MatDestroy(&Evecs_im);
  MatCreate(PETSC_COMM_WORLD, &Evecs_re);
  MatSetSizes(Evecs_re, PETSC_DECIDE, PETSC_DECIDE, dim, dim);
  MatSetFromOptions(Evecs_re);
  MatSetUp(Evecs_re);
  MatCreate(PETSC_COMM_WORLD, &Evecs_im);
  MatSetSizes(Evecs_im, PETSC_DECIDE, PETSC_DECIDE, dim, dim);
  MatSetFromOptions(Evecs_im);
  MatSetUp(Evecs_im);

  for (PetscInt i=0; i<dim; i++){
    for (PetscInt j=0; j<dim; j++){
      double val_re, val_im;
      MatGetValue(A_re, i, j, &val_re);
      MatGetValue(A_im, i, j, &val_im);
      // Set A_ij
      MatSetValue(A, i, j, val_re, INSERT_VALUES);               // Top-left
      MatSetValue(A, i, j + dim, -val_im, INSERT_VALUES);       // Top-right
      MatSetValue(A, i + dim, j, val_im, INSERT_VALUES);        // Bottom-left
      MatSetValue(A, i + dim, j + dim, val_re, INSERT_VALUES);  // Bottom-right
    }
  }
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); 
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

  /* Create Slepc's eigensolver */
  EPS eigensolver;       
  EPSCreate(PETSC_COMM_WORLD, &eigensolver);
  EPSSetOperators(eigensolver, A, NULL);
  EPSSetProblemType(eigensolver, EPS_NHEP);
  EPSSetFromOptions(eigensolver);
  EPSSetDimensions(eigensolver,neigvals,PETSC_DEFAULT,PETSC_DEFAULT);

  // Solve eigenvalue problem
  EPSSolve(eigensolver);
  // EPSView(eigensolver, PETSC_VIEWER_STDOUT_WORLD);

  /* Get information about convergence */
  // int its, nev, maxit;
  // EPSType type;
  // double tol;
  // EPSGetIterationNumber(eigensolver,&its);
  // EPSGetType(eigensolver,&type);
  // EPSGetDimensions(eigensolver,&nev,NULL,NULL);
  // EPSGetTolerances(eigensolver,&tol,&maxit);
  int nconv = 0;
  EPSGetConverged(eigensolver, &nconv );
  if (nconv < neigvals){
    printf("ERROR: Only %d eigenvalues converged, but %d requested!\n", nconv, neigvals);
    exit(1);
  }
  // PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n",type);
  // PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenpairs: %D\n",nev);
  // PetscPrintf(PETSC_COMM_WORLD," Number of iterations taken: %D / %D\n",its, maxit);
  // PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g\n",(double)tol);
  // PetscPrintf(PETSC_COMM_WORLD," Number of converged eigenpairs: %D\n\n",nconv);

  // Get the result
  eigvals_re->clear();
  eigvals_im->clear();
  eigvals_re->resize(nconv);
  eigvals_im->resize(nconv);
  MatZeroEntries(Evecs_re);
  MatZeroEntries(Evecs_im);
  Vec eigvecLarge_re;
  Vec eigvecLarge_im;
  MatCreateVecs(A, &eigvecLarge_re, NULL);
  MatCreateVecs(A, &eigvecLarge_im, NULL);
  Vec v_re, v_im;
  MatCreateVecs(A_re, &v_re, NULL);
  MatCreateVecs(A_re, &v_im, NULL);
  double kr, ki;
  int found_eigenpairs = 0;
  for (int j=0; j<nconv; j++) {
    // Get the j-th eigenvalue
    EPSGetEigenvalue(eigensolver, j, &kr, &ki);
    // Get the eigenvector of the larger matrix
    EPSGetEigenvector( eigensolver, j, eigvecLarge_re, eigvecLarge_im);
    for (PetscInt i=0; i<2*dim; i++){
      double val_re, val_im;
      VecGetValues(eigvecLarge_re, 1, &i, &val_re);
      VecGetValues(eigvecLarge_im, 1, &i, &val_im);
    }
    // Error of the eigensolver for this eigenvalue
    // double error;
    // EPSComputeError( eigensolver, j, EPS_ERROR_RELATIVE, &error );
    // printf("%1.8e + i%1.8e (err %f)\n", kr, ki, error);

    // TODO: Check for multiplicity of an eigenvalue here? 

    // Set up the smaller vector v = evec[:dim] + 1j evec[dim:] and test whether this is an eigenvector of the complex matrix A_re + i A_im.
    // v_re = eigvecLarge_re[:dim] - eigvecLarge_im[dim:]
    // v_im = eigvecLarge_im[:dim] + eigvecLarge_re[dim:]
    VecZeroEntries(v_re);
    VecZeroEntries(v_im);
    for (PetscInt i=0; i<dim; i++){
      double val_re, val_im;
      PetscInt idx = i;
      VecGetValues(eigvecLarge_re, 1, &idx, &val_re);
      VecGetValues(eigvecLarge_im, 1, &idx, &val_im);
      VecSetValue(v_re, i, val_re, INSERT_VALUES);
      VecSetValue(v_im, i, val_im, INSERT_VALUES);
      idx += dim;
      VecGetValues(eigvecLarge_re, 1, &idx, &val_re);
      VecGetValues(eigvecLarge_im, 1, &idx, &val_im);
      VecSetValue(v_re, i, -val_im, ADD_VALUES);
      VecSetValue(v_im, i, val_re, ADD_VALUES);
    }
    VecAssemblyBegin(v_re); VecAssemblyEnd(v_re);
    VecAssemblyBegin(v_im); VecAssemblyEnd(v_im);

    // Test the norm of v, should be zero for half of them (unless multiplicity):
    double norm_v_re, norm_v_im;
    VecNorm(v_re, NORM_2, &norm_v_re);
    VecNorm(v_im, NORM_2, &norm_v_im);
    double norm_v = sqrt(norm_v_re*norm_v_re + norm_v_im*norm_v_im);
    // printf("Norm of eigenvector %d: %1.14e \n", j, norm_v);

    // Discard if this is a zero vector 
    if (!(norm_v > 1e-10)) {
      // printf("Discarding eigenvalue %1.8e + i%1.8e with zero eigenvector.\n", kr, ki);
      continue;
    }

    // Else, normalize v and store the eigenpair
    VecScale(v_re, 1.0/norm_v);
    VecScale(v_im, 1.0/norm_v);
    // printf("Setting eigenvalue %d: %1.8e + i%1.8e, vnorm=%1.8e\n", found_eigenpairs, kr, ki, norm_v);
    for (PetscInt i=0; i<dim; i++){
      double val_re, val_im;
      VecGetValues(v_re, 1, &i, &val_re);
      VecGetValues(v_im, 1, &i, &val_im);
      MatSetValue(Evecs_re, i, found_eigenpairs, val_re, INSERT_VALUES);
      MatSetValue(Evecs_im, i, found_eigenpairs, val_im, INSERT_VALUES);
      // printf("Setting eigenvector %d at row %d col %d\n", j, i, found_eigenpairs);
    }

    /* Store the eigenvalue */
    eigvals_re->at(found_eigenpairs) = kr;
    eigvals_im->at(found_eigenpairs) = ki;
    found_eigenpairs++;
  }
  
  eigvals_re->resize(found_eigenpairs);
  eigvals_im->resize(found_eigenpairs);

  MatAssemblyBegin(Evecs_re, MAT_FINAL_ASSEMBLY); 
  MatAssemblyEnd(Evecs_re, MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(Evecs_im, MAT_FINAL_ASSEMBLY); 
  MatAssemblyEnd(Evecs_im, MAT_FINAL_ASSEMBLY);

  if (print) {
    // Print each eigenvalue
    printf("Eigen decomposition of U^dagger V done.\n");
    for (size_t i=0; i<eigvals_re->size(); i++){
      double phase = atan2(eigvals_im->at(i), eigvals_re->at(i));
      double abs = sqrt( eigvals_re->at(i)*eigvals_re->at(i) + eigvals_im->at(i)*eigvals_im->at(i) );
      printf("Eigenvalue %d: %f + i*%f, abs=%f, phase = %f \n", i, eigvals_re->at(i), eigvals_im->at(i), abs,  phase*180.0/M_PI);
    }
    MatView(Evecs_re, NULL);
    MatView(Evecs_im, NULL);

  }

  /* Clean up*/
  EPSDestroy(&eigensolver);
  MatDestroy(&A);
  VecDestroy(&eigvecLarge_re);
  VecDestroy(&eigvecLarge_im);
  VecDestroy(&v_re);
  VecDestroy(&v_im);
}


int testEigenComplex(const Mat A_re, const Mat A_im, const std::unique_ptr<std::vector<double>>& eigvals_re, const std::unique_ptr<std::vector<double>>& eigvals_im, const Mat& Evecs_re, const Mat& Evecs_im) {
  PetscInt dim;
  MatGetSize(A_re, &dim, NULL);
  Vec v_re, v_im;
  MatCreateVecs(Evecs_re, &v_re, NULL);
  MatCreateVecs(Evecs_im, &v_im, NULL);
    Vec Av_re, Av_im;
    MatCreateVecs(A_re, &Av_re, NULL);
    MatCreateVecs(A_re, &Av_im, NULL);
 
  // Test total number of eigenvalues = dim
  if (eigvals_re->size()!= dim){
    printf("Error in RiemannianDistance: Wrong number of eigenvalues received: %d instead of %d\n", (int)eigvals_re->size(), dim);
    exit(1);
  }
  for (size_t i=0; i<eigvals_re->size(); i++){
    // Get eigenvector i
    MatGetColumnVector(Evecs_re, v_re, i);
    MatGetColumnVector(Evecs_im, v_im, i);

    // Test norm one of each eigenvector
    double norm_re, norm_im;
    VecNorm(v_re, NORM_2, &norm_re);
    VecNorm(v_im, NORM_2, &norm_im);
    double norm = sqrt(norm_re*norm_re + norm_im*norm_im);
    if (fabs(norm - 1.0) > 1e-12){
      printf("Error in RiemannianDistance: Eigenvector %d has norm %f != 1!\n", (int)i, norm);
      exit(1);
    }

    // Test that A v = lambda v
    // Av_re = A_re v_re - A_im v_im
    MatMult(A_re, v_re, Av_re);
    MatMult(A_im, v_im, Av_im);
    VecAXPY(Av_re, -1.0, Av_im); // Av_re = A_re v_re - A_im v_im
    // Av_im = A_re v_im + A_im v_re
    MatMult(A_re, v_im, Av_im);
    MatMultAdd(A_im, v_re, Av_im, Av_im); // Av_im = A_re v_im + A_im v_re
    // Test difference elementwise
    for (size_t j=0; j<dim; j++){
      double val_re, val_im, valAv_re, valAv_im;
      VecGetValues(v_re, 1, (PetscInt*)&j, &val_re);
      VecGetValues(v_im, 1, (PetscInt*)&j, &val_im);
      VecGetValues(Av_re, 1, (PetscInt*)&j, &valAv_re);
      VecGetValues(Av_im, 1, (PetscInt*)&j, &valAv_im);
      double diff_re = valAv_re - (eigvals_re->at(i)*val_re - eigvals_im->at(i)*val_im);
      double diff_im = valAv_im - (eigvals_re->at(i)*val_im + eigvals_im->at(i)*val_re);
      if (fabs(diff_re) > 1e-10 || fabs(diff_im) > 1e-10){
        printf("Error in RiemannianDistance: Vector %d is not an eigenvector! Diff = %f + i*%f\n", (int)i, diff_re, diff_im);
        exit(1);
      }
    }

    // Test orthogonality of v to all other eigenvectors
    for (size_t j=0; j<eigvals_re->size(); j++){
      if (j==i) continue;
      Vec v2_re, v2_im;
      MatCreateVecs(Evecs_re, &v2_re, NULL);
      MatCreateVecs(Evecs_im, &v2_im, NULL);
      MatGetColumnVector(Evecs_re, v2_re, j);
      MatGetColumnVector(Evecs_im, v2_im, j);
      double dot_re1, dot_im1, dot_re2, dot_im2;
      VecDot(v_re, v2_re, &dot_re1);
      VecDot(v_im, v2_im, &dot_re2);
      VecDot(v_re, v2_im, &dot_im1);
      VecDot(v_im, v2_re, &dot_im2);
      double dot_re = dot_re1 + dot_re2;
      double dot_im = dot_im1 - dot_im2;
      if (fabs(dot_re) > 1e-5 || fabs(dot_im) > 1e-5){
        printf("WARNING in RiemannianDistance: Eigenvector %d is not orthogonal to eigenvector %d! Dot = %1.14e + i*%1.14e\n", (int)i, (int)j, dot_re, dot_im);
        // exit(1);
      }
      VecDestroy(&v2_re);
      VecDestroy(&v2_im);
    }

  }

  // Cleanup
  VecDestroy(&v_re);
  VecDestroy(&v_im);
  VecDestroy(&Av_re);
  VecDestroy(&Av_im);

  return 1;
}


void reconstructMatrixFromEigenComplex(const std::unique_ptr<std::vector<double>>& eigvals_re, const std::unique_ptr<std::vector<double>>& eigvals_im, const Mat& Evecs_re, const Mat& Evecs_im, Mat& Aout_re, Mat& Aout_im, const bool do_log, const Mat& Atest_re, const Mat& Atest_im){
// Compute A_out = Evecs * diag(evals) * Evecs^dagger, or the log thereof 
// A_out_re = V_re * (D_re V_re^T + D_im V_im^T) - V_im * (D_im V_re^T - D_re V_im^T)
// A_out_im = V_im * (D_re V_re^T + D_im V_im^T) + V_re * (D_im V_re^T - D_re V_im^T)

  PetscInt dim;
  MatGetSize(Evecs_re, &dim, NULL);

  MatDuplicate(Evecs_re, MAT_DO_NOT_COPY_VALUES, &Aout_re);
  MatDuplicate(Evecs_im, MAT_DO_NOT_COPY_VALUES, &Aout_im);
  MatZeroEntries(Aout_re);
  MatZeroEntries(Aout_im);

  // First compute D*V^dagger, for V=Evecs, and D=diag(evals) or D=log(diag(evals))
  // LVT_1 = D_re V_re^T + D_im V_im^T
  // LVT_2 = D_im V_re^T - D_re V_im^T
  Mat VreT, VimT;
  MatTranspose(Evecs_re, MAT_INITIAL_MATRIX, &VreT);
  MatTranspose(Evecs_im, MAT_INITIAL_MATRIX, &VimT);
  Mat LVT_1, LVT_2;
  MatDuplicate(VreT, MAT_DO_NOT_COPY_VALUES, &LVT_1);
  MatDuplicate(VreT, MAT_DO_NOT_COPY_VALUES, &LVT_2);
  for (size_t row=0; row<dim; row++){
    // scale each row by the diagonal matrix D
    double diag_re = eigvals_re->at(row); // for reconstruction of A
    double diag_im = eigvals_im->at(row);
    if (do_log) { // for reconstruction of log(A), compute log of d = diag_re + i diag_im. Know that d is on unit circle, so log(d) = i*phi with phi = atan2(diag_im, diag_re)
      double phi = atan2(diag_im, diag_re);
      diag_re = 0.0;
      diag_im = phi;
    }
    for (size_t col=0; col<dim; col++){
      double val_re, val_im;
      MatGetValue(VreT, row, col, &val_re);
      MatGetValue(VimT, row, col, &val_im);

      MatSetValue(LVT_1, row, col, diag_re * val_re + diag_im * val_im, INSERT_VALUES);
      MatSetValue(LVT_2, row, col, diag_im * val_re - diag_re * val_im, INSERT_VALUES);
    }
  }
  MatAssemblyBegin(LVT_1, MAT_FINAL_ASSEMBLY); 
  MatAssemblyEnd(LVT_1, MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(LVT_2, MAT_FINAL_ASSEMBLY); 
  MatAssemblyEnd(LVT_2, MAT_FINAL_ASSEMBLY);

  // printf("log(diag) @ Evecs: \n");
  // MatView(LVT_1, NULL);
  // MatView(LVT_2, NULL);

  // Now compute Aout_re and Aout_im
  // A_test_re = V_re * LVT_1 - V_im * LVT_2
  // A_test_im = V_im * LVT_1 + V_re * LVT_2
  Mat temp_re, temp_im;
  MatMatMult(Evecs_re, LVT_1, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp_re);
  MatMatMult(Evecs_im, LVT_2, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp_im);
  MatAXPY(Aout_re, 1.0, temp_re, DIFFERENT_NONZERO_PATTERN);
  MatAXPY(Aout_re, -1.0, temp_im, DIFFERENT_NONZERO_PATTERN);
  MatDestroy(&temp_re);
  MatDestroy(&temp_im);
  MatMatMult(Evecs_im, LVT_1, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp_re);
  MatMatMult(Evecs_re, LVT_2, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp_im);
  MatAXPY(Aout_im, 1.0, temp_re, DIFFERENT_NONZERO_PATTERN);
  MatAXPY(Aout_im, 1.0, temp_im, DIFFERENT_NONZERO_PATTERN);
  MatDestroy(&temp_re);
  MatDestroy(&temp_im);
  
  // Test the reconstruction if requested
  if (Atest_re != NULL && Atest_im != NULL){
    Mat A_diff_re, A_diff_im;
    MatDuplicate(Atest_re, MAT_COPY_VALUES, &A_diff_re);
    MatDuplicate(Atest_im, MAT_COPY_VALUES, &A_diff_im);
    MatAXPY(A_diff_re, -1.0, Aout_re, SAME_NONZERO_PATTERN);
    MatAXPY(A_diff_im, -1.0, Aout_im, SAME_NONZERO_PATTERN);
    double norm_re, norm_im;
    MatNorm(A_diff_re, NORM_FROBENIUS, &norm_re);
    MatNorm(A_diff_im, NORM_FROBENIUS, &norm_im);
    // if (do_log)
    //   printf("Reconstruction test of log(A): ||log(A) - log(A)_test||_F = %1.14e + i*%1.14e\n", norm_re, norm_im);
    // else
    //   printf("Reconstruction test: ||A - A_test||_F = %1.14e + i*%1.14e\n", norm_re, norm_im);
    if (norm_re > 1e-6 || norm_im > 1e-6){
      printf("Error: Reconstruction from eigen decomposition failed!\n");
      exit(1);
    }
    MatDestroy(&A_diff_re);
    MatDestroy(&A_diff_im);
  }

  // Cleanup
  MatDestroy(&VreT);
  MatDestroy(&VimT);
  MatDestroy(&LVT_1);
  MatDestroy(&LVT_2);
}

// computeMatrixLogTaylor(UdagV_re, UdagV_im, logUdagV_re, logUdagV_im);
double computeMatrixLogTaylor(const Mat& A_re, const Mat& A_im, Mat& logA_re, Mat& logA_im){
  PetscInt dim;
  MatGetSize(A_re, &dim, NULL);

  // Initialize logA_re and logA_im to zero
  MatDuplicate(A_re, MAT_DO_NOT_COPY_VALUES, &logA_re);
  MatDuplicate(A_im, MAT_DO_NOT_COPY_VALUES, &logA_im);
  MatZeroEntries(logA_re);
  MatZeroEntries(logA_im);

  // Compute log(A) = - sum_{n=1}^\infty (I - A)^n / n
  Mat I_minus_A_re, I_minus_A_im;
  MatDuplicate(A_re, MAT_COPY_VALUES, &I_minus_A_re);
  MatDuplicate(A_im, MAT_COPY_VALUES, &I_minus_A_im);
  MatScale(I_minus_A_re, -1.0); // -A_re
  MatScale(I_minus_A_im, -1.0); // -A_im
  MatShift(I_minus_A_re, 1.0); // I - A_re

  // Test norm of I_minus_A
  double Anorm_re, Anorm_im;
  MatNorm(I_minus_A_re, NORM_FROBENIUS, &Anorm_re);
  MatNorm(I_minus_A_im, NORM_FROBENIUS, &Anorm_im);
  if (Anorm_re >= 1.0 || Anorm_im >= 1.0){
    printf("ERROR for Log with Taylor: Norm of I - A is >= 1 (%1.14e + i*%1.14e). Taylor expansion may not converge!\n", Anorm_re, Anorm_im);
    exit(1);
  }

  // Now do the Taylor expansion
  Mat term_re, term_im;
  // First term: I-A
  MatDuplicate(I_minus_A_re, MAT_COPY_VALUES, &term_re); // = I-A_re
  MatDuplicate(I_minus_A_im, MAT_COPY_VALUES, &term_im); // =  -A_im
  MatAXPY(logA_re, -1.0, term_re, DIFFERENT_NONZERO_PATTERN);
  MatAXPY(logA_im, -1.0, term_im, DIFFERENT_NONZERO_PATTERN);
  // Higher order terms
  int max_terms = 1000;
  double tol = 1e-12;
  double err_est = 0.0;
  for (int n=2; n<=max_terms; n++){
    // Update log_A += -term_new / n where term_new = term * (I - A)
    //  real part of term * (I-A) is  term_re*I_minus_A_re - term_im*I_minus_A_im
    // imaginary part of term * (I-A) is  term_im*I_minus_A_re + term_re*I_minus_A_im  
    Mat term_new_re, term_new_im, tmp;
    MatMatMult(term_re, I_minus_A_re, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &term_new_re);
    MatMatMult(term_im, I_minus_A_im, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &tmp);
    MatAXPY(term_new_re, -1.0, tmp, DIFFERENT_NONZERO_PATTERN); // term_new_re = term_re*I_minus_A_re - term_im*I_minus_A_im
    MatDestroy(&tmp);
    // Update logA_re 
    MatAXPY(logA_re, -1.0/n, term_new_re, DIFFERENT_NONZERO_PATTERN);

    MatMatMult(term_im, I_minus_A_re, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &term_new_im);
    MatMatMult(term_re, I_minus_A_im, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &tmp);
    MatAXPY(term_new_im, 1.0, tmp, DIFFERENT_NONZERO_PATTERN); // term_new_im = term_im*I_minus_A_re + term_re*I_minus_A_im
    MatDestroy(&tmp);
    // Update logA_im
    MatAXPY(logA_im, -1.0/n, term_new_im, DIFFERENT_NONZERO_PATTERN);

    // Update term 
    MatCopy(term_new_re, term_re, DIFFERENT_NONZERO_PATTERN);
    MatCopy(term_new_im, term_im, DIFFERENT_NONZERO_PATTERN);
    MatDestroy(&term_new_re);
    MatDestroy(&term_new_im);

    // Check convergence by norm of term / n
    double norm_re, norm_im;
    MatNorm(term_re, NORM_FROBENIUS, &norm_re);
    MatNorm(term_im, NORM_FROBENIUS, &norm_im);
    err_est = sqrt( (norm_re/n)*(norm_re/n) + (norm_im/n)*(norm_im/n) );
    if (norm_re/n < tol && norm_im/n < tol){
      printf("Matrix log converged after %d terms.\n", n);
      break;  
    }
  }

  // Ceanup
  MatDestroy(&I_minus_A_re);
  MatDestroy(&I_minus_A_im);
  MatDestroy(&term_re);
  MatDestroy(&term_im);


  return err_est;
}