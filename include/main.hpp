#pragma once

#include "config.hpp"

/**
 * @brief Run a Quandary simulation or optimization.
 *
 * This function handles MPI/PETSc initialization (if not already done),
 * and runs the simulation, optimization, or gradient evaluation based on
 * the config settings.
 *
 * @param config The loaded configuration
 * @param quietmode If true, suppress most console output
 * @param argc Argument count for MPI_Init (0 for defaults)
 * @param argv Argument vector for MPI_Init (nullptr for defaults)
 * @param petsc_argc Number of PETSc command-line arguments (0 for defaults)
 * @param petsc_argv PETSc command-line arguments (nullptr for defaults)
 * @return 0 on success, non-zero on error
 */
int runQuandary(const Config& config, bool quietmode = false,
                int argc = 0, char** argv = nullptr,
                int petsc_argc = 0, char** petsc_argv = nullptr);
