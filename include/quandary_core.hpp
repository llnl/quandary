#pragma once

/**
 * @brief Main entry point for Quandary simulations and optimizations.
 *
 * This function handles MPI/PETSc initialization, config loading, and runs
 * the simulation, optimization, or gradient evaluation based on config settings.
 *
 * @param argc Argument count (same as main)
 * @param argv Argument vector (same as main)
 * @return 0 on success, non-zero on error
 */
int runQuandary(int argc, char** argv);
