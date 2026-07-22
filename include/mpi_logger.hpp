#pragma once

#include <cstdlib>
#include <iostream>
#include <mpi.h>

/**
 * @brief MPI-aware logger that handles rank filtering and quiet mode.
 *
 * Encapsulates MPI rank and quiet mode to simplify logging calls throughout
 * the codebase. Only rank 0 outputs messages, and quiet mode can suppress output.
 */
class MPILogger {
 private:
  int mpi_rank;
  bool quiet_mode;

  /**
   * @brief Get MPI rank, handling the case where MPI may not be initialized.
   *
   * If MPI is initialized, returns the rank from MPI_COMM_WORLD.
   * If MPI is not initialized (e.g., in unit tests), returns 0.
   */
  static int getMpiRank() {
    int rank = 0;
    int mpi_initialized = 0;
    MPI_Initialized(&mpi_initialized);
    if (mpi_initialized) {
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    }
    return rank;
  }

 public:
  /**
  * @brief Constructs an MPI logger with automatic rank detection.
  *
  * Gets the MPI rank from MPI_COMM_WORLD if MPI is initialized,
  * otherwise uses rank 0 (for unit tests).
  *
  * @param quiet If true, suppresses non-error output
  */
  explicit MPILogger(bool quiet) : mpi_rank(getMpiRank()), quiet_mode(quiet) {}

  /**
  * @brief Constructs an MPI logger with explicit rank.
  *
  * @param rank MPI rank of the current process
  * @param quiet If true, suppresses non-error output
  */
  MPILogger(int rank, bool quiet = false) : mpi_rank(rank), quiet_mode(quiet) {}

  /**
  * @brief Logs a message to stdout (rank 0 only, respects quiet mode).
  *
  * @param message Message to log
  */
  void log(const std::string& message) const {
    if (!quiet_mode && mpi_rank == 0) {
      std::cout << message << std::endl;
    }
  }

  /**
  * @brief Logs an error message to stderr (rank 0 only).
  *
  * @param message Error message to log
  */
  void error(const std::string& message) const {
    if (mpi_rank == 0) {
      std::cerr << "ERROR: " << message << std::endl;
    }
  }

  /**
  * @brief Logs an error and terminates the program.
  *
  * @param message Error message before exit
  */
  [[noreturn]]
  void exitWithError(const std::string& message) const {
    error(message);
    exit(1);
  }
};
