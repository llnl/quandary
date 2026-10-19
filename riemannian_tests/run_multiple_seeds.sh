#!/bin/bash

# Script to run Quandary with N different random seeds
# Usage: ./run_multiple_seeds.sh <config.toml> <N> [MPI_RANKS] [TIMEOUT]
#
# Arguments:
#   config.toml - Path to the TOML configuration file
#   N           - Number of random seeds to run
#   MPI_RANKS   - Number of MPI ranks (optional, default=1 for serial)
#   TIMEOUT     - Timeout in seconds per seed (optional, default=600)
#
# Output:
#   Creates directories {config_basename}_seed_0001, {config_basename}_seed_0002, etc.
#   Each directory contains:
#     - Modified config file with unique random seed
#     - All Quandary output files (optim_history.dat, etc.)
#     - quandary_output.log (stdout/stderr capture)
#     - run_status.txt (SUCCESS, FAILED, or TIMEOUT)
#
# Timeout:
#   Each seed run has a configurable timeout to prevent hanging.
#   Default is 600 seconds (10 minutes). Pass 0 to disable timeout.
#   Exit code 124 indicates timeout was reached.

set -e  # Exit on error

# Check arguments
if [ "$#" -lt 2 ] || [ "$#" -gt 4 ]; then
    echo "Error: Wrong number of arguments"
    echo "Usage: $0 <config.toml> <N> [MPI_RANKS] [TIMEOUT]"
    echo ""
    echo "Arguments:"
    echo "  config.toml - Path to TOML configuration file"
    echo "  N           - Number of random seeds to run"
    echo "  MPI_RANKS   - Number of MPI ranks (optional, default=1)"
    echo "  TIMEOUT     - Timeout in seconds per seed (optional, default=600, 0=disabled)"
    echo ""
    echo "Examples:"
    echo "  $0 config_QFT2_dual.toml 10              # Serial, 10 min timeout"
    echo "  $0 config_QFT2_dual.toml 10 4            # 4 ranks, 10 min timeout"
    echo "  $0 config_QFT2_dual.toml 10 4 1200       # 4 ranks, 20 min timeout"
    echo "  $0 config_QFT2_dual.toml 10 1 0          # Serial, no timeout"
    exit 1
fi

CONFIG_FILE="$1"
N_SEEDS="$2"
MPI_RANKS="${3:-1}"       # Default to 1 rank (serial) if not specified
TIMEOUT_SECS="${4:-600}"  # Default to 600 seconds (10 minutes) if not specified

# Validate inputs
if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: Config file '$CONFIG_FILE' not found"
    exit 1
fi

if ! [[ "$N_SEEDS" =~ ^[0-9]+$ ]] || [ "$N_SEEDS" -le 0 ]; then
    echo "Error: N must be a positive integer"
    exit 1
fi

if ! [[ "$MPI_RANKS" =~ ^[0-9]+$ ]] || [ "$MPI_RANKS" -le 0 ]; then
    echo "Error: MPI_RANKS must be a positive integer"
    exit 1
fi

if ! [[ "$TIMEOUT_SECS" =~ ^[0-9]+$ ]]; then
    echo "Error: TIMEOUT must be a non-negative integer (0 to disable)"
    exit 1
fi

# Find Quandary executable
QUANDARY_EXEC=""
if [ -f "../../build/quandary" ]; then
    QUANDARY_EXEC="../../build/quandary"
elif [ -f "../build/quandary" ]; then
    QUANDARY_EXEC="../build/quandary"
elif command -v quandary &> /dev/null; then
    QUANDARY_EXEC="quandary"
else
    echo "Error: Cannot find Quandary executable"
    echo "Searched: ../../build/quandary, ../build/quandary, and \$PATH"
    exit 1
fi

# Determine MPI command if using multiple ranks
MPI_CMD=""
if [ "$MPI_RANKS" -gt 1 ]; then
    if command -v mpirun &> /dev/null; then
        MPI_CMD="mpirun -np $MPI_RANKS"
    elif command -v mpiexec &> /dev/null; then
        MPI_CMD="mpiexec -n $MPI_RANKS"
    elif command -v srun &> /dev/null; then
        MPI_CMD="srun -n $MPI_RANKS"
    else
        echo "Error: MPI_RANKS=$MPI_RANKS but no MPI launcher found"
        echo "Searched for: mpirun, mpiexec, srun"
        exit 1
    fi
fi

# Check for timeout command availability
TIMEOUT_CMD=""
USE_MANUAL_WATCHDOG=false
if [ "$TIMEOUT_SECS" -gt 0 ]; then
    if command -v timeout &> /dev/null; then
        # Use SIGKILL after 10 seconds if SIGTERM doesn't work
        TIMEOUT_CMD="timeout --kill-after=10 $TIMEOUT_SECS"
    elif command -v gtimeout &> /dev/null; then
        # macOS with GNU coreutils installed
        TIMEOUT_CMD="gtimeout --kill-after=10 $TIMEOUT_SECS"
    else
        # Use manual watchdog implementation
        USE_MANUAL_WATCHDOG=true
    fi
fi

echo "Using Quandary executable: $QUANDARY_EXEC"
if [ -n "$MPI_CMD" ]; then
    echo "MPI execution: $MPI_CMD"
else
    echo "MPI execution: Serial (1 rank)"
fi
echo "Config file: $CONFIG_FILE"
echo "Number of seeds: $N_SEEDS"
if [ "$TIMEOUT_SECS" -gt 0 ]; then
    echo "Timeout per seed: ${TIMEOUT_SECS}s"
else
    echo "Timeout per seed: disabled"
fi
echo ""

# Determine number of digits for directory names
N_DIGITS=${#N_SEEDS}
if [ "$N_DIGITS" -lt 4 ]; then
    N_DIGITS=4
fi

# Get base name of config file (without path and without .toml extension)
CONFIG_BASENAME=$(basename "$CONFIG_FILE")
CONFIG_PREFIX="${CONFIG_BASENAME%.toml}"

# Run for each seed
for i in $(seq 1 $N_SEEDS); do
    # Format seed number with leading zeros
    SEED_NUM=$(printf "%0${N_DIGITS}d" $i)
    DIR_NAME="${CONFIG_PREFIX}_seed_${SEED_NUM}"

    echo "====================================="
    echo "Running seed $i of $N_SEEDS (directory: $DIR_NAME)"
    echo "====================================="

    # Create directory
    mkdir -p "$DIR_NAME"

    # Copy config file and modify random seed
    # If the config has a rand_seed parameter, modify it; otherwise append it
    if grep -q "^rand_seed" "$CONFIG_FILE"; then
        # Replace existing rand_seed
        sed "s/^rand_seed.*$/rand_seed = $i/" "$CONFIG_FILE" > "$DIR_NAME/$CONFIG_BASENAME"
    else
        # Append rand_seed to [master] section if it exists, otherwise to end of file
        if grep -q "^\[master\]" "$CONFIG_FILE"; then
            awk -v seed="$i" '/^\[master\]/ {print; print "rand_seed = " seed; next} 1' "$CONFIG_FILE" > "$DIR_NAME/$CONFIG_BASENAME"
        else
            cp "$CONFIG_FILE" "$DIR_NAME/$CONFIG_BASENAME"
            echo "" >> "$DIR_NAME/$CONFIG_BASENAME"
            echo "# Random seed added by run_multiple_seeds.sh" >> "$DIR_NAME/$CONFIG_BASENAME"
            echo "rand_seed = $i" >> "$DIR_NAME/$CONFIG_BASENAME"
        fi
    fi

    # Run Quandary in the seed directory
    cd "$DIR_NAME"

    # Build full command
    if [ -n "$MPI_CMD" ]; then
        FULL_CMD="$MPI_CMD $QUANDARY_EXEC $CONFIG_BASENAME"
    else
        FULL_CMD="$QUANDARY_EXEC $CONFIG_BASENAME"
    fi

    # Run with output to both console and log file
    echo "Starting Quandary (seed=$i)..."
    echo "Command: $FULL_CMD"
    START_TIME=$(date +%s)

    # Close stdin and redirect to /dev/null to prevent MPI from waiting for input
    # Use timeout to prevent infinite hangs (configurable via TIMEOUT_SECS)
    RUN_STATUS="UNKNOWN"
    if [ -n "$TIMEOUT_CMD" ]; then
        # Run with timeout (sends SIGTERM, then SIGKILL after 10s)
        if $TIMEOUT_CMD $FULL_CMD </dev/null 2>&1 | tee quandary_output.log; then
            END_TIME=$(date +%s)
            ELAPSED=$((END_TIME - START_TIME))
            echo "Seed $i completed successfully in ${ELAPSED}s"
            echo "SUCCESS" > run_status.txt
            RUN_STATUS="SUCCESS"
        else
            EXIT_CODE=$?
            END_TIME=$(date +%s)
            ELAPSED=$((END_TIME - START_TIME))
            if [ $EXIT_CODE -eq 124 ] || [ $EXIT_CODE -eq 137 ]; then
                # 124 = timeout, 137 = SIGKILL (128+9)
                echo "ERROR: Seed $i timed out after ${ELAPSED}s (exceeded ${TIMEOUT_SECS}s limit)"
                echo "TIMEOUT" > run_status.txt
                RUN_STATUS="TIMEOUT"
                # Clean up any lingering MPI processes
                pkill -9 -f "quandary.*config.*seed_$(printf "%0${N_DIGITS}d" $i)" 2>/dev/null || true
            else
                echo "ERROR: Seed $i failed after ${ELAPSED}s (exit code: $EXIT_CODE)"
                echo "FAILED" > run_status.txt
                RUN_STATUS="FAILED"
            fi
        fi
    elif [ "$USE_MANUAL_WATCHDOG" = true ]; then
        # Run with manual background watchdog (for macOS without timeout command)
        $FULL_CMD </dev/null 2>&1 | tee quandary_output.log &
        CMD_PID=$!

        # Wait for process with manual timeout
        WAIT_COUNT=0
        TIMED_OUT=false
        while kill -0 $CMD_PID 2>/dev/null; do
            sleep 1
            WAIT_COUNT=$((WAIT_COUNT + 1))
            if [ $WAIT_COUNT -ge $TIMEOUT_SECS ]; then
                echo "Manual timeout reached after ${TIMEOUT_SECS}s, killing process..."
                kill -9 $CMD_PID 2>/dev/null || true
                pkill -9 -P $CMD_PID 2>/dev/null || true
                # Also kill by name pattern
                pkill -9 -f "$QUANDARY_EXEC" 2>/dev/null || true
                sleep 1
                TIMED_OUT=true
                break
            fi
        done

        END_TIME=$(date +%s)
        ELAPSED=$((END_TIME - START_TIME))

        if [ "$TIMED_OUT" = true ]; then
            echo "ERROR: Seed $i timed out after ${ELAPSED}s (exceeded ${TIMEOUT_SECS}s limit)"
            echo "TIMEOUT" > run_status.txt
            RUN_STATUS="TIMEOUT"
        else
            # Get exit code
            wait $CMD_PID 2>/dev/null
            EXIT_CODE=$?
            if [ $EXIT_CODE -eq 0 ]; then
                echo "Seed $i completed successfully in ${ELAPSED}s"
                echo "SUCCESS" > run_status.txt
                RUN_STATUS="SUCCESS"
            else
                echo "ERROR: Seed $i failed after ${ELAPSED}s (exit code: $EXIT_CODE)"
                echo "FAILED" > run_status.txt
                RUN_STATUS="FAILED"
            fi
        fi
    else
        # Run without any timeout protection
        if $FULL_CMD </dev/null 2>&1 | tee quandary_output.log; then
            END_TIME=$(date +%s)
            ELAPSED=$((END_TIME - START_TIME))
            echo "Seed $i completed successfully in ${ELAPSED}s"
            echo "SUCCESS" > run_status.txt
            RUN_STATUS="SUCCESS"
        else
            EXIT_CODE=$?
            END_TIME=$(date +%s)
            ELAPSED=$((END_TIME - START_TIME))
            echo "ERROR: Seed $i failed after ${ELAPSED}s (exit code: $EXIT_CODE)"
            echo "FAILED" > run_status.txt
            RUN_STATUS="FAILED"
        fi
    fi

    cd ..

    # Rename directory based on status
    if [ "$RUN_STATUS" = "FAILED" ]; then
        NEW_DIR_NAME="${DIR_NAME}_failed"
        if [ -d "$NEW_DIR_NAME" ]; then
            rm -rf "$NEW_DIR_NAME"
        fi
        mv "$DIR_NAME" "$NEW_DIR_NAME"
        echo "Renamed directory to: $NEW_DIR_NAME"
    elif [ "$RUN_STATUS" = "TIMEOUT" ]; then
        NEW_DIR_NAME="${DIR_NAME}_timeout"
        if [ -d "$NEW_DIR_NAME" ]; then
            rm -rf "$NEW_DIR_NAME"
        fi
        mv "$DIR_NAME" "$NEW_DIR_NAME"
        echo "Renamed directory to: $NEW_DIR_NAME"
    fi

    echo ""
done

echo "====================================="
echo "All runs completed!"
echo "====================================="

# Summary
echo ""
echo "Summary:"
SUCCESS_COUNT=$(find ${CONFIG_PREFIX}_seed_* -maxdepth 1 -name "run_status.txt" -exec grep -l "SUCCESS" {} \; 2>/dev/null | wc -l)
FAILED_COUNT=$(find ${CONFIG_PREFIX}_seed_*_failed -maxdepth 1 -name "run_status.txt" 2>/dev/null | wc -l)
TIMEOUT_COUNT=$(find ${CONFIG_PREFIX}_seed_*_timeout -maxdepth 1 -name "run_status.txt" 2>/dev/null | wc -l)
echo "  Successful: $SUCCESS_COUNT / $N_SEEDS"
echo "  Failed:     $FAILED_COUNT / $N_SEEDS"
echo "  Timeout:    $TIMEOUT_COUNT / $N_SEEDS"
echo ""
echo "Results saved in directories:"
echo "  Success: ${CONFIG_PREFIX}_seed_XXXX"
echo "  Failed:  ${CONFIG_PREFIX}_seed_XXXX_failed"
echo "  Timeout: ${CONFIG_PREFIX}_seed_XXXX_timeout"
