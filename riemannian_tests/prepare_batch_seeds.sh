#!/bin/bash

# Script to prepare batch submission directories for multiple random seeds
# Usage: ./prepare_batch_seeds.sh <config.toml> <N> <MPI_RANKS> <FLUX_TIME> [FLUX_QUEUE]
#
# This script creates directory structure and submission scripts but does NOT
# start any jobs. Submit jobs manually with: cd <dir> && flux batch run.sh
#
# Arguments:
#   config.toml - Path to the TOML configuration file
#   N           - Number of random seeds to prepare
#   MPI_RANKS   - Number of MPI ranks per job
#   FLUX_TIME   - Time limit per job (e.g., 30m, 1.5h, 2h)
#   FLUX_QUEUE  - Flux queue name (optional, default=pbatch)
#
# Output:
#   Creates directories {config_basename}_seed_0001, {config_basename}_seed_0002, etc.
#   Each directory contains:
#     - Modified config file with unique random seed
#     - run.sh script for Flux batch submission

set -e  # Exit on error

# Check arguments
if [ "$#" -lt 4 ] || [ "$#" -gt 5 ]; then
    echo "Error: Wrong number of arguments"
    echo "Usage: $0 <config.toml> <N> <MPI_RANKS> <FLUX_TIME> [FLUX_QUEUE]"
    echo ""
    echo "Arguments:"
    echo "  config.toml - Path to TOML configuration file"
    echo "  N           - Number of random seeds to prepare"
    echo "  MPI_RANKS   - Number of MPI ranks per job"
    echo "  FLUX_TIME   - Time limit per job (e.g., 30m, 1.5h, 2h)"
    echo "  FLUX_QUEUE  - Flux queue name (optional, default=pbatch)"
    echo ""
    echo "Examples:"
    echo "  $0 config_QFT2_dual.toml 10 4 1h"
    echo "  $0 config_QFT2_dual.toml 10 16 30m pdebug"
    echo "  $0 config_QFT3_lbfgs.toml 50 8 2.5h pbatch"
    echo ""
    echo "Note: This script creates directories and scripts but does NOT submit jobs."
    echo "      To submit: cd <directory> && flux batch run.sh"
    exit 1
fi

CONFIG_FILE="$1"
N_SEEDS="$2"
MPI_RANKS="$3"
FLUX_TIME="$4"
FLUX_QUEUE="${5:-pbatch}"     # Default to pbatch queue

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

echo "Using Quandary executable: $QUANDARY_EXEC"
echo "Config file: $CONFIG_FILE"
echo "Number of seeds: $N_SEEDS"
echo "MPI ranks per job: $MPI_RANKS"
echo "Flux queue: $FLUX_QUEUE"
echo "Time limit: $FLUX_TIME"
echo ""

# Determine number of digits for directory names
N_DIGITS=${#N_SEEDS}
if [ "$N_DIGITS" -lt 4 ]; then
    N_DIGITS=4
fi

# Get base name of config file (without path and without .toml extension)
CONFIG_BASENAME=$(basename "$CONFIG_FILE")
CONFIG_PREFIX="${CONFIG_BASENAME%.toml}"

# Calculate number of nodes (assuming you want 1 node per job for now)
# Adjust this logic if you need multiple nodes
N_NODES=1

# Prepare for each seed
for i in $(seq 1 $N_SEEDS); do
    # Format seed number with leading zeros
    SEED_NUM=$(printf "%0${N_DIGITS}d" $i)
    DIR_NAME="${CONFIG_PREFIX}_seed_${SEED_NUM}"

    echo "====================================="
    echo "Preparing seed $i of $N_SEEDS (directory: $DIR_NAME)"
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
            echo "# Random seed added by prepare_batch_seeds.sh" >> "$DIR_NAME/$CONFIG_BASENAME"
            echo "rand_seed = $i" >> "$DIR_NAME/$CONFIG_BASENAME"
        fi
    fi

    # Create Flux batch submission script
    cat > "$DIR_NAME/run.sh" <<EOF
#!/bin/sh
# Flux batch submission script for seed $i
# Submit with: flux batch run.sh
#
#flux: -N $N_NODES
#flux: -q $FLUX_QUEUE
#flux: -B rqspam
#flux: -t $FLUX_TIME

# Optional: Set environment variables if needed
# export HSA_XNACK=1

# Run Quandary with Flux
flux run -N $N_NODES -n $MPI_RANKS $QUANDARY_EXEC $CONFIG_BASENAME
EOF

    chmod +x "$DIR_NAME/run.sh"

    echo "Created: $DIR_NAME/"
    echo "  - $CONFIG_BASENAME (with rand_seed = $i)"
    echo "  - run.sh (Flux submission script)"
    echo ""
done

echo "====================================="
echo "All seed directories prepared!"
echo "====================================="
echo ""
echo "Directories created: ${CONFIG_PREFIX}_seed_0001, ${CONFIG_PREFIX}_seed_0002, ..., ${CONFIG_PREFIX}_seed_$(printf "%0${N_DIGITS}d" $N_SEEDS)"
echo ""
echo "To submit jobs:"
echo "  cd ${CONFIG_PREFIX}_seed_0001 && flux batch run.sh"
echo "  cd ${CONFIG_PREFIX}_seed_0002 && flux batch run.sh"
echo "  ..."
echo ""
echo "Or submit all at once with:"
echo "  for dir in ${CONFIG_PREFIX}_seed_*; do (cd \$dir && flux batch run.sh); done"
