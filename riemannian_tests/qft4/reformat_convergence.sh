#!/bin/bash
# Reformat convergence_raw.dat to match ksp_convergence format

input_file="convergence_raw.dat"
output_file="ksp_convergence_formatted.dat"

# Write header (omit unknown metadata - plotting script will handle gracefully)
cat > "$output_file" << 'HEADER'
# KSP convergence history
# KSP_TYPE: formatted
# iteration  residual_norm  solution_norm
HEADER

# Parse and reformat data
sed 's/KSP it //' "$input_file" | \
  sed 's/: residual norm = / /' | \
  sed 's/, solution norm = / /' | \
  awk '{printf "%d %s %s\n", $1, $2, $3}' >> "$output_file"

echo "Reformatted data written to: $output_file"
wc -l "$output_file"
head -10 "$output_file"
