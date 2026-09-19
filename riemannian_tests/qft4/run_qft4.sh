#!/bin/sh
# Submit with the command "flux batch run_qft4.sh"
#flux: -N 4
#flux: -q pdebug
#flux: -B rqspam
#flux: -t 15m

export HSA_XNACK=1
export MPICH_GPU_SUPPORT_ENABLED=1

flux run -N 4 -x --gpus-per-node=4 --tasks-per-node=4 quandary ./config_qft4.toml --petsc-options "-vec_type kokkos -mat_type aijkokkos"
