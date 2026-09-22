#!/bin/sh
# Submit with the command "flux batch run-cpu.sh" 
#flux: -N 1
#flux: -q pbatch
#flux: -B rqspam
#flux: -t 1h

export HSA_XNACK=1

# flux run -N 1 -n 16 -x quandary ./config_qft4.toml
flux run -N 1 -n 16 -x quandary ./config_qft4_Jtrace.toml
