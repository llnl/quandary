#!/bin/sh
# Submit with the command "flux batch run-cpu.sh" 
#flux: -N 1
#flux: -q pbatch
#flux: -B rqspam
#flux: -t 1.3h

export HSA_XNACK=1

flux run -N 1 -n 16 -x quandary ./config_QFT4_lbfgs.toml
