# Make sure you have the location of quandary.py in your PYTHONPATH. E.g. with
#   > export PYTHONPATH=/path/to/quandary/:$PYTHONPATH
# Further, make sure that your quandary executable is in your $PATH variable. E.g. with
#   > export PATH=/path/to/quandary/:$PATH
from quandary import * 
import numpy as np

## Two qubit test case: CNOT gate, two levels each, no guard levels, dipole-dipole coupling 5KHz ##

# 01 transition frequencies [GHz] per oscillator
freq01 = [4.80595, 4.8601] 

# Coupling strength [GHz] (Format [0<->1, 0<->2, ..., 1<->2, ... ])
Jkl = [0.005]  # Dipole-Dipole coupling of qubit 0<->1

# Frequency of rotations for computational frame [GHz] per oscillator
favg = sum(freq01)/len(freq01)
rotfreq = favg*np.ones(len(freq01))

# Set the pulse duration (ns)
T = 200.0

# Set up the CNOT target gate
unitary = np.identity(4)
unitary[2,2] = 0.0
unitary[3,3] = 0.0
unitary[2,3] = 1.0
unitary[3,2] = 1.0
# print("Target gate: ", unitary)

# Flag for printing out more information to screen
verbose = False

# For reproducability: Random number generator seed
rand_seed=1234

# Set up the Quandary configuration for this test case. Make sure to pass all of the above to the corresponding fields, compare help(Quandary)!
quandary = Quandary(freq01=freq01, Jkl=Jkl, rotfreq=rotfreq, T=T, targetgate=unitary, verbose=verbose, rand_seed=rand_seed) 


# Submit jobs on increasing numbers of cores
ncores_all = [1, 2, 4]
for i in range(len(ncores_all)):

  # Number of cores
  ncores = ncores_all[i]

  # Specify job name, this will also be the output directory for this job
  jobname = 'ncores_'+str(ncores)

  # Define slurm parameters
  maxtime = "00:02:00"
  account = "sqoc"
  nodes = int(np.ceil(ncores / 36)) # for quartz: 36 cores per node
  batchargs = [maxtime, account, nodes]

  # Submit quandary jobs to batch system. 
  # Currently, this will generate scripts to submit SLURM batch jobs. The batch system can be changed through the definitions at the end of quandary.py
  t, pt, qt, infidelity, expectedEnergy, population = quandary.optimize(datadir=jobname, maxcores=ncores, batchargs=batchargs)


  # AFTER THE JOBS FINISHED, load results into python if needed. E.g. with 
  # time, pt, qt, uT, expectedEnergy, population, popt, infidelity, optim_hist = quandary.get_results(datadir=jobname)
  # print("Fidelity = ", 1.0 - infidelity)



## Another example: Run quandary on varying parameters 
#mynsplines = [10,20,30]
#for i in range(len(mynsplines)):
#  nsplines = mynsplines[i]
#
#  # Change the quandary configurations.
#  quandary.nsplines = nsplines
#  quandary.update()
#  # Alternatively, you could also just set up a new Quandary instance with the new configuration. In this case it would be 
#  # quandary = Quandary(nsplines=nsplines, freq01=freq01, Jkl=Jkl, rotfreq=rotfreq, T=T, targetgate=unitary, verbose=verbose, rand_seed=rand_seed) 
#
#  # Execute quandary locally in a specified data directory
#  dir="nsplines_"+str(nsplines)
#  t, pt, qt, infidelity, expectedEnergy, population = quandary.optimize(datadir=dir)
