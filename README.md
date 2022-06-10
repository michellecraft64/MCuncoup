# MCuncoup
Program files referenced in "The effects of background noise on a biophysical model of olfactory bulb mitral cells" by M. Craft &amp; C. Ly, #### (journal info. later)

MCsigFast.m -- runs the MC model, saves spike times in flName.dat (flName is string input), spike times are appended to the flat .dat file, separated by \n; assuming sub-directory ~/1dats/ exists.

drv_varySig.m -- driver file the calls McsigFast.m, specifies applied current (Iapp), flName (in flNameB variable), input noise (sigm vector), and required number of spikes (ReqNums variable).

XPP files to calc bifurcation diagram &amp; PRC in Appendix:
MCunc.ode

Sub-directories: 
~/1dats/  Contains .dat files from long Monte Carlo runs of biopysical MC model (spike times)

~/phenomModel/ Contains Matlab scripts to plot the phenomenological ISI density model. These are analytic solutions to the steady-state density & derived entities.

scpt_simpSDE.m -- the model with additive (constant sigma) noise; on mon. incr. CV/std b/c Gaussian

scpt_SDEmulti.m -- the model with S dependent noise, simple ramp to get non-mon. CV/std

scpt_mixGauss.m -- mixture of 2 Gaussians (way more params) that can get non-mon. CV/std
