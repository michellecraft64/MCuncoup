# MCuncoup
Program files referenced in "The effects of background noise on a biophysical model of olfactory bulb mitral cells." Bulletin of Mathematical Biology (107) Vol. 84 pp:1--20. https://link.springer.com/article/10.1007/s11538-022-01066-8

MCsigFast.m -- runs the MC model, saves spike times in flName.dat (flName is string input), spike times are appended to the flat .dat file, separated by \n; assuming sub-directory ~/1dats/ exists.

drv_varySig.m -- driver file the calls McsigFast.m, specifies applied current (Iapp), flName (in flNameB variable), input noise (sigm vector), and required number of spikes (ReqNums variable).

XPP files to calc bifurcation diagram and PRC in Appendix:
MCunc.ode

MAT files: 
dFI.mat -- the FI curve, firing rate (frate) in Hz and input current (Iap_v) in muA/cm^2, with no noise

dFIsgmPt5.mat -- same as dFI.mat but with input noise sigma=0.5

dFIsgm1.mat -- same other dFI[].mat but with input noise sigma=1

dFIsgm1pt5.mat -- same other dFI[].mat but with input noise sigma=1.5

Sub-directories: 
~/1dats/  Contains .dat files from long Monte Carlo runs of biopysical MC model (spike times)

~/phenomModel/ Contains Matlab scripts to plot the phenomenological ISI density model. These are analytic solutions to the steady-state density & derived entities.

scpt_simpSDE.m -- the model with additive (constant sigma) noise; on mon. incr. CV/std b/c Gaussian

scpt_SDEmulti.m -- the model with S dependent noise, simple ramp to get non-mon. CV/std

scpt_mixGauss.m -- mixture of 2 Gaussians (way more params) that can get non-mon. CV/std
