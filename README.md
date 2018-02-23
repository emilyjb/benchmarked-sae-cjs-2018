# benchmarked-sae-cjs-2018
Code for Benchmarked Small Area Prediction (Canadian Journal of Statistics, 2018)
## Additive error models -- code for simulation 1 and NRI application

* augfuns.R: Functions for augmented model benchmarking 
* estfuns.R: Functions for estimation of additive error model
* genfuns.R: Functions to generate data for additive error simulation
* fixedpars.R: Establish fixed parameters for simulation 1
* storeoutput.R: Empty matrices to store output. 
* runsim1.R : Runs simulation 1
* outputforsupplementtables.R: Obtain output for tables in supplement
* NRIEstimationFromFakeData.R: Runs code for data analysis. Due to confidentiality of the real data, we use a simulated data set with properties similar to properties of the NRI data. Estimates (and SEs) for the simulated data set are as follows:
> - beta0: -4.02 (0.13)
> - beta1: 1.03 (0.03)
> - sigma2u: 0.0068 (0.0025)
* Generated data sets:
> - Sig2u005datlistgen.Rdata: Data to reproduce simulation 1 with sigma2u = 0.005
> - Sig2u06datlistgen.Rdata: Data to reproduce simulation 1 with sigma2u = 0.06
> - FakeDataNRI.Rdata: Data simulated to have properties similar to NRI data

## Nonadditive error models -- code for simulation 2
