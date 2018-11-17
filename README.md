# Quick guide (with HDF5)

This is the brief version of the usual workflow of fastGenoGAM. It involves:

- Reading in data through `GenoGAMDataSet()` to get a `GenoGAMDataSet` object. This is done with the HDF5 backend.
  
- Computing size factors with `computeSizeFactors()`
  
- Compute the model with `genogam()` to get the final `GenoGAM` object

- compute position-wise log p-values with `computeSignificance()`

## Getting started
The example dataset that is part of the package is restricted to the first 100kb of the first yeast chromosome. 

1. We set some meta variables

```r
library(fastGenoGAM)

## A.
## specify folder and experiment design path
wd <- system.file("extdata/Set1", package='fastGenoGAM')
folder <- file.path(wd, "bam")
expDesign <- file.path(wd, "experimentDesign.txt")

## B.
## set HDF5 folder where the data will be stored
## Note, don't usually use /tmp because otherwise all your data and results get deleted later
hdf5_folder <- tempdir()
settings <- GenoGAMSettings(hdf5Control = list(dir = hdf5_folder))

## C. 
## register parallel backend (default is "the number of cores" - 2)
## Here we also use the SnowParam backend which is advised for larger data
## For yeast MulticoreParam should also do fine
BiocParallel::register(BiocParallel::SnowParam(worker = 2))

## D.
## specify chunk and overhang size. It is possible to skip this,
## but the automatic specification would prevent us from using
## multiple tiles in such a small example.
chunkSize <- 50000
overhangSize <- 1000

## E.
## the experiment design reflecting the underlying GAM
## x = position
design <- ~ s(x) + s(x, by = genotype)
```

2. Read in the data to obtain a `GenoGAMDataSet` object.
Warnings about out-of-bound ranges can be ignored, as they occur during the 
shifting process when the genome gets tiled. It is trimmed accordingly. I have not yet figured out how to silent those warnings.
The usual `suppressWarnings` does not seem to work in the parallel Snow environment.

```r
## build the GenoGAMDataSet with HDF5 backend
ggd <- GenoGAMDataSet(
  expDesign, directory = folder,
  chunkSize = chunkSize, overhangSize = overhangSize,
  design = design,
  settings = settings, hdf5 = TRUE
)
```

3. Compute Size factors

```r
## compute size factors
ggd <- computeSizeFactors(ggd)
```

## Fitting the model

1. Compute model with fixed hyperparameters

```r
## compute model without parameter estimation to save time in vignette
result <- genogam(ggd, lambda = 4601, theta = 4.51)

## Check the fit and standard error
fits(result)
se(result)
```

2. Compute log p-values
```r
computeSignificance(result)

## check p-values
pval(result)
````

Plot results of the center 10kb, where both tiles are joined together. A proper plotting functions is in development.
```r
idx <- 45000:55000
upper_input <- fits(result)$chrI[,"s(x)"][idx] + 2*se(result)$chrI[,"s(x)"][idx]
lower_input <- fits(result)$chrI[,"s(x)"][idx] - 2*se(result)$chrI[,"s(x)"][idx]

upper_ip <- fits(result)$chrI[,"s(x):genotype"][idx] + 2*se(result)$chrI[,"s(x):genotype"][idx]
lower_ip <- fits(result)$chrI[,"s(x):genotype"][idx] - 2*se(result)$chrI[,"s(x):genotype"][idx]

par(mfrow = c(2,1))
plot(fits(result)$chrI[,"s(x)"][idx], type = "l")
lines(upper_input, lty = 2)
lines(lower_input, lty = 2)

plot(fits(result)$chrI[,"s(x):genotype"][idx], type = "l")
lines(upper_ip, lty = 2)
lines(lower_ip, lty = 2)
```

# FAQ

**Problem:** An error occurs during model fitting along those lines:

> error: matrix multiplication: incompatible matrix dimensions: 22333830147200x5360120267008000 and 4294972496x1

**Solution:** First, make sure you have all Armadillo dependencies installed correctly. See [here](http://arma.sourceforge.net/download.html)

Second, the error is most likely related to the fact, that Armadillo is using 32bit matrices, thus causing problems for large matrices fastGenoGAM is using. The solution is to enable `ARMA_64BIT_WORD`, which is not enabled in RcppArmadillo by default. This should have been done during compilation, but if it fails for some reason you can do it manually with `#define ARMA_64BIT_WORD 1` in `my_R_Directory/lib/R/library/RcppArmadillo/include/RcppArmadilloConfig.h`.
See [here](https://stackoverflow.com/questions/40592054/large-matrices-in-rcpparmadillo-via-the-arma-64bit-word-define)
