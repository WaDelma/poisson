# poisson    
**poisson** is library for generating n-dimensional poisson-disk distributions.    

[![Build Status](https://travis-ci.org/WaDelma/poisson.svg?branch=master)](https://travis-ci.org/WaDelma/poisson)   

## Using **poisson**    
**poisson** is used through struct PoissonDisk.    
PoissonDisk is generic for Rng supplied and VecLike.    
VecLike describes what traits from std, nalgebra and num crates needs to be implemented for samples so the algorithmn can function.    
When constructing a PoissonDisk boolean for periodicity can be used to make the distribution to wrap around.    
When PoissonDisk is created using with_samples the algorithm calculates approximately the radius that is needed to create distribution with that many samples.    
Construction of PoissonDisk with method with_samples for non-perioditic distributions is currently restricted to 2, 3 and 4 dimensions (due need of contanstants for each dimension).    
Samples can be given to the method for creating distribution to make it take those account when creating the distribution.    

Documentation is [here](https://WaDelma.github.io/poisson/)    

# TODO
   * Can quadtree be compressed? Skip Quadtree?    
   * Benchmarks!    
