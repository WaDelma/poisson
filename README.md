# poisson    
**poisson** is library for generating n-dimensional poisson-disk distributions.    

[![Build Status](https://travis-ci.org/WaDelma/poisson.svg?branch=master)](https://travis-ci.org/WaDelma/poisson)   

## Using **poisson**    
**poisson** is used through struct PoissonGen which is built using struct PoissonDisk.    

PoissonDisk is generic to Rng supplied and VecLike which describes generation dimension.    
VecLike describes what traits from std, nalgebra and num crates are needed to be implemented for samples so the algorithmn can function.    
When building a PoissonDisk periodicity can be set from builder to make the distribution to wrap around.    

When PoissonGen is created using build_samples the algorithm calculates approximately the radius that is needed to create distribution with that many samples.    
Building of PoissonDisk with method build_samples for non-perioditic distributions is currently restricted to 2, 3 and 4 dimensions (due need of contanstants for each dimension).  

Samples can be given to the method for creating distribution to make it take those account when creating the distribution.    

Documentation is [here](https://WaDelma.github.io/poisson/)    

# TODO
   * Can quadtree be compressed? Skip Quadtree?    
   * Benchmarks!    
