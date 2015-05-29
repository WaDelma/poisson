# poisson 
**poisson** is library for generating 2d poisson-disk distributions.

[![Build Status](https://travis-ci.org/WaDelma/poisson.svg?branch=master)](https://travis-ci.org/WaDelma/poisson)   

## Using **poisson**
**poisson** is used through struct PoissonDisk.   
PoissonDisk can be created with perioditic constructor for distribution to wrap around.    
Function calc_radius can be used to calculate radius based on amount of samples wanted and relative radius.    
Samples can be given to the method for creating distribution to make it take those account when creating the distribution.    

Documentation is [here](https://WaDelma.github.io/poisson/)  

# TODO
   * Generalise to arbitary dimensions
     * Requires to be generic over vector dimensions
     * Requires hyperoctree
   * Can quadtree be compressed? Skip Quadtree?
   * Benchmarks!
