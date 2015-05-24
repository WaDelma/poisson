# poisson
**poisson** is library for generating 2d poisson-disk distributions.

## Using **poisson**
**poisson** is used through struct PoissonDisk. 
Function calc_radius can be used to calculate radius based on amount of samples wanted and relative radius.


# TODO
   * Generalise to arbitary dimensions
     * Requires to be generic over vector dimensions
     * Requires hyperoctree
   * Generating takes consideration of already existing points
     * Multiple different radii?
   * Can quadtree be compressed? Skip Quadtree?
   * Benchmarks!