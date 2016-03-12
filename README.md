# poisson    
**poisson** is library for generating n-dimensional [poisson-disk distributions](http://devmag.org.za/2009/05/03/poisson-disk-sampling/).    

[![Build Status](https://travis-ci.org/WaDelma/poisson.svg?branch=master)](https://travis-ci.org/WaDelma/poisson)   

## Using **poisson**    
**poisson** is used through struct PoissonGen which is built using PoissonDisk.    

There is currently two different algorithms in this library:
* Ebeida's algorithm that generates uniform maximum poisson-disk distributions.
* Bridson's algorithm that generates non-uniform non-maximal poisson-disk distributions.

Building of PoissonDisk with method build_samples for non-perioditic distributions is currently restricted to 2, 3 and 4 dimensions (because there is needed constants for only those dimensions).  

Documentation is [here](https://WaDelma.github.io/poisson/)    
