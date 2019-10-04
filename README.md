# poisson

[![Documentation][di]][dl] [![Crates.io][ri]][rl] [![License: MIT][li]][ll] [![Build Status][ti]][tl] [![Coverage Status][ci]][cl]   

[di]: https://docs.rs/poisson/badge.svg
[dl]: https://docs.rs/poisson

[ri]: https://img.shields.io/crates/v/poisson.svg
[rl]: https://crates.io/crates/poisson/

[li]: https://img.shields.io/badge/License-MIT-blue.svg
[ll]: https://opensource.org/licenses/MIT

[ti]: https://travis-ci.org/WaDelma/poisson.svg?branch=master
[tl]: https://travis-ci.org/WaDelma/poisson

[ci]: https://coveralls.io/repos/github/WaDelma/poisson/badge.svg?branch=master
[cl]: https://coveralls.io/github/WaDelma/poisson?branch=master

This is a library for generating n-dimensional [poisson-disk distributions](http://mollyrocket.com/casey/stream_0014.html).    

It generates distribution of points in [0, 1)<sup>d</sup> where:

 * For each point there is disk of certain radius which doesn't intersect
 with any other disk of other points
 * Samples fill the space uniformly

Due it's blue noise properties poisson-disk distribution
can be used for object placement in procedural texture/world generation,
as source distribution for digital stipling,
as distribution for sampling in rendering or for (re)meshing.

# Usage

Works with nalgebra 0.16 and rand 0.5

```rust
extern crate nalgebra as na;

use rand::FromEntropy;
use rand::rngs::SmallRng;

use poisson::{Builder, Type, algorithm};

fn main() {
    let poisson =
        Builder::<_, na::Vector2<f64>>::with_radius(0.1, Type::Normal)
            .build(SmallRng::from_entropy(), algorithm::Ebeida);
    println!("{:?}", poisson.generate());
}
```
