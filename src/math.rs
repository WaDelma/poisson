use ::VecLike;

use rand::{Rand, Rng};

use num::Zero;

use na::Dim;

use std::fmt::{Debug, Formatter};
use std::fmt::Result as FmtResult;
use std::f64::consts::PI;

#[derive(Clone, Copy)]
pub struct Hypercube<T>  {
    pub min: T,
    pub max: T,
}

impl<T: VecLike<T>> Hypercube<T> {
    #[inline]
    pub fn new(min: T, max: T) -> Self {
        Hypercube{min: min, max: max}
    }

    #[inline]
    pub fn volume(&self) -> f64 {
        self.edge().powi(T::dim(None) as i32)
    }

    #[inline]
    pub fn edge(&self) -> f64 {
        self.max[0] - self.min[0]
    }

    #[inline]
    pub fn center(&self) -> T {
        (self.max + self.min) / 2.0
    }

    #[inline]
    pub fn random_point_inside<F>(&self, rand: &mut F) -> T where F: Rng {
        let mut t = T::zero();
        let dim = T::dim(None);
        for n in 0..dim {
            t[n] = f64::rand(rand).mul_add(self.edge(), self.min[n]);
        }
        t
    }
}

impl<T: VecLike<T>> Debug for Hypercube<T> {
    #[inline]
    fn fmt(&self, f: &mut Formatter) -> FmtResult {
        write!(f, "Hypercube({:?}, {:?})", self.min, self.max)
    }
}

impl<T: VecLike<T>> PartialEq for Hypercube<T> {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.min == other.min && self.max == other.max
    }
}

pub enum Intersection {
    Out,
    Over,
    In,
}

#[inline]
pub fn test_intersection<T: VecLike<T>>(rect: Hypercube<T>, sample: T, radius: f64) -> Intersection {
    let radius2 = radius * radius;
    let dims = T::dim(None);
    let mut min = 0.0;
    let mut max = 0.0;
    for i in 0usize..dims {
        let cur_min = rect.min[i];
        let cur_max = rect.max[i];
        let cur_sample = sample[i];
        let dmin = cur_min - cur_sample;
        if dmin > 0.0 {
            if dmin > radius {
                return Intersection::Out;
            }
            min += dmin * dmin;
            if min > radius2 {
                return Intersection::Out;
            }
            let temp_max = cur_max - cur_sample;
            max += temp_max * temp_max;
            continue;
        }
        let dmin = cur_sample - cur_max;
        if dmin > 0.0 {
            if dmin > radius {
                return Intersection::Out;
            }
            min += dmin * dmin;
            if min > radius2 {
                return Intersection::Out;
            }
            let temp_max = cur_sample - cur_min;
            max += temp_max * temp_max;
            continue;
        }
        let temp_max = (cur_sample - cur_min).max(cur_max - cur_sample);
        max += temp_max * temp_max;
    }
    if max > radius2 {
        return Intersection::Over;
    }
    return Intersection::In;
}

lazy_static! {
    static ref MAX_PACKING_DENSITIES: [f64; 7] = [
        1.0 / 6.0 * PI * 3f64.sqrt(),
        1.0 / 6.0 * PI * 2f64.sqrt(),
        1.0 / 16.0 * PI.powi(2),
        1.0 / 30.0 * PI.powi(2) * 2f64.sqrt(),
        1.0 / 144.0 * PI.powi(3) * 3f64.sqrt(),
        1.0 / 105.0 * PI.powi(3),
        1.0 / 384.0 * PI.powi(4),
        ];
    // gamma((index + 2) / 2 + 1)
    static ref GAMMA: [f64; 7] = [
        1.0,
        (3.0 * PI.sqrt()) / 4.0,
        2.0,
        (15.0 * PI.sqrt()) / 8.0,
        6.0,
        (105.0 * PI.sqrt()) / 16.0,
        24.0,
        ];
    static ref MAX_RADII: [f64; 7] = [
            precalc(2),
            precalc(3),
            precalc(4),
            precalc(5),
            precalc(6),
            precalc(7),
            precalc(8),
        ];
    //TODO: Paper provides needed constants only for 2, 3 and 4 dimensions.
    static ref ALPHA: [f64; 3] = [
            1.0997,
            2.2119,
            4.1114,
        ];
    static ref BETA: [f64; 3] = [
            -0.4999,
            -0.3538,
            -0.3056,
        ];
}

fn precalc(dim: usize) -> f64 {
    let index = dim - 2;
    (MAX_PACKING_DENSITIES[index] * GAMMA[index]) / PI.powf(dim as f64 / 2.0)
}

fn newton(samples: u32, dim: usize) -> u32 {
    let mut n = 1f64;
    let alpha = ALPHA[dim - 2];
    let beta = BETA[dim - 2];
    for _ in 0..5 {
        n = n - (n + alpha * n.powf(beta + 1.0) - samples as f64) / (1.0 + alpha * (beta + 1.0) * n.powf(beta));
        if n < 1.0 {
            return 1;
        }
    }
    n as u32
}

pub fn calc_radius<T: VecLike<T>>(samples: u32, relative_radius: f64, periodicity: bool) -> f64 {
    let dim = T::dim(None) as usize;
    let samples = if periodicity { samples } else {
        newton(samples, dim)
    };
    (MAX_RADII[dim - 2] / samples as f64).powf(1.0 / dim as f64) * relative_radius
}
