use ::VecLike;

use rand::{Rand, Rng};

use num::Zero;

use na::Dim;

use std::fmt::{Debug, Formatter};
use std::fmt::Result as FmtResult;

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
        let mut result = 1.0;
        let dim = T::dim(None);
        for _ in 0..dim {
            result *= self.edge();
        }
        result
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
    let dims = T::dim(None);//2usize;
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
