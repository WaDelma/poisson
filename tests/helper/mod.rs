#![allow(unused)]
use poisson::{Ebeida, PoissonType, PoissonIter, PoissonDisk, VecLike, FloatLike};

use rand::{SeedableRng, XorShiftRng};

use std::fmt::Debug;

extern crate num;
use self::num::Float;

use na::Norm;

pub fn print_v<F: FloatLike + Debug, V: VecLike<F>>(v: V) -> String {
    let mut result = "(".to_owned();
    for i in v.iter() {
        result.push_str(&format!("{:?}, ", i));
    }
    if V::dim(None) != 0 {
        result.pop();
    }
    result.push(')');
    result
}


pub enum When {
    Always,
    Sometimes,
    Never,
}

pub fn test_with_samples<T>(samples: usize, relative_radius: f64, seeds: u32, ptype: PoissonType)
    where T: Debug + VecLike<f64> + Copy
{
    test_with_samples_prefilled(samples, relative_radius, seeds, ptype, |_| |_| None::<T>, When::Always);
}

pub fn test_with_samples_prefilled<'r, T, F, I>(samples: usize, relative_radius: f64, seeds: u32, ptype: PoissonType, mut prefiller: F, valid: When)
    where T: 'r + Debug + VecLike<f64> + Copy, F: FnMut(f64) -> I, I: FnMut(Option<T>) -> Option<T>
{
    use self::When::*;
    for i in 0..seeds {
        unsafe{::poisson::SEED = i as usize}
        let mut prefilled = vec![];
        let rand = XorShiftRng::from_seed([i + 1, seeds - i + 1, (i + 1) * (i + 1), 1]);
        let mut poisson = PoissonDisk::with_samples(samples, relative_radius, ptype);//new(rand);
        let mut poisson_iter = poisson.build(rand, Ebeida).into_iter();
        let mut poisson = vec![];
        let mut prefil = (prefiller)(poisson_iter.radius());
        let mut last = None;
        loop {
            while let Some(p) = (prefil)(last) {
                match valid {
                    Always => assert!(poisson_iter.stays_legal(p), "All prefilled should be accepted by the algorithm. \
                                    {} was rejected.", print_v(p)),
                    Never => assert!(!poisson_iter.stays_legal(p), "All prefilled should be rejected by the algorithm. \
                                    {} was allowed even though {:?} was last to be generated.", print_v(p), last.map(print_v)),
                    _ => {},
                }
                prefilled.push(p);
                poisson_iter.restrict(p);
            }
            if let Some(pp) = poisson_iter.next() {
                last = Some(pp);
                poisson.push(pp);
            } else {
                break;
            }
        }
        let radius = poisson_iter.radius();
        let poisson_type = poisson_iter.poisson_type();
        let poisson = poisson
            .into_iter()
            .chain(if let Always = valid {
                prefilled
            } else {
                vec![]
            }.into_iter());
        test_poisson(poisson, radius, poisson_type);
        // break;
    }
}

pub fn test_poisson<I, T>(poisson: I, radius: f64, poisson_type: PoissonType)
    where I: Iterator<Item=T>, T: Debug + VecLike<f64> + Copy
{
    use poisson::PoissonType::*;
    let dim = T::dim(None);
    let mut vecs = vec![];
    let mut hints = vec![];
    {
        let mut iter = poisson.into_iter();
        while let Some(v) = iter.next() {
            if let (low, Some(high)) = iter.size_hint() {
                hints.push((low, high));
            } else {
                panic!("There wasn't hint for {}th iteration.", hints.len());
            }
            vecs.push(v);
        }
    }
    let len = hints.len();
    for (n, (l, h)) in hints.into_iter().enumerate() {
        let remaining = len - (n + 1);
        assert!(l <= remaining, "Lower bound of hint should be smaller than or equal to actual: {} <= {}", l, remaining);
        assert!(h >= remaining, "Upper bound of hint should be larger than or equal to actual: {} >= {}", h, remaining);
    }
    //TODO: Figure out how to check if distribution is maximal.
    // let packing_density = vecs.len() as f64 * ::sphere::sphere_volume(poisson.radius(), dim as u64);
    // println!("{}", packing_density);
    // panic!();
    let vecs = match poisson_type {
        Perioditic => {
            let mut vecs2 = vec![];
            for n in 0..3i64.pow(dim as u32) {
                let mut t = T::zero();
                let mut div = n;
                for i in t.iter_mut() {
                    let rem = div % 3;
                    div /= 3;
                    *i = (rem - 1) as f64;
                }
                for v in &vecs {
                    vecs2.push(*v + t);
                }
            }
            vecs2
        },
        Normal => vecs,
    };
    assert_legal_poisson(&vecs, radius);
}

pub fn assert_legal_poisson<T>(vecs: &Vec<T>, radius: f64)
    where T: Debug + VecLike<f64> + Copy
{
    for &v1 in vecs {
        for &v2 in vecs {
            if v1 == v2 {
                continue;
            }
            let dist = (v1 - v2).norm();
            assert!(dist >= radius * 2.,
                    "Poisson-disk distribution requirement not met: There exists 2 vectors with \
                     distance to each other of {} which is smaller than smallest allowed one {}. \
                     The samples: [{:?}, {:?}]",
                    dist,
                    radius * 2.,
                    v1,
                    v2);
        }
    }
}
