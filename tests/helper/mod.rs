#![allow(unused)]
use poisson::{Type, Builder, Vector, Float, algorithm};

use rand::{SeedableRng, XorShiftRng};

extern crate num_traits as num;
use self::num::NumCast;

use na::Norm;

use std::fmt::Debug;

pub fn print_v<F: Float, V: Vector<F>>(v: V) -> String {
    let mut result = "(".to_owned();
    for i in v.iter() {
        result.push_str(&format!("{}, ", i.to_f64().unwrap()));
    }
    if V::dimension(None) != 0 {
        result.pop();
    }
    result.push(')');
    result
}

#[derive(Clone, Copy)]
pub enum When {
    Always,
    Sometimes,
    Never,
}

pub fn test_with_samples<T>(samples: usize, relative_radius: f64, seeds: u32, ptype: Type)
    where T: Debug + Vector<f64> + Copy
{
    test_with_samples_prefilled(samples, relative_radius, seeds, ptype, |_| |_| None::<T>, When::Always);
}

pub fn test_with_samples_prefilled<'r, T, F, I>(samples: usize, relative_radius: f64, seeds: u32, ptype: Type, mut prefiller: F, valid: When)
    where T: 'r + Debug + Vector<f64> + Copy, F: FnMut(f64) -> I, I: FnMut(Option<T>) -> Option<T>
{
    test_algo(samples, relative_radius, seeds, ptype, &mut prefiller, valid, algorithm::Ebeida);
    test_algo(samples, relative_radius, seeds, ptype, &mut prefiller, valid, algorithm::Bridson);
}

fn test_algo<'r, T, F, I, A>(samples: usize, relative_radius: f64, seeds: u32, ptype: Type, prefiller: &mut F, valid: When, algo: A)
    where T: 'r + Debug + Vector<f64> + Copy, F: FnMut(f64) -> I, I: FnMut(Option<T>) -> Option<T>, A: algorithm::Creator<f64, T>
{
    use self::When::*;
    for i in 0..seeds {
        let mut prefilled = vec![];
        let rand = XorShiftRng::from_seed([i + 1, seeds - i + 1, (i + 1) * (i + 1), 1]);
        let mut poisson_iter = Builder::with_samples(samples, relative_radius, ptype).build(rand, algo).into_iter();
        let mut poisson = vec![];
        let mut prefill = (prefiller)(poisson_iter.radius());
        let mut last = None;
        loop {
            while let Some(p) = (prefill)(last) {
                match valid {
                    Always => assert!(poisson_iter.stays_legal(p), "All prefilled should be accepted by the '{:?}' algorithm. \
                                    {} was rejected.", algo, print_v(p)),
                    Never => assert!(!poisson_iter.stays_legal(p), "All prefilled should be rejected by the '{:?}' algorithm. \
                                    {} was allowed even though {:?} was last to be generated.", algo, print_v(p), last.map(print_v)),
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
        test_poisson(poisson, radius, poisson_type, algo);
    }
}

pub fn test_poisson<F, I, T, A>(poisson: I, radius: F, poisson_type: Type, algo: A)
    where I: Iterator<Item=T>, F: Float, T: Debug + Vector<F> + Copy, A: algorithm::Creator<F, T>
{
    use poisson::Type::*;
    let dim = T::dimension(None);
    let mut vecs = vec![];
    let mut hints = vec![];
    {
        let mut iter = poisson.into_iter();
        while let Some(v) = iter.next() {
            if let (low, Some(high)) = iter.size_hint() {
                hints.push((low, high));
            } else {
                panic!("There wasn't hint for {}th iteration for the '{:?}' algorithm.", hints.len(), algo);
            }
            vecs.push(v);
        }
    }
    let len = hints.len();
    for (n, (l, h)) in hints.into_iter().enumerate() {
        let remaining = len - (n + 1);
        assert!(l <= remaining, "For the '{:?}' algorithm the lower bound of hint should be smaller than or equal to actual: {} <= {}", algo, l, remaining);
        assert!(h >= remaining, "For the '{:?}' algorithm the upper bound of hint should be larger than or equal to actual: {} >= {}", algo, h, remaining);
    }

    let vecs = match poisson_type {
        Perioditic => {
            let mut vecs2 = vec![];
            for n in 0..3i64.pow(dim as u32) {
                let mut t = T::zero();
                let mut div = n;
                for i in t.iter_mut() {
                    let rem = div % 3;
                    div /= 3;
                    *i = NumCast::from(rem - 1).unwrap();
                }
                for v in &vecs {
                    vecs2.push(*v + t);
                }
            }
            vecs2
        },
        Normal => vecs,
    };

    //TODO: Figure out how to check if distribution is maximal.
    assert_legal_poisson(&vecs, radius, algo);
}

pub fn assert_legal_poisson<F, T, A>(vecs: &Vec<T>, radius: F, algo: A)
    where F: Float, T: Debug + Vector<F> + Copy, A: algorithm::Creator<F, T>
{
    for &v1 in vecs {
        for &v2 in vecs {
            if v1 == v2 {
                continue;
            }
            let dist = (v1 - v2).norm();
            assert!(dist > radius * F::cast(2),
                    "Poisson-disk distribution requirement not met while generating using the '{:?}' algorithm: There exists 2 vectors with \
                     distance to each other of {} which is smaller than smallest allowed one {}. \
                     The samples: [{:?}, {:?}]",
                    algo,
                    dist.to_f64().unwrap(),
                    radius.to_f64().unwrap() * 2.,
                    v1,
                    v2);
        }
    }
}
