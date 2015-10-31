use poisson::{PoissonDisk, VecLike};

use rand::{SeedableRng, XorShiftRng};

use std::fmt::Debug;

use na::Norm;

pub fn test_with_samples<T: Debug + VecLike>(samples: u32,
                                         relative_radius: f64,
                                         seeds: u32,
                                         periodicity: bool) {
    for i in 0..seeds {
        let rand = XorShiftRng::from_seed([i + 1, seeds - i + 1, (i + 1) * (i + 1), 1]);
        let mut poisson = PoissonDisk::new(rand);
        if periodicity {
            poisson = poisson.perioditic();
        }
        let mut poisson = poisson.build_samples::<T>(samples, relative_radius);
        let vecs = poisson.generate();
        let vecs = if periodicity {
            let mut vecs2 = vec![];
            let dim = T::dim(None);
            for n in 0..3i64.pow(dim as u32) {
                let mut t = T::zero();
                let mut div = n;
                for i in 0..dim {
                    let rem = div % 3;
                    div /= 3;
                    t[i] = (rem - 1) as f64;
                }
                for v in &vecs {
                    vecs2.push(*v + t);
                }
            }
            vecs2
        } else {
            vecs
        };
        assert_legal_poisson(&vecs, poisson.radius());
    }
}

pub fn assert_legal_poisson<T: Debug + VecLike>(vecs: &Vec<T>, radius: f64) {
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
