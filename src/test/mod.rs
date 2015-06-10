use ::{PoissonDisk, Sample, VecLike};

use rand::{SeedableRng, XorShiftRng};

use na::Norm;

mod common;
mod dim2;
mod dim3;
mod dim4;

fn test_with_samples<T: VecLike<T>>(samples: u32, relative_radius: f64, seeds: u32, periodicity: bool) {
    for i in 0..seeds {
        let rand = XorShiftRng::from_seed([i + 1, seeds - i + 1, (i + 1) * (i + 1), 1]);
        let mut poisson = PoissonDisk::<_, T>::with_samples(rand, samples, relative_radius, periodicity);
        let mut vecs = vec![];
        poisson.create(&mut vecs);
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
                    vecs2.push(Sample::new(v.pos + t, v.radius()));
                }
            }
            vecs2
        } else {
            vecs
        };
        assert_legal_poisson(&vecs);
    }
}

// fn test_with_seeds<T: VecLike<T>>(radius: f64, seeds: u32, periodicity: bool) {
//     test_with_seeds_prefill::<T, _>(radius, seeds, periodicity, &mut |_, _|{});
// }

fn test_with_seeds_prefill<T: VecLike<T>, F>(radius: f64, seeds: u32, periodicity: bool, filler: &mut F) where F: FnMut(&mut Vec<Sample<T>>, u32) {
    for i in 0..seeds {
        let rand = XorShiftRng::from_seed([i + 1, seeds - i + 1, (i + 1) * (i + 1), 1]);
        let mut poisson = PoissonDisk::with_radius(rand, radius, periodicity);
        let mut vecs = vec![];
        filler(&mut vecs, i);
        poisson.create(&mut vecs);
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
                    vecs2.push(Sample::new(v.pos + t, v.radius()));
                }
            }
            vecs2
        } else {
            vecs
        };
        assert_legal_poisson(&vecs);
    }
}

pub fn assert_legal_poisson<T: VecLike<T>>(vecs: &Vec<Sample<T>>) {
    for &v1 in vecs {
        for &v2 in vecs {
            if v1.pos == v2.pos {
                continue;
            }
            let diff = v1.pos - v2.pos;
            let dist = diff.norm();
            let allowed_dist = v1.radius() + v2.radius();
            assert!(dist >= allowed_dist, "Poisson-disk distribution requirement not met: There exists 2 vectors with distance to each other of {} which is smaller than smallest allowed one {}. The samples: [{:?}, {:?}]", dist, allowed_dist, v1, v2);
        }
    }
}
