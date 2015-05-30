use ::{PoissonDisk, Sample, VecLike};

use rand::{SeedableRng, XorShiftRng};

use na::Norm;

mod common;
#[cfg(not(feature = "no2d"))]
mod dim2;
#[cfg(not(feature = "no3d"))]
mod dim3;

fn test_with_seeds<T: VecLike<T>>(radius: f64, seeds: u32, periodicity: bool) {
    test_with_seeds_prefill::<T, _>(radius, seeds, periodicity, &mut |_, _|{});
}

fn test_with_seeds_prefill<T: VecLike<T>, F>(radius: f64, seeds: u32, periodicity: bool, filler: &mut F) where F: FnMut(&mut Vec<Sample<T>>, u32) {
    for i in 0..seeds {
        let rand = XorShiftRng::from_seed([i + 1, seeds - i + 1, (i + 1) * (i + 1), 1]);
        let mut poisson = if periodicity {
            PoissonDisk::perioditic(rand, radius)
        } else {
            PoissonDisk::new(rand, radius)
        };
        let mut vecs = vec![];
        filler(&mut vecs, i);
        poisson.create(&mut vecs);
        let vecs = if periodicity {
            let mut vecs2 = vec![];
            let dim = T::dim(None);
            for n in 0..3u64.pow(dim as u32) {
                let mut t = T::zero();
                for i in 0..dim {
                    t[i] = match (n / (1 + i as u64 * 2)) % 3 {
                        0 => -1.0,
                        1 => 0.0,
                        2 => 1.0,
                        j @ _ => unreachable!("This shouldn't be possible, but Rust cannot figure it out. Failed with: {}", j)
                    }
                }
                for v in &vecs {
                    vecs2.push(*v + t);
                }
            }
            vecs2
        } else {
            vecs
        };
        assert_legal_poisson(&vecs);
    }
}

fn assert_legal_poisson<T: VecLike<T>>(vecs: &Vec<Sample<T>>) {
    for &v1 in vecs {
        for &v2 in vecs {
            if v1.pos == v2.pos {
                continue;
            }
            let diff = v1.pos - v2.pos;
            let dist = diff.norm();
            let allowed_dist = v1.get_radius() + v2.get_radius();
            assert!(dist >= allowed_dist, "Poisson-disk distribution requirement not met: There exists 2 vectors with distance to each other of {} which is smaller than smallest allowed one {}. The samples: [{:?}, {:?}]", dist, allowed_dist, v1, v2);
        }
    }
}
