use ::{PoissonType, VecLike, FloatLike};

use num::{NumCast, Float};

use na::Dim;

// const TAU: f64
//     = 6.283185307179586476925286766559005768394338798750211641949;
const HALF_TAU: f64
    = 3.141592653589793238462643383279502884197169399375105820974;

lazy_static! {
    static ref MAX_PACKING_DENSITIES: [f64; 7] = [
        1. / 6. * HALF_TAU * 3f64.sqrt(),
        1. / 6. * HALF_TAU * 2f64.sqrt(),
        1. / 16. * HALF_TAU.powi(2),
        1. / 30. * HALF_TAU.powi(2) * 2f64.sqrt(),
        1. / 144. * HALF_TAU.powi(3) * 3f64.sqrt(),
        1. / 105. * HALF_TAU.powi(3),
        1. / 384. * HALF_TAU.powi(4),
        ];
    // gamma((index + 2) / 2 + 1)
    static ref GAMMA: [f64; 7] = [
        1.,
        (3. * HALF_TAU.sqrt()) / 4.,
        2.,
        (15. * HALF_TAU.sqrt()) / 8.,
        6.,
        (105. * HALF_TAU.sqrt()) / 16.,
        24.,
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
    (MAX_PACKING_DENSITIES[index] * GAMMA[index]) / HALF_TAU.powf(dim as f64 / 2.)
}

fn newton(samples: usize, dim: usize) -> usize {
    let mut n = 1f64;
    let alpha = ALPHA[dim - 2];
    let beta = BETA[dim - 2];
    for _ in 0..5 {
        n = n -
            (n + alpha * n.powf(beta + 1.) - samples as f64) /
            (1. + alpha * (beta + 1.) * n.powf(beta));
        if n < 1. {
            return 1;
        }
    }
    n as usize
}

/// Calculates radius from approximate samples and relative radius.
/// The amount of samples should be larger than 0.
/// The relative radius should be [0, 1].
/// For non-perioditic this is supported only for 2, 3 and 4 dimensional generation.
pub fn calc_radius<F, V>(samples: usize, relative: F, poisson_type: PoissonType) -> F
    where F: FloatLike,
          V: VecLike<F>
{
    use PoissonType::*;
    assert!(PoissonType::Perioditic == poisson_type || V::dim(None) < 5);
    assert!(samples > 0);
    assert!(relative >= F::f(0));
    assert!(relative <= F::f(1));
    let dim = V::dim(None);
    let samples = match poisson_type {
        Perioditic => samples,
        Normal => newton(samples, dim),
    };
    let max_radii: F = NumCast::from(MAX_RADII[dim - 2]).unwrap();
    (max_radii / F::f(samples)).powf(F::f(1) / F::f(dim)) * relative
}