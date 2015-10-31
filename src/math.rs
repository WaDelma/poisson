use VecLike;
use na::Dim;
use std::f64::consts::PI;

lazy_static! {
    static ref MAX_PACKING_DENSITIES: [f64; 7] = [
        1. / 6. * PI * 3f64.sqrt(),
        1. / 6. * PI * 2f64.sqrt(),
        1. / 16. * PI.powi(2),
        1. / 30. * PI.powi(2) * 2f64.sqrt(),
        1. / 144. * PI.powi(3) * 3f64.sqrt(),
        1. / 105. * PI.powi(3),
        1. / 384. * PI.powi(4),
        ];
    // gamma((index + 2) / 2 + 1)
    static ref GAMMA: [f64; 7] = [
        1.,
        (3. * PI.sqrt()) / 4.,
        2.,
        (15. * PI.sqrt()) / 8.,
        6.,
        (105. * PI.sqrt()) / 16.,
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
    (MAX_PACKING_DENSITIES[index] * GAMMA[index]) / PI.powf(dim as f64 / 2.)
}

fn newton(samples: u32, dim: usize) -> u32 {
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
    n as u32
}

pub fn calc_radius<T>(samples: u32, relative_radius: f64, periodicity: bool) -> f64
    where T: VecLike
{
    let dim = T::dim(None) as usize;
    let samples = if periodicity {
        samples
    } else {
        newton(samples, dim)
    };
    (MAX_RADII[dim - 2] / samples as f64).powf(1. / dim as f64) * relative_radius
}
