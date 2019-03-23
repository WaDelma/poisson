#![feature(test)]

extern crate test;
use test::{Bencher, black_box};

extern crate poisson;
use poisson::{Builder, Type, algorithm};

extern crate rand;
use rand::{SeedableRng, XorShiftRng};

extern crate nalgebra as na;
pub type Vect = na::Vector3<f64>;

#[bench]
fn bench_ebeida_3d_1_80_normal(b: &mut Bencher) {
    let rand = XorShiftRng::from_seed([1, 2, 3, 4]);
    let poisson =
        Builder::<_, Vect>::with_samples(1, 0.8, Type::Normal)
            .build(rand, algorithm::Ebeida);
    b.iter(|| black_box(poisson.generate()));
}

#[bench]
fn bench_ebeida_3d_10_80_normal(b: &mut Bencher) {
    let rand = XorShiftRng::from_seed([1, 2, 3, 4]);
    let poisson =
        Builder::<_, Vect>::with_samples(10, 0.8, Type::Normal)
            .build(rand, algorithm::Ebeida);
    b.iter(|| black_box(poisson.generate()));
}

#[bench]
fn bench_ebeida_3d_100_80_normal(b: &mut Bencher) {
    let rand = XorShiftRng::from_seed([1, 2, 3, 4]);
    let poisson =
        Builder::<_, Vect>::with_samples(100, 0.8, Type::Normal)
            .build(rand, algorithm::Ebeida);
    b.iter(|| black_box(poisson.generate()));
}

#[bench]
fn bench_bridson_3d_1_80_normal(b: &mut Bencher) {
    let rand = XorShiftRng::from_seed([1, 2, 3, 4]);
    let poisson =
        Builder::<_, Vect>::with_samples(1, 0.8, Type::Normal)
            .build(rand, algorithm::Bridson);
    b.iter(|| black_box(poisson.generate()));
}

#[bench]
fn bench_bridson_3d_10_80_normal(b: &mut Bencher) {
    let rand = XorShiftRng::from_seed([1, 2, 3, 4]);
    let poisson =
        Builder::<_, Vect>::with_samples(10, 0.8, Type::Normal)
            .build(rand, algorithm::Bridson);
    b.iter(|| black_box(poisson.generate()));
}

#[bench]
fn bench_bridson_3d_100_80_normal(b: &mut Bencher) {
    let rand = XorShiftRng::from_seed([1, 2, 3, 4]);
    let poisson =
        Builder::<_, Vect>::with_samples(100, 0.8, Type::Normal)
            .build(rand, algorithm::Bridson);
    b.iter(|| black_box(poisson.generate()));
}

#[bench]
fn bench_ebeida_3d_1_80_perioditic(b: &mut Bencher) {
    let rand = XorShiftRng::from_seed([1, 2, 3, 4]);
    let poisson =
        Builder::<_, Vect>::with_samples(1, 0.8, Type::Perioditic)
            .build(rand, algorithm::Ebeida);
    b.iter(|| black_box(poisson.generate()));
}

#[bench]
fn bench_ebeida_3d_10_80_perioditic(b: &mut Bencher) {
    let rand = XorShiftRng::from_seed([1, 2, 3, 4]);
    let poisson =
        Builder::<_, Vect>::with_samples(10, 0.8, Type::Perioditic)
            .build(rand, algorithm::Ebeida);
    b.iter(|| black_box(poisson.generate()));
}

#[bench]
fn bench_ebeida_3d_100_80_perioditic(b: &mut Bencher) {
    let rand = XorShiftRng::from_seed([1, 2, 3, 4]);
    let poisson =
        Builder::<_, Vect>::with_samples(100, 0.8, Type::Perioditic)
            .build(rand, algorithm::Ebeida);
    b.iter(|| black_box(poisson.generate()));
}

#[bench]
fn bench_bridson_3d_1_80_perioditic(b: &mut Bencher) {
    let rand = XorShiftRng::from_seed([1, 2, 3, 4]);
    let poisson =
        Builder::<_, Vect>::with_samples(1, 0.8, Type::Perioditic)
            .build(rand, algorithm::Bridson);
    b.iter(|| black_box(poisson.generate()));
}

#[bench]
fn bench_bridson_3d_10_80_perioditic(b: &mut Bencher) {
    let rand = XorShiftRng::from_seed([1, 2, 3, 4]);
    let poisson =
        Builder::<_, Vect>::with_samples(10, 0.8, Type::Perioditic)
            .build(rand, algorithm::Bridson);
    b.iter(|| black_box(poisson.generate()));
}

#[bench]
fn bench_bridson_3d_100_80_perioditic(b: &mut Bencher) {
    let rand = XorShiftRng::from_seed([1, 2, 3, 4]);
    let poisson =
        Builder::<_, Vect>::with_samples(100, 0.8, Type::Perioditic)
            .build(rand, algorithm::Bridson);
    b.iter(|| black_box(poisson.generate()));
}
