#![feature(test)]

extern crate test;
use test::{black_box, Bencher};

use poisson::{algorithm, Builder, Type};

use rand::{rngs::SmallRng, SeedableRng};

extern crate nalgebra as na;
pub type Vect = na::Vector2<f64>;

const SEED: [u8; 16] = [
    (3 + 2741) as u8,
    (7 + 2729) as u8,
    (13 + 2713) as u8,
    (19 + 2707) as u8,
    (29 + 2693) as u8,
    (37 + 2687) as u8,
    (43 + 2677) as u8,
    (53 + 2663) as u8,
    (61 + 2657) as u8,
    (71 + 2633) as u8,
    (79 + 2609) as u8,
    (89 + 2591) as u8,
    (101 + 2557) as u8,
    (107 + 2549) as u8,
    (113 + 2539) as u8,
    (131 + 2521) as u8,
];

#[bench]
fn bench_ebeida_2d_1_80_normal(b: &mut Bencher) {
    let rand = SmallRng::from_seed(SEED);
    let poisson =
        Builder::<_, Vect>::with_samples(1, 0.8, Type::Normal).build(rand, algorithm::Ebeida);
    b.iter(|| black_box(poisson.generate()));
}

#[bench]
fn bench_ebeida_2d_10_80_normal(b: &mut Bencher) {
    let rand = SmallRng::from_seed(SEED);
    let poisson =
        Builder::<_, Vect>::with_samples(10, 0.8, Type::Normal).build(rand, algorithm::Ebeida);
    b.iter(|| black_box(poisson.generate()));
}

#[bench]
fn bench_ebeida_2d_100_80_normal(b: &mut Bencher) {
    let rand = SmallRng::from_seed(SEED);
    let poisson =
        Builder::<_, Vect>::with_samples(100, 0.8, Type::Normal).build(rand, algorithm::Ebeida);
    b.iter(|| black_box(poisson.generate()));
}

#[bench]
fn bench_bridson_2d_1_80_normal(b: &mut Bencher) {
    let rand = SmallRng::from_seed(SEED);
    let poisson =
        Builder::<_, Vect>::with_samples(1, 0.8, Type::Normal).build(rand, algorithm::Bridson);
    b.iter(|| black_box(poisson.generate()));
}

#[bench]
fn bench_bridson_2d_10_80_normal(b: &mut Bencher) {
    let rand = SmallRng::from_seed(SEED);
    let poisson =
        Builder::<_, Vect>::with_samples(10, 0.8, Type::Normal).build(rand, algorithm::Bridson);
    b.iter(|| black_box(poisson.generate()));
}

#[bench]
fn bench_bridson_2d_100_80_normal(b: &mut Bencher) {
    let rand = SmallRng::from_seed(SEED);
    let poisson =
        Builder::<_, Vect>::with_samples(100, 0.8, Type::Normal).build(rand, algorithm::Bridson);
    b.iter(|| black_box(poisson.generate()));
}

#[bench]
fn bench_ebeida_2d_1_80_perioditic(b: &mut Bencher) {
    let rand = SmallRng::from_seed(SEED);
    let poisson =
        Builder::<_, Vect>::with_samples(1, 0.8, Type::Perioditic).build(rand, algorithm::Ebeida);
    b.iter(|| black_box(poisson.generate()));
}

#[bench]
fn bench_ebeida_2d_10_80_perioditic(b: &mut Bencher) {
    let rand = SmallRng::from_seed(SEED);
    let poisson =
        Builder::<_, Vect>::with_samples(10, 0.8, Type::Perioditic).build(rand, algorithm::Ebeida);
    b.iter(|| black_box(poisson.generate()));
}

#[bench]
fn bench_ebeida_2d_100_80_perioditic(b: &mut Bencher) {
    let rand = SmallRng::from_seed(SEED);
    let poisson =
        Builder::<_, Vect>::with_samples(100, 0.8, Type::Perioditic).build(rand, algorithm::Ebeida);
    b.iter(|| black_box(poisson.generate()));
}

#[bench]
fn bench_bridson_2d_1_80_perioditic(b: &mut Bencher) {
    let rand = SmallRng::from_seed(SEED);
    let poisson =
        Builder::<_, Vect>::with_samples(1, 0.8, Type::Perioditic).build(rand, algorithm::Bridson);
    b.iter(|| black_box(poisson.generate()));
}

#[bench]
fn bench_bridson_2d_10_80_perioditic(b: &mut Bencher) {
    let rand = SmallRng::from_seed(SEED);
    let poisson =
        Builder::<_, Vect>::with_samples(10, 0.8, Type::Perioditic).build(rand, algorithm::Bridson);
    b.iter(|| black_box(poisson.generate()));
}

#[bench]
fn bench_bridson_2d_100_80_perioditic(b: &mut Bencher) {
    let rand = SmallRng::from_seed(SEED);
    let poisson = Builder::<_, Vect>::with_samples(100, 0.8, Type::Perioditic)
        .build(rand, algorithm::Bridson);
    b.iter(|| black_box(poisson.generate()));
}
