use ::{math, PoissonDisk, Sample};

extern crate image;

extern crate rand;

use na::Vec2 as naVec2;
pub type Vec2 = naVec2<f64>;

use std::fs::File;
use std::path::Path;

#[test]
fn visualise_64th_radius() {
    let radius = 2f64.sqrt() / 2f64 / 64f64;
    let size = 2 << 9;
    let mut poisson = PoissonDisk::new(rand::weak_rng(), radius);
    let mut samples = vec![];
    poisson.create(&mut samples);
    visualise(&samples, size, false, &Path::new("poisson_64th.png"));
}

#[test]
fn visualise_32th_of_max_radius() {
    let radius = 2f64.sqrt() / 2f64 / 32f64;
    let size = 2 << 9;
    let mut poisson = PoissonDisk::new(rand::weak_rng(), radius);
    let mut samples = vec![];
    poisson.create(&mut samples);
    visualise(&samples, size, false, &Path::new("poisson_32th.png"));
}

#[test]
fn visualise_16th_of_max_radius() {
    let radius = 2f64.sqrt() / 2f64 / 16f64;
    let size = 2 << 9;
    let mut poisson = PoissonDisk::new(rand::weak_rng(), radius);
    let mut samples = vec![];
    poisson.create(&mut samples);
    visualise(&samples, size, false, &Path::new("poisson_16th.png"));
}

#[test]
fn visualise_8th_of_max_radius() {
    let radius = 2f64.sqrt() / 2f64 / 8f64;
    let size = 2 << 9;
    let mut poisson = PoissonDisk::new(rand::weak_rng(), radius);
    let mut samples = vec![];
    poisson.create(&mut samples);
    visualise(&samples, size, false, &Path::new("poisson_8th.png"));
}

#[test]
fn visualise_4th_of_max_radius() {
    let radius = 2f64.sqrt() / 2f64 / 4f64;
    let size = 2 << 9;
    let mut poisson = PoissonDisk::new(rand::weak_rng(), radius);
    let mut samples = vec![];
    poisson.create(&mut samples);
    visualise(&samples, size, false, &Path::new("poisson_4th.png"));
}

#[test]
fn visualise_2th_of_max_radius() {
    let radius = 2f64.sqrt() / 2f64 / 2f64;
    let size = 2 << 9;
    let mut poisson = PoissonDisk::new(rand::weak_rng(), radius);
    let mut samples = vec![];
    poisson.create(&mut samples);
    visualise(&samples, size, false, &Path::new("poisson_2th.png"));
}

#[test]
fn visualise_max_radius() {
    let radius = 2f64.sqrt() / 2f64;
    let size = 2 << 9;
    let mut poisson = PoissonDisk::new(rand::weak_rng(), radius);
    let mut samples = vec![];
    poisson.create(&mut samples);
    visualise(&samples, size, false, &Path::new("poisson.png"));
}

#[test]
fn visualise_64th_radius_periodic() {
    let radius = 2f64.sqrt() / 2f64 / 64f64;
    let size = 2 << 9;
    let mut poisson = PoissonDisk::perioditic(rand::weak_rng(), radius);
    let mut samples = vec![];
    poisson.create(&mut samples);
    visualise(&samples, size, true, &Path::new("periodic_poisson_64th.png"));
}

#[test]
fn visualise_32th_of_max_radius_periodic() {
    let radius = 2f64.sqrt() / 2f64 / 32f64;
    let size = 2 << 9;
    let mut poisson = PoissonDisk::perioditic(rand::weak_rng(), radius);
    let mut samples = vec![];
    poisson.create(&mut samples);
    visualise(&samples, size, true, &Path::new("periodic_poisson_32th.png"));
}

#[test]
fn visualise_16th_of_max_radius_periodic() {
    let radius = 2f64.sqrt() / 2f64 / 16f64;
    let size = 2 << 9;
    let mut poisson = PoissonDisk::perioditic(rand::weak_rng(), radius);
    let mut samples = vec![];
    poisson.create(&mut samples);
    visualise(&samples, size, true, &Path::new("periodic_poisson_16th.png"));
}

#[test]
fn visualise_8th_of_max_radius_periodic() {
    let radius = 2f64.sqrt() / 2f64 / 8f64;
    let size = 2 << 9;
    let mut poisson = PoissonDisk::perioditic(rand::weak_rng(), radius);
    let mut samples = vec![];
    poisson.create(&mut samples);
    visualise(&samples, size, true, &Path::new("periodic_poisson_8th.png"));
}

#[test]
fn visualise_4th_of_max_radius_periodic() {
    let radius = 2f64.sqrt() / 2f64 / 4f64;
    let size = 2 << 9;
    let mut poisson = PoissonDisk::perioditic(rand::weak_rng(), radius);
    let mut samples = vec![];
    poisson.create(&mut samples);
    visualise(&samples, size, true, &Path::new("periodic_poisson_4th.png"));
}

#[test]
fn visualise_2th_of_max_radius_periodic() {
    let radius = 2f64.sqrt() / 2f64 / 2f64;
    let size = 2 << 9;
    let mut poisson = PoissonDisk::perioditic(rand::weak_rng(), radius);
    let mut samples = vec![];
    poisson.create(&mut samples);
    visualise(&samples, size, true, &Path::new("periodic_poisson_2th.png"));
}

#[test]
fn visualise_max_radius_periodic() {
    let radius = 2f64.sqrt() / 2f64;
    let size = 2 << 9;
    let mut poisson = PoissonDisk::perioditic(rand::weak_rng(), radius);
    let mut samples = vec![];
    poisson.create(&mut samples);
    visualise(&samples, size, true, &Path::new("periodic_poisson.png"));
}

#[test]
fn visualise_64th_of_max_radius_prefilled_with_32th_of_max_radius() {
    let radius = 2f64.sqrt() / 2f64;
    let size = 2 << 9;
    let mut samples = vec![];

    let mut poisson = PoissonDisk::new(rand::weak_rng(), radius / 32f64);
    poisson.create(&mut samples);

    let mut poisson = PoissonDisk::new(rand::weak_rng(), radius / 64f64);
    poisson.create(&mut samples);

    visualise(&samples, size, false, &Path::new("prefilled_poisson_64th_32th.png"));
}

#[test]
fn visualise_32th_of_max_radius_prefilled_with_16th_of_max_radius() {
    let radius = 2f64.sqrt() / 2f64;
    let size = 2 << 9;
    let mut samples = vec![];

    let mut poisson = PoissonDisk::new(rand::weak_rng(), radius / 16f64);
    poisson.create(&mut samples);

    let mut poisson = PoissonDisk::new(rand::weak_rng(), radius / 32f64);
    poisson.create(&mut samples);

    visualise(&samples, size, false, &Path::new("prefilled_poisson_32th_16th.png"));
}

#[test]
fn visualise_16th_of_max_radius_prefilled_with_8th_of_max_radius() {
    let radius = 2f64.sqrt() / 2f64;
    let size = 2 << 9;
    let mut samples = vec![];

    let mut poisson = PoissonDisk::new(rand::weak_rng(), radius / 8f64);
    poisson.create(&mut samples);

    let mut poisson = PoissonDisk::new(rand::weak_rng(), radius / 16f64);
    poisson.create(&mut samples);

    visualise(&samples, size, false, &Path::new("prefilled_poisson_16th_8th.png"));
}

#[test]
fn visualise_8th_of_max_radius_prefilled_with_4th_of_max_radius() {
    let radius = 2f64.sqrt() / 2f64;
    let size = 2 << 9;
    let mut samples = vec![];

    let mut poisson = PoissonDisk::new(rand::weak_rng(), radius / 4f64);
    poisson.create(&mut samples);

    let mut poisson = PoissonDisk::new(rand::weak_rng(), radius / 8f64);
    poisson.create(&mut samples);

    visualise(&samples, size, false, &Path::new("prefilled_poisson_8th_4th.png"));
}

#[test]
fn visualise_4th_of_max_radius_prefilled_with_2th_of_max_radius() {
    let radius = 2f64.sqrt() / 2f64;
    let size = 2 << 9;
    let mut samples = vec![];

    let mut poisson = PoissonDisk::new(rand::weak_rng(), radius / 2f64);
    poisson.create(&mut samples);

    let mut poisson = PoissonDisk::new(rand::weak_rng(), radius / 4f64);
    poisson.create(&mut samples);

    visualise(&samples, size, false, &Path::new("prefilled_poisson_4th_2nd.png"));
}

#[test]
fn visualise_2th_of_max_radius_prefilled_with_max_radius() {
    let radius = 2f64.sqrt() / 2f64;
    let size = 2 << 9;
    let mut samples = vec![];

    let mut poisson = PoissonDisk::new(rand::weak_rng(), radius);
    poisson.create(&mut samples);

    let mut poisson = PoissonDisk::new(rand::weak_rng(), radius / 2f64);
    poisson.create(&mut samples);

    visualise(&samples, size, false, &Path::new("prefilled_poisson_2nd_max.png"));
}

#[test]
fn visualise_64th_32th_16th_8th_4th() {
    let radius = 2f64.sqrt() / 2f64;
    let size = 2 << 9;
    let mut samples = vec![];

    let mut poisson = PoissonDisk::new(rand::weak_rng(), radius / 4f64);
    poisson.create(&mut samples);

    let mut poisson = PoissonDisk::new(rand::weak_rng(), radius / 8f64);
    poisson.create(&mut samples);

    let mut poisson = PoissonDisk::new(rand::weak_rng(), radius / 16f64);
    poisson.create(&mut samples);

    let mut poisson = PoissonDisk::new(rand::weak_rng(), radius / 32f64);
    poisson.create(&mut samples);

    let mut poisson = PoissonDisk::new(rand::weak_rng(), radius / 64f64);
    poisson.create(&mut samples);

    visualise(&samples, size, false, &Path::new("prefilled_poisson_64th_32th_16th_8th_4th.png"));
}


#[test]
fn visualise_perioditic_64th_of_max_radius_prefilled_with_32th_of_max_radius() {
    let radius = 2f64.sqrt() / 2f64;
    let size = 2 << 9;
    let mut samples = vec![];

    let mut poisson = PoissonDisk::perioditic(rand::weak_rng(), radius / 32f64);
    poisson.create(&mut samples);

    let mut poisson = PoissonDisk::perioditic(rand::weak_rng(), radius / 64f64);
    poisson.create(&mut samples);

    visualise(&samples, size, true, &Path::new("prefilled_perioditic_poisson_64th_32th.png"));
}

#[test]
fn visualise_perioditic_32th_of_max_radius_prefilled_with_16th_of_max_radius() {
    let radius = 2f64.sqrt() / 2f64;
    let size = 2 << 9;
    let mut samples = vec![];

    let mut poisson = PoissonDisk::perioditic(rand::weak_rng(), radius / 16f64);
    poisson.create(&mut samples);

    let mut poisson = PoissonDisk::perioditic(rand::weak_rng(), radius / 32f64);
    poisson.create(&mut samples);

    visualise(&samples, size, true, &Path::new("prefilled_perioditic_poisson_32th_16th.png"));
}

#[test]
fn visualise_perioditic_16th_of_max_radius_prefilled_with_8th_of_max_radius() {
    let radius = 2f64.sqrt() / 2f64;
    let size = 2 << 9;
    let mut samples = vec![];

    let mut poisson = PoissonDisk::perioditic(rand::weak_rng(), radius / 8f64);
    poisson.create(&mut samples);

    let mut poisson = PoissonDisk::perioditic(rand::weak_rng(), radius / 16f64);
    poisson.create(&mut samples);

    visualise(&samples, size, true, &Path::new("prefilled_perioditic_poisson_16th_8th.png"));
}

#[test]
fn visualise_perioditic_8th_of_max_radius_prefilled_with_4th_of_max_radius() {
    let radius = 2f64.sqrt() / 2f64;
    let size = 2 << 9;
    let mut samples = vec![];

    let mut poisson = PoissonDisk::perioditic(rand::weak_rng(), radius / 4f64);
    poisson.create(&mut samples);

    let mut poisson = PoissonDisk::perioditic(rand::weak_rng(), radius / 8f64);
    poisson.create(&mut samples);

    visualise(&samples, size, true, &Path::new("prefilled_perioditic_poisson_8th_4th.png"));
}

#[test]
fn visualise_perioditic_4th_of_max_radius_prefilled_with_2th_of_max_radius() {
    let radius = 2f64.sqrt() / 2f64;
    let size = 2 << 9;
    let mut samples = vec![];

    let mut poisson = PoissonDisk::perioditic(rand::weak_rng(), radius / 2f64);
    poisson.create(&mut samples);

    let mut poisson = PoissonDisk::perioditic(rand::weak_rng(), radius / 4f64);
    poisson.create(&mut samples);

    visualise(&samples, size, true, &Path::new("prefilled_perioditic_poisson_4th_2nd.png"));
}

#[test]
fn visualise_perioditic_2th_of_max_radius_prefilled_with_max_radius() {
    let radius = 2f64.sqrt() / 2f64;
    let size = 2 << 9;
    let mut samples = vec![];

    let mut poisson = PoissonDisk::perioditic(rand::weak_rng(), radius);
    poisson.create(&mut samples);

    let mut poisson = PoissonDisk::perioditic(rand::weak_rng(), radius / 2f64);
    poisson.create(&mut samples);

    visualise(&samples, size, true, &Path::new("prefilled_perioditic_poisson_2nd_max.png"));
}

#[test]
fn visualise_perioditic_64th_32th_16th_8th_4th() {
    let radius = 2f64.sqrt() / 2f64;
    let size = 2 << 9;
    let mut samples = vec![];

    let mut poisson = PoissonDisk::perioditic(rand::weak_rng(), radius / 4f64);
    poisson.create(&mut samples);

    let mut poisson = PoissonDisk::perioditic(rand::weak_rng(), radius / 8f64);
    poisson.create(&mut samples);

    let mut poisson = PoissonDisk::perioditic(rand::weak_rng(), radius / 16f64);
    poisson.create(&mut samples);

    let mut poisson = PoissonDisk::perioditic(rand::weak_rng(), radius / 32f64);
    poisson.create(&mut samples);

    let mut poisson = PoissonDisk::perioditic(rand::weak_rng(), radius / 64f64);
    poisson.create(&mut samples);

    visualise(&samples, size, true, &Path::new("prefilled_perioditic_poisson_64th_32th_16th_8th_4th.png"));
}

fn visualise(samples: &Vec<Sample>, size: u32, periodicity: bool, path: &Path) {
    let mut imgbuf = image::ImageBuffer::new(size, size);

    let color = image::Rgb([255 as u8, 255 as u8, 255 as u8]);
    let middle = image::Rgb([255 as u8, 0 as u8, 0 as u8]);
        for sample in samples {
            let xx = (sample.pos.x * size as f64) as i32;
            let yy = (sample.pos.y * size as f64) as i32;
            let radius = (sample.get_radius() * size as f64) as i32;
            for x in -radius..radius {
                for y in -radius..radius {
                    if x * x + y * y < radius * radius {
                        let color = if x == 0 || y == 0 {
                            middle
                        } else {
                            color
                        };
                        let cur_x = if periodicity {
                            math::modulo(xx + x, size as i32)
                        } else {
                            let cur_x = xx + x;
                            if cur_x < 0 || cur_x >= size as i32 {
                                continue;
                            }
                            cur_x
                        };
                        let cur_y = if periodicity {
                            math::modulo(yy + y, size as i32)
                        } else {
                            let cur_y = yy + y;
                            if cur_y < 0 || cur_y >= size as i32 {
                                continue;
                            }
                            cur_y
                        };
                        imgbuf.put_pixel(cur_x as u32, cur_y as u32, color);
                    }
                }
            }
        }
    let ref mut fout = File::create(path).unwrap();
    let _ = image::ImageRgb8(imgbuf).save(fout, image::PNG);
}
