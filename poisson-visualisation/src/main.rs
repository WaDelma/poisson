use clap::{App, Arg, ArgMatches, arg_enum, _clap_count_exprs, value_t};

use poisson::{Builder, Type, algorithm::{Bridson, Ebeida}};

use rand::{rngs::SmallRng, FromEntropy, Rng, seq::SliceRandom, SeedableRng};

use nalgebra::Vector2;

use image::{ImageBuffer, Rgb};

use lab::Lab;

use fnv::FnvHasher;

use std::hash::Hasher;

arg_enum! {
    #[derive(PartialEq, Debug)]
    pub enum Algo {
        Ebeida,
        Bridson
    }
}

arg_enum! {
    #[derive(PartialEq, Debug)]
    pub enum Style {
        Plain,
        Colorful,
        Dot
    }
}

fn main() {
    let app = App::with_defaults("Poisson visualisation")
        .about("Visualisation for poisson library")
        .arg(
            Arg::with_name("OUTPUT")
                .help("Output file that's generated")
                .required(true)
                .index(1)
        )
        .arg(
            Arg::with_name("SEED")
                .help("Seed for the generation")
                .index(2)
        )
        .arg(
            Arg::with_name("radius")
                .short("r")
                .takes_value(true)
                .help("Radius of the disks")
        )
        .arg(
            Arg::with_name("width")
                .short("w")
                .takes_value(true)
                .help("Width of the generated image")
        )
        .arg(
            Arg::with_name("height")
                .short("h")
                .takes_value(true)
                .help("Height of the generated image")
        )
        .arg(
            Arg::with_name("style")
                .short("s")
                .takes_value(true)
                .help("Style for the disks")
                .possible_values(&Style::variants())
        )
        .arg(
            Arg::with_name("algo")
                .short("a")
                .help("Algorithm that's used to generate image")
                .takes_value(true)
                .possible_values(&Algo::variants())
        );
    visualise(app.get_matches());
}

fn visualise(m: ArgMatches) {
    let width = value_t!(m, "width", u32).unwrap_or(1024);
    let height = value_t!(m, "height", u32).unwrap_or(1024);
    let radius = value_t!(m, "radius", f32).unwrap_or(0.02);
    let algo = value_t!(m, "algo", Algo).unwrap_or(Algo::Ebeida);
    let style = value_t!(m, "style", Style).unwrap_or(Style::Plain);
    let name = m.value_of("OUTPUT").unwrap();
    let master_rng = m.value_of("SEED").map(|s| {
        let mut fnv = FnvHasher::with_key(0);
        for b in s.bytes() {
            fnv.write_u8(b);
        }
        SmallRng::seed_from_u64(fnv.finish())
    }).unwrap_or_else(SmallRng::from_entropy);

    let mut style_rng = master_rng.clone();

    let builder = Builder::<_, Vector2<f32>>::with_radius(radius, Type::Normal);
    let points = if algo == Algo::Ebeida {
        builder.build(master_rng, Ebeida).generate()
    } else {
        builder.build(master_rng, Bridson).generate()
    };

    let mut ps = points.clone();
    ps.shuffle(&mut style_rng);

    let mut image = ImageBuffer::new(width, height);
    for p in points {
        let pp = ps.pop().unwrap();
        let col = Rgb {
            data: Lab {
                l: style_rng.gen::<f32>() * 80. + 10.,
                a: pp.x * 256. - 128.,
                b: pp.y * 256. - 128.
            }.to_rgb()
        };

        let x = p.x * width as f32;
        let y = p.y * height as f32;
        let (rx, ry) = if style == Style::Dot {
            (0.2 * radius * width as f32, 0.2 * radius * height as f32)
        } else {
            (radius * width as f32, radius * height as f32)
        };
        for xx in -rx as i32..rx as i32 {
            for yy in -ry as i32..ry as i32 {
                let xx = xx as f32;
                let yy = yy as f32;
                let xxx = (x + xx) as i32;
                let yyy = height as i32 - (y + yy) as i32;
                if xxx < 0 || xxx >= width as i32 {
                    // Outside of the picture horizontally
                    continue;
                }
                if yyy < 0 || yyy >= height as i32 {
                    // Outside of the picture vertically
                    continue;
                }
                if xx * xx / (rx * rx) + yy * yy / (ry * ry) > 1. {
                    // Outside of the disk
                    continue;
                }
                let xxx = xxx as u32;
                let yyy = yyy as u32;
                if style == Style::Colorful {
                    image[(xxx, yyy)] = col;
                } else {
                    image[(xxx, yyy)] = Rgb { data: [255, 255, 255] };
                }
                if style == Style::Plain && (xx == 0. || yy == 0.) {
                    image[(xxx, yyy)] = Rgb { data: [255, 0, 0] };
                }
            }
        }
    }
    image.save(name).unwrap();
}
