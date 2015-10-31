extern crate image;
use {Grid, VecLike};

use modulo::Mod;

use std::fs::{self, File};
use std::path::{PathBuf};
use std::ops::{Deref, DerefMut};

pub fn print_v<V: VecLike>(v: V) -> String {
    let mut result = "(".to_owned();
    for i in 0..V::dim(None) {
        result.push_str(&format!("{}, ", v[i]));
    }
    if V::dim(None) != 0 {
        result.pop();
    }
    result.push(')');
    result
}

pub fn visualise<V: VecLike>(level: usize,
                             grid: &Grid<V>,
                             indices: &Vec<V>,
                             r: f64,
                             periodicity: bool) {
    let size = 2u32.pow(9);//512;
    let samples: Vec<V> = grid.data.iter().filter_map(|v| *v).collect();
    let mut imgbuf = image::ImageBuffer::new(size, size);

    let color = image::Rgb([255 as u8, 255 as u8, 255 as u8]);
    let half = image::Rgb([127 as u8, 127 as u8, 127 as u8]);
    let green = image::Rgb([0 as u8, 255 as u8, 0 as u8]);
    let middle = image::Rgb([255 as u8, 0 as u8, 0 as u8]);
    let cells_per_cell = 2usize.pow(level as u32);

    let grid_w = (grid.cell / cells_per_cell as f64) * size as f64;
    for sample in indices {
        let x_start = (sample[0] * grid_w) as i32;
        let y_start = (sample[1] * grid_w) as i32;
        let x_end = ((sample[0] + 1.) * grid_w) as i32;
        let y_end = ((sample[1] + 1.) * grid_w) as i32;
        for x in x_start..x_end {
            for y in y_start..y_end {
                draw_pixel(&mut imgbuf, size, x, y, half, periodicity);
            }
        }
    }

    let grid_w = (grid.cell * size as f64) as u32;
    for x in 0..size {
        for y in 0..size {
            if x % grid_w == 0 || y % grid_w == 0 {
                imgbuf.put_pixel(x as u32, y as u32, green);
            }
        }
    }
    let r = 0.5 * r;
    for sample in samples {
        let xx = (sample[0] * size as f64) as i32;
        let yy = (sample[1] * size as f64) as i32;
        let radius = (r * size as f64) as i32;
        draw_pixel(&mut imgbuf, size, xx, yy, middle, periodicity);
        draw_circle(&mut imgbuf, size, xx, yy, radius, color, periodicity);
    }
    let mut p = PathBuf::new();
    p.push("visualise");
    let _ = fs::create_dir(p.clone());
    p.push(format!("{}", level));
    let ref mut fout = File::create(p.with_extension("png")).unwrap();
    let _ = image::ImageRgb8(imgbuf).save(fout, image::PNG);
}

fn draw_circle<C>(imgbuf: &mut image::ImageBuffer<image::Rgb<u8>, C>,
                  size: u32,
                  x0: i32,
                  y0: i32,
                  radius: i32,
                  color: image::Rgb<u8>,
                  periodicity: bool)
    where C: Deref<Target = [u8]> + DerefMut
{
    let mut x = radius;
    let mut y = 0;
    let mut decision_over_2 = 1 - x;

    while y <= x {
        draw_pixel(imgbuf, size, x + x0, y + y0, color, periodicity);
        draw_pixel(imgbuf, size, y + x0, x + y0, color, periodicity);
        draw_pixel(imgbuf, size, -x + x0, y + y0, color, periodicity);
        draw_pixel(imgbuf, size, -y + x0, x + y0, color, periodicity);
        draw_pixel(imgbuf, size, -x + x0, -y + y0, color, periodicity);
        draw_pixel(imgbuf, size, -y + x0, -x + y0, color, periodicity);
        draw_pixel(imgbuf, size, x + x0, -y + y0, color, periodicity);
        draw_pixel(imgbuf, size, y + x0, -x + y0, color, periodicity);
        y += 1;
        if decision_over_2 <= 0 {
            decision_over_2 += 2 * y + 1;
        } else {
            x -= 1;
            decision_over_2 += 2 * (y - x) + 1;
        }
    }
}

fn draw_pixel<C>(imgbuf: &mut image::ImageBuffer<image::Rgb<u8>, C>,
                 size: u32,
                 mut x: i32,
                 mut y: i32,
                 color: image::Rgb<u8>,
                 periodicity: bool)
    where C: Deref<Target = [u8]> + DerefMut
{
    if periodicity {
        x = x.modulo(size as i32);
        y = y.modulo(size as i32);
    } else {
        if x < 0 || x >= size as i32 {
            return;
        }
        if y < 0 || y >= size as i32 {
            return;
        }
    }
    imgbuf.put_pixel(x as u32, y as u32, color);
}
