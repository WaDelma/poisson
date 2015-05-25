use rand::Rng;

use na::Vec2 as naVec2;
pub type Vec2 = naVec2<f64>;

use std::fmt::{Debug, Formatter};
use std::fmt::Result as FmtResult;

#[derive(Clone, Copy)]
pub struct Rect {
    pub min: Vec2,
    pub max: Vec2,
}

impl Rect {
    #[inline]
    pub fn new(min: Vec2, max: Vec2) -> Self {
        Rect{min: min, max: max}
    }

    #[inline]
    pub fn area(&self) -> f64 {
        self.width() * self.height()
    }

    #[inline]
    pub fn width(&self) -> f64 {
        self.max.x - self.min.x
    }

    #[inline]
    pub fn height(&self) -> f64 {
        self.max.y - self.min.y
    }

    #[inline]
    pub fn center_x(&self) -> f64 {
        (self.max.x + self.min.x) / 2f64
    }

    #[inline]
    pub fn center_y(&self) -> f64 {
        (self.max.y + self.min.y) / 2f64
    }

    #[inline]
    pub fn random_point_inside<F>(&self, rand: &mut F) -> Vec2 where F: Rng {
        Vec2::new(
            self.min.x + rand.gen::<f64>() * self.width(),
            self.min.y + rand.gen::<f64>() * self.height()
        )
    }
}

impl Debug for Rect {
    #[inline]
    fn fmt(&self, f: &mut Formatter) -> FmtResult {
        write!(f, "Rect({:?}, {:?})", self.min, self.max)
    }
}

impl PartialEq for Rect {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.min == other.min && self.max == other.max
    }
}

pub enum Intersection {
    Out,
    Over,
    In,
}

#[inline]
pub fn test_intersection(rect: Rect, sample: Vec2, radius: f64) -> Intersection {
    let radius2 = radius * radius;
    let dims = 2usize;
    let mut min = 0f64;
    let mut max = 0f64;
    for i in 0usize..dims {
        let cur_min = rect.min[i];
        let cur_max = rect.max[i];
        let cur_sample = sample[i];
        let dmin = cur_min - cur_sample;
        if dmin > 0f64 {
            if dmin > radius {
                return Intersection::Out;
            }
            min += dmin * dmin;
            if min > radius2 {
                return Intersection::Out;
            }
            let temp_max = cur_max - cur_sample;
            max += temp_max * temp_max;
            continue;
        }
        let dmin = cur_sample - cur_max;
        if dmin > 0f64 {
            if dmin > radius {
                return Intersection::Out;
            }
            min += dmin * dmin;
            if min > radius2 {
                return Intersection::Out;
            }
            let temp_max = cur_sample - cur_min;
            max += temp_max * temp_max;
            continue;
        }
        let temp_max = (cur_sample - cur_min).max(cur_max - cur_sample);
        max += temp_max * temp_max;
    }
    if max > radius2 {
        return Intersection::Over;
    }
    return Intersection::In;
}
