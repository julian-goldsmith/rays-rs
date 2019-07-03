extern crate cgmath;
extern crate png;

use std::path::Path;
use std::fs::File;
use std::io::BufWriter;
use cgmath::prelude::*;
use cgmath::{Matrix2, Vector2, Vector3};
use png::HasParameters;

// https://nelari.us/post/raytracer_with_rust_and_zig/
// https://www.cl.cam.ac.uk/teaching/1999/AGraphHCI/SMAG/node2.html
// http://www.realtimerendering.com/raytracing/Ray%20Tracing%20in%20a%20Weekend.pdf

// P(t) = E + tD, t >= 0
#[derive(Copy, Clone)]
pub struct Ray {
    pub origin: Vector3<f32>,
    pub direction: Vector3<f32>,
}

impl Ray {
    pub fn point_at_distance(&self, dist: f32) -> Vector3<f32> {
        self.origin + self.direction * dist
    }
}

#[derive(Copy, Clone)]
pub struct Intersection {
    pub pos: Vector3<f32>,
    pub dist: f32,
    pub normal: Vector3<f32>,
}

pub trait Intersect {
    fn intersect(&self, r: &Ray) -> Option<Intersection>;
}

#[derive(Copy, Clone)]
pub struct Sphere {
    pub center: Vector3<f32>,
    pub radius: f32,
}

impl Intersect for Sphere {
    fn intersect(&self, r: &Ray) -> Option<Intersection> {
        let r_origin_local = r.origin - self.center;
        let a = r.direction.dot(r.direction);
        let b = 2.0 * r_origin_local.dot(r.direction);
        let c = r_origin_local.dot(r_origin_local) - self.radius * self.radius;
        let discriminant = b * b - 4.0 * a * c;

        if discriminant > 0.0 {
            let dist = if a < 0.0 {
                -b + discriminant.sqrt()
            } else {
                -b - discriminant.sqrt()
            } / (2.0 * a);

            let pos = r.point_at_distance(dist);
            let intersection = Intersection {
                pos,
                dist,
                normal: (pos - self.center) / self.radius,
            };

            Some(intersection)
        } else {
            None
        }
    }
}

pub struct World {
    pub spheres: Vec<Sphere>,
}

impl World {
    fn render(&self, width: usize, height: usize) -> Vec<u8> {
        let mut pixels = Vec::new();
        let aspect = (height as f32) / (width as f32);
        let origin = Vector3::zero();
        let screen_space = Vector2::new(4.0, 4.0 * aspect);
        let lower_left_corner = (screen_space / -2.0).extend(-1.0);     // FIXME: Go from top left.

        let screen_space_transform = Matrix2::new(
            screen_space.x / (width as f32), 0.0,
            0.0, screen_space.y / (height as f32));

        for y in 0..height {
            for x in 0..width {
                let uv = screen_space_transform * Vector2::new(x as f32, y as f32);

                let r = Ray {
                    origin,
                    direction: lower_left_corner + uv.extend(0.0),
                };

                let curr_ixn = self.spheres.iter().
                    map(|s| s.intersect(&r)).
                    filter(|ixn| ixn.is_some()).
                    map(|ixn| ixn.unwrap()).
                    min_by(|a, b| a.dist.partial_cmp(&b.dist).unwrap());

                let color = 255.0 * match curr_ixn {
                    Some(ixn) => {
                        0.5 * (ixn.normal + Vector3::new(1.0, 1.0, 1.0))
                    },
                    None => Vector3::new(0.1, 0.4, 0.1),
                };

                pixels.push(color.x as u8);
                pixels.push(color.y as u8);
                pixels.push(color.z as u8);
            };
        };

        pixels
    }
}

fn write_png(path: &Path, data: &[u8], width: usize, height: usize) {
    let file = File::create(path).unwrap();
    let ref mut w = BufWriter::new(file);

    let mut encoder = png::Encoder::new(w, width as u32, height as u32);
    encoder.set(png::ColorType::RGB).set(png::BitDepth::Eight);

    let mut writer = encoder.write_header().unwrap();
    writer.write_image_data(data).unwrap();
}

fn main() {
    let width = 900;
    let height = 700;
    let path = Path::new("/var/www/rays.png");

    let world = World {
        spheres: vec![
            Sphere {
                center: Vector3::new(0.0, 0.0, -2.0),
                radius: 0.5,
            },
            Sphere {
                center: Vector3::new(1.0, 0.0, -1.5),
                radius: 0.5,
            },
            Sphere {
                center: Vector3::new(-1.0, 0.0, -1.5),
                radius: 0.5,
            },
        ],
    };
    let data = world.render(width, height);

    write_png(&path, &data, width, height);
}
