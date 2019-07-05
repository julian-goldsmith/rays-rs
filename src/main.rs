extern crate cgmath;
extern crate png;
extern crate rand;

use std::path::Path;
use std::fs::File;
use std::io::BufWriter;
use cgmath::prelude::*;
use cgmath::{Matrix2, Matrix4, Point3, Vector2, Vector3, Vector4};
use png::HasParameters;
use rand::Rng;

// https://nelari.us/post/raytracer_with_rust_and_zig/
// https://www.cl.cam.ac.uk/teaching/1999/AGraphHCI/SMAG/node2.html
// http://www.realtimerendering.com/raytracing/Ray%20Tracing%20in%20a%20Weekend.pdf

// P(t) = E + tD, t >= 0
#[derive(Copy, Clone)]
pub struct Ray {
    pub origin: Point3<f32>,
    pub direction: Vector4<f32>,
}

impl Ray {
    pub fn point_at_distance(&self, dist: f32) -> Point3<f32> {
        self.origin + self.direction.truncate() * dist
    }
}

#[derive(Copy, Clone)]
pub struct Intersection {
    pub pos: Point3<f32>,
    pub dist: f32,
    pub normal: Vector4<f32>,
}

pub trait Intersect {
    fn intersect(&self, r: &Ray) -> Option<Intersection>;
}

#[derive(Copy, Clone)]
pub struct Sphere {
    pub center: Point3<f32>,
    pub radius: f32,
}

impl Intersect for Sphere {
    fn intersect(&self, r: &Ray) -> Option<Intersection> {
        let r_proj = self.center.to_vec().project_on(r.origin.to_vec() + r.direction.truncate()).extend(1.0);
        let discriminant = r_proj.magnitude2() + self.radius * self.radius - self.center.distance2(r.origin);

        if discriminant > 0.0 {
            let dist = r_proj.magnitude() - discriminant.sqrt();
            let pos = r.point_at_distance(dist);
            
            Some(Intersection {
                pos,
                dist,
                normal: (pos - self.center).normalize().extend(1.0),
            })
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
        let num_samples = 4;
        let mut pixels = Vec::new();
        let mut rng = rand::thread_rng();

        let aspect = (height as f32) / (width as f32);
        let screen_space = Vector2::new(4.0, 4.0 * aspect);
        let lower_left_corner = (screen_space / -2.0).extend(-1.0);

        let screen_space_transform = Matrix4::new(
            screen_space.x / (width as f32), 0.0, 0.0, 0.0,
            0.0, screen_space.y / (height as f32), 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0,
            lower_left_corner.x, lower_left_corner.y, lower_left_corner.z, 1.0);

        for y in 0..height {
            for x in 0..width {
                let mut color = Vector3::zero();

                for _ in 0..num_samples {
                    let ofs_x = rng.gen_range(-0.5, 0.5);
                    let ofs_y = rng.gen_range(-0.5, 0.5);
                    let uv = screen_space_transform *
                        Vector4::new(x as f32 + ofs_x, y as f32 + ofs_y, 0.0, 1.0);
                    let direction = uv.truncate().normalize().extend(1.0);

                    let r = Ray {
                        origin: Point3::new(0.0, 0.0, 0.0),
                        direction,
                    };

                    let curr_ixn = self.spheres.iter().
                        filter_map(|s| s.intersect(&r)).
                        min_by(|a, b| a.dist.partial_cmp(&b.dist).unwrap());

                    color += match curr_ixn {
                        Some(ixn) => {
                            0.5 * (ixn.normal.truncate() + Vector3::new(1.0, 1.0, 1.0))
                        },
                        None => Vector3::new(0.1, 0.4, 0.1),
                    };
                };

                color *= 255.0 / (num_samples as f32);
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
    let width = 800;
    let height = 800;
    let path = Path::new("/var/www/rays.png");

    let world = World {
        spheres: vec![
            Sphere {
                center: Point3::new(0.0, 0.0, -2.0),
                radius: 0.5,
            },
            Sphere {
                center: Point3::new(1.0, 0.0, -1.5),
                radius: 0.5,
            },
            Sphere {
                center: Point3::new(-1.0, 0.0, -1.5),
                radius: 0.5,
            },
        ],
    };
    let data = world.render(width, height);

    write_png(&path, &data, width, height);
}
