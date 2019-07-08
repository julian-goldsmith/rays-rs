extern crate cgmath;
extern crate png;
extern crate rand;
extern crate rand_distr;

use std::path::Path;
use std::fs::File;
use std::io::BufWriter;
use cgmath::prelude::*;
use cgmath::{Matrix4, Point3, Vector2, Vector3, Vector4};
use png::HasParameters;
use rand::rngs::ThreadRng;
use rand_distr::{Distribution, Normal, UnitSphere};

// https://nelari.us/post/raytracer_with_rust_and_zig/
// https://www.cl.cam.ac.uk/teaching/1999/AGraphHCI/SMAG/node2.html
// http://www.realtimerendering.com/raytracing/Ray%20Tracing%20in%20a%20Weekend.pdf
// https://www.tutorialspoint.com/computer_graphics/3d_transformation.htm

// P(t) = E + tD, t >= 0
#[derive(Copy, Clone)]
pub struct Ray {
    pub origin: Point3<f32>,
    pub direction: Vector3<f32>,
}

impl Ray {
    pub fn point_at_distance(&self, dist: f32) -> Point3<f32> {
        self.origin + self.direction * dist
    }
}

#[derive(Copy, Clone)]
pub struct Intersection {
    pub pos: Point3<f32>,
    pub dist: f32,
    pub normal: Vector3<f32>,
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
        let r_proj = self.center.to_vec().project_on(r.origin.to_vec() + r.direction).extend(1.0);
        let discriminant = r_proj.magnitude2() + self.radius * self.radius - self.center.distance2(r.origin);

        if discriminant > 0.0 {
            let dist = r_proj.magnitude() - discriminant.sqrt();
            let pos = r.point_at_distance(dist);
            
            Some(Intersection {
                pos,
                dist,
                normal: (pos - self.center).normalize(),
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
        let ndist = Normal::new(0.0, 0.25).unwrap();                   // TODO: Adjust.

        let origin = Point3::new(0.0, 0.0, 0.0);

        let aspect = (height as f32) / (width as f32);
        let sst = Matrix4::new(
            1.0 / width as f32, 0.0, 0.0, 0.0,
            0.0, 1.0 / height as f32, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0,
            -0.5, -0.5, 0.0, 1.0);

        let persp = cgmath::frustum(aspect * -0.5, aspect * 0.5, -0.5, 0.5, 0.1, 10.0);
        let view = sst * Matrix4::look_at(origin, Point3::new(0.0, 0.0, -1.0), Vector3::new(0.0, 1.0, 0.0));

        for y in (0..height).rev() {
            for x in 0..width {
                let uv = Vector2::new(x as f32, y as f32);

                let color: Vector3<f32> = (0..num_samples).
                    fold(Vector3::zero(),
                         |acc, _| acc + self.sample(&mut rng, &ndist, origin, uv, &persp, &view)) *
                    255.0 / (num_samples as f32);

                pixels.push(color.x as u8);
                pixels.push(color.y as u8);
                pixels.push(color.z as u8);
            };
        };

        pixels
    }

    fn sample(&self, rng: &mut ThreadRng, ndist: &Normal<f32>, origin: Point3<f32>, uv: Vector2<f32>,
              persp: &Matrix4<f32>, view: &Matrix4<f32>) -> Vector3<f32> {
        let jitter = Vector2::new(ndist.sample(rng), ndist.sample(rng));
        let sample_uv = Vector4::new(uv.x + jitter.x, uv.y + jitter.y, 0.0, 1.0);
        let direction = (persp * view * sample_uv).truncate().normalize();

        let r = Ray { origin, direction };

        let curr_ixn = self.spheres.iter().
            filter_map(|s| s.intersect(&r)).
            min_by(|a, b| a.dist.partial_cmp(&b.dist).unwrap());

        match curr_ixn {
            Some(ixn) => {
                let bounce: Vector3<f32> = UnitSphere.sample(rng).into();
                let target = ixn.pos + ixn.normal + bounce;

                0.5 * ((ixn.pos - target).normalize() + Vector3::new(1.0, 1.0, 1.0))
            },
            None => Vector3::new(0.1, 0.4, 0.1),
        }
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
    let width = 1920;
    let height = 1080;
    let path = Path::new("/var/www/rays.png");

    let world = World {
        spheres: vec![
            Sphere {
                center: Point3::new(0.0, 0.0, -5.0),
                radius: 0.5,
            },
            Sphere {
                center: Point3::new(0.0, -100.5, -105.0),
                radius: 100.0,
            },
        ],
    };
    let data = world.render(width, height);

    write_png(&path, &data, width, height);
}
