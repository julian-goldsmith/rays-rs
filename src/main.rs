extern crate cgmath;
extern crate png;
extern crate rand;
extern crate rand_distr;
extern crate stl;

use std::borrow::Borrow;
use std::fs::File;
use std::io::BufWriter;
use std::path::Path;
use std::sync::Arc;
use std::thread;
use cgmath::prelude::*;
use cgmath::{Matrix3, Matrix4, Point2, Point3, Vector2, Vector3};
use cgmath::num_traits::NumCast;
use png::HasParameters;
use rand_distr::{Distribution, UnitSphere};

// https://nelari.us/post/raytracer_with_rust_and_zig/
// https://www.cl.cam.ac.uk/teaching/1999/AGraphHCI/SMAG/node2.html
// http://www.realtimerendering.com/raytracing/Ray%20Tracing%20in%20a%20Weekend.pdf
// https://www.tutorialspoint.com/computer_graphics/3d_transformation.htm

// P(t) = E + tD, t >= 0
#[derive(Copy, Clone, Debug)]
pub struct Ray {
    pub origin: Point3<f32>,
    pub direction: Vector3<f32>,
}

impl Ray {
    pub fn new(origin: Point3<f32>, direction: Vector3<f32>) -> Ray {
        Ray {
            origin,
            direction: direction.normalize(),
        }
    }

    pub fn point_at_distance(&self, dist: f32) -> Point3<f32> {
        self.origin + self.direction * dist
    }
}

#[derive(Copy, Clone, Debug)]
pub struct Intersection {
    pub pos: Point3<f32>,
    pub dist: f32,
    pub normal: Vector3<f32>,
    pub color: Vector3<f32>,
}

pub trait Intersect {
    fn intersect(&self, r: &Ray) -> Option<Intersection>;
}

#[derive(Copy, Clone, Debug)]
pub struct Triangle {
    pub v0: Point3<f32>,
    pub v1: Point3<f32>,
    pub v2: Point3<f32>,
    pub edge0: Vector3<f32>,
    pub edge1: Vector3<f32>,
    pub edge2: Vector3<f32>,
    pub normal: Vector3<f32>,
    pub color: Vector3<f32>,
}

impl Triangle {
    pub fn new(v0: Point3<f32>, v1: Point3<f32>, v2: Point3<f32>, normal: Vector3<f32>, color: Vector3<f32>) -> Triangle {
        Triangle {
            v0,
            v1,
            v2,
            edge0: v1 - v0,
            edge1: v2 - v0,
            edge2: v2 - v1,
            normal,
            color,
        }
    }
}

impl Intersect for Triangle {
    // http://www.lighthouse3d.com/tutorials/maths/ray-triangle-intersection/
    fn intersect(&self, r: &Ray) -> Option<Intersection> {
        let h = r.direction.cross(self.edge1);
        let a = self.edge0.dot(h);

        let epsilon = 0.0000001;
        if a > -epsilon && a < epsilon {
            // This ray is parallel to the triangle.
            return None;
        }

        let f = 1.0 / a;
        let s = r.origin - self.v0;
        let u = f * s.dot(h);

        if u < 0.0 || u > 1.0 {
            return None;
        }

        let q = s.cross(self.edge0);
        let v = f * r.direction.dot(q);

        if v < 0.0 || u + v > 1.0 {
            return None;
        }

        let dist = f * self.edge1.dot(q);
        if dist > epsilon {
            // Ray intersection
            Some(Intersection {
                pos: r.origin + r.direction * dist,
                dist,
                normal: self.normal,
                color: self.color,
            })
        } else {
            // Line intersection only
            None
        }
    }
}

#[derive(Copy, Clone, Debug)]
pub struct Sphere {
    pub center: Point3<f32>,
    pub radius: f32,
    pub color: Vector3<f32>,
}

impl Intersect for Sphere {
    fn intersect(&self, r: &Ray) -> Option<Intersection> {
        let oc = r.origin - self.center;
        let b = -oc.dot(r.direction);
        let discriminant = b * b + self.radius * self.radius - oc.magnitude2();

        if discriminant > 0.0 {
            let dist = if b < discriminant.sqrt() {
                b + discriminant.sqrt()
            } else {
                b - discriminant.sqrt()
            };

            if dist > 0.0001 {
                let pos = r.point_at_distance(dist);

                Some(Intersection {
                    pos,
                    dist,
                    normal: (pos - self.center).normalize(),
                    color: self.color,
                })
            } else {
                None
            }
        } else {
            None
        }
    }
}

#[derive(Clone, Debug)]
pub struct World {
    pub origin: Point3<f32>,
    pub look_at: Point3<f32>,

    pub spheres: Vec<Sphere>,
    pub triangles: Vec<Triangle>,
}

impl World {
    fn sample(&self, r: Ray, depth: usize) -> Vector3<f32> {
        let mut rng = rand::thread_rng();

        if depth > 4 {
            return Vector3::new(1.0, 0.0, 0.0);
        };

        let curr_ixn = self.spheres.iter().map(|s| s as &dyn Intersect).
            chain(self.triangles.iter().map(|t| t as &dyn Intersect)).
            filter_map(|i| i.intersect(&r)).
            min_by(|a, b| a.dist.partial_cmp(&b.dist).unwrap());

        match curr_ixn {
            Some(ixn) => {
                let num_bounces = 4;
                let mut color = Vector3::zero();

                for _ in 0..num_bounces {
                    let jitter = (Vector3::from(UnitSphere.sample(&mut rng)) - Vector3::new(0.5, 0.5, 0.5)) * 0.5;
                    let tr = Ray::new(ixn.pos, ixn.normal + jitter);

                    let bounce_sample = self.sample(tr, depth + 1);
                    color += ixn.color.mul_element_wise(bounce_sample);
                };

                color / num_bounces as f32
            },
            _ => {
                Vector3::new(0.6, 0.6, 0.6)
            },
        }
    }
}

pub struct Tile {
    pub canvas_base: Point2<f32>,
    pub pos: Vector2<f32>,
    pub pixels: [[Vector3<f32>; Tile::SIZE]; Tile::SIZE],
}

impl Tile {
    pub const SIZE: usize = 16;

    pub fn new(x: usize, y: usize) -> Tile {
        Tile {
            pos: Vector2::new(x as f32, y as f32),
            pixels: [[Vector3::zero(); Tile::SIZE]; Tile::SIZE],
            canvas_base: Point2::new((x * Tile::SIZE) as f32, (y * Tile::SIZE) as f32),
        }
    }

    pub fn get_canvas_pos<TI: Copy + NumCast, TO: NumCast>(&self, pos: Vector2<TI>) -> Point2<TO> {
        (self.canvas_base + pos.cast().unwrap()).cast().unwrap()
    }
}

pub struct RenderCanvas {
    pub size: Vector2<usize>,
    pub pixels: Box<[u8]>,
}

impl RenderCanvas {
    pub fn new(width: usize, height: usize) -> RenderCanvas {
        RenderCanvas {
            size: Vector2::new(width, height),
            pixels: vec![0; 3 * width * height].into_boxed_slice(),
        }
    }

    #[inline]
    pub fn set_pixel(&mut self, x: usize, y: usize, color: Vector3<f32>) {
        let value = Vector3::new(color.x.sqrt(), color.y.sqrt(), color.z.sqrt());                               // Apply gamma.
        let pixel_base = 3 * (y * self.size.x + x);
        self.pixels[pixel_base + 0] = (value.x * 255.999) as u8;
        self.pixels[pixel_base + 1] = (value.y * 255.999) as u8;
        self.pixels[pixel_base + 2] = (value.z * 255.999) as u8;
    }

    pub fn fill_tile(&mut self, tile: &Tile) {
        for y in 0..Tile::SIZE {
            let row = tile.pixels[y];

            for x in 0..Tile::SIZE {
                let tile_pos = Vector2::new(x, y);
                let canvas_pos = tile.get_canvas_pos(tile_pos);
                let canvas_y = self.size.y - canvas_pos.y;                                                      // Flip upside-down.

                // We need this in case the canvas width or height isn't divisible by the tile size.
                if canvas_pos.x >= self.size.x || canvas_y >= self.size.y {
                    continue;
                };

                self.set_pixel(canvas_pos.x, canvas_y, row[x]);
            };
        };
    }

    fn write_png(&self, path: &Path) {
        let file = File::create(path).unwrap();
        let ref mut w = BufWriter::new(file);

        let size = self.size;
        let mut encoder = png::Encoder::new(w, size.x as u32, size.y as u32);
        encoder.set(png::ColorType::RGB).set(png::BitDepth::Eight);

        let mut writer = encoder.write_header().unwrap();
        writer.write_image_data(&self.pixels).unwrap();
    }
}

pub struct Renderer {
    pub canvas: RenderCanvas,
    pub aspect: f32,
    pub num_samples: usize,

    pub world: Arc<World>,

    view: Matrix4<f32>,
    pv: Matrix4<f32>,
}

impl Renderer {
    fn new(width: usize, height: usize, num_samples: usize, world: World) -> Renderer {
        let aspect = (width as f32) / (height as f32);

        let persp = cgmath::frustum(-0.5, 0.5, aspect * -0.5, aspect * 0.5, 0.1, 100.0);
        let view = Matrix4::look_at(world.origin, world.look_at, Vector3::new(0.0, 1.0, 0.0));

        Renderer {
            canvas: RenderCanvas::new(width, height),
            aspect,
            num_samples,
            world: Arc::new(world),
            view,
            pv: persp * view,
        }
    }

    fn render(&mut self, path: &Path) {
        let origin = self.view.transform_point(Point3::origin());
        let canvas_size = self.canvas.size.cast::<f32>().unwrap();
        let pv = self.pv;
        let num_samples = self.num_samples;
        let size = self.canvas.size;
        let num_tiles_x = (size.x + Tile::SIZE - 1) / Tile::SIZE;
        let num_tiles_y = (size.y + Tile::SIZE - 1) / Tile::SIZE;
        let mut handles = Vec::new();

        for tile_y in 0..num_tiles_y {
            for tile_x in 0..num_tiles_x {
                let world = self.world.clone();
                let handle = thread::spawn(move || {
                    Renderer::render_tile(world.clone(), origin, canvas_size, pv,
                                          num_samples, tile_x, tile_y)
                });
                handles.push(handle);
            };
        };

        for handle in handles {
            let tile = handle.join().unwrap();
            self.canvas.fill_tile(&tile);
        };

        self.canvas.write_png(&path);
    }

    fn render_tile(world: Arc<World>, origin: Point3<f32>, canvas_size: Vector2<f32>, pv: Matrix4<f32>,
                   num_samples: usize, tile_x: usize, tile_y: usize) -> Tile {
        let mut tile = Tile::new(tile_x, tile_y);
        let inv_size = 1.0 / canvas_size;
        let world: &World = world.borrow();

        let uvt = Matrix3::new(
            inv_size.x, 0.0, 0.0,
            0.0, inv_size.y, 0.0,
            tile.canvas_base.x * inv_size.x - 0.5,
            tile.canvas_base.y * inv_size.y - 0.5, 1.0);

        for y in 0..Tile::SIZE {
            for x in 0..Tile::SIZE {
                let tile_pos = Point2::new(x, y);
                let uv = tile_pos.cast().unwrap();
                let color = &mut tile.pixels[y][x];

                for _ in 0..num_samples {
                    let jitter = Vector2::new(rand::random::<f32>() - 0.5, rand::random::<f32>() - 0.5);
                    let sample_uv = uvt.transform_point(uv + jitter) - Point2::origin();
                    let direction = pv.transform_vector(sample_uv.extend(0.1));                                 // TODO: Don't hardcode near frustum distance.
                    let r = Ray::new(origin, direction);

                    *color += world.sample(r, 0);
                };

                *color /= num_samples as f32;
            };
        };

        tile
    }
}

fn main() {
    let mut stl_file = File::open("monkey.stl").expect("Couldn't open STL file");
    let stl = stl::read_stl(&mut stl_file).expect("Couldn't read STL file.");
    let triangles = stl.triangles.iter().
        map(|t| Triangle::new(t.v1.into(), t.v2.into(), t.v3.into(), t.normal.into(), Vector3::new(0.7, 0.3, 0.6))).
        collect::<Vec<Triangle>>();

    let world = World {
        origin: Point3::new(0.0, 0.0, 3.0),                                                                     // FIXME: Origin doesn't work properly (currently, it's flipped in Z).
        look_at: Point3::new(0.0, 0.0, -5.0),

        spheres: vec![
            Sphere {
                center: Point3::new(0.0, -100.5, -1.0),
                radius: 100.0,
                color: Vector3::new(0.1, 0.1, 0.9),
            },
        ],
        triangles,
    };

    let mut renderer = Renderer::new(800, 450, 1, world);
    renderer.render(&Path::new("rays.png"));
}
