use std::env;
use std::hint::black_box;
use std::time::Instant;
#[inline(never)]
fn dot(a: &[f64], b: &[f64]) -> f64 {
    let mut sum = 0.0;
    for i in 0..a.len() { sum += a[i] * b[i]; }
    sum
}
fn main() {
    let args: Vec<usize> = env::args().skip(1).map(|x| x.parse().unwrap()).collect();
    assert_eq!(args.len(), 3);
    let (n, iterations, seed) = (args[0], args[1], args[2]);
    assert!(n <= 1048576 && iterations > 0 && iterations <= 100000000 && seed <= 1000000);
    let count = n.max(1) * 32;
    let a: Vec<f64> = (0..count).map(|i| (((i * 17 + seed * 13) % 101) as i32 - 50) as f64 / 16.0).collect();
    let b: Vec<f64> = (0..count).map(|i| (((i * 29 + (seed + 1) * 7) % 103) as i32 - 51) as f64 / 16.0).collect();
    for i in 0..iterations.min(128).max(8) { let off = (i % 32) * n; black_box(dot(&a[off..off+n], &b[off..off+n])); }
    let mut checksum = 0.0;
    let start = Instant::now();
    for i in 0..iterations { let off = (i % 32) * n; checksum += dot(&a[off..off+n], &b[off..off+n]); }
    let elapsed = start.elapsed().as_nanos();
    assert!(checksum.is_finite());
    println!("{{\"elapsed_ns\":{},\"iterations\":{},\"checksum\":{}}}", elapsed, iterations, checksum);
}
