//! RiemannRho: Library for approximating nontrivial zeros of the Riemann zeta function.
//!
//! This module provides functions to compute Hardy's Z-function using the Riemann-Siegel formula,
//! find zeros via bisection, and estimate the nth zero. It supports optional higher-order terms
//! for precision (High-Order Correction Mode). No visualization or I/O is included here; that's for binaries.

use std::f64::consts::PI;

/// Computes the Riemann-Siegel theta function approximation.
///
/// # Arguments
/// * `t` - The imaginary part (must be positive).
///
/// # Returns
/// The approximated theta(t).
pub fn theta(t: f64) -> f64 {
    let log_term = (t / 2.0) * (t / (2.0 * PI)).ln();
    log_term - t / 2.0 - PI / 8.0 + 1.0 / (48.0 * t) + 7.0 / (5760.0 * t.powi(3)) - 31.0 / (80640.0 * t.powi(5))
}

/// Numerical derivative helper for higher-order terms.
///
/// # Arguments
/// * `f` - The function to differentiate.
/// * `x` - Point of evaluation.
/// * `n` - Order of derivative.
/// * `h` - Step size for finite differences.
///
/// # Returns
/// The nth derivative approximation.
fn derivative<F: Fn(f64) -> f64 + Copy>(f: F, x: f64, n: u32, h: f64) -> f64 {
    if n == 0 {
        f(x)
    } else {
        (derivative(f, x + h, n - 1, h) - derivative(f, x - h, n - 1, h)) / (2.0 * h)
    }
}

/// Computes Hardy's Z-function using Riemann-Siegel with optional higher terms.
///
/// # Arguments
/// * `t` - Imaginary part.
/// * `terms` - Number of remainder terms (0: basic, 1: +C1, 2: +C1+C2).
///
/// # Returns
/// Z(t) value.
pub fn z_func(t: f64, terms: u32) -> f64 {
    let sqrt_t_over_2pi = (t / (2.0 * PI)).sqrt();
    let nu = sqrt_t_over_2pi.floor() as i64;
    let p = sqrt_t_over_2pi - nu as f64;

    let mut sum = 0.0;
    for k in 1..=(nu as usize) {
        let sqrt_k = (k as f64).sqrt();
        let arg = theta(t) - t * (k as f64).ln();
        sum += arg.cos() / sqrt_k;
    }
    sum *= 2.0;

    let sign = if nu % 2 == 0 { -1.0 } else { 1.0 };
    let a = t / (2.0 * PI);
    let scale = a.powf(-0.25);

    let psi = move |pp: f64| (2.0 * PI * (pp * pp - pp - 1.0 / 16.0)).cos() / (2.0 * PI * pp).cos();

    let c0 = psi(p);
    let mut r = sign * scale * c0;

    if terms >= 1 {
        let psi3 = derivative(psi, p, 3, 1e-5);
        let c1 = -1.0 / (96.0 * PI * PI) * psi3;
        r += sign * scale * c1 * a.powf(-0.5);
    }

    if terms >= 2 {
        let psi2 = derivative(psi, p, 2, 1e-5);
        let psi6 = derivative(psi, p, 6, 1e-5);
        let c2 = 1.0 / (18432.0 * PI.powi(4) as f64) * psi6 + 1.0 / (64.0 * PI.powi(2) as f64) * psi2;
        r += sign * scale * c2 * a.powf(-1.0);
    }

    sum + r
}

/// Finds a zero of Z(t) in [a, b] using bisection.
///
/// # Arguments
/// * `a` - Lower bound.
/// * `b` - Upper bound.
/// * `tol` - Tolerance.
/// * `terms` - Remainder terms count.
///
/// # Returns
/// Some(t) if zero found, None if no sign change.
#[allow(unused_assignments)]
pub fn find_zero(mut a: f64, mut b: f64, tol: f64, terms: u32) -> Option<f64> {
    let mut za = z_func(a, terms);
    let mut zb = z_func(b, terms);
    if za * zb > 0.0 {
        return None;
    }
    loop {
        let mid = (a + b) / 2.0;
        let zm = z_func(mid, terms);
        if (b - a) < tol {
            return Some(mid);
        }
        if za * zm > 0.0 {
            a = mid;
            za = zm;
        } else {
            b = mid;
            zb = zm;
        }
    }
}

/// Estimates the imaginary part of the nth zero.
///
/// # Arguments
/// * `n` - Zero index (starting from 1).
///
/// # Returns
/// Approximated t_n.
pub fn estimate_t(n: f64) -> f64 {
    if n < 1.0 {
        return 0.0;
    }
    let mut t = 2.0 * PI * n * n.ln();
    for _ in 0..20 {
        let a = t / (2.0 * PI);
        let log_a = a.ln();
        let nt = a * log_a - a + 7.0 / 8.0;
        let dnt = log_a;
        t -= (nt - n) * (2.0 * PI / dnt);
    }
    t
}