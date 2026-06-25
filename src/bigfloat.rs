//! Arbitrary-precision evaluation of Hardy's Z-function (optional `bigfloat` feature).
//!
//! The `f64` implementation in the crate root is limited at large `t` by catastrophic
//! cancellation: the main-sum argument `theta(t) - t*ln(k)` grows like `t*ln(t)`, so once
//! it exceeds ~1e15 the `f64` value no longer pins the fractional part that `cos` needs.
//! Computing the same Riemann-Siegel formula in arbitrary precision keeps that argument
//! accurate, breaking the ceiling. (At small `t` the dominant error is the asymptotic
//! truncation of the remainder series, which extra precision cannot fix — so this module
//! is most useful for large `t`.)

use crate::Precision;
use dashu_float::DBig;
use std::str::FromStr;

/// ~99 significant digits of pi; bounds the usable precision of this module.
const PI_DIGITS: &str = "3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211707";

/// Largest precision (decimal digits) supported, limited by [`PI_DIGITS`].
pub const MAX_DIGITS: usize = 90;

fn pi(prec: usize) -> DBig {
    DBig::from_str(PI_DIGITS)
        .unwrap()
        .with_precision(prec)
        .value()
}
fn int(n: i64, prec: usize) -> DBig {
    DBig::from(n).with_precision(prec).value()
}
fn from_f64(x: f64, prec: usize) -> DBig {
    // f64 is binary; round-trip through a decimal string to seed the base-10 DBig.
    DBig::from_str(&format!("{x:.17e}"))
        .unwrap()
        .with_precision(prec)
        .value()
}
fn to_f64(x: &DBig) -> f64 {
    x.to_f64().value()
}
fn sqrt(x: &DBig, prec: usize) -> DBig {
    x.clone().powf(&(int(1, prec) / int(2, prec)))
}

/// Riemann-Siegel theta in arbitrary precision.
fn theta(t: &DBig, prec: usize) -> DBig {
    let pi = pi(prec);
    let two_pi = int(2, prec) * &pi;
    let half_t = t / int(2, prec);
    let ln_term = (t / &two_pi).ln();
    let t2 = t * t;
    let t3 = &t2 * t;
    let t5 = &t3 * &t2;
    &half_t * &ln_term - &half_t - &pi / int(8, prec)
        + int(1, prec) / (int(48, prec) * t)
        + int(7, prec) / (int(5760, prec) * &t3)
        - int(31, prec) / (int(80640, prec) * &t5)
}

/// The Riemann-Siegel `Psi` function `cos(2*pi*(p^2-p-1/16)) / cos(2*pi*p)`.
fn psi(p: &DBig, pi: &DBig, prec: usize) -> DBig {
    let two_pi = int(2, prec) * pi;
    let inner = p * p - p - int(1, prec) / int(16, prec);
    (&two_pi * &inner).cos() / (&two_pi * p).cos()
}

/// Hardy's Z-function at `t`, evaluated with `prec` decimal digits of working precision.
///
/// Mirrors [`crate::z_func`] term-for-term (main sum plus `C0`/`C1`/`C2`), but in
/// arbitrary precision. `prec` is clamped to [`MAX_DIGITS`].
pub fn z_func(t: f64, prec: usize, precision: Precision) -> f64 {
    if !t.is_finite() || t <= 0.0 {
        return f64::NAN;
    }
    let prec = prec.clamp(16, MAX_DIGITS);
    let t = from_f64(t, prec);
    let pi = pi(prec);
    let two_pi = int(2, prec) * &pi;
    let a = &t / &two_pi;
    let sqrt_a = sqrt(&a, prec);
    let nu = to_f64(&sqrt_a).floor() as i64;
    let p = &sqrt_a - int(nu, prec);

    let theta_t = theta(&t, prec);
    let mut sum = int(0, prec);
    for k in 1..=nu.max(0) {
        let kd = int(k, prec);
        let arg = &theta_t - &t * kd.clone().ln();
        sum += arg.cos() / sqrt(&kd, prec);
    }
    sum *= int(2, prec);

    let sign = if nu % 2 == 0 {
        int(-1, prec)
    } else {
        int(1, prec)
    };
    let scale = a.clone().powf(&(int(-1, prec) / int(4, prec)));
    let mut r = &sign * &scale * psi(&p, &pi, prec);

    let terms = precision.correction_terms();
    if terms >= 1 {
        // C1 = -psi'''(p) / (96 pi^2); finite differences are stable in high precision.
        let h = int(1, prec) / int(1000, prec);
        let d3 = (psi(&(&p + int(2, prec) * &h), &pi, prec)
            - int(2, prec) * psi(&(&p + &h), &pi, prec)
            + int(2, prec) * psi(&(&p - &h), &pi, prec)
            - psi(&(&p - int(2, prec) * &h), &pi, prec))
            / (int(2, prec) * &h * &h * &h);
        let c1 = int(-1, prec) / (int(96, prec) * &pi * &pi) * d3;
        r += &sign * &scale * &c1 * a.clone().powf(&(int(-1, prec) / int(2, prec)));
    }
    if terms >= 2 {
        let h = int(1, prec) / int(100, prec);
        let two = int(2, prec);
        let d2 = (psi(&(&p + &h), &pi, prec) - &two * psi(&p, &pi, prec)
            + psi(&(&p - &h), &pi, prec))
            / (&h * &h);
        // 6th derivative via the standard central stencil.
        let d6 = (psi(&(&p - int(3, prec) * &h), &pi, prec)
            - int(6, prec) * psi(&(&p - &two * &h), &pi, prec)
            + int(15, prec) * psi(&(&p - &h), &pi, prec)
            - int(20, prec) * psi(&p, &pi, prec)
            + int(15, prec) * psi(&(&p + &h), &pi, prec)
            - int(6, prec) * psi(&(&p + &two * &h), &pi, prec)
            + psi(&(&p + int(3, prec) * &h), &pi, prec))
            / (&h * &h * &h * &h * &h * &h);
        let pi2 = &pi * &pi;
        let pi4 = &pi2 * &pi2;
        let c2 = int(1, prec) / (int(18432, prec) * &pi4) * d6
            + int(1, prec) / (int(64, prec) * &pi2) * d2;
        r += &sign * &scale * &c2 * a.clone().powf(&(int(-1, prec)));
    }

    to_f64(&(sum + r))
}

/// Finds a zero of `Z(t)` in `[a, b]` using bisection on the arbitrary-precision
/// [`z_func`]. Returns `None` for invalid input or no bracketed sign change.
pub fn find_zero(a: f64, b: f64, tol: f64, prec: usize, precision: Precision) -> Option<f64> {
    if !a.is_finite() || !b.is_finite() || a >= b {
        return None;
    }
    let mut lo = a;
    let mut hi = b;
    let z_lo = z_func(lo, prec, precision);
    let z_hi = z_func(hi, prec, precision);
    if !z_lo.is_finite() || !z_hi.is_finite() || z_lo * z_hi > 0.0 {
        return None;
    }
    let mut s_lo = z_lo.signum();
    for _ in 0..200 {
        let mid = lo + (hi - lo) / 2.0;
        if (hi - lo) < tol {
            return Some(mid);
        }
        let z_mid = z_func(mid, prec, precision);
        if z_mid == 0.0 {
            return Some(mid);
        }
        if s_lo * z_mid.signum() > 0.0 {
            lo = mid;
            s_lo = z_mid.signum();
        } else {
            hi = mid;
        }
    }
    Some(lo + (hi - lo) / 2.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn agrees_with_f64_at_moderate_t() {
        // Same formula, so the big version must match the f64 version where f64 is good.
        for &t in &[15.0, 30.0, 50.0, 123.4] {
            let big = z_func(t, 50, Precision::Order2);
            let small = crate::z_func(t, Precision::Order2);
            assert!((big - small).abs() < 1e-6, "t={t}: big={big}, f64={small}");
        }
    }

    #[test]
    fn recovers_first_zero() {
        let z = find_zero(14.0, 15.0, 1e-12, 50, Precision::Order2).unwrap();
        assert!((z - 14.134725141734693).abs() < 1e-3, "got {z}");
    }

    #[test]
    fn large_t_value_is_precision_stable() {
        // Near the millionth zero the value is tiny and must not drift with precision,
        // unlike f64 whose main-sum argument has lost absolute accuracy by this height.
        let t = 600269.6770190396;
        let lo = z_func(t, 25, Precision::Order1);
        let hi = z_func(t, 75, Precision::Order1);
        assert!((lo - hi).abs() < 1e-9, "not stable: {lo} vs {hi}");
    }
}
