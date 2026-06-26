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

/// Precomputed per-term constants for the main sum: `ln(k)` and `1/sqrt(k)` for
/// `k = 1..=nu_max`, indexed by `k - 1`.
///
/// These depend only on `k` (not `t`), so a single table is reused across every
/// evaluation of a [`find_zero`] bisection — turning the dominant `ln`/`exp` cost from
/// per-iteration into one-time. `1/sqrt(k)` is derived as `exp(-ln(k)/2)`, reusing the
/// logarithm rather than paying a second one inside `powf`.
struct SumTables {
    ln_k: Vec<DBig>,
    inv_sqrt_k: Vec<DBig>,
}

impl SumTables {
    fn build(nu_max: usize, prec: usize) -> Self {
        if nu_max == 0 {
            return SumTables {
                ln_k: Vec::new(),
                inv_sqrt_k: Vec::new(),
            };
        }
        let neg_half = int(-1, prec) / int(2, prec);
        let one = |k: usize| -> (DBig, DBig) {
            let l = int(k as i64, prec).ln();
            let inv = (&l * &neg_half).exp();
            (l, inv)
        };

        let threads = worker_threads();
        let parts: Vec<Vec<(DBig, DBig)>> = if threads <= 1 || nu_max < 64 {
            vec![(1..=nu_max).map(&one).collect()]
        } else {
            // Building the table is itself nu_max logarithms/exponentials; parallelize it
            // so it does not dominate once the secant root finder needs few iterations.
            let chunk = nu_max.div_ceil(threads);
            std::thread::scope(|s| {
                let mut handles = Vec::new();
                let mut start = 1;
                while start <= nu_max {
                    let end = (start + chunk - 1).min(nu_max);
                    let one = &one;
                    handles.push(s.spawn(move || (start..=end).map(one).collect::<Vec<_>>()));
                    start = end + 1;
                }
                handles.into_iter().map(|h| h.join().unwrap()).collect()
            })
        };

        let mut ln_k = Vec::with_capacity(nu_max);
        let mut inv_sqrt_k = Vec::with_capacity(nu_max);
        for part in parts {
            for (l, inv) in part {
                ln_k.push(l);
                inv_sqrt_k.push(inv);
            }
        }
        SumTables { ln_k, inv_sqrt_k }
    }
}

/// Number of worker threads to use for the parallel table build and main sum.
fn worker_threads() -> usize {
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
        .min(8)
}

/// `2 * sum_{k=1}^{nu} cos(theta_t - t*ln(k)) / sqrt(k)`, the Riemann-Siegel main sum.
///
/// The per-term `cos` (with range reduction of an argument that grows like `t*ln t`) is
/// the irreducible cost at large `t`; it is spread across worker threads. The `ln`/`sqrt`
/// factors come from the precomputed `tables`, so the inner loop is one `cos` and a
/// couple of big multiplies per term.
fn main_sum(t: &DBig, theta_t: &DBig, nu: usize, prec: usize, tables: &SumTables) -> DBig {
    let two = int(2, prec);
    if nu == 0 {
        return int(0, prec);
    }

    let term = |k: usize| -> DBig {
        let arg = theta_t - t * &tables.ln_k[k - 1];
        arg.cos() * &tables.inv_sqrt_k[k - 1]
    };

    let threads = worker_threads();
    if threads <= 1 || nu < 64 {
        let mut acc = int(0, prec);
        for k in 1..=nu {
            acc += term(k);
        }
        return acc * two;
    }

    let chunk = nu.div_ceil(threads);
    let partials: Vec<DBig> = std::thread::scope(|s| {
        let mut handles = Vec::new();
        let mut start = 1;
        while start <= nu {
            let end = (start + chunk - 1).min(nu);
            let term = &term;
            handles.push(s.spawn(move || {
                let mut acc = int(0, prec);
                for k in start..=end {
                    acc += term(k);
                }
                acc
            }));
            start = end + 1;
        }
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });

    let mut acc = int(0, prec);
    for p in partials {
        acc += p;
    }
    acc * two
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
    let nu = nu_of(t, prec);
    let tables = SumTables::build(nu, prec);
    z_eval(t, prec, precision, &tables)
}

/// `nu = floor(sqrt(t / (2*pi)))`, the number of terms in the main sum.
fn nu_of(t: f64, prec: usize) -> usize {
    let a = from_f64(t, prec) / (int(2, prec) * pi(prec));
    to_f64(&sqrt(&a, prec)).floor().max(0.0) as usize
}

/// Evaluates `Z(t)` using a prebuilt `tables` (which must cover `nu(t)`).
fn z_eval(t: f64, prec: usize, precision: Precision, tables: &SumTables) -> f64 {
    if !t.is_finite() || t <= 0.0 {
        return f64::NAN;
    }
    let t = from_f64(t, prec);
    let pi = pi(prec);
    let two_pi = int(2, prec) * &pi;
    let a = &t / &two_pi;
    let sqrt_a = sqrt(&a, prec);
    let nu = to_f64(&sqrt_a).floor() as i64;
    let p = &sqrt_a - int(nu, prec);

    let theta_t = theta(&t, prec);
    let sum = main_sum(&t, &theta_t, nu.max(0) as usize, prec, tables);

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
    let prec = prec.clamp(16, MAX_DIGITS);
    // nu is monotonic in t, so the table built for the upper bound covers every t in
    // [a, b]; build it once and reuse it across all bisection iterations.
    let tables = SumTables::build(nu_of(b, prec), prec);

    let mut lo = a;
    let mut hi = b;
    let mut f_lo = z_eval(lo, prec, precision, &tables);
    let mut f_hi = z_eval(hi, prec, precision, &tables);
    if !f_lo.is_finite() || !f_hi.is_finite() || f_lo * f_hi > 0.0 {
        return None;
    }

    // Safeguarded secant: a false-position step when it lands comfortably inside the
    // bracket, otherwise bisection. Each big-precision evaluation is expensive, so the
    // ~6-10 secant steps this needs are a large saving over ~33 bisection steps.
    let mut force_bisect = false;
    for _ in 0..200 {
        let width = hi - lo;
        if width < tol {
            return Some(lo + width / 2.0);
        }

        let mid = if force_bisect {
            lo + width / 2.0
        } else {
            let s = lo - f_lo * (hi - lo) / (f_hi - f_lo);
            if s.is_finite() && s > lo + 0.05 * width && s < hi - 0.05 * width {
                s
            } else {
                lo + width / 2.0
            }
        };

        let f_mid = z_eval(mid, prec, precision, &tables);
        if f_mid == 0.0 {
            return Some(mid);
        }
        if f_lo * f_mid < 0.0 {
            hi = mid;
            f_hi = f_mid;
        } else {
            lo = mid;
            f_lo = f_mid;
        }
        force_bisect = (hi - lo) > 0.5 * width;
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
