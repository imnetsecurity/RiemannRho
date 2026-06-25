//! RiemannRho: Library for approximating nontrivial zeros of the Riemann zeta function.
//!
//! This module provides functions to compute Hardy's Z-function using the Riemann-Siegel formula,
//! find zeros via a bracketed root finder, and estimate the nth zero. It supports optional
//! higher-order remainder terms for precision (High-Order Correction Mode). No visualization or
//! I/O is included here; that belongs to the binaries.

use std::f64::consts::PI;

/// Maximum number of root-finding iterations before [`find_zero`] gives up.
///
/// With `f64` a bracket can be halved at most ~60 times before reaching the limit of
/// representable precision, so this also guards against an infinite loop when `tol` is
/// set to 0 (or smaller than the floating-point resolution). The cap is generous to
/// leave room for the secant-accelerated steps.
const MAX_ROOT_ITERS: u32 = 200;

/// Number of remainder correction terms to include when evaluating [`z_func`].
///
/// The corrections follow the Riemann-Siegel expansion; higher variants reduce the
/// asymptotic error at the cost of a few extra evaluations.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Precision {
    /// Main sum plus the leading `C0` remainder term only.
    Base,
    /// Adds the `C1` correction.
    Order1,
    /// Adds the `C1` and `C2` corrections (the "high-order" mode).
    Order2,
}

impl Precision {
    /// Number of remainder correction terms beyond `C0`.
    fn correction_terms(self) -> u32 {
        match self {
            Precision::Base => 0,
            Precision::Order1 => 1,
            Precision::Order2 => 2,
        }
    }
}

/// Computes the Riemann-Siegel theta function approximation.
///
/// The asymptotic expansion is accurate for large `t`; for small `t` (roughly `t < 10`)
/// the error grows because the omitted terms are no longer negligible.
///
/// # Arguments
/// * `t` - The imaginary part (must be positive).
///
/// # Returns
/// The approximated theta(t).
pub fn theta(t: f64) -> f64 {
    let log_term = (t / 2.0) * (t / (2.0 * PI)).ln();
    log_term - t / 2.0 - PI / 8.0 + 1.0 / (48.0 * t) + 7.0 / (5760.0 * t.powi(3))
        - 31.0 / (80640.0 * t.powi(5))
}

/// Approximates the `n`th derivative of `f` at `x` using a central-difference stencil
/// with one level of Richardson extrapolation.
///
/// A fixed tiny step would be catastrophically unstable for high orders: dividing by
/// `(2h)^n` amplifies the ~1e-16 rounding noise of `f` by many orders of magnitude
/// (for `n = 6` by roughly 1e13). Here we use explicit stencils with a step appropriate
/// to the order and Richardson-extrapolate `D(h)` and `D(h/2)` to cancel the leading
/// `O(h^2)` truncation error.
///
/// Accuracy still degrades near the poles of the Riemann-Siegel `Psi` function
/// (fractional part `p ≈ 0.25` or `0.75`).
fn nth_derivative<F: Fn(f64) -> f64 + Copy>(f: F, x: f64, n: u32, h: f64) -> f64 {
    let stencil = |h: f64| -> f64 {
        match n {
            0 => f(x),
            2 => (f(x - h) - 2.0 * f(x) + f(x + h)) / (h * h),
            3 => {
                (-f(x - 2.0 * h) + 2.0 * f(x - h) - 2.0 * f(x + h) + f(x + 2.0 * h))
                    / (2.0 * h.powi(3))
            }
            6 => {
                (f(x - 3.0 * h) - 6.0 * f(x - 2.0 * h) + 15.0 * f(x - h) - 20.0 * f(x)
                    + 15.0 * f(x + h)
                    - 6.0 * f(x + 2.0 * h)
                    + f(x + 3.0 * h))
                    / h.powi(6)
            }
            // Fallback for any other order: recursive central difference.
            _ => {
                fn rec<G: Fn(f64) -> f64 + Copy>(g: G, x: f64, n: u32, h: f64) -> f64 {
                    if n == 0 {
                        g(x)
                    } else {
                        (rec(g, x + h, n - 1, h) - rec(g, x - h, n - 1, h)) / (2.0 * h)
                    }
                }
                rec(f, x, n, h)
            }
        }
    };

    let d_h = stencil(h);
    let d_h2 = stencil(h / 2.0);
    // Richardson extrapolation for an O(h^2) central stencil: (4*D(h/2) - D(h)) / 3.
    (4.0 * d_h2 - d_h) / 3.0
}

/// Computes Hardy's Z-function using the Riemann-Siegel formula with optional
/// higher-order remainder terms.
///
/// The main sum runs over `floor(sqrt(t / 2pi))` terms, giving the `O(sqrt(t))`
/// cost that makes large-`t` evaluation feasible. The remainder correction constants
/// follow Edwards: `C1 = -psi'''/(96 pi^2)` and
/// `C2 = psi''''''/(18432 pi^4) + psi''/(64 pi^2)`.
///
/// # Arguments
/// * `t` - Imaginary part (must be > 0).
/// * `precision` - Which remainder correction terms to include.
///
/// # Returns
/// Z(t) value. Returns `NaN` for non-positive or non-finite `t`.
pub fn z_func(t: f64, precision: Precision) -> f64 {
    if !t.is_finite() || t <= 0.0 {
        return f64::NAN;
    }

    let sqrt_t_over_2pi = (t / (2.0 * PI)).sqrt();
    let nu = sqrt_t_over_2pi.floor() as i64;
    let p = sqrt_t_over_2pi - nu as f64;

    // Main sum: 2 * sum_{k=1}^{nu} cos(theta(t) - t*ln(k)) / sqrt(k).
    // theta(t) is constant across the sum, so it is computed once here rather than
    // re-evaluated inside the loop.
    let theta_t = theta(t);
    let mut sum = 0.0;
    for k in 1..=nu.max(0) {
        let kf = k as f64;
        sum += (theta_t - t * kf.ln()).cos() / kf.sqrt();
    }
    sum *= 2.0;

    // Remainder: (-1)^{nu-1} * a^{-1/4} * [C0 + C1*a^{-1/2} + C2*a^{-1} + ...].
    let sign = if nu % 2 == 0 { -1.0 } else { 1.0 };
    let a = t / (2.0 * PI);
    let scale = a.powf(-0.25);

    let psi = move |pp: f64| (2.0 * PI * (pp * pp - pp - 1.0 / 16.0)).cos() / (2.0 * PI * pp).cos();

    let c0 = psi(p);
    let mut r = sign * scale * c0;

    let terms = precision.correction_terms();
    if terms >= 1 {
        let psi3 = nth_derivative(psi, p, 3, 0.05);
        let c1 = -1.0 / (96.0 * PI * PI) * psi3;
        r += sign * scale * c1 * a.powf(-0.5);
    }

    if terms >= 2 {
        let psi2 = nth_derivative(psi, p, 2, 0.05);
        let psi6 = nth_derivative(psi, p, 6, 0.05);
        let c2 = 1.0 / (18432.0 * PI.powi(4)) * psi6 + 1.0 / (64.0 * PI.powi(2)) * psi2;
        r += sign * scale * c2 * a.powf(-1.0);
    }

    sum + r
}

/// Finds a zero of Z(t) in `[a, b]`.
///
/// Uses a bracketed root finder: a secant (false-position) step when it lands safely
/// inside the bracket, otherwise a bisection step. A bisection step is also forced
/// whenever the bracket fails to shrink by at least half, which guarantees convergence
/// at least as fast as plain bisection while typically being faster.
///
/// **Note:** the bracket must be evaluated with the same `precision` it will be refined
/// with. Because `Z(t)` differs slightly between precisions, a sign change that brackets
/// a zero at one precision may not bracket it at another.
///
/// # Arguments
/// * `a` - Lower bound.
/// * `b` - Upper bound.
/// * `tol` - Tolerance on the bracket width.
/// * `precision` - Remainder correction terms to use.
///
/// # Returns
/// `Some(t)` if a sign change is bracketed, `None` if the inputs are invalid or there
/// is no sign change in the interval.
pub fn find_zero(a: f64, b: f64, tol: f64, precision: Precision) -> Option<f64> {
    // Reject invalid input outright. Without this an NaN bound (e.g. from a bad nth
    // estimate) made an unguarded loop spin forever, because `(b - a) < tol` is never
    // true when `b - a` is NaN.
    if !a.is_finite() || !b.is_finite() || a >= b {
        return None;
    }

    let mut lo = a;
    let mut hi = b;
    let mut f_lo = z_func(lo, precision);
    let mut f_hi = z_func(hi, precision);
    if !f_lo.is_finite() || !f_hi.is_finite() || f_lo * f_hi > 0.0 {
        return None;
    }

    let mut force_bisect = false;
    for _ in 0..MAX_ROOT_ITERS {
        let width = hi - lo;
        if width < tol {
            return Some(lo + width / 2.0);
        }

        // Secant / false-position candidate, accepted only when it lies comfortably
        // inside the bracket; otherwise (or when progress stalled) bisect.
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

        let f_mid = z_func(mid, precision);
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
        // If the bracket barely shrank, force a bisection next iteration so width is
        // guaranteed to keep falling toward `tol`.
        force_bisect = (hi - lo) > 0.5 * width;
    }
    Some(lo + (hi - lo) / 2.0)
}

/// Estimates the imaginary part of the `n`th nontrivial zero.
///
/// Solves `N(t) = n` (where `N(t)` is the Riemann-von Mangoldt zero-counting estimate)
/// with Newton's method. For very small `n` the asymptotic counting formula is poor, so
/// the first handful of zeros are returned from a small table of known values to provide
/// a robust bracket for [`find_zero`].
///
/// `n` is taken as `f64` so that very large ordinals can be requested in scientific
/// notation; note that beyond `2^53` the integer value can no longer be represented
/// exactly.
///
/// # Arguments
/// * `n` - Zero index (1-based).
///
/// # Returns
/// Approximated `t_n`, or `NaN` for `n < 1`.
pub fn estimate_t(n: f64) -> f64 {
    if !n.is_finite() || n < 1.0 {
        return f64::NAN;
    }

    // Known imaginary parts of the first few zeros. The asymptotic formula below is
    // unreliable for small n (and would produce NaN for n = 1, where the old initial
    // guess 2*pi*n*ln(n) collapsed to 0).
    const KNOWN_ZEROS: [f64; 5] = [
        14.134725141734693,
        21.022039638771555,
        25.01085758014569,
        30.424876125859513,
        32.93506158773919,
    ];
    let idx = n as usize;
    if (n.fract() == 0.0) && idx >= 1 && idx <= KNOWN_ZEROS.len() {
        return KNOWN_ZEROS[idx - 1];
    }

    // Initial guess from the asymptotic gram-point spacing, guarded so the argument of
    // the logarithms stays comfortably positive even for small n.
    let mut t = (2.0 * PI * n).max(2.0 * PI * n * n.ln()).max(10.0);
    for _ in 0..50 {
        let a = t / (2.0 * PI);
        let log_a = a.ln();
        // N(t) ~ a*ln(a) - a + 7/8, with dN/dt = ln(a) / (2*pi).
        let n_t = a * log_a - a + 7.0 / 8.0;
        let dn_dt = log_a / (2.0 * PI);
        let step = (n_t - n) / dn_dt;
        t -= step;
        if step.abs() < 1e-9 {
            break;
        }
    }
    t
}

/// Smooth part of the Riemann-von Mangoldt zero-counting function:
/// `N(t) = theta(t)/pi + 1 + S(t)`, returning `theta(t)/pi + 1` (i.e. neglecting the
/// oscillating term `S(t)`, which averages to 0).
///
/// Rounding this to the nearest integer predicts how many nontrivial zeros have
/// imaginary part in `(0, t]`, and comparing it with the number actually found is a
/// (Turing-flavored) consistency check on the Riemann hypothesis up to height `t`.
pub fn expected_zero_count(t: f64) -> f64 {
    if !t.is_finite() || t <= 0.0 {
        return 0.0;
    }
    theta(t) / PI + 1.0
}

/// Default number of scan samples per average zero spacing.
const DEFAULT_SCAN_RESOLUTION: f64 = 10.0;

/// Walks the sign changes of `Z(t)` over `(t_start, t_max]`, calling `on_zero` with each
/// refined zero in increasing order. Iteration stops early if `on_zero` returns `false`.
///
/// The scan step is the local average zero spacing `2*pi/ln(t/2pi)` divided by
/// `resolution`. A higher `resolution` is more likely to separate very close pairs of
/// zeros at the cost of more evaluations.
fn scan_zeros<F: FnMut(f64) -> bool>(
    t_start: f64,
    t_max: f64,
    precision: Precision,
    resolution: f64,
    mut on_zero: F,
) {
    let mut t = t_start.max(1.0);
    let mut z = z_func(t, precision);
    while t < t_max {
        let a = (t / (2.0 * PI)).max(1.0);
        let spacing = (2.0 * PI / a.ln().max(0.05)).max(0.05);
        let step = (spacing / resolution).min(t_max - t);
        let t_next = t + step;
        let z_next = z_func(t_next, precision);
        if z.is_finite() && z_next.is_finite() && z != 0.0 && z * z_next < 0.0 {
            if let Some(root) = find_zero(t, t_next, 1e-9, precision) {
                if root <= t_max && !on_zero(root) {
                    return;
                }
            }
        }
        t = t_next;
        z = z_next;
    }
}

/// Locates every zero of `Z(t)` with imaginary part in `(0, t_max]`, in increasing order.
///
/// Cost is `O(t_max * sqrt(t_max))`, so this is intended for exploration over moderate
/// heights rather than astronomically large `t_max`.
pub fn zeros_below(t_max: f64, precision: Precision) -> Vec<f64> {
    let mut zeros = Vec::new();
    scan_zeros(1.0, t_max, precision, DEFAULT_SCAN_RESOLUTION, |root| {
        zeros.push(root);
        true
    });
    zeros
}

/// Counts the nontrivial zeros with imaginary part in `(0, t_max]` found on the
/// critical line. Compare with `round(`[`expected_zero_count`]`(t_max))`.
pub fn count_zeros_below(t_max: f64, precision: Precision) -> usize {
    let mut count = 0usize;
    scan_zeros(1.0, t_max, precision, DEFAULT_SCAN_RESOLUTION, |_| {
        count += 1;
        true
    });
    count
}

/// Largest plausible magnitude of `S(t)` at the heights this tool explores. `S(t)` grows
/// only like `O(log t)` and stays below ~1 for small `t`, so an implied `|S|` beyond this
/// signals a miscount (a missed or spurious zero) rather than a genuine fluctuation.
const MAX_PLAUSIBLE_S: f64 = 2.5;

/// Outcome of a Turing-flavored zero-count verification (see [`verify_zero_count`]).
#[derive(Clone, Copy, Debug)]
pub struct CountReport {
    /// Number of zeros found on the critical line in `(0, t_max]`.
    pub found: usize,
    /// Theoretical smooth count `theta(t_max)/pi + 1`.
    pub expected: f64,
    /// Implied `S(t_max) = found - expected`. The true zero count is
    /// `expected + S`, so `S` should be a small `O(1)` fluctuation around 0.
    pub s: f64,
    /// Whether the implied `S` is within the plausible range (i.e. no zero appears to be
    /// missing or doubled).
    pub consistent: bool,
    /// Scan samples per average spacing that produced this result (raised automatically
    /// when the initial scan appeared to miss zeros).
    pub resolution: f64,
}

/// Counts zeros up to `t_max` and checks the tally against the theoretical count.
///
/// The Riemann-von Mangoldt formula is `N(t) = theta(t)/pi + 1 + S(t)`, where `S(t)`
/// oscillates around 0. So the found count should *not* equal `round(theta/pi + 1)`
/// exactly — it equals that smooth value plus the `O(1)` term `S(t)` (for example at
/// `t = 50` there are 10 zeros while `theta/pi + 1 = 9.42`, i.e. `S = 0.58`). We instead
/// check that the *implied* `S = found - expected` is a plausible small fluctuation.
///
/// A coarse scan can step over an unusually close pair of zeros, which shows up as a
/// suspiciously *negative* `S`; in that case the scan resolution is raised and the count
/// retried, so the check is self-healing rather than merely reporting the discrepancy.
pub fn verify_zero_count(t_max: f64, precision: Precision) -> CountReport {
    let expected = expected_zero_count(t_max);

    let mut resolution = DEFAULT_SCAN_RESOLUTION;
    let mut found = count_zeros_below(t_max, precision);
    // Refine only while the count looks too low (S < -1), the signature of a missed
    // zero. A finer grid cannot remove a genuine positive fluctuation, so don't bother.
    for _ in 0..4 {
        if (found as f64) - expected >= -1.0 {
            break;
        }
        resolution *= 4.0;
        let mut c = 0usize;
        scan_zeros(1.0, t_max, precision, resolution, |_| {
            c += 1;
            true
        });
        if c <= found {
            break; // refinement found nothing new; the count has stabilized
        }
        found = c;
    }

    let s = found as f64 - expected;
    CountReport {
        found,
        expected,
        s,
        consistent: s.abs() <= MAX_PLAUSIBLE_S,
        resolution,
    }
}

/// Finds the exact `n`th nontrivial zero (1-based) by scanning zeros sequentially from
/// the bottom, so the result is guaranteed to be the `n`th — not merely a zero near an
/// asymptotic estimate (which [`estimate_t`] alone can miss).
///
/// Because it scans the whole range `(0, t_n]`, this is practical only for moderate `n`;
/// for very large `n` use [`estimate_t`] to bracket a single zero instead.
///
/// # Returns
/// `Some(t_n)`, or `None` if `n == 0` or the location could not be estimated.
pub fn nth_zero(n: u64, precision: Precision) -> Option<f64> {
    if n == 0 {
        return None;
    }
    // Scan a little past the asymptotic estimate to be sure the nth zero is included.
    let est = estimate_t(n as f64);
    if !est.is_finite() {
        return None;
    }
    let t_max = est + 50.0;

    let mut count = 0u64;
    let mut result = None;
    scan_zeros(1.0, t_max, precision, DEFAULT_SCAN_RESOLUTION, |root| {
        count += 1;
        if count == n {
            result = Some(root);
            false
        } else {
            true
        }
    });
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Reference imaginary parts of the first nontrivial zeros (Odlyzko's tables).
    const REFERENCE_ZEROS: [f64; 5] = [
        14.134725141734693,
        21.022039638771555,
        25.01085758014569,
        30.424876125859513,
        32.93506158773919,
    ];

    #[test]
    fn finds_first_zero_within_tolerance() {
        let z = find_zero(14.0, 15.0, 1e-10, Precision::Order2)
            .expect("first zero should be bracketed");
        assert!(
            (z - REFERENCE_ZEROS[0]).abs() < 1e-3,
            "got {z}, expected ~{}",
            REFERENCE_ZEROS[0]
        );
    }

    #[test]
    fn high_order_is_at_least_as_accurate_as_base() {
        let base = find_zero(14.0, 15.0, 1e-12, Precision::Base).unwrap();
        let high = find_zero(14.0, 15.0, 1e-12, Precision::Order2).unwrap();
        let err_base = (base - REFERENCE_ZEROS[0]).abs();
        let err_high = (high - REFERENCE_ZEROS[0]).abs();
        assert!(
            err_high <= err_base + 1e-9,
            "high-order error {err_high} should not exceed base error {err_base}"
        );
    }

    #[test]
    fn brackets_several_known_zeros() {
        // Search a window around each reference zero and confirm we recover it.
        for &z0 in REFERENCE_ZEROS.iter() {
            let z = find_zero(z0 - 0.5, z0 + 0.5, 1e-10, Precision::Order2)
                .unwrap_or_else(|| panic!("no zero bracketed near {z0}"));
            assert!((z - z0).abs() < 1e-2, "got {z}, expected ~{z0}");
        }
    }

    #[test]
    fn secant_and_bisection_agree() {
        // The accelerated finder must land on the same root as plain bisection would.
        for &z0 in REFERENCE_ZEROS.iter() {
            let z = find_zero(z0 - 0.4, z0 + 0.4, 1e-12, Precision::Base).unwrap();
            assert!((z - z0).abs() < 1e-2, "got {z}, expected ~{z0}");
            // The reported point must actually be (nearly) a zero of Z.
            assert!(z_func(z, Precision::Base).abs() < 1e-6);
        }
    }

    #[test]
    fn find_zero_rejects_invalid_input() {
        assert_eq!(find_zero(f64::NAN, 15.0, 1e-10, Precision::Base), None);
        assert_eq!(find_zero(15.0, 14.0, 1e-10, Precision::Base), None); // a >= b
        assert_eq!(find_zero(5.0, 6.0, 1e-10, Precision::Base), None); // no sign change
    }

    #[test]
    fn estimate_t_handles_small_n() {
        // n = 1 used to produce NaN (and an infinite root-finding loop downstream).
        assert!((estimate_t(1.0) - REFERENCE_ZEROS[0]).abs() < 1e-6);
        assert!(estimate_t(0.0).is_nan());
        assert!(estimate_t(-5.0).is_nan());
    }

    #[test]
    fn estimate_t_is_finite_for_large_n() {
        let t = estimate_t(1_000_000.0);
        assert!(t.is_finite() && t > 0.0);
    }

    #[test]
    fn z_func_rejects_nonpositive_t() {
        assert!(z_func(0.0, Precision::Base).is_nan());
        assert!(z_func(-1.0, Precision::Base).is_nan());
        assert!(z_func(f64::NAN, Precision::Base).is_nan());
    }

    #[test]
    fn counts_known_number_of_zeros() {
        // There are 10 nontrivial zeros with 0 < gamma <= 50 (the next is at ~52.97).
        assert_eq!(count_zeros_below(50.0, Precision::Order2), 10);
        // The first five recovered zeros must match the reference values in order.
        let zeros = zeros_below(35.0, Precision::Order2);
        assert_eq!(zeros.len(), 5);
        for (got, &want) in zeros.iter().zip(REFERENCE_ZEROS.iter()) {
            assert!((got - want).abs() < 1e-3, "got {got}, expected ~{want}");
        }
    }

    #[test]
    fn counted_zeros_match_expected_count() {
        // The Turing-flavored check: the count equals theta(T)/pi+1 plus the small S(T)
        // term, so found - expected (= S) must be a modest O(1) fluctuation, not zero.
        for &t in &[50.0, 100.0, 200.0, 300.0] {
            let found = count_zeros_below(t, Precision::Order2) as f64;
            let expected = expected_zero_count(t);
            assert!(
                (found - expected).abs() < 1.5,
                "implausible S at T={t}: found {found}, expected {expected}"
            );
        }
    }

    #[test]
    fn nth_zero_is_exact() {
        // Sequential scanning must return the true nth zero, in order.
        for (i, &want) in REFERENCE_ZEROS.iter().enumerate() {
            let n = (i + 1) as u64;
            let got = nth_zero(n, Precision::Order2).unwrap();
            assert!(
                (got - want).abs() < 1e-3,
                "zero #{n}: got {got}, expected ~{want}"
            );
        }
        assert_eq!(nth_zero(0, Precision::Base), None);
    }

    #[test]
    fn verify_zero_count_is_consistent() {
        for &t in &[50.0, 100.0, 300.0] {
            let report = verify_zero_count(t, Precision::Order2);
            assert!(
                report.consistent,
                "T={t}: found {} vs expected {:.3} at resolution {}",
                report.found, report.expected, report.resolution
            );
            // A clean range should not need refinement beyond the default resolution.
            assert_eq!(report.resolution, 10.0);
        }
    }
}
