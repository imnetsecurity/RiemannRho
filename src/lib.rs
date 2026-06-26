//! RiemannRho: Library for approximating nontrivial zeros of the Riemann zeta function.
//!
//! This module provides functions to compute Hardy's Z-function using the Riemann-Siegel formula,
//! find zeros via a bracketed root finder, and estimate the nth zero. It supports optional
//! higher-order remainder terms for precision (High-Order Correction Mode). No visualization or
//! I/O is included here; that belongs to the binaries.

use std::f64::consts::PI;

/// Optional arbitrary-precision evaluation (enable with the `bigfloat` feature).
#[cfg(feature = "bigfloat")]
pub mod bigfloat;

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
///
/// # Examples
/// ```
/// use riemannrho::{z_func, Precision};
/// // Z(t) is real and changes sign across the first zero near t = 14.1347.
/// let left = z_func(14.0, Precision::Order2);
/// let right = z_func(14.3, Precision::Order2);
/// assert!(left * right < 0.0);
/// ```
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
///
/// # Examples
/// ```
/// use riemannrho::{find_zero, Precision};
/// let z = find_zero(14.0, 15.0, 1e-10, Precision::Order2).unwrap();
/// assert!((z - 14.134725).abs() < 1e-3);
/// // No sign change in [5, 6], so nothing is bracketed.
/// assert!(find_zero(5.0, 6.0, 1e-10, Precision::Order2).is_none());
/// ```
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

/// The Chebyshev function `psi(x) = sum_{p^k <= x} ln(p)`, computed directly by sieving.
///
/// This is the "true" prime-power count that the explicit formula reconstructs from the
/// zeros; the prime number theorem says `psi(x) ~ x`. Intended for moderate `x`, as it
/// sieves all primes up to `x`.
pub fn chebyshev_psi(x: f64) -> f64 {
    if !x.is_finite() || x < 2.0 {
        return 0.0;
    }
    let n = x.floor() as usize;
    let mut is_composite = vec![false; n + 1];
    let mut sum = 0.0;
    for p in 2..=n {
        if is_composite[p] {
            continue;
        }
        // p is prime: add ln(p) once for every prime power p^k <= x.
        let ln_p = (p as f64).ln();
        let mut power = p as u64;
        while power <= n as u64 {
            sum += ln_p;
            match power.checked_mul(p as u64) {
                Some(next) => power = next,
                None => break,
            }
        }
        let mut m = p * p;
        while m <= n {
            is_composite[m] = true;
            m += p;
        }
    }
    sum
}

/// Approximates the Chebyshev function `psi(x)` from the nontrivial zeros via Riemann's
/// explicit (von Mangoldt) formula — reconstructing the primes *from the zeros*.
///
/// With `rho = 1/2 + i*gamma` paired with its conjugate, the formula is
///
/// `psi(x) = x - sum_{gamma>0} (cos(g) + 2*gamma*sin(g)) * sqrt(x) / (1/4 + gamma^2)`
/// `          - ln(2*pi) - (1/2) ln(1 - x^-2)`,  where `g = gamma * ln(x)`.
///
/// `gammas` are the positive imaginary parts of the zeros (e.g. from [`zeros_below`]); more
/// zeros sharpen the reconstruction of the prime-power staircase. The series converges
/// slowly and oscillates (Gibbs-like) near prime powers, so a finite sum is an
/// approximation, sharpest away from the jumps.
///
/// Returns `NaN` for `x <= 1`.
pub fn psi_from_zeros(x: f64, gammas: &[f64]) -> f64 {
    if !x.is_finite() || x <= 1.0 {
        return f64::NAN;
    }
    let ln_x = x.ln();
    let sqrt_x = x.sqrt();
    let oscillation: f64 = gammas
        .iter()
        .map(|&g| {
            let arg = g * ln_x;
            (arg.cos() + 2.0 * g * arg.sin()) / (0.25 + g * g)
        })
        .sum();
    x - sqrt_x * oscillation - (2.0 * PI).ln() - 0.5 * (1.0 - x.powi(-2)).ln()
}

/// Normalized ("unfolded") gaps between consecutive zeros, rescaled so the mean spacing
/// is 1.
///
/// Raw gaps shrink with height because the zero density grows like `ln(t)/2pi`; dividing
/// out that local density makes spacings at different heights comparable. Unfolding via the
/// smooth count `theta/pi`, the gap between consecutive zeros `a < b` becomes
/// `(theta(b) - theta(a)) / pi`, whose mean over many zeros is 1.
///
/// `zeros` must be sorted ascending (as returned by [`zeros_below`]); the result has one
/// fewer element.
pub fn normalized_spacings(zeros: &[f64]) -> Vec<f64> {
    zeros
        .windows(2)
        .map(|w| (theta(w[1]) - theta(w[0])) / PI)
        .collect()
}

/// The Wigner surmise for GUE nearest-neighbor spacings:
/// `P(s) = (32/pi^2) s^2 exp(-4 s^2/pi)`.
///
/// By the Montgomery-Odlyzko law the normalized spacings of the zeta zeros follow the
/// eigenvalue statistics of the Gaussian Unitary Ensemble of random matrices, of which
/// this is the standard closed-form approximation. The defining feature is *level
/// repulsion*: `P(0) = 0`, in sharp contrast to a random (Poisson) sequence's
/// `P(s) = e^{-s}` with `P(0) = 1`.
pub fn wigner_surmise(s: f64) -> f64 {
    (32.0 / (PI * PI)) * s * s * (-4.0 * s * s / PI).exp()
}

/// The `n`th Gram point: the solution of `theta(g) = n * pi`.
///
/// Gram points are defined and strictly increasing for `n >= -1` (where `theta` is past
/// its minimum at `t = 2*pi` and monotonically increasing). They underpin the classical
/// method for isolating zeros: by *Gram's law* the sign of `(-1)^n Z(g_n)` is usually
/// positive, and each Gram interval `[g_{n-1}, g_n)` usually contains exactly one zero —
/// which is the basis of Turing's method for verifying zero counts.
///
/// # Returns
/// `g_n`, or `NaN` for `n < -1`.
///
/// # Examples
/// ```
/// use riemannrho::gram_point;
/// assert!((gram_point(0) - 17.845600).abs() < 1e-4);
/// assert!(gram_point(-2).is_nan());
/// ```
pub fn gram_point(n: i64) -> f64 {
    if n < -1 {
        return f64::NAN;
    }
    let target = n as f64 * PI;
    // N(g_n) = theta(g_n)/pi + 1 + S = n + 1 + S, so the height where N = n+1 is an
    // excellent starting guess; fall back to just above the theta minimum for n = -1.
    let mut t = if n >= 0 {
        let est = estimate_t((n + 1) as f64);
        if est.is_finite() {
            est
        } else {
            10.0
        }
    } else {
        9.6
    };
    for _ in 0..60 {
        let d_theta = 0.5 * (t / (2.0 * PI)).ln(); // theta'(t)
        if d_theta.abs() < 1e-12 {
            break;
        }
        let step = (theta(t) - target) / d_theta;
        t -= step;
        // Stay in the monotonically increasing region t > 2*pi.
        if t <= 2.0 * PI {
            t = 2.0 * PI + 0.1;
        }
        if step.abs() < 1e-10 {
            break;
        }
    }
    t
}

/// A Gram point is "good" when `(-1)^n Z(g_n) > 0` (Gram's law). Good Gram points are the
/// anchors of Turing's method; runs of "bad" ones form Gram blocks.
pub fn gram_point_is_good(n: i64, precision: Precision) -> bool {
    let z = z_func(gram_point(n), precision);
    let sign = if n.rem_euclid(2) == 0 { 1.0 } else { -1.0 };
    sign * z > 0.0
}

/// Outcome of a Gram-block / Turing-method zero count (see [`count_zeros_gram`]).
#[derive(Clone, Copy, Debug)]
pub struct GramCount {
    /// First good Gram point used as the lower anchor (`g_{lower_index}`).
    pub lower: f64,
    /// Last good Gram point not exceeding `t_max` (`g_{upper_index}`).
    pub upper: f64,
    pub lower_index: i64,
    pub upper_index: i64,
    /// Zeros actually found in `(lower, upper]`.
    pub count: usize,
    /// Zeros predicted by the Gram indices: `upper_index - lower_index`.
    pub expected: usize,
    /// Gram points in `[lower, upper]` that violate Gram's law.
    pub gram_law_failures: usize,
    /// Non-trivial Gram blocks (spanning two or more Gram intervals).
    pub gram_blocks: usize,
    /// Gram blocks whose zero count differs from their length (Rosser's-rule violations).
    pub rosser_violations: usize,
    /// `true` when every block satisfied Rosser's rule and the total matched `expected`.
    pub verified: bool,
}

/// Counts nontrivial zeros between consecutive good Gram points using the classical
/// Gram-block / Turing method, rather than a blind scan.
///
/// Gram's law — `(-1)^n Z(g_n) > 0`, with one zero per Gram interval — usually holds, but
/// fails periodically (first near `n = 126`). Turing's method copes by working with *Gram
/// blocks*: a maximal run `[g_a, g_b]` bounded by good Gram points with only bad ones
/// inside. By *Rosser's rule* a block of `b - a` intervals contains exactly `b - a` zeros.
/// This routine forms those blocks, counts the sign changes of `Z` in each, checks them
/// against Rosser's rule, and confirms the total equals `upper_index - lower_index` — the
/// rigorous prediction `N(g_b) - N(g_a)` since `S` vanishes at good Gram points.
///
/// Returns `None` when `t_max` is below the first Gram interval. Note this is the *method*,
/// not a certified proof: rigorous bounds would require interval arithmetic, and Rosser's
/// rule itself has (much higher) exceptions this does not special-case.
///
/// # Examples
/// ```
/// use riemannrho::{count_zeros_gram, Precision};
/// let report = count_zeros_gram(100.0, Precision::Order2).unwrap();
/// assert!(report.verified);
/// assert_eq!(report.count, report.expected);
/// ```
pub fn count_zeros_gram(t_max: f64, precision: Precision) -> Option<GramCount> {
    if !t_max.is_finite() || t_max < gram_point(0) {
        return None;
    }

    // Collect Gram nodes (index, g_n, is_good) with g_n <= t_max, starting at n = -1.
    let mut nodes: Vec<(i64, f64, bool)> = Vec::new();
    let mut n: i64 = -1;
    loop {
        let g = gram_point(n);
        if g > t_max {
            break;
        }
        nodes.push((n, g, gram_point_is_good(n, precision)));
        n += 1;
    }

    // Turing's method anchors counts at good Gram points: trim bad ones off both ends.
    while nodes.last().is_some_and(|&(_, _, good)| !good) {
        nodes.pop();
    }
    while nodes.first().is_some_and(|&(_, _, good)| !good) {
        nodes.remove(0);
    }
    if nodes.len() < 2 {
        return None;
    }

    let (lower_index, lower, _) = nodes[0];
    let (upper_index, upper, _) = *nodes.last().unwrap();

    let mut count = 0usize;
    let mut gram_blocks = 0usize;
    let mut rosser_violations = 0usize;

    // Walk the good Gram points; each consecutive pair bounds one Gram block.
    let good: Vec<(i64, f64)> = nodes
        .iter()
        .filter(|&&(_, _, g)| g)
        .map(|&(idx, gp, _)| (idx, gp))
        .collect();
    for pair in good.windows(2) {
        let (a_idx, a_g) = pair[0];
        let (b_idx, b_g) = pair[1];
        let length = (b_idx - a_idx) as usize;

        // Count sign changes of Z in the block, at a finer resolution since anomalies
        // (close zero pairs) cluster where Gram's law fails.
        let mut found = 0usize;
        scan_zeros(a_g, b_g, precision, 2.0 * DEFAULT_SCAN_RESOLUTION, |_| {
            found += 1;
            true
        });

        count += found;
        if length >= 2 {
            gram_blocks += 1;
        }
        if found != length {
            rosser_violations += 1;
        }
    }

    let gram_law_failures = nodes.iter().filter(|&&(_, _, g)| !g).count();
    let expected = (upper_index - lower_index) as usize;

    Some(GramCount {
        lower,
        upper,
        lower_index,
        upper_index,
        count,
        expected,
        gram_law_failures,
        gram_blocks,
        rosser_violations,
        verified: count == expected && rosser_violations == 0,
    })
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

/// Collects every zero of `Z(t)` in `(start, end]`, in increasing order, spreading the
/// scan across worker threads for large ranges.
///
/// The range is tiled into contiguous chunks `(a_i, b_i]` with `a_{i+1} = b_i`, so the
/// per-chunk results partition the zeros with no overlap. Small ranges (or single-core
/// systems) fall back to a serial scan to avoid thread overhead.
fn collect_zeros(start: f64, end: f64, precision: Precision, resolution: f64) -> Vec<f64> {
    let span = end - start;
    let threads = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
        .min(8);

    if threads <= 1 || span < 200.0 {
        let mut zeros = Vec::new();
        scan_zeros(start, end, precision, resolution, |root| {
            zeros.push(root);
            true
        });
        return zeros;
    }

    let chunk = span / threads as f64;
    let handles: Vec<_> = (0..threads)
        .map(|i| {
            let a = start + i as f64 * chunk;
            let b = if i == threads - 1 {
                end
            } else {
                start + (i + 1) as f64 * chunk
            };
            std::thread::spawn(move || {
                let mut zeros = Vec::new();
                scan_zeros(a, b, precision, resolution, |root| {
                    zeros.push(root);
                    true
                });
                zeros
            })
        })
        .collect();

    let mut all = Vec::new();
    for handle in handles {
        if let Ok(part) = handle.join() {
            all.extend(part);
        }
    }
    all.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
    all
}

/// Locates every zero of `Z(t)` with imaginary part in `(0, t_max]`, in increasing order.
///
/// Cost is `O(t_max * sqrt(t_max))`; the scan is parallelized across CPU cores for large
/// `t_max`, but this is still intended for exploration over moderate heights rather than
/// astronomically large `t_max`.
///
/// # Examples
/// ```
/// use riemannrho::{zeros_below, Precision};
/// // Five zeros lie below t = 35.
/// assert_eq!(zeros_below(35.0, Precision::Order2).len(), 5);
/// ```
pub fn zeros_below(t_max: f64, precision: Precision) -> Vec<f64> {
    collect_zeros(1.0, t_max, precision, DEFAULT_SCAN_RESOLUTION)
}

/// Counts the nontrivial zeros with imaginary part in `(0, t_max]` found on the
/// critical line. Compare with `round(`[`expected_zero_count`]`(t_max))`.
///
/// # Examples
/// ```
/// use riemannrho::{count_zeros_below, Precision};
/// // There are 10 nontrivial zeros with 0 < t <= 50.
/// assert_eq!(count_zeros_below(50.0, Precision::Order2), 10);
/// ```
pub fn count_zeros_below(t_max: f64, precision: Precision) -> usize {
    collect_zeros(1.0, t_max, precision, DEFAULT_SCAN_RESOLUTION).len()
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
        let c = collect_zeros(1.0, t_max, precision, resolution).len();
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
///
/// # Examples
/// ```
/// use riemannrho::{nth_zero, Precision};
/// assert!((nth_zero(1, Precision::Order2).unwrap() - 14.134725).abs() < 1e-3);
/// assert_eq!(nth_zero(0, Precision::Order2), None);
/// ```
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
    fn nth_zero_matches_higher_reference_zeros() {
        // Validate the whole pipeline well above the first handful of zeros, where the
        // main sum has several terms (nu > 1). Ordinates from Odlyzko's tables.
        let references = [
            (10u64, 49.7738324777),
            (50, 143.1118458076),
            (100, 236.5242296658),
        ];
        for (n, want) in references {
            let got = nth_zero(n, Precision::Order2).unwrap();
            assert!(
                (got - want).abs() < 1e-5,
                "zero #{n}: got {got}, expected ~{want}"
            );
        }
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

    #[test]
    fn parallel_scan_is_ordered_and_complete() {
        // Exercises the multi-threaded path (span > 200) and the ordered partition.
        let zeros = zeros_below(400.0, Precision::Order2);
        for w in zeros.windows(2) {
            assert!(w[0] < w[1], "zeros not strictly increasing: {w:?}");
        }
        // Implied S must be a small fluctuation, i.e. nothing dropped or duplicated.
        let s = zeros.len() as f64 - expected_zero_count(400.0);
        assert!(s.abs() < 1.5, "implausible S over parallel range: {s}");
    }

    #[test]
    fn gram_points_match_known_values() {
        // Known Gram points g_0..g_5 (Haselgrove / standard tables).
        let known = [
            (0i64, 17.8455995405),
            (1, 23.1702827012),
            (2, 27.6701822178),
            (3, 31.7179799547),
            (4, 35.4671842971),
            (5, 38.9992099640),
        ];
        for (n, want) in known {
            let got = gram_point(n);
            assert!(
                (got - want).abs() < 1e-4,
                "g_{n}: got {got}, expected ~{want}"
            );
        }
        assert!(gram_point(-2).is_nan());
    }

    #[test]
    fn grams_law_mostly_holds() {
        // For low n, (-1)^n Z(g_n) > 0 (Gram's law) and each interval holds one zero.
        for n in 0..=20i64 {
            let g = gram_point(n);
            let sign = if n % 2 == 0 { 1.0 } else { -1.0 };
            assert!(
                sign * z_func(g, Precision::Order2) > 0.0,
                "Gram's law violated at n={n}"
            );
        }
    }

    #[test]
    fn gram_count_is_verified_and_matches_a_blind_scan() {
        for &t in &[50.0, 100.0, 300.0] {
            let gc = count_zeros_gram(t, Precision::Order2).expect("range above first Gram point");
            assert!(
                gc.verified,
                "T={t}: count {} vs expected {} ({} Rosser violations)",
                gc.count, gc.expected, gc.rosser_violations
            );
            // The Gram-block count must agree with an independent parallel scan over the
            // same (lower, upper] range.
            let scanned = collect_zeros(
                gc.lower,
                gc.upper,
                Precision::Order2,
                DEFAULT_SCAN_RESOLUTION,
            )
            .len();
            assert_eq!(
                scanned, gc.count,
                "T={t}: blind scan {scanned} vs gram {}",
                gc.count
            );
        }
    }

    #[test]
    fn gram_count_resolves_gram_blocks_via_rosser() {
        // Gram's law first fails near n = 126 (t ~ 282), so scanning to 300 must form a
        // Gram block that Rosser's rule still resolves to the correct, verified count.
        let gc = count_zeros_gram(300.0, Precision::Order2).unwrap();
        assert!(
            gc.gram_law_failures > 0,
            "expected a Gram's-law failure below 300"
        );
        assert!(gc.gram_blocks > 0, "expected a non-trivial Gram block");
        assert!(gc.rosser_violations == 0 && gc.verified);
    }

    #[test]
    fn gram_count_rejects_tiny_range() {
        assert!(count_zeros_gram(5.0, Precision::Base).is_none());
    }

    #[test]
    fn chebyshev_psi_matches_closed_form() {
        assert_eq!(chebyshev_psi(1.0), 0.0);
        assert!((chebyshev_psi(2.0) - 2f64.ln()).abs() < 1e-12);
        // psi(10) sums ln p over prime powers <= 10: 2,4,8 (3*ln2), 3,9 (2*ln3), 5, 7.
        let want = 3.0 * 2f64.ln() + 2.0 * 3f64.ln() + 5f64.ln() + 7f64.ln();
        assert!((chebyshev_psi(10.0) - want).abs() < 1e-12);
    }

    #[test]
    fn explicit_formula_reconstructs_psi_from_zeros() {
        // Collect the first ~300 positive zeros and rebuild psi(x) from them. Away from
        // prime-power jumps the truncated series is accurate to well under 1.
        let height = estimate_t(300.0) * 1.05 + 10.0;
        let mut gammas = zeros_below(height, Precision::Order2);
        gammas.truncate(300);
        assert!(
            gammas.len() >= 300,
            "needed 300 zeros, got {}",
            gammas.len()
        );

        for &x in &[20.5, 30.5, 50.5] {
            let approx = psi_from_zeros(x, &gammas);
            let actual = chebyshev_psi(x);
            assert!(
                (approx - actual).abs() < 0.5,
                "psi({x}): zeros gave {approx}, actual {actual}"
            );
        }
        assert!(psi_from_zeros(1.0, &gammas).is_nan());
    }

    #[test]
    fn explicit_formula_sharpens_with_more_zeros() {
        // At a fixed flat point, more zeros reduce the reconstruction error on average.
        let height = estimate_t(400.0) * 1.05 + 10.0;
        let all = zeros_below(height, Precision::Order2);
        let x = 50.5;
        let actual = chebyshev_psi(x);
        let err = |n: usize| (psi_from_zeros(x, &all[..n.min(all.len())]) - actual).abs();
        assert!(
            err(300) < err(20),
            "more zeros should help: {} vs {}",
            err(300),
            err(20)
        );
    }

    #[test]
    fn normalized_spacings_have_unit_mean() {
        let zeros = zeros_below(1500.0, Precision::Order2);
        let spac = normalized_spacings(&zeros);
        assert!(
            spac.len() > 1000,
            "expected many spacings, got {}",
            spac.len()
        );
        let mean = spac.iter().sum::<f64>() / spac.len() as f64;
        assert!((mean - 1.0).abs() < 0.05, "mean normalized spacing {mean}");
    }

    #[test]
    fn zeros_exhibit_gue_level_repulsion() {
        // The Montgomery-Odlyzko law: spacings follow GUE statistics, which repel — far
        // fewer tiny gaps than a Poisson process (where 1 - e^-0.5 = 39% fall below 0.5).
        let zeros = zeros_below(1500.0, Precision::Order2);
        let spac = normalized_spacings(&zeros);
        let n = spac.len() as f64;
        let frac_below = |t: f64| spac.iter().filter(|&&x| x < t).count() as f64 / n;
        assert!(
            frac_below(0.5) < 0.25,
            "too many small gaps: {}",
            frac_below(0.5)
        );
        assert!(
            frac_below(0.1) < 0.03,
            "unexpected near-coincidences: {}",
            frac_below(0.1)
        );
    }

    #[test]
    fn wigner_surmise_is_normalized_and_repels() {
        // P(0) = 0 (level repulsion), and the distribution integrates to 1 with mean 1.
        assert_eq!(wigner_surmise(0.0), 0.0);
        let (mut area, mut mean) = (0.0, 0.0);
        let h = 1e-3;
        let mut s = h / 2.0;
        while s < 6.0 {
            let p = wigner_surmise(s);
            area += p * h;
            mean += s * p * h;
            s += h;
        }
        assert!((area - 1.0).abs() < 1e-2, "integral {area}");
        assert!((mean - 1.0).abs() < 1e-2, "mean {mean}");
    }
}
