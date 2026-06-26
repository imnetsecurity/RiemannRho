//! Random Matrix Theory (RMT) diagnostics for an *arbitrary* real spectrum.
//!
//! The Montgomery-Odlyzko law makes the zeta zeros a textbook example of GUE statistics,
//! but the same diagnostics apply to any sorted list of levels — nuclear energy levels,
//! quantum-billiard spectra, MIMO channel eigenvalues, the singular values of a neural
//! network's weight matrix, a cleaned financial covariance spectrum, and so on. This module
//! takes a generic spectrum and answers: does it look uncorrelated (Poisson / integrable)
//! or correlated/repelling (GOE/GUE/GSE / chaotic)?
//!
//! Tools provided:
//! * the **spacing-ratio statistic** `<r>` (Atas et al.) — needs no unfolding, so it is the
//!   robust first thing to reach for;
//! * **polynomial unfolding** to a mean spacing of 1;
//! * the **number variance** `Sigma^2(L)`;
//! * reference nearest-neighbor distributions and `Sigma^2` curves for each ensemble.

use std::fmt;

/// A random-matrix universality class (plus Poisson for uncorrelated spectra).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Ensemble {
    /// Uncorrelated levels (integrable systems): no level repulsion.
    Poisson,
    /// Gaussian Orthogonal Ensemble (time-reversal-symmetric chaotic systems).
    Goe,
    /// Gaussian Unitary Ensemble (broken time-reversal symmetry) — the zeta zeros.
    Gue,
    /// Gaussian Symplectic Ensemble (spin-1/2, time-reversal symmetric).
    Gse,
}

impl fmt::Display for Ensemble {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            Ensemble::Poisson => "Poisson (uncorrelated)",
            Ensemble::Goe => "GOE",
            Ensemble::Gue => "GUE",
            Ensemble::Gse => "GSE",
        };
        f.write_str(s)
    }
}

/// Mean of the consecutive spacing ratio `<r>` for each ensemble (Atas et al. 2013).
/// These are the reference values [`ensemble_from_ratio`] classifies against.
pub const RATIO_POISSON: f64 = 0.386_294; // 2 ln 2 - 1
pub const RATIO_GOE: f64 = 0.535_898;
pub const RATIO_GUE: f64 = 0.602_658;
pub const RATIO_GSE: f64 = 0.676_169;

/// Nearest-neighbor spacing distribution `P(s)` for the Poisson case: `e^{-s}` (no
/// repulsion, `P(0) = 1`).
pub fn nnsd_poisson(s: f64) -> f64 {
    (-s).exp()
}

/// Wigner surmise for the GOE: `(pi/2) s exp(-pi s^2/4)`.
pub fn nnsd_goe(s: f64) -> f64 {
    use std::f64::consts::PI;
    (PI / 2.0) * s * (-PI * s * s / 4.0).exp()
}

/// Wigner surmise for the GUE: `(32/pi^2) s^2 exp(-4 s^2/pi)`.
pub fn nnsd_gue(s: f64) -> f64 {
    use std::f64::consts::PI;
    (32.0 / (PI * PI)) * s * s * (-4.0 * s * s / PI).exp()
}

/// Wigner surmise for the GSE: `(2^18/(3^6 pi^3)) s^4 exp(-(64/9pi) s^2)`.
pub fn nnsd_gse(s: f64) -> f64 {
    use std::f64::consts::PI;
    let coeff = 262_144.0 / (729.0 * PI * PI * PI);
    coeff * s.powi(4) * (-(64.0 / (9.0 * PI)) * s * s).exp()
}

/// Reference nearest-neighbor distribution for `ensemble`, evaluated at `s`.
pub fn nnsd(ensemble: Ensemble, s: f64) -> f64 {
    match ensemble {
        Ensemble::Poisson => nnsd_poisson(s),
        Ensemble::Goe => nnsd_goe(s),
        Ensemble::Gue => nnsd_gue(s),
        Ensemble::Gse => nnsd_gse(s),
    }
}

/// Classifies a spectrum by the nearest reference value of the mean spacing ratio `<r>`.
///
/// # Examples
/// ```
/// use riemannrho::rmt::{ensemble_from_ratio, Ensemble};
/// assert_eq!(ensemble_from_ratio(0.60), Ensemble::Gue);
/// assert_eq!(ensemble_from_ratio(0.39), Ensemble::Poisson);
/// ```
pub fn ensemble_from_ratio(r: f64) -> Ensemble {
    let candidates = [
        (Ensemble::Poisson, RATIO_POISSON),
        (Ensemble::Goe, RATIO_GOE),
        (Ensemble::Gue, RATIO_GUE),
        (Ensemble::Gse, RATIO_GSE),
    ];
    candidates
        .into_iter()
        .min_by(|a, b| {
            (a.1 - r)
                .abs()
                .partial_cmp(&(b.1 - r).abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(e, _)| e)
        .unwrap()
}

/// Consecutive spacing ratios `r_i = min(s_i, s_{i+1}) / max(s_i, s_{i+1})` of a sorted
/// spectrum, where `s_i = levels[i+1] - levels[i]`.
///
/// The ratio is invariant to the local density, so — unlike the spacings themselves — it
/// needs **no unfolding**, which makes it the most robust ensemble discriminator.
pub fn spacing_ratios(levels: &[f64]) -> Vec<f64> {
    let gaps: Vec<f64> = levels.windows(2).map(|w| w[1] - w[0]).collect();
    gaps.windows(2)
        .filter_map(|g| {
            let (a, b) = (g[0], g[1]);
            let max = a.max(b);
            if max > 0.0 {
                Some(a.min(b) / max)
            } else {
                None
            }
        })
        .collect()
}

/// Mean consecutive spacing ratio `<r>` (see [`spacing_ratios`]). `NaN` if undefined.
pub fn mean_spacing_ratio(levels: &[f64]) -> f64 {
    let r = spacing_ratios(levels);
    if r.is_empty() {
        return f64::NAN;
    }
    r.iter().sum::<f64>() / r.len() as f64
}

/// Solves the linear system `m c = b` (square, row-major) by Gaussian elimination with
/// partial pivoting. Returns `None` if singular.
fn solve(mut m: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = b.len();
    for col in 0..n {
        let pivot = (col..n).max_by(|&a, &c| {
            m[a][col]
                .abs()
                .partial_cmp(&m[c][col].abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        })?;
        if m[pivot][col].abs() < 1e-12 {
            return None;
        }
        m.swap(col, pivot);
        b.swap(col, pivot);
        let pivot_row = m[col].clone();
        let pivot_b = b[col];
        for row in (col + 1)..n {
            let factor = m[row][col] / pivot_row[col];
            for (mk, &pk) in m[row].iter_mut().zip(pivot_row.iter()).skip(col) {
                *mk -= factor * pk;
            }
            b[row] -= factor * pivot_b;
        }
    }
    let mut x = vec![0.0; n];
    for row in (0..n).rev() {
        let mut s = b[row];
        for k in (row + 1)..n {
            s -= m[row][k] * x[k];
        }
        x[row] = s / m[row][row];
    }
    Some(x)
}

/// Unfolds a sorted spectrum to a mean spacing of 1 by fitting the integrated density
/// (the staircase `N(E)`) with a degree-`degree` polynomial and mapping each level through
/// it.
///
/// Unfolding removes the smooth, system-specific variation in level density so that the
/// remaining fluctuations can be compared with universal RMT predictions. The fit is done
/// in a centered/scaled coordinate for numerical stability. Returns the unfolded levels
/// (same length), or an empty vector if the fit is degenerate.
pub fn unfold(levels: &[f64], degree: usize) -> Vec<f64> {
    let n = levels.len();
    if n < degree + 1 {
        return Vec::new();
    }
    let mean = levels.iter().sum::<f64>() / n as f64;
    let var = levels.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / n as f64;
    let scale = var.sqrt();
    if scale <= 0.0 {
        return Vec::new();
    }
    let u: Vec<f64> = levels.iter().map(|&x| (x - mean) / scale).collect();
    // Staircase value at each level (number of levels at or below it).
    let y: Vec<f64> = (0..n).map(|i| i as f64).collect();

    // Normal equations for least-squares polynomial fit: M[j][k] = sum u^(j+k).
    let d = degree;
    let mut power_sums = vec![0.0; 2 * d + 1];
    for &ui in &u {
        let mut p = 1.0;
        for ps in power_sums.iter_mut() {
            *ps += p;
            p *= ui;
        }
    }
    let mut rhs = vec![0.0; d + 1];
    for (yi, &ui) in y.iter().zip(&u) {
        let mut p = 1.0;
        for r in rhs.iter_mut() {
            *r += yi * p;
            p *= ui;
        }
    }
    let matrix: Vec<Vec<f64>> = (0..=d)
        .map(|j| (0..=d).map(|k| power_sums[j + k]).collect())
        .collect();
    let coeffs = match solve(matrix, rhs) {
        Some(c) => c,
        None => return Vec::new(),
    };

    u.iter()
        .map(|&ui| {
            let mut p = 1.0;
            let mut acc = 0.0;
            for &c in &coeffs {
                acc += c * p;
                p *= ui;
            }
            acc
        })
        .collect()
}

/// The number variance `Sigma^2(L)`: the variance of the count of unfolded levels in a
/// window of length `L`, averaged over window position.
///
/// For uncorrelated levels `Sigma^2(L) = L`; correlated (RMT) spectra are far more rigid,
/// growing only logarithmically. `unfolded` must be sorted ascending and have unit mean
/// spacing (e.g. from [`unfold`]). Returns `NaN` if the range is shorter than `L`.
pub fn number_variance(unfolded: &[f64], l: f64) -> f64 {
    if unfolded.len() < 2 || l <= 0.0 || l.is_nan() {
        return f64::NAN;
    }
    let lo = unfolded[0];
    let hi = unfolded[unfolded.len() - 1];
    if hi - lo <= l {
        return f64::NAN;
    }
    let samples = 1000usize;
    let start_max = hi - l;
    let mut sum = 0.0;
    let mut sumsq = 0.0;
    for k in 0..samples {
        let start = lo + (start_max - lo) * k as f64 / (samples - 1) as f64;
        let upper = unfolded.partition_point(|&x| x < start + l);
        let lower = unfolded.partition_point(|&x| x < start);
        let c = (upper - lower) as f64;
        sum += c;
        sumsq += c * c;
    }
    let m = samples as f64;
    let mean = sum / m;
    sumsq / m - mean * mean
}

/// `Sigma^2(L)` for the Poisson case: `L`.
pub fn number_variance_poisson(l: f64) -> f64 {
    l
}

/// Large-`L` `Sigma^2(L)` for the GUE: `(1/pi^2)(ln(2 pi L) + gamma + 1)`.
pub fn number_variance_gue(l: f64) -> f64 {
    use std::f64::consts::PI;
    const EULER_GAMMA: f64 = 0.577_215_664_901_532_9;
    (1.0 / (PI * PI)) * ((2.0 * PI * l).ln() + EULER_GAMMA + 1.0)
}

/// Large-`L` `Sigma^2(L)` for the GOE: `(2/pi^2)(ln(2 pi L) + gamma + 1 - pi^2/8)`.
pub fn number_variance_goe(l: f64) -> f64 {
    use std::f64::consts::PI;
    const EULER_GAMMA: f64 = 0.577_215_664_901_532_9;
    (2.0 / (PI * PI)) * ((2.0 * PI * l).ln() + EULER_GAMMA + 1.0 - PI * PI / 8.0)
}

/// Summary of a spectrum's RMT diagnostics (see [`analyze`]).
#[derive(Clone, Copy, Debug)]
pub struct RmtReport {
    /// Number of levels analyzed.
    pub levels: usize,
    /// Mean consecutive spacing ratio `<r>`.
    pub mean_ratio: f64,
    /// Ensemble classified from `<r>`.
    pub ensemble: Ensemble,
}

/// Runs the unfolding-free spacing-ratio diagnostic and classifies the spectrum.
///
/// # Examples
/// ```
/// use riemannrho::rmt::{analyze, Ensemble};
/// // An equally spaced ("picket-fence") spectrum is perfectly rigid: <r> = 1.
/// let levels: Vec<f64> = (0..100).map(|i| i as f64).collect();
/// let report = analyze(&levels).unwrap();
/// assert!((report.mean_ratio - 1.0).abs() < 1e-9);
/// ```
pub fn analyze(levels: &[f64]) -> Option<RmtReport> {
    if levels.len() < 3 {
        return None;
    }
    let mean_ratio = mean_spacing_ratio(levels);
    if !mean_ratio.is_finite() {
        return None;
    }
    Some(RmtReport {
        levels: levels.len(),
        mean_ratio,
        ensemble: ensemble_from_ratio(mean_ratio),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{zeros_below, Precision};

    #[test]
    fn ratio_classifies_reference_values() {
        assert_eq!(ensemble_from_ratio(RATIO_POISSON), Ensemble::Poisson);
        assert_eq!(ensemble_from_ratio(RATIO_GOE), Ensemble::Goe);
        assert_eq!(ensemble_from_ratio(RATIO_GUE), Ensemble::Gue);
        assert_eq!(ensemble_from_ratio(RATIO_GSE), Ensemble::Gse);
    }

    #[test]
    fn equally_spaced_spectrum_is_rigid() {
        let levels: Vec<f64> = (0..200).map(|i| i as f64 * 2.5).collect();
        // All gaps equal, so every ratio is 1.
        assert!((mean_spacing_ratio(&levels) - 1.0).abs() < 1e-12);
        // And the number variance of the unit-spacing picket fence stays tiny.
        let unfolded: Vec<f64> = (0..200).map(|i| i as f64).collect();
        assert!(number_variance(&unfolded, 5.0) < 0.5);
    }

    #[test]
    fn zeta_zeros_classify_as_gue() {
        // The Montgomery-Odlyzko law: the zeros' spacings are GUE.
        let zeros = zeros_below(2000.0, Precision::Order2);
        let r = mean_spacing_ratio(&zeros);
        assert!((0.56..0.64).contains(&r), "ratio {r} not GUE-like");
        assert_eq!(ensemble_from_ratio(r), Ensemble::Gue);
    }

    #[test]
    fn zeta_zeros_are_more_rigid_than_poisson() {
        // Unfold the zeros exactly via the smooth count, then Sigma^2 should be far below
        // the Poisson value L (level rigidity), but positive.
        let zeros = zeros_below(2000.0, Precision::Order2);
        let mut unfolded = Vec::with_capacity(zeros.len());
        let mut acc = 0.0;
        unfolded.push(0.0);
        for w in zeros.windows(2) {
            acc += (crate::theta(w[1]) - crate::theta(w[0])) / std::f64::consts::PI;
            unfolded.push(acc);
        }
        let s2 = number_variance(&unfolded, 3.0);
        assert!(s2 > 0.1 && s2 < 3.0, "Sigma^2(3) = {s2} not sub-Poisson");
    }

    #[test]
    fn unfolding_removes_density_trend() {
        // Start rigid (spacing 1 in u), warp by a smooth nonlinear density, then unfold:
        // the recovered spectrum should be rigid again (mean spacing 1, small variance).
        let warped: Vec<f64> = (0..300)
            .map(|i| {
                let u = i as f64;
                u + 0.0005 * u * u // smooth, monotonically increasing density change
            })
            .collect();
        let unfolded = unfold(&warped, 4);
        assert_eq!(unfolded.len(), warped.len());
        let gaps: Vec<f64> = unfolded.windows(2).map(|w| w[1] - w[0]).collect();
        let mean = gaps.iter().sum::<f64>() / gaps.len() as f64;
        let var = gaps.iter().map(|&g| (g - mean).powi(2)).sum::<f64>() / gaps.len() as f64;
        assert!((mean - 1.0).abs() < 0.05, "unfolded mean spacing {mean}");
        assert!(var < 0.05, "unfolded spacing variance {var} too large");
    }
}
