//! Prime-distribution and primality tooling — the honest, practical bridge from the zeta
//! zeros to cryptography.
//!
//! The Riemann hypothesis (and its generalization, GRH) is *not* a cryptographic primitive
//! — the zeros are not keys. Its real relevance is twofold, both captured here:
//!
//! * **Deterministic primality.** Under GRH a composite `n` has a Miller-Rabin witness
//!   below `2 (ln n)^2` (Bach's bound), turning Miller-Rabin into a *deterministic* test —
//!   the foundation of prime generation for RSA/DH keys. (For `u64` a fixed witness set is
//!   provably sufficient *unconditionally*, which is what [`is_prime`] uses in practice.)
//! * **Prime distribution.** RH pins down how tightly the primes follow `li(x)`: the error
//!   `|pi(x) - li(x)|` is bounded by `(1/8pi) sqrt(x) ln x` (Schoenfeld), which underlies
//!   reasoning about prime density and key sizes.

/// `(a * b) mod m`, computed via `u128` to avoid overflow for `u64` operands.
fn mulmod(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 * b as u128) % m as u128) as u64
}

/// `base^exp mod m` by square-and-multiply.
fn powmod(mut base: u64, mut exp: u64, m: u64) -> u64 {
    if m == 1 {
        return 0;
    }
    let mut result = 1u64;
    base %= m;
    while exp > 0 {
        if exp & 1 == 1 {
            result = mulmod(result, base, m);
        }
        base = mulmod(base, base, m);
        exp >>= 1;
    }
    result
}

/// Whether `n` passes the Miller-Rabin test for the single base `a` (i.e. `a` is *not* a
/// witness to `n`'s compositeness). Assumes `n` is odd and `2 <= a < n`.
fn miller_rabin_passes(n: u64, a: u64) -> bool {
    let mut d = n - 1;
    let mut s = 0u32;
    while d.is_multiple_of(2) {
        d >>= 1;
        s += 1;
    }
    let mut x = powmod(a, d, n);
    if x == 1 || x == n - 1 {
        return true;
    }
    for _ in 0..s.saturating_sub(1) {
        x = mulmod(x, x, n);
        if x == n - 1 {
            return true;
        }
    }
    false
}

/// Witness set that makes Miller-Rabin deterministic for every `u64` (in fact for all
/// `n < 3.3e24`).
const U64_WITNESSES: [u64; 12] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];

/// Deterministic primality test for `u64`, unconditional (no reliance on GRH).
///
/// Uses Miller-Rabin with a fixed witness set proven sufficient for all 64-bit integers,
/// so the answer is exact — the workhorse for generating cryptographic primes.
///
/// # Examples
/// ```
/// use riemannrho::primes::is_prime;
/// assert!(is_prime(2_147_483_647));   // the Mersenne prime 2^31 - 1
/// assert!(!is_prime(561));            // a Carmichael number: fools the Fermat test
/// ```
pub fn is_prime(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for &p in &U64_WITNESSES {
        if n == p {
            return true;
        }
        if n.is_multiple_of(p) {
            return false;
        }
    }
    U64_WITNESSES
        .iter()
        .all(|&a| a >= n || miller_rabin_passes(n, a))
}

/// Result of the GRH-conditional deterministic primality test (see [`is_prime_grh`]).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct GrhPrimality {
    /// Verdict (conditional on the Generalized Riemann Hypothesis).
    pub is_prime: bool,
    /// Bach's bound `floor(2 (ln n)^2)`: the largest base that needs testing.
    pub bound: u64,
    /// Number of bases actually tested.
    pub bases_tested: u64,
}

/// Deterministic primality test conditional on the **Generalized Riemann Hypothesis**.
///
/// Under GRH, a composite `n` has a Miller-Rabin witness `a` with `2 <= a <= 2 (ln n)^2`
/// (Bach), so testing every base up to that bound decides primality. This is the direct
/// cryptographic payoff of the Riemann hypothesis; compare it with the unconditional
/// [`is_prime`] (which they always agree with for `u64`).
///
/// # Examples
/// ```
/// use riemannrho::primes::is_prime_grh;
/// let r = is_prime_grh(1_000_003);
/// assert!(r.is_prime);
/// assert!(r.bases_tested <= r.bound);
/// ```
pub fn is_prime_grh(n: u64) -> GrhPrimality {
    if n < 2 {
        return GrhPrimality {
            is_prime: false,
            bound: 0,
            bases_tested: 0,
        };
    }
    if n == 2 || n == 3 {
        return GrhPrimality {
            is_prime: true,
            bound: 0,
            bases_tested: 0,
        };
    }
    if n.is_multiple_of(2) {
        return GrhPrimality {
            is_prime: false,
            bound: 0,
            bases_tested: 1,
        };
    }
    let ln_n = (n as f64).ln();
    let bound = (2.0 * ln_n * ln_n).floor() as u64;
    let bound = bound.min(n - 2).max(2);

    let mut is_prime = true;
    let mut bases_tested = 0u64;
    for a in 2..=bound {
        bases_tested += 1;
        if !miller_rabin_passes(n, a) {
            is_prime = false;
            break;
        }
    }
    GrhPrimality {
        is_prime,
        bound,
        bases_tested,
    }
}

/// The prime-counting function `pi(x)`: the number of primes `<= x`, by sieving.
///
/// Intended for moderate `x` (it allocates a sieve of size `x`).
///
/// # Examples
/// ```
/// use riemannrho::primes::prime_pi;
/// assert_eq!(prime_pi(10), 4);    // 2, 3, 5, 7
/// assert_eq!(prime_pi(100), 25);
/// ```
pub fn prime_pi(x: u64) -> u64 {
    if x < 2 {
        return 0;
    }
    let n = x as usize;
    let mut is_composite = vec![false; n + 1];
    let mut count = 0u64;
    for p in 2..=n {
        if !is_composite[p] {
            count += 1;
            let mut m = p * p;
            while m <= n {
                is_composite[m] = true;
                m += p;
            }
        }
    }
    count
}

/// The logarithmic integral `li(x)`, the Riemann-hypothesis main term for `pi(x)`.
///
/// Evaluated via the convergent series `li(x) = gamma + ln(ln x) + sum_{k>=1} (ln x)^k /
/// (k * k!)`. Returns `NaN` for `x <= 1`.
pub fn logarithmic_integral(x: f64) -> f64 {
    if !x.is_finite() || x <= 1.0 {
        return f64::NAN;
    }
    const EULER_GAMMA: f64 = 0.577_215_664_901_532_9;
    let ln_x = x.ln();
    let mut sum = 0.0;
    let mut term = 1.0; // (ln x)^k / k!  for k starting at 1
    for k in 1..1000 {
        term *= ln_x / k as f64; // now (ln x)^k / k!
        let contribution = term / k as f64;
        sum += contribution;
        if k as f64 > ln_x && contribution < 1e-15 * sum.abs() {
            break;
        }
    }
    EULER_GAMMA + ln_x.ln() + sum
}

/// Schoenfeld's RH-conditional bound on the prime-counting error:
/// `|pi(x) - li(x)| < (1/(8 pi)) sqrt(x) ln x`.
///
/// Proven for `x >= 2657` (the formula is returned for any `x > 1`). This is the precise
/// sense in which the Riemann hypothesis controls the distribution of the primes.
pub fn rh_prime_count_bound(x: f64) -> f64 {
    use std::f64::consts::PI;
    if !x.is_finite() || x <= 1.0 {
        return f64::NAN;
    }
    (1.0 / (8.0 * PI)) * x.sqrt() * x.ln()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn is_prime_handles_edges_carmichaels_and_large() {
        for n in [0u64, 1, 4, 100, 561, 1105, 1729, 41041, 4_294_967_295] {
            assert!(!is_prime(n), "{n} should be composite");
        }
        for n in [2u64, 3, 97, 7919, 1_000_003, 2_147_483_647, 4_294_967_291] {
            assert!(is_prime(n), "{n} should be prime");
        }
        // Stress the u64 range and overflow path: 2^64 - 59 is prime, 2^64 - 1 is not.
        assert!(is_prime(18_446_744_073_709_551_557));
        assert!(!is_prime(18_446_744_073_709_551_615));
    }

    #[test]
    fn grh_test_agrees_with_unconditional() {
        for n in [2u64, 17, 561, 7919, 104_729, 1_000_003, 1_000_000] {
            assert_eq!(is_prime_grh(n).is_prime, is_prime(n), "disagreement at {n}");
        }
        // Bach's bound stays tiny even for large n (2 (ln n)^2 grows very slowly).
        let r = is_prime_grh(1_000_000_007);
        assert!(r.is_prime && r.bound < 1000, "bound {}", r.bound);
    }

    #[test]
    fn prime_pi_matches_known_values() {
        assert_eq!(prime_pi(1), 0);
        assert_eq!(prime_pi(1000), 168);
        assert_eq!(prime_pi(1_000_000), 78498);
    }

    #[test]
    fn rh_bound_contains_the_actual_li_error() {
        // |pi(x) - li(x)| must sit comfortably inside the RH bound.
        let x = 1_000_000.0;
        let pi = prime_pi(x as u64) as f64;
        let li = logarithmic_integral(x);
        let err = (pi - li).abs();
        let bound = rh_prime_count_bound(x);
        assert!(err < bound, "error {err} exceeded RH bound {bound}");
        // li is a far better estimate than x/ln x.
        assert!(err < (pi - x / x.ln()).abs());
    }
}
