# RiemannRho: Rust Library for Riemann Zeta Nontrivial Zeros Approximation

## About

RiemannRho is a high-performance Rust library and command-line tool dedicated to the numerical approximation of nontrivial zeros of the Riemann zeta function \(\zeta(s)\) on the critical line \(\operatorname{Re}(s) = 1/2\). Utilizing the Riemann-Siegel asymptotic formula, RiemannRho computes Hardy's Z-function \(Z(t)\), facilitating precise location of zeros corresponding to \(\zeta(1/2 + it) = 0\). Engineered for accuracy, efficiency, and extensibility, this tool serves researchers in analytic number theory, supporting investigations into the Riemann hypothesis through scalable computations of high-order zeros.

RiemannRho approximates nontrivial zeros of the Riemann zeta function $\zeta(s)$ on the critical line. These zeros, particularly under the unproven Riemann hypothesis (which posits that all nontrivial zeros have real part 1/2), provide insights into prime number distribution and oscillatory phenomena. Some representative use cases:

1. **Analytic Number Theory and Prime Distribution.** The zeros control the oscillations of prime numbers around their expected positions via explicit formulas connected to the prime number theorem ($\pi(x) \approx \frac{x}{\ln x}$). Computing them helps explore bounds on the error term in prime-counting functions.

2. **Testing the Riemann Hypothesis.** Extensive numerical computation of zeros has so far found no counterexamples. RiemannRho can locate individual zeros and, via `--count T`, verify that the number of zeros found on the critical line up to height `T` matches the theoretical count — the numerical core of testing the hypothesis over a range.

3. **Physics (Quantum Chaos and Random Matrix Theory).** Zero spacings statistically resemble eigenvalue spacings of random matrices, with connections to energy levels in chaotic quantum systems.

4. **Education and Exploration.** The tool computes individual zeros (the first at $\approx 14.1347$) with optional visualizations of $Z(t)$, useful for teaching complex analysis, number theory, and asymptotic methods.

For large-scale, record-breaking computations, dedicated tools such as Odlyzko's remain the reference; RiemannRho targets accessible, dependency-free exploration on ordinary hardware.

Core strengths include:
- **Riemann-Siegel Evaluation**: Computes Hardy's $Z(t)$ with the main sum plus optional remainder correction terms (up to $C_2$) for reduced error.
- **High-Order Correction Mode**: Adds the $C_1$ and $C_2$ terms to noticeably reduce the asymptotic error at moderate $t$ (see [Accuracy](#accuracy)).
- **Visualization**: Optional export of a D3.js plot to HTML, rendering the $Z(t)$ curve with a zero marker.
- **$O(\sqrt{t})$ Main Sum**: Each $Z(t)$ evaluation costs $O(\sqrt{t})$ thanks to the Riemann-Siegel formula. Note that with 64-bit floats, accuracy and runtime degrade well before astronomically large ordinals.

Licensed under MIT, RiemannRho promotes open collaboration in mathematical software, ensuring robust, reproducible results without external dependencies.

## Features

- **Range-Based Zero Detection**: Define intervals \([low, high]\) for bisection-based zero isolation via \(Z(t)\).
- **Ordinal Zero Approximation**: The `--nth` option estimates the \(n\)th zero via the Riemann-von Mangoldt counting formula with Newton's iteration, then scans for and brackets the actual sign change.
- **High-Order Correction Mode**: Invoke `--high-order` to include the $C_1$ and $C_2$ correction terms, reducing the asymptotic error.
- **Browser-Based Plots**: Post-calculation prompt generates an HTML file (default `zeta_plot.html`, configurable with `--out`) with a D3.js visualization: line, axes, and a red zero marker.
- **Flexible Interfaces**: Command-line parameters or interactive prompts, with customizable tolerance for convergence control.
- **Dependency-Free by Default**: The default build uses only Rust's standard library. Arbitrary precision is an opt-in feature (`bigfloat`) that pulls in a single pure-Rust crate.

## Installation

RiemannRho requires Rust (1.70+). Installation steps:

1. Install Rust via [rustup](https://rustup.rs/).
2. Clone the source:
git clone https://github.com/imnetsecurity/riemannrho.git
cd riemannrho
3. Compile:
cargo build --release
Produces an optimized executable at `target/release/riemannrho`.

The tool operates independently, with no runtime dependencies.

## Usage

Execute with arguments or via interactive mode.

### Command Syntax

```
./target/release/riemannrho [low] [high] [tol] [--high-order] [--nth N] [--out FILE]
./target/release/riemannrho --count T [--high-order]
```

- `low`, `high`: Search interval bounds (omit when using `--nth`).
- `tol`: Root-finder bracket-width tolerance (default: 1e-10).
- `--high-order`: Include the $C_1$ and $C_2$ correction terms.
- `--nth N`: Target the Nth zero (`N >= 1`). For `N <= 100000` the zero is found
  *exactly* by scanning sequentially; beyond that an asymptotic estimate of a single
  nearby zero is used.
- `--count T`: Count the zeros with $0 < t \le T$ and compare with the theoretical
  count $\theta(T)/\pi + 1$ (a Turing-flavored consistency check — see below).
- `--list T`: Print every zero with $0 < t \le T$, one per line.
- `--gram N`: Print the first $N$ Gram points and check Gram's law at each.
- `--turing T`: Verify the zero count up to $T$ via Gram blocks and Rosser's rule (see
  [Gram Points & Turing's Method](#gram-points--turings-method)).
- `--digits D`: Locate the zero in arbitrary precision with $D$ decimal digits (requires
  the optional `bigfloat` feature; most useful at large $t$ — see [Arbitrary Precision](#arbitrary-precision)).
- `--out FILE`: Output path for the generated plot (default: `zeta_plot.html`).
- `-h`, `--help`: Print usage.

Upon completion, the tool prints the zero approximation and asks: "Do you want a D3.js visualization? (yes/no)". Answering `yes` writes the plot HTML for browser viewing.

> **Note:** `tol` controls the *root-finding* precision, not the accuracy of the
> underlying asymptotic model. For small `t` the Riemann-Siegel approximation
> itself limits how close the result can get to the true zero.

### Usage Examples

1. **First zero with high-order correction**:
   ```
   ./target/release/riemannrho 14 15 1e-10 --high-order
   ```
   Sample output:
   ```
   Approximate imaginary part of the nontrivial zero: 14.1348228651
   Do you want a D3.js visualization? (yes/no)
   ```
   (True value: 14.1347251417; the residual error is the limit of the asymptotic series at this small `t`.)

2. **Nth zero directly**:
   ```
   ./target/release/riemannrho --nth 1000000 --high-order
   ```
   Yields t ≈ 600269.6770 (the millionth zero).

3. **Interactive mode** (no arguments):
   ```
   ./target/release/riemannrho
   ```
   Prompts for `low`, `high`, and `tol`.

4. **Count and verify zeros up to a height**:
   ```
   ./target/release/riemannrho --count 100 --high-order
   ```
   Sample output:
   ```
   Zeros found with 0 < t <= 100: 29
   Smooth estimate theta(T)/pi + 1: 29.0024
   Implied S(T) = found - estimate: -0.0024  (true count = estimate + S)
   Consistent: S(T) is a normal small fluctuation, so every zero up to T was
   found on the critical line (Turing-flavored check passes).
   ```

5. **List every zero up to a height**:
   ```
   ./target/release/riemannrho --list 35 --high-order
   ```
   Prints each zero with its index, then the total.

6. **Gram points and Gram's law**:
   ```
   ./target/release/riemannrho --gram 10 --high-order
   ```
   Prints the first 10 Gram points $g_n$ and confirms $(-1)^n Z(g_n) > 0$ at each.

## Arbitrary Precision

By default RiemannRho is fully dependency-free and uses 64-bit floating point. At large
$t$ this hits a hard ceiling: the main-sum argument $\theta(t) - t\ln k$ grows like
$t\ln t$, so once it exceeds ~$10^{15}$, `f64` can no longer pin the fractional part that
$\cos$ needs (catastrophic cancellation).

The optional **`bigfloat`** feature breaks this ceiling by evaluating the same
Riemann-Siegel formula in arbitrary precision (via the pure-Rust `dashu-float` crate). Build
with it and pass `--digits D`:

```
cargo build --release --features bigfloat
./target/release/riemannrho --nth 10000 --high-order --digits 25
```

For the 10000th zero this extends the result from the `f64` value `9881.1023608652` to
`9881.102360865374067`. (At small $t$ the dominant error is the asymptotic truncation of
the remainder series, which extra precision cannot fix — so `--digits` matters most at large
$t$.) The default build needs neither the feature nor the dependency.

## Gram Points & Turing's Method

A Gram point $g_n$ solves $\theta(g_n) = n\pi$. By *Gram's law*, $(-1)^n Z(g_n)$ is
usually positive and each interval $(g_{n-1}, g_n]$ usually contains exactly one zero — the
classical foundation for isolating zeros. RiemannRho computes $g_n$ by Newton's method on
$\theta$ and exposes them through `--gram`.

Gram's law is not exact: it first fails near $n = 126$ ($t \approx 282$) and periodically
thereafter. **Turing's method** copes by working with *Gram blocks* — a maximal run
$[g_a, g_b]$ bounded by good Gram points (where Gram's law holds) with only bad ones
inside. By *Rosser's rule*, a block spanning $b - a$ intervals contains exactly $b - a$
zeros. Since $S$ vanishes at good Gram points, $N(g_b) - N(g_a) = b - a$ exactly, so the
total over all blocks is the rigorous prediction.

`--turing T` runs this: it forms the Gram blocks up to $T$, counts the sign changes of
$Z$ in each, checks them against Rosser's rule, and confirms the grand total matches the
Gram-index prediction.

```
./target/release/riemannrho --turing 300
```
```
Gram-block zero count over (g_-1, g_136] = (9.6669, 298.8423]:
  zeros found:        137
  expected (Turing):  137
  Gram's-law failures: 2
  Gram blocks (len>1): 2
  Rosser violations:   0
  verified: yes - every Gram block resolved and the total matches N(g) exactly
```

This is the *method*, not a certified proof: rigorous bounds would need interval
arithmetic, and Rosser's rule itself has (much higher) exceptions this does not
special-case.

## Counting and Verification

`--count T` locates every zero of $Z(t)$ on the critical line with $0 < t \le T$ and
compares the tally with the smooth part of the Riemann-von Mangoldt formula,
$N(T) = \theta(T)/\pi + 1 + S(T)$.

Crucially, the count does **not** equal $\theta(T)/\pi + 1$ — it equals that smooth value
plus the oscillating term $S(T)$. For example at $T = 50$ there are 10 zeros while
$\theta(T)/\pi + 1 = 9.42$, i.e. $S(50) = 0.58$. The tool reports the *implied*
$S(T) = \text{found} - \text{estimate}$ and checks that it is a plausible small
fluctuation. A normal $|S(T)|$ means every zero up to $T$ was found on the critical line —
a practical (Turing-flavored) check, which is the numerical heart of testing the Riemann
hypothesis over a range.

If the scan appears to have missed an unusually close pair of zeros (a suspiciously
negative implied $S$), the resolution is automatically increased and the count retried, so
the check is self-healing.

This is a heuristic, not a rigorous proof: full rigor requires bounding $S(T)$ via Gram
points and Turing's method proper.

## Visualization Details

Generated plots feature:
- **Horizontal Axis**: Imaginary part t.
- **Vertical Axis**: Z(t) values.
- **Line Graph**: Interpolated curve from sampled data.
- **Zero Indicator**: Prominent red line if detected.
- **Zero Indicator**: A red vertical line at the located zero (if one was found).

Loads D3.js v7 from its public CDN, so plot viewing requires an internet connection.

## Accuracy

The Riemann-Siegel formula is asymptotic, so its error shrinks as `t` grows. The
correction terms reduce the error at moderate `t`. For the first zero (`t ≈ 14.13`):

| Mode                     | Result            | Absolute error |
|--------------------------|-------------------|----------------|
| Base (`C0`)              | 14.1371961027     | ~2.5e-3        |
| `--high-order` (`+C1+C2`)| 14.1348228651     | ~9.8e-5        |

For larger `t` the agreement improves substantially. Unit tests in `src/lib.rs`
check several known zeros (14.1347, 21.0220, 25.0109, 30.4249, 32.9351).

## Limitations

- The asymptotic approximation is least accurate for small `t` (roughly `t < 10`).
- Computations use 64-bit floating point, so both accuracy and runtime degrade for
  very large ordinals long before any theoretical limit. `--nth` for an extremely
  large `N` will be slow and increasingly imprecise.
- The correction terms compute derivatives of the Riemann-Siegel `Psi` function
  numerically; accuracy is reduced near its poles (fractional part `p ≈ 0.25` or `0.75`).
- Plot viewing requires a browser with internet access (D3.js CDN).

## Contributing

We encourage contributions! Fork, branch, and pull request with enhancements. Adhere to Rust standards and test additions. Discuss major features via issues.

## License

MIT License—refer to [LICENSE](LICENSE).

## Acknowledgments

Inspired by Riemann's seminal work and contemporary efforts (e.g., Odlyzko's computations). Gratitude to the Rust ecosystem for enabling precise numerical tools.

Contact: [iman.akbari@imnetsecurity.com] for inquiries.
