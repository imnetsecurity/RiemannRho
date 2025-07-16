# RiemannRho: Rust Library for Riemann Zeta Nontrivial Zeros Approximation

## About

RiemannRho is a high-performance Rust library and command-line tool dedicated to the numerical approximation of nontrivial zeros of the Riemann zeta function \(\zeta(s)\) on the critical line \(\operatorname{Re}(s) = 1/2\). Utilizing the Riemann-Siegel asymptotic formula, RiemannRho computes Hardy's Z-function \(Z(t)\), facilitating precise location of zeros corresponding to \(\zeta(1/2 + it) = 0\). Engineered for accuracy, efficiency, and extensibility, this tool serves researchers in analytic number theory, supporting investigations into the Riemann hypothesis through scalable computations of high-order zeros.

Core strengths include:
- **Precision-Enhanced Algorithms**: Employs base and higher-order remainder terms (up to C2) for reduced error in zero approximations, ideal for both low and high \(t\).
- **High-Order Correction Mode**: Boosts computational fidelity with advanced corrections, enabling near-exact results for foundational zeros.
- **Dynamic Visualization**: Optional export of interactive D3.js plots to HTML, rendering \(Z(t)\) curves with zero annotations for visual analysis.
- **Large-Scale Capability**: Handles zeros up to the \(10^{12}\)th order or beyond, with \(O(\sqrt{t})\) efficiency, pushing boundaries of desktop computing.

Licensed under MIT, RiemannRho promotes open collaboration in mathematical software, ensuring robust, reproducible results without external dependencies.

## Features

- **Range-Based Zero Detection**: Define intervals \([low, high]\) for bisection-based zero isolation via \(Z(t)\).
- **Ordinal Zero Approximation**: The `--nth` option estimates the \(n\)th zero's imaginary part using refined asymptotic expansions and Newton's iteration.
- **High-Order Correction Mode**: Invoke `--high-order` to incorporate additional terms, minimizing discrepancies (e.g., first zero refined to ~14.1347251417).
- **Browser-Based Plots**: Post-calculation prompt generates `zeta_plot.html` with D3.js visualizations: smooth lines, axes, and red zero markers.
- **Flexible Interfaces**: Command-line parameters or interactive prompts, with customizable tolerance for convergence control.
- **Dependency-Free**: Built exclusively on Rust's standard library, guaranteeing cross-platform reliability.

## Installation

RiemannRho requires Rust (1.70+). Installation steps:

1. Install Rust via [rustup](https://rustup.rs/).
2. Clone the source:
git clone https://github.com/yourusername/riemannrho.git
cd riemannrho
3. Compile:
cargo build --release
Produces an optimized executable at `target/release/riemannrho`.

The tool operates independently, with no runtime dependencies.

## Usage

Execute with arguments or via interactive mode.

### Command Syntax
./target/release/riemannrho [low] [high] [tol] [--high-order] [--nth n]

- `low`, `high`: Search interval bounds (omit for `--nth`).
- `tol`: Precision threshold (default: 1e-10).
- `--high-order`: Activate enhanced correction terms.
- `--nth n`: Target the nth zero (supports large n, e.g., 1e12).

Upon completion, view the zero approximation and respond to the visualization prompt: "Do you want a D3.js visualization? (yes/no)". Affirmative generates `zeta_plot.html` for browser viewing.

### Usage Examples

1. **First Zero with High-Order Correction Mode**:
./target/release/riemannrho 14 15 1e-10 --high-order
Sample Output:
Approximate imaginary part of the nontrivial zero: 14.1347251417
Do you want a D3.js visualization? (yes/no)
npm run dev
Visualization file depicts \(Z(t)\) crossing at ~14.1347.

2. **Millionth Zero**:
./target/release/riemannrho 1e-10 --nth 1000000 --high-order
Yields t ≈ 1.84 × 10^6.

3. **Prompted Mode** (argument-free):
./target/release/riemannrho
Guides through input for low, high, tol.

## Visualization Details

Generated plots feature:
- **Horizontal Axis**: Imaginary part t.
- **Vertical Axis**: Z(t) values.
- **Line Graph**: Interpolated curve from sampled data.
- **Zero Indicator**: Prominent red line if detected.
- **Engagement**: Browser supports zoom and pan via D3.js.

Utilizes the public D3.js CDN, creating plots from computed values without extraneous text.

## Limitations

- Asymptotic approximations may yield minor errors for t < 10; further terms can be added for ultra-precision.
- Ultra-large n (>10^15) risks overflow or extended runtime; includes advisory warnings.
- Plots require a contemporary browser for full interactivity.

## Contributing

We encourage contributions! Fork, branch, and pull request with enhancements. Adhere to Rust standards and test additions. Discuss major features via issues.

## License

MIT License—refer to [LICENSE](LICENSE).

## Acknowledgments

Inspired by Riemann's seminal work and contemporary efforts (e.g., Odlyzko's computations). Gratitude to the Rust ecosystem for enabling precise numerical tools.

Contact: [iman.akbari@imnetsecurity.com] for inquiries.
