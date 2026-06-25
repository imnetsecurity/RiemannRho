//! Binary executable for RiemannRho: command-line interface, computation, and visualization.

use riemannrho::{
    estimate_t, find_zero, gram_point, nth_zero, verify_zero_count, z_func, zeros_below, Precision,
};
use std::env;
use std::fs::File;
use std::io::{self, Write};
use std::process::ExitCode;

const DEFAULT_TOL: f64 = 1e-10;
const DEFAULT_PLOT_PATH: &str = "zeta_plot.html";
/// Largest ordinal for which `--nth` scans sequentially for the *exact* nth zero;
/// beyond this it falls back to the asymptotic estimate of a single nearby zero.
const NTH_EXACT_MAX: f64 = 100_000.0;

fn print_usage(program: &str) {
    eprintln!(
        "RiemannRho - approximate nontrivial zeros of the Riemann zeta function.

USAGE:
    {program} [low] [high] [tol] [--high-order] [--nth N] [--out FILE]
    {program} --count T [--high-order]
    {program} --list T [--high-order]
    {program} --gram N [--high-order]
    {program}                      (interactive mode)

ARGUMENTS:
    low, high    Search interval bounds for the root finder (omit when using --nth).
    tol          Bracket-width tolerance for the root finder (default: {DEFAULT_TOL}).

OPTIONS:
    --high-order     Include higher-order Riemann-Siegel correction terms.
    --nth N          Target the Nth zero (N >= 1). Exact for N <= {NTH_EXACT_MAX:.0},
                     asymptotic estimate beyond that.
    --count T        Count zeros with 0 < t <= T and compare with the theoretical
                     count theta(T)/pi + 1 (a Turing-flavored consistency check).
    --list T         Print every zero with 0 < t <= T, one per line.
    --gram N         Print the first N Gram points and check Gram's law at each.
    --out FILE       Path for the generated plot (default: {DEFAULT_PLOT_PATH}).
    -h, --help       Print this help.

EXAMPLES:
    {program} 14 15 1e-10 --high-order
    {program} --nth 1 --high-order
    {program} --count 100 --high-order
    {program} --list 50
    {program} --gram 10 --high-order
    {program}"
    );
}

/// Parses `f64` arguments with a descriptive error instead of panicking.
fn parse_f64(value: &str, what: &str) -> Result<f64, String> {
    value
        .trim()
        .parse::<f64>()
        .map_err(|_| format!("invalid value for {what}: {value:?}"))
}

struct CliArgs {
    low: Option<f64>,
    high: Option<f64>,
    tol: f64,
    precision: Precision,
    nth: Option<f64>,
    count: Option<f64>,
    list: Option<f64>,
    gram: Option<u64>,
    out: String,
}

fn parse_args(args: &[String]) -> Result<CliArgs, String> {
    let mut cli = CliArgs {
        low: None,
        high: None,
        tol: DEFAULT_TOL,
        precision: Precision::Base,
        nth: None,
        count: None,
        list: None,
        gram: None,
        out: DEFAULT_PLOT_PATH.to_string(),
    };

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--high-order" => cli.precision = Precision::Order2,
            "--nth" => {
                i += 1;
                let v = args
                    .get(i)
                    .ok_or_else(|| "--nth requires a value".to_string())?;
                let n = parse_f64(v, "--nth")?;
                if n < 1.0 {
                    return Err(format!("--nth must be >= 1 (got {n})"));
                }
                cli.nth = Some(n);
            }
            "--count" => {
                i += 1;
                let v = args
                    .get(i)
                    .ok_or_else(|| "--count requires a value".to_string())?;
                let t = parse_f64(v, "--count")?;
                if t <= 0.0 {
                    return Err(format!("--count must be > 0 (got {t})"));
                }
                cli.count = Some(t);
            }
            "--list" => {
                i += 1;
                let v = args
                    .get(i)
                    .ok_or_else(|| "--list requires a value".to_string())?;
                let t = parse_f64(v, "--list")?;
                if t <= 0.0 {
                    return Err(format!("--list must be > 0 (got {t})"));
                }
                cli.list = Some(t);
            }
            "--gram" => {
                i += 1;
                let v = args
                    .get(i)
                    .ok_or_else(|| "--gram requires a value".to_string())?;
                let n: u64 = v
                    .trim()
                    .parse()
                    .map_err(|_| format!("invalid value for --gram: {v:?}"))?;
                cli.gram = Some(n);
            }
            "--out" => {
                i += 1;
                cli.out = args
                    .get(i)
                    .ok_or_else(|| "--out requires a value".to_string())?
                    .clone();
            }
            other => {
                if cli.low.is_none() {
                    cli.low = Some(parse_f64(other, "low")?);
                } else if cli.high.is_none() {
                    cli.high = Some(parse_f64(other, "high")?);
                } else {
                    cli.tol = parse_f64(other, "tol")?;
                }
            }
        }
        i += 1;
    }

    if cli.tol <= 0.0 {
        return Err(format!("tol must be positive (got {})", cli.tol));
    }
    Ok(cli)
}

/// Reads and parses one line from stdin in interactive mode.
fn read_f64(
    prompt: &str,
    lines: &mut impl Iterator<Item = io::Result<String>>,
) -> Result<f64, String> {
    println!("{prompt}");
    let line = lines
        .next()
        .ok_or_else(|| "unexpected end of input".to_string())?
        .map_err(|e| format!("failed to read input: {e}"))?;
    parse_f64(&line, "input")
}

/// Scans `[center - radius, center + radius]` for the sign-change bracket of `Z(t)`
/// whose midpoint is closest to `center`, returning it as `(lo, hi)`.
fn bracket_near(center: f64, radius: f64, precision: Precision) -> Option<(f64, f64)> {
    const SAMPLES: usize = 400;
    let lo = center - radius;
    let step = 2.0 * radius / SAMPLES as f64;

    let mut best: Option<(f64, f64)> = None;
    let mut best_dist = f64::INFINITY;
    let mut t_prev = lo;
    let mut z_prev = z_func(lo, precision);
    for i in 1..=SAMPLES {
        let t = lo + i as f64 * step;
        let z = z_func(t, precision);
        if z_prev.is_finite() && z.is_finite() && z_prev * z <= 0.0 {
            let dist = ((t_prev + t) / 2.0 - center).abs();
            if dist < best_dist {
                best_dist = dist;
                best = Some((t_prev, t));
            }
        }
        t_prev = t;
        z_prev = z;
    }
    best
}

/// Resolves the search interval and tolerance, prompting interactively if needed.
fn resolve_interval(cli: &CliArgs) -> Result<(f64, f64, f64), String> {
    if let Some(n) = cli.nth {
        // For moderate, integer n, scan sequentially for the *exact* nth zero; otherwise
        // fall back to the asymptotic estimate of a single nearby zero.
        let center = if n.fract() == 0.0 && (1.0..=NTH_EXACT_MAX).contains(&n) {
            nth_zero(n as u64, cli.precision)
                .ok_or_else(|| format!("could not locate zero #{n}"))?
        } else {
            let est = estimate_t(n);
            if !est.is_finite() {
                return Err(format!("could not estimate the {n}th zero"));
            }
            if est > 1e30 {
                eprintln!("Warning: n={n} is extremely large; computation may be very slow.");
            }
            est
        };
        // Search a window a few average spacings wide and bracket the sign change
        // nearest the center. A naive [c - spacing, c + spacing] window can straddle
        // two zeros, in which case the endpoints share a sign and none is detected.
        let spacing = (2.0 * std::f64::consts::PI / (center / (2.0 * std::f64::consts::PI)).ln())
            .abs()
            .max(1.0);
        let (low, high) = bracket_near(center, 3.0 * spacing, cli.precision).ok_or_else(|| {
            format!("no zero found near the located position t ~= {center:.4} for n = {n}")
        })?;
        Ok((low, high, cli.tol))
    } else if let (Some(low), Some(high)) = (cli.low, cli.high) {
        if low >= high {
            return Err(format!("low ({low}) must be less than high ({high})"));
        }
        Ok((low, high, cli.tol))
    } else {
        let stdin = io::stdin();
        let mut lines = stdin.lines();
        let low = read_f64(
            "Enter low bound (or use --nth for ordinal-zero mode):",
            &mut lines,
        )?;
        let high = read_f64("Enter high bound:", &mut lines)?;
        let tol = read_f64("Enter tolerance (e.g., 1e-10):", &mut lines)?;
        if low >= high {
            return Err(format!("low ({low}) must be less than high ({high})"));
        }
        if tol <= 0.0 {
            return Err(format!("tol must be positive (got {tol})"));
        }
        Ok((low, high, tol))
    }
}

fn run() -> Result<(), String> {
    let args: Vec<String> = env::args().collect();
    let program = args
        .first()
        .map(String::as_str)
        .unwrap_or("riemannrho")
        .to_string();

    if args.iter().any(|a| a == "-h" || a == "--help") {
        print_usage(&program);
        return Ok(());
    }

    let cli = parse_args(&args)?;

    // Gram-point mode: list the first N Gram points and illustrate Gram's law.
    if let Some(n) = cli.gram {
        println!(
            "{:>6}  {:>14}  {:>14}  Gram's law",
            "n", "g_n", "(-1)^n Z(g_n)"
        );
        for k in 0..n as i64 {
            let g = gram_point(k);
            let sign = if k % 2 == 0 { 1.0 } else { -1.0 };
            let val = sign * z_func(g, cli.precision);
            let ok = if val > 0.0 { "holds" } else { "VIOLATED" };
            println!("{k:>6}  {g:>14.6}  {val:>14.6}  {ok}");
        }
        return Ok(());
    }

    // Listing mode: print every zero up to T, one per line.
    if let Some(t_max) = cli.list {
        let zeros = zeros_below(t_max, cli.precision);
        for (i, z) in zeros.iter().enumerate() {
            println!("{:>8}  {:.10}", i + 1, z);
        }
        println!("{} zeros with 0 < t <= {t_max}", zeros.len());
        return Ok(());
    }

    // Counting mode is a standalone report; it does not produce a single zero or plot.
    if let Some(t_max) = cli.count {
        let report = verify_zero_count(t_max, cli.precision);
        println!("Zeros found with 0 < t <= {t_max}: {}", report.found);
        println!("Smooth estimate theta(T)/pi + 1: {:.4}", report.expected);
        println!(
            "Implied S(T) = found - estimate: {:+.4}  (true count = estimate + S)",
            report.s
        );
        if report.resolution > 10.0 {
            println!(
                "(scan refined to {} samples per spacing to resolve close zeros)",
                report.resolution
            );
        }
        if report.consistent {
            println!(
                "Consistent: S(T) is a normal small fluctuation, so every zero up to T was \
                 found on the critical line (Turing-flavored check passes)."
            );
        } else {
            println!(
                "Suspicious: |S(T)| is implausibly large, suggesting a missed or spurious \
                 zero. Try --high-order, or a slightly different T (T may be near a zero)."
            );
        }
        return Ok(());
    }

    let (low, high, tol) = resolve_interval(&cli)?;

    let zero = find_zero(low, high, tol, cli.precision);
    match zero {
        Some(zero_t) => {
            println!("Approximate imaginary part of the nontrivial zero: {zero_t:.10}");
        }
        None => {
            println!("No sign change detected in [{low}, {high}]. Adjust the interval or try a different n.");
        }
    }

    // Visualization prompt.
    println!("Do you want a D3.js visualization? (yes/no)");
    let mut input = String::new();
    if io::stdin()
        .read_line(&mut input)
        .map_err(|e| e.to_string())?
        == 0
    {
        return Ok(()); // EOF: treat as "no".
    }
    if input.trim().eq_ignore_ascii_case("yes") {
        generate_d3_plot(&cli.out, low, high, zero, cli.precision)
            .map_err(|e| format!("error generating plot: {e}"))?;
    }
    Ok(())
}

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

/// Generates a D3.js HTML plot of Z(t) over the search interval (binary-only).
fn generate_d3_plot(
    path: &str,
    low: f64,
    high: f64,
    zero: Option<f64>,
    precision: Precision,
) -> std::io::Result<()> {
    let mut file = File::create(path)?;
    let num_points = 200;
    let step = (high - low) / (num_points as f64 - 1.0);

    let mut data_str = String::from("const data = [\n");
    let mut min_z = f64::INFINITY;
    let mut max_z = f64::NEG_INFINITY;
    for i in 0..num_points {
        let t = low + (i as f64) * step;
        let z = z_func(t, precision);
        min_z = min_z.min(z);
        max_z = max_z.max(z);
        data_str.push_str(&format!("  {{ t: {t}, z: {z} }},\n"));
    }
    data_str.push_str("];\n");

    let zero_line = if let Some(zero_t) = zero {
        format!(
            r#"g.append("line")
            .attr("x1", x({zero_t}))
            .attr("y1", 0)
            .attr("x2", x({zero_t}))
            .attr("y2", height)
            .attr("stroke", "red")
            .attr("stroke-width", 2);"#
        )
    } else {
        String::new()
    };

    let y_min = min_z - 0.1 * (max_z - min_z);
    let y_max = max_z + 0.1 * (max_z - min_z);

    let html = format!(
        r#"<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Riemann Zeta Z(t) Plot</title>
    <script src="https://d3js.org/d3.v7.min.js"></script>
    <style>
        .chart {{
            margin: 20px;
        }}
        .axis path,
        .axis line {{
            stroke: #000;
            shape-rendering: crispEdges;
        }}
        .line {{
            fill: none;
            stroke: steelblue;
            stroke-width: 1.5px;
        }}
    </style>
</head>
<body>
    <div class="chart">
        <svg width="800" height="500"></svg>
    </div>

    <script>
        {data_str}

        const svg = d3.select("svg");
        const margin = {{ top: 20, right: 20, bottom: 30, left: 50 }};
        const width = +svg.attr("width") - margin.left - margin.right;
        const height = +svg.attr("height") - margin.top - margin.bottom;

        const g = svg.append("g")
            .attr("transform", `translate(${{margin.left}},${{margin.top}})`);

        const x = d3.scaleLinear()
            .domain([{low}, {high}])
            .range([0, width]);

        const y = d3.scaleLinear()
            .domain([{y_min}, {y_max}])
            .range([height, 0]);

        g.append("g")
            .attr("class", "axis axis--x")
            .attr("transform", `translate(0,${{height}})`)
            .call(d3.axisBottom(x));

        g.append("g")
            .attr("class", "axis axis--y")
            .call(d3.axisLeft(y));

        const line = d3.line()
            .x(d => x(d.t))
            .y(d => y(d.z));

        g.append("path")
            .datum(data)
            .attr("class", "line")
            .attr("d", line);

        {zero_line}
    </script>
</body>
</html>"#
    );

    file.write_all(html.as_bytes())?;
    println!("Visualization generated in {path}. Open it in a web browser.");
    Ok(())
}
