//! Binary executable for RiemannRho: Handles CLI, computation, and visualization.

use riemannrho::{estimate_t, find_zero, z_func}; // Import from lib
use std::env;
use std::fs::File;
use std::io::{self, BufRead, Write};

fn main() {
    // Parse command-line arguments
    let args: Vec<String> = env::args().collect();
    let mut i = 1;
    let mut low: Option<f64> = None;
    let mut high: Option<f64> = None;
    let mut tol: f64 = 1e-10;
    let mut terms: u32 = 0;
    let mut nth: Option<f64> = None;

    while i < args.len() {
        match args[i].as_str() {
            "--high-order" => terms = 2,
            "--nth" => {
                i += 1;
                nth = Some(args[i].parse().expect("Invalid nth value"));
            }
            _ => {
                if low.is_none() {
                    low = Some(args[i].parse().expect("Invalid low"));
                } else if high.is_none() {
                    high = Some(args[i].parse().expect("Invalid high"));
                } else {
                    tol = args[i].parse().expect("Invalid tol");
                }
            }
        }
        i += 1;
    }

    // Determine interval
    let (l, h) = if let Some(n) = nth {
        let est = estimate_t(n);
        if est > 1e30 {
            println!("Warning: n={} is extremely large; computation may take forever or overflow.", n);
        }
        (est - 1.0, est + 1.0)
    } else if let (Some(ll), Some(hh)) = (low, high) {
        (ll, hh)
    } else {
        // Interactive mode
        println!("Enter low bound (or use --nth for high-order correction mode):");
        let stdin = io::stdin();
        let mut lines = stdin.lines();

        let low_str = lines.next().unwrap().unwrap();
        let ll: f64 = low_str.trim().parse().expect("Invalid low");

        println!("Enter high bound:");
        let high_str = lines.next().unwrap().unwrap();
        let hh: f64 = high_str.trim().parse().expect("Invalid high");

        println!("Enter tolerance (e.g., 1e-10):");
        let tol_str = lines.next().unwrap().unwrap();
        tol = tol_str.trim().parse().expect("Invalid tol");

        (ll, hh)
    };

    // Compute zero
    let zero = find_zero(l, h, tol, terms);
    match zero {
        Some(zero_t) => {
            println!("Approximate imaginary part of the nontrivial zero: {:.10}", zero_t);
        }
        None => {
            println!("No sign change detected in [{}, {}]. Adjust interval or try smaller n.", l, h);
        }
    }

    // Prompt for visualization
    println!("Do you want a D3.js visualization? (yes/no)");
    let mut input = String::new();
    io::stdin().read_line(&mut input).unwrap();
    if input.trim().to_lowercase() == "yes" {
        if let Err(e) = generate_d3_plot(l, h, zero, terms) {
            println!("Error generating plot: {}", e);
        }
    }
}

/// Generates D3.js HTML plot (binary-only; not in lib).
fn generate_d3_plot(low: f64, high: f64, zero: Option<f64>, terms: u32) -> std::io::Result<()> {
    let mut file = File::create("zeta_plot.html")?;
    let num_points = 200;
    let step = (high - low) / (num_points as f64 - 1.0);

    let mut data_str = String::new();
    data_str.push_str("const data = [\n");
    let mut min_z = f64::INFINITY;
    let mut max_z = f64::NEG_INFINITY;
    for i in 0..num_points {
        let t = low + (i as f64) * step;
        let z = z_func(t, terms);
        min_z = min_z.min(z);
        max_z = max_z.max(z);
        data_str.push_str(&format!("  {{ t: {}, z: {} }},\n", t, z));
    }
    data_str.push_str("];\n");

    let zero_line = if let Some(zero_t) = zero {
        format!(
            r#"g.append("line")
            .attr("x1", x({}))
            .attr("y1", 0)
            .attr("x2", x({}))
            .attr("y2", height)
            .attr("stroke", "red")
            .attr("stroke-width", 2);"#,
            zero_t, zero_t
        )
    } else {
        String::new()
    };

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
            .domain([{}, {}])
            .range([0, width]);

        const y = d3.scaleLinear()
            .domain([{}, {}])
            .range([height, 0]);

        g.append("g")
            .attr("class", "axis axis--x")
            .attr("transform", `translate(0,${{height}}`)
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
</html>"#,
        low, high, min_z - 0.1 * (max_z - min_z), max_z + 0.1 * (max_z - min_z)
    );

    file.write_all(html.as_bytes())?;
    println!("Visualization generated in zeta_plot.html. Open it in a web browser.");
    Ok(())
}