// Released under MIT License.
// Copyright (c) 2026 Ladislav Bartos

use std::f32;

use clap::Parser;

// Calculate number of lipid flip-flops in a simulation trajectory.
#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about,
    long_about = "Calculate number of flip-flop events."
)]
pub struct Args {
    #[arg(
        short = 's',
        long = "structure",
        help = "Input structure file",
        long_help = "Path to a gro, pdb, or tpr file containing the system structure."
    )]
    pub structure: String,

    #[arg(
        short = 'f',
        long = "trajectory",
        help = "Input trajectory file",
        long_help = "Path to the xtc file(s) containing the trajectory to be analyzed.",
        num_args = 1..,
    )]
    pub trajectories: Vec<String>,

    #[arg(
        short = 'n',
        long = "index",
        help = "Input index file",
        long_help = "Path to an ndx file containing groups associated with the system."
    )]
    pub index: Option<String>,

    #[arg(
        long = "heads",
        default_value_t = String::from("name PO4 P"),
        help = "Selection of atoms representing lipid heads. Only one atom per lipid molecule!",
        long_help = "Groan selection language query selecting atoms representing lipid heads. There should only be one selected atom per lipid molecule."
    )]
    pub heads: String,

    #[arg(
        long = "tails",
        help = "Selection of atoms representing the ends of lipid tails",
        long_help = "Groan selection language query selecting atoms representing ends of lipid tails."
    )]
    pub tails: String,

    #[arg(
        short = 'b',
        long = "begin",
        help = "Time of the first frame to read (in ps) [default: 0.0]",
        long_help = "Time of the first frame to read from the trajectory (in ps). All previous frames will be skipped.\n\n[default: 0.0]"
    )]
    pub begin: Option<f32>,

    #[arg(
        short = 'e',
        long = "end",
        help = "Time of the last frame to read (in ps) [default: NaN]",
        long_help = "Time of the last frame to read from the trajectory (in ps). All following frames will be skipped.\n\n[default: NaN]"
    )]
    pub end: Option<f32>,

    #[arg(
        long = "transition",
        value_parser = parse_transition,
        default_value = "55.0 125.0",
        help = "Angle range (in degrees) defining the transition zone between leaflets",
        long_help = "Angle range that defines the transition zone between the upper and lower leaflets. \
                    Lipids whose orientation angle falls within this zone are treated as belonging to neither leaflet.",
    )]
    pub transition: [f32; 2],

    #[arg(
        long = "window",
        default_value_t = 0,
        help = "Number of neighboring frames on each side of the sliding window",
        long_help = "Number of neighboring frames on each side of the current frame to include \
                 in the sliding window for running average of the lipid angle. The window \
                 spans 2N+1 frames in total. Set to 0 to disable smoothing."
    )]
    pub window: usize,
}

fn parse_transition(s: &str) -> Result<[f32; 2], String> {
    let parts: Vec<&str> = s.split_whitespace().collect();

    if parts.len() != 2 {
        return Err(format!(
            "Expected two space-separated values in the form \"MIN MAX\", got {} value(s)",
            parts.len()
        ));
    }

    let min: f32 = parts[0]
        .parse()
        .map_err(|e| format!("Invalid MIN value '{}': {}", parts[0], e))?;
    let max: f32 = parts[1]
        .parse()
        .map_err(|e| format!("Invalid MAX value '{}': {}", parts[1], e))?;

    if !min.is_finite() || !max.is_finite() {
        return Err("MIN and MAX must be finite numbers".to_string());
    }

    if min >= max {
        return Err(format!(
            "MIN ({}) must be strictly less than MAX ({})",
            min, max
        ));
    }

    Ok([min, max])
}
