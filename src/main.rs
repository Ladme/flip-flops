// Released under MIT License.
// Copyright (c) 2026 Ladislav Bartos

use std::{cmp::Ordering, process};

use clap::Parser;
use groan_rs::{
    prelude::{ProgressPrinter, TrajMasterRead, Vector3D},
    system::System,
};

use crate::{structures::Leaflet, structures::Lipid};

mod argparse;
mod structures;

fn construct_lipids(system: &System) -> Vec<Lipid> {
    let mut lipids = Vec::new();

    for head in system.group_iter("FlipFlopsReserved-Heads").unwrap() {
        let residue = head.get_residue_number();
        let mut matching_tails = Vec::new();

        for tail in system.group_iter("FlipFlopsReserved-Tails").unwrap() {
            let tail_res = tail.get_residue_number();
            match tail_res.cmp(&residue) {
                Ordering::Equal => matching_tails.push(tail.get_index()),
                Ordering::Greater => continue, // assuming the atoms are ordered
                Ordering::Less => (),
            }
        }

        lipids.push(Lipid {
            residue,
            head: head.get_index(),
            tails: matching_tails,
            angles: Vec::new(),
            leaflets: Vec::new(),
        })
    }

    lipids
}

fn analyze_frame(
    frame: &System,
    lipids: &mut [Lipid],
) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    for lipid in lipids.iter_mut() {
        let head_pos = frame.get_atom(lipid.head).unwrap().get_position().unwrap();
        let mut angles = Vec::new();

        for tail_index in lipid.tails.iter() {
            let tail_pos = frame.get_atom(*tail_index).unwrap().get_position().unwrap();
            let vector = head_pos.vector_to(tail_pos, frame.get_box().unwrap());
            angles.push(vector.angle(&Vector3D::new(0.0, 0.0, 1.0)).to_degrees());
        }

        let average = angles.iter().sum::<f32>() / angles.len() as f32;
        lipid.angles.push(average);
    }

    Ok(())
}

fn running_average(values: &[f32], window: usize) -> Vec<f32> {
    if window == 0 || values.len() <= 1 {
        return values.to_vec();
    }

    let n = values.len();
    let mut result = Vec::with_capacity(n);

    for i in 0..n {
        let start = i.saturating_sub(window);
        let end = (i + window + 1).min(n);
        let window = &values[start..end];
        let sum: f32 = window.iter().sum();
        result.push(sum / window.len() as f32);
    }

    result
}

fn assign_leaflets(lipids: &mut [Lipid], window: usize, transition: [f32; 2]) {
    let transition_min = transition[0];
    let transition_max = transition[1];

    for lipid in lipids.iter_mut() {
        let running_angles = running_average(&lipid.angles, window);

        for angle in running_angles {
            if angle < transition_min {
                lipid.leaflets.push(Leaflet::LOWER);
            } else if angle > transition_max {
                lipid.leaflets.push(Leaflet::UPPER);
            } else {
                lipid.leaflets.push(Leaflet::TRANSITIONING);
            }
        }
    }
}

fn detect_flipflops(lipids: &[Lipid]) -> (usize, usize) {
    let mut total_flipflops_down = 0;
    let mut total_flipflops_up = 0;

    for lipid in lipids.iter() {
        let mut lipid_flipflops_down = 0;
        let mut lipid_flipflops_up = 0;
        let mut last: Option<&Leaflet> = None;

        for (frame, leaflet) in lipid.leaflets.iter().enumerate() {
            match leaflet {
                Leaflet::TRANSITIONING => continue,
                current => {
                    if let Some(prev) = last {
                        match (prev, current) {
                            (Leaflet::UPPER, Leaflet::LOWER) => {
                                println!(
                                    "Lipid {} flipped from the UPPER to the LOWER leaflet in frame {}.",
                                    lipid.residue, frame
                                );
                                lipid_flipflops_down += 1;
                            }
                            (Leaflet::LOWER, Leaflet::UPPER) => {
                                println!(
                                    "Lipid {} flipped from the LOWER to the UPPER leaflet in frame {}.",
                                    lipid.residue, frame
                                );
                                lipid_flipflops_up += 1;
                            }
                            _ => (),
                        }
                    }

                    last = Some(current);
                }
            }
        }

        let lipid_flipflops = lipid_flipflops_down + lipid_flipflops_up;
        if lipid_flipflops > 0 {
            println!("Lipid {}: {} flip-flops", lipid.residue, lipid_flipflops);
            println!("  > UPPER->LOWER: {}", lipid_flipflops_down);
            println!("  > LOWER->UPPER: {}", lipid_flipflops_up);
            println!("\n");
        }

        total_flipflops_down += lipid_flipflops_down;
        total_flipflops_up += lipid_flipflops_up;
    }

    (total_flipflops_down, total_flipflops_up)
}

fn run() -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let args = argparse::Args::parse();

    let mut system = System::from_file(&args.structure)?;

    if let Some(ndx) = args.index {
        system.read_ndx(ndx)?;
    }

    system.group_create("FlipFlopsReserved-Heads", &args.heads)?;
    system.group_create("FlipFlopsReserved-Tails", &args.tails)?;
    system.group_create(
        "FlipFlopsReserved-Relevant",
        "FlipFlopsReserved-Heads FlipFlopsReserved-Tails",
    )?;

    println!("Analyzing system topology...");
    let mut lipids = construct_lipids(&system);

    println!("Reading trajectory...");
    let reader = system
        .group_xtc_iter(&args.trajectory, "FlipFlopsReserved-Relevant")?
        .print_progress(ProgressPrinter::default());

    let reader = match (args.begin, args.end) {
        (None, None) => reader.with_range(0.0, f32::MAX)?,
        (Some(start), None) => reader.with_range(start, f32::MAX)?,
        (None, Some(end)) => reader.with_range(0.0, end)?,
        (Some(start), Some(end)) => reader.with_range(start, end)?,
    };

    for frame in reader {
        let frame = frame?;
        analyze_frame(&frame, &mut lipids)?;
    }

    println!("Assigning lipids to leaflets...");
    assign_leaflets(&mut lipids, args.window, args.transition);

    println!("Detecting flip-flops...");
    let (flipflops_down, flipflops_up) = detect_flipflops(&lipids);

    println!(
        "\nTOTAL NUMBER OF FLIP-FLOPS: {}",
        flipflops_down + flipflops_up
    );
    println!("  > UPPER->LOWER: {}", flipflops_down);
    println!("  > LOWER->UPPER: {}", flipflops_up);

    Ok(())
}

fn main() {
    if let Err(e) = run() {
        eprintln!("{}", e);
        process::exit(1);
    }

    process::exit(0);
}
