# flip-flops

Simple program for calculating the number of flip-flop events in Gromacs simulations.

Fast implementation of the lipid orientation method from [Li et al., 2024](https://doi.org/10.1073/pnas.2319476121).

## Installation
0) [Install Rust](https://rust-lang.org/tools/install/).

1) Clone this repository:
```
$ git clone https://github.com/Ladme/flip-flops.git
```

2) Navigate to the repository.

3) Build the program using `cargo`:
```
$ cargo build --release
```

A binary file named `flip-flops` will be built inside `target/release/`.


## Usage

```
Calculate number of flip-flop events.

Usage: flip-flops [OPTIONS] --structure <STRUCTURE> --trajectory <TRAJECTORY> --tails <TAILS>

Options:
  -s, --structure <STRUCTURE>
          Path to a gro, pdb, or tpr file containing the system structure.

  -f, --trajectory <TRAJECTORY>
          Path to the xtc file containing the trajectory to be analyzed.

  -n, --index <INDEX>
          Path to an ndx file containing groups associated with the system.

      --heads <HEADS>
          Groan selection language query selecting atoms representing lipid heads. There should only be one selected atom per lipid molecule.
          
          [default: "name PO4 P"]

      --tails <TAILS>
          Groan selection language query selecting atoms representing ends of lipid tails.

  -b, --begin <BEGIN>
          Time of the first frame to read from the trajectory (in ps). All previous frames will be skipped.
          
          [default: 0.0]

  -e, --end <END>
          Time of the last frame to read from the trajectory (in ps). All following frames will be skipped.
          
          [default: NaN]

      --transition <TRANSITION>
          Angle range that defines the transition zone between the upper and lower leaflets. Lipids whose orientation angle falls within this zone are treated as belonging to neither leaflet.
          
          [default: "55.0 125.0"]

      --window <WINDOW>
          Number of neighboring frames on each side of the current frame to include in the sliding window for running average of the lipid angle. The window spans 2N+1 frames in total. Set to 0 to disable smoothing.
          
          [default: 0]

  -h, --help
          Print help (see a summary with '-h')

  -V, --version
          Print version
```

## Limitations

This program has a lot of limitations. The most important ones include:
- the membrane is assumed to be planar and built in the `xy` plane,
- each lipid molecule is assumed to be composed of a single residue,
- simulation box is expected to be orthogonal.