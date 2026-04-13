// Released under MIT License.
// Copyright (c) 2026 Ladislav Bartos

#[derive(Debug, Clone)]
pub(crate) enum Leaflet {
    UPPER,
    LOWER,
    TRANSITIONING,
}

#[derive(Debug, Clone)]
pub(crate) struct Lipid {
    pub(crate) residue: usize,
    pub(crate) head: usize,
    pub(crate) tails: Vec<usize>,
    pub(crate) angles: Vec<f32>,
    pub(crate) leaflets: Vec<Leaflet>,
}
