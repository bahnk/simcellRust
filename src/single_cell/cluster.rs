use crate::single_cell::cell::Cell;

pub struct Cluster {
    name: String,
    cells: Vec<Cell>,
}

impl Cluster {
    pub fn new(name: String, cells: Vec<Cell>) -> Self {
        Cluster {
            name: name,
            cells: cells,
        }
    }
}
