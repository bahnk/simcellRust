use crate::single_cell::gene::Gene;

pub struct Cell {
    pub name: String,
    pub cluster: String,
    pub genes: Vec<Gene>,
}

impl Cell {
    pub fn new(name: String, cluster: String, genes: Vec<Gene>) -> Self {
        Cell {
            name: name,
            cluster: cluster,
            genes: genes,
        }
    }
}

impl std::fmt::Display for Cell {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "({}, {})\n", self.name, self.cluster);
        for gene in &self.genes {
            write!(f, "{}\n", gene);
        }
        Ok(())
    }
}
