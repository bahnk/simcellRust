use crate::distributions::negative_binomial::ZINB;
/// There are variable genes and non variable genes.
/// Non variable genes have the same parameters for every cells.
/// Variable genes have different parameters between clusters.
/// Cells from the same cluster have the same parameters.

/// G: number of genes
/// N: number of cells
/// V: number of variable genes
/// D: distance matrix
///
/// d = 1 means two clusters have the same parameters for all the variable genes
/// d = 0 means two clusters have different parameters for all the variable genes
use crate::single_cell::cell::Cell;
use crate::single_cell::gene::Gene;

use rand::distributions::Distribution; // for sample
use rand::rngs::StdRng;
use statrs::distribution::Normal;
use statrs::distribution::Uniform;

use std::convert::TryFrom;

use std::io::Write;
use std::path::Path;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::fs::File;
use std::fs::create_dir;

pub struct Dataset {
    n_genes: u32,
    cells: Vec<Cell>,
}

impl Dataset {
    pub fn new(n_genes: u32, n_cells: u32, n_clusters: u32, rng: &mut StdRng) -> Self {
        /* The constant genes have the same parameters for every cells, while
         * variable genes have the same parameters for the cells in the same
         * cluster. */
        let pct_var_genes: f64 = 0.3;
        let n_var_genes: u32 = (n_genes as f64 * pct_var_genes).ceil() as u32;

        // The constant genes have mean and variance following a normal
        // distribution and dropout rate is uniform.
        let mean_dist = Normal::new(100.0, 20.0).unwrap();
        let var_dist = Normal::new(200.0, 50.0).unwrap();
        let dropout_dist = Uniform::new(0.0, 1.0).unwrap();

        // We set all the genes as constant and will replace with the variable
        // genes later.
        let mut genes: Vec<Gene> = Vec::<Gene>::new();
        for g in 1..n_genes + 1 {
            // Gene ID and symbol.
            let fill: usize = n_genes.to_string().len();
            let symbol: String = format!("Gene{:0fill$}", g, fill = fill);
            let id: String = format!("ENS{:0fill$}", g, fill = fill);

            // ZINB distribution parameters.
            let mean: f64 = mean_dist.sample(rng);
            let variance: f64 = var_dist.sample(rng);
            let dropout: f64 = dropout_dist.sample(rng);
            let zinb: ZINB = ZINB::new(mean, variance, dropout);

            // We create the gene.
            let gene: Gene = Gene::new(id, symbol, "Constant".to_string(), zinb);
            genes.push(gene);
        }

        // We create all the cells with the same parameters.
        let mut cells: Vec<Cell> = Vec::<Cell>::new();
        for n in 1..n_cells + 1 {
            let fill: usize = n_cells.to_string().len();
            let cell_name: String = format!("Cell{:0fill$}", n, fill = fill);
            let cell: Cell = Cell::new(cell_name, "NA".to_string(), genes.clone());
            cells.push(cell);
        }

        // Keeps tracks for the cells and genes that were created.
        let mut accum_genes: u32 = 0;
        let mut accum_cells: u32 = 0;

        // We create the cells for each cluster.
        for cluster in 1..n_clusters + 1 {
            // Cluster name.
            let fill: usize = n_clusters.to_string().len();
            let cluster_name: String = format!("Cluster{:0fill$}", cluster, fill = fill);

            // Number of cells
            let n: u32;
            if cluster <= n_cells % n_clusters {
                n = n_cells / n_clusters + 1;
            } else {
                n = n_cells / n_clusters;
            }

            // Number of genes.
            let g: u32;
            if cluster <= n_var_genes % n_clusters {
                g = n_var_genes / n_clusters + 1;
            } else {
                g = n_var_genes / n_clusters;
            }

            // Cells and genes for the current cluster.
            let lower_cell: usize = usize::try_from(accum_cells).unwrap();
            let upper_cell: usize = usize::try_from(accum_cells + n).unwrap();
            let lower_gene: usize = usize::try_from(accum_genes).unwrap();
            let upper_gene: usize = usize::try_from(accum_genes + g).unwrap();

            for cell in cells[lower_cell..upper_cell].iter_mut() {
                cell.cluster = cluster_name.clone();

                for gene in cell.genes[lower_gene..upper_gene].iter_mut() {
                    // ZINB distribution parameters.
                    let mean: f64 = mean_dist.sample(rng);
                    let variance: f64 = var_dist.sample(rng);
                    let dropout: f64 = dropout_dist.sample(rng);
                    let zinb: ZINB = ZINB::new(mean, variance, dropout);

                    *gene = Gene::new(
                        gene.id.clone(),
                        gene.symbol.clone(),
                        cluster_name.clone(),
                        zinb,
                    );
                }
            }

            accum_genes += g;
            accum_cells += n;
        }

        Dataset {
            n_genes: n_genes,
            cells: cells,
        }
    }

    pub fn write(&self, path: &Path, rng: &mut StdRng) -> () {
        if ! path.exists() {
            create_dir(path);
        }

        
        ///////////////////////////////////////////////////////////////////////
        // MATRIX

        // entries
        let n_genes: usize = self.n_genes as usize;
        let n_cells: usize = self.cells.len();
        let mut entries: Vec<(usize, usize, u64)> = Vec::<(usize, usize, u64)>::new();
        for i in 1..n_cells+1 {
            for j in 1..n_genes+1 {
                let count: u64 = self.cells[i-1].genes[j-1].distribution.sample(rng);
                if count > 0 {
                    entries.push((j, i, count));
                }

            }
        }

        entries.sort();

        let n_entries: usize = entries.len();
 
        let mut matrix_path = path.join("matrix.mtx.gz");
        let matrix = File::create(matrix_path).unwrap();
        let mut encoder = GzEncoder::new(matrix, Compression::default());

        // header
        encoder.write(b"%%MatrixMarket matrix coordinate integer general\n").unwrap();
        encoder.write(b"%metadata_json: {\"program\": \"simcell\"}\n").unwrap();
        encoder.write_fmt(format_args!("{} {} {}\n", n_genes, n_cells, n_entries)).unwrap();

        // entries
        for entry in entries {
            encoder.write((format!("{} {} {}\n", entry.0, entry.1, entry.2)).as_bytes());
        }

        encoder.finish().unwrap();

        ///////////////////////////////////////////////////////////////////////
        // CELLS

        let mut barcodes_path = path.join("barcodes.tsv.gz");
        let barcodes = File::create(barcodes_path).unwrap();
        let mut encoder = GzEncoder::new(barcodes, Compression::default());
        for cell in &self.cells {
            encoder.write((cell.name.clone()+"\n").as_bytes());
        }
        encoder.finish().unwrap();
 
        let mut obs_path = path.join("obs.tsv.gz");
        let obs = File::create(obs_path).unwrap();
        let mut encoder = GzEncoder::new(obs, Compression::default());
        for cell in &self.cells {
            encoder.write(
                format!("{}\t{}\n", cell.name.clone(), cell.cluster.clone()).as_bytes()
            );
        }
        encoder.finish().unwrap();
 
        ///////////////////////////////////////////////////////////////////////
        // GENES
        let mut features_path = path.join("features.tsv.gz");
        let features = File::create(features_path).unwrap();
        let mut encoder = GzEncoder::new(features, Compression::default());
        let genes = self.cells[0].genes.clone();
        for gene in genes {
            encoder.write(
                format!("{}\t{}\tGene Expression\n", gene.id, gene.symbol).as_bytes()
            );
        }
        encoder.finish().unwrap();

        let mut var_path = path.join("var.tsv.gz");
        let var = File::create(var_path).unwrap();
        let mut encoder = GzEncoder::new(var, Compression::default());
        let genes = self.cells[0].genes.clone();
        for gene in genes {
            encoder.write(
                format!("{}\t{}\t{}\n", gene.id, gene.symbol, gene.status).as_bytes()
            );
        }
        encoder.finish().unwrap();
 
    } 
}
