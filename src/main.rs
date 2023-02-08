pub mod distributions;
pub mod single_cell;

use rand::rngs::StdRng;
use rand::SeedableRng;

use single_cell::dataset::Dataset;

use std::path::Path;

/*
fn write_anndata() -> Result<(), E> {
    let file = File::create("anndata.h5").unwrap();
    //let group = file.create_group("X")?;
    //Ok(())
}
// */

fn main() {
    let mut rng = StdRng::seed_from_u64(123);
    let dataset: Dataset = Dataset::new(50, 100, 3, &mut rng);
    let outdir: &Path = Path::new("tenx");
    dataset.write(outdir, &mut rng);
}

