pub mod distributions;
pub mod single_cell;

use rand::rngs::StdRng;
use rand::SeedableRng;

use single_cell::dataset::Dataset;

extern crate hdf5;
use hdf5::File;

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
    //write_anndata()?;
    let file = File::create("anndata.h5").unwrap();
    let group = file.create_group("X").unwrap();
    let builder = group.new_dataset_builder();
}

