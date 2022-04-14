use std::path::Path;
use tantivy::directory::MmapDirectory;
use tantivy::schema::*;
use tantivy::{Index, IndexBuilder};

pub use tantivy::doc;

pub fn schema() -> Schema {
    let mut builder = SchemaBuilder::new();

    let score_options = NumericOptions::default()
        .set_indexed()
        .set_fieldnorm()
        .set_fast(Cardinality::SingleValue);
    builder.add_text_field("smile", TEXT | STORED);
    builder.add_json_field("descriptors", TEXT | STORED);
    builder.add_text_field("fragment_parent", TEXT | STORED);
    builder.add_text_field("canonical_tautomer", TEXT | STORED);
    builder.add_text_field("smarts", TEXT | STORED);
    builder.add_text_field("cxsmiles", TEXT | STORED);
    // builder.add_text_field("coords", TEXT | STORED);
    builder.add_text_field("numatoms", TEXT | STORED);
    builder.add_text_field("numbonds", TEXT | STORED);
    builder.add_text_field("inchikey", TEXT | STORED);
    // builder.add_text_field("jsonMol", TEXT | STORED);
    // builder.add_text_field("CrippenClogP", TEXT | STORED);
    // // builder.add_f64_field("CrippenClogP", FAST | STORED);
    // // builder.add_f64_field("CrippenClogP", score_options);
    // builder.add_text_field("CrippenMR", TEXT | STORED);
    // builder.add_text_field("FractionCSP3", TEXT | STORED);
    // builder.add_f64_field("pca_smthg", score_options);
    builder.add_f64_field("chi0n", score_options);
    builder.add_f64_field("kappa1", INDEXED);
    builder.add_f64_field("kappa2", INDEXED);
    builder.add_f64_field("kappa3", INDEXED);
    builder.add_f64_field("labute_asa", INDEXED);

    builder.build()
}

pub fn create_index(p: impl AsRef<Path>) -> eyre::Result<(Schema, Index)> {
    let schema = schema();

    let builder = IndexBuilder::new().schema(schema.clone());

    let index = builder.create_in_dir(p)?;
    println!("{:?} schema {:?} index", schema, index);
    Ok((schema, index))
}

pub fn open_index(p: impl AsRef<Path>) -> eyre::Result<Index> {
    let directory = MmapDirectory::open(p)?;

    let index = Index::open(directory)?;

    Ok(index)
}
