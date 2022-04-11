use serde_json::Value;
use super::prelude::*;

pub const NAME: &'static str = "index-pubchem-sdf";

pub fn command() -> Command<'static> {
    Command::new(NAME)
        .arg(Arg::new("sdf").required(true).long("sdf").takes_value(true))
        .arg(
            Arg::new("index")
                .required(true)
                .long("index")
                .short('i')
                .takes_value(true),
        )
        .arg(
            Arg::new("limit")
                .required(false)
                .long("limit")
                .takes_value(true),
        )
}

pub fn action(matches: &ArgMatches) -> eyre::Result<()> {
    let sdf_path = matches.value_of("sdf").unwrap();
    let index_dir = matches.value_of("index").unwrap();
    let limit = matches.value_of("limit");

    let mol_iter = MolBlockIter::from_gz_file(sdf_path)
        .map_err(|e| eyre::eyre!("could not read gz file: {:?}", e))?;

    let mol_iter: Box<dyn Iterator<Item = _>> = if let Some(limit) = limit {
        Box::new(mol_iter.take(limit.parse()?))
    } else {
        Box::new(mol_iter)
    };

    let (schema, index) = create_index(index_dir)?;

    let mut index_writer = index.writer_with_num_threads(1, 50 * 1024 * 1024)?;

        for mol_block in mol_iter {
            let mol = match Molecule::new(&mol_block, "") {
                Some(m) => m,
                None => continue,
            };

            let smile = schema.get_field("smile").unwrap();
            let descriptors = schema.get_field("descriptors").unwrap();
            let smarts = schema.get_field("smarts").unwrap();
            let cxsmiles = schema.get_field("cxsmiles").unwrap();
            let coords = schema.get_field("coords").unwrap();
            let numbonds = schema.get_field("numbonds").unwrap();
            let numatoms = schema.get_field("numatoms").unwrap();
            let inchikey = schema.get_field("inchikey").unwrap();
            let jsonMol = schema.get_field("jsonMol").unwrap();
            let CrippenClogP = schema.get_field("CrippenClogP").unwrap();
            let CrippenMR = schema.get_field("CrippenMR").unwrap();
            let FractionCSP3 = schema.get_field("FractionCSP3").unwrap();
            let chi0n = schema.get_field("chi0n").unwrap();
            let kappa1 = schema.get_field("kappa1").unwrap();
            let kappa2 = schema.get_field("kappa2").unwrap();
            let kappa3 = schema.get_field("kappa3").unwrap();
            let labuteASA = schema.get_field("labuteASA").unwrap();

            // println!("{:?} descriptors", descriptors);
            fn flatten<T>(nested: Vec<Vec<T>>) -> Vec<T> {
                nested.into_iter().flatten().collect()
            }

            let json: serde_json::Value = serde_json::from_str(&mol.get_descriptors())?;
            println!("{:?} json", json);
            let jsonMolBonds = (&mol.get_JsonMolecule()).bonds.len();
            println!("{:?} jsonMolBonds", jsonMolBonds);
            let jsonMolAtoms = (&mol.get_JsonMolecule()).atoms.len();
            println!("{:?} jsonMolAtoms", jsonMolAtoms);
            // let coord = (mol.get_coords().iter().map(|c| c.clone()).collect::<Vec<f32>>());
            let coord = flatten(mol.get_coords()).clone();
            println!("{:?} coord", coord);
            // let clog = json.get("CrippenClogP").unwrap().clone();
            let mut clog = json.as_object().unwrap()["CrippenClogP"].clone();
            println!("{:?} clog ", clog);
            for (name, obj) in json.as_object().unwrap().iter() {
                println!("{} is {:?}", name, obj);
            }

            let crip = json.as_object().unwrap()["CrippenMR"].clone();
            println!("{:?} crip ", crip);
            let fract = json.as_object().unwrap()["FractionCSP3"].clone();
            let chi = json.as_object().unwrap()["chi0n"].clone();
            let k1 = json.as_object().unwrap()["kappa1"].clone();
            let k2 = json.as_object().unwrap()["kappa2"].clone();
            let k3 = json.as_object().unwrap()["kappa2"].clone();
            let lab = json.as_object().unwrap()["labuteASA"].clone();
            println!("{:?} chi {:?} kappa1 {:} kappa2 {:?} lab" , chi, k1, k2, lab);

            println!("{:?} smarts",mol.get_smarts("") );
            println!("{:?} get_cxsmiles",mol.get_cxsmiles("") );
            let doc = doc!(
                smile => mol.get_smiles(""),
                descriptors => json,
                smarts => mol.get_smarts(""),
                cxsmiles => mol.get_cxsmiles(""),
                coords => format!("{:?}", coord).as_str(), //.flatten().join(", "),
                numatoms => jsonMolAtoms.to_string().as_str(),
                numbonds => jsonMolBonds.to_string().as_str(),
                inchikey => mol.get_inchikey(""),
                jsonMol => &*mol.get_JsonMolecule().name.as_str(),
                CrippenClogP => clog.to_string(),
                CrippenMR => crip.to_string(),
                FractionCSP3 => fract.to_string(),
                chi0n => chi.to_string(),
                kappa1 => k1.to_string(),
                kappa2 => k2.to_string(),
                kappa3 => k3.to_string(),
                labuteASA => lab.to_string(),
        );
            println!("{:?} doc ", doc);

            index_writer.add_document(doc)?;
        }

    index_writer.commit()?;

    Ok(())
}
