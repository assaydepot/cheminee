use super::prelude::*;
use std::env;
pub const NAME: &str = "index-pubchem-sdf";
// pub const NAME: &'static str = "index-pubchem-sdf";

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
        let mut mol = match Molecule::new(&mol_block, "") {
            Some(m) => m,
            None => continue,
        };
        println!("mol {:?}", mol);
        mol.cleanup("");
        mol.neutralize("");
        mol.reionize("");
        println!("mol  cleanup {:?}", mol);

        mol.fragment_parent("");
        let fr_parent = mol.get_smiles("");
        mol.canonical_tautomer("");
        let canonical = mol.get_smiles("");
        println!("{:?} mol after canonical", mol);

        let smile = schema.get_field("smile").unwrap();
        let descriptors = schema.get_field("descriptors").unwrap();
        let fragment_parent = schema.get_field("fragment_parent").unwrap();
        let canonical_tautomer = schema.get_field("canonical_tautomer").unwrap();
        let smarts = schema.get_field("smarts").unwrap();
        let cxsmiles = schema.get_field("cxsmiles").unwrap();
        // let coords = schema.get_field("coords").unwrap();
        let numbonds = schema.get_field("numbonds").unwrap();
        let numatoms = schema.get_field("numatoms").unwrap();
        let inchikey = schema.get_field("inchikey").unwrap();
        // let jsonMol = schema.get_field("jsonMol").unwrap();
        // let CrippenClogP = schema.get_field("CrippenClogP").unwrap();
        // let CrippenMR = schema.get_field("CrippenMR").unwrap();
        // let FractionCSP3 = schema.get_field("FractionCSP3").unwrap();
        let chi0n = schema.get_field("chi0n").unwrap();
        let kappa1 = schema.get_field("kappa1").unwrap();
        let kappa2 = schema.get_field("kappa2").unwrap();
        let kappa3 = schema.get_field("kappa3").unwrap();
        let labute_asa = schema.get_field("labute_asa").unwrap();

        // println!("{:?} descriptors", descriptors);
        // fn flatten<T>(nested: Vec<Vec<T>>) -> Vec<T> {
        //     nested.into_iter().flatten().collect()
        // }
        //
        let json: serde_json::Value = serde_json::from_str(&mol.get_descriptors())?;
        // println!("{:?} json", json);
        let json_mol_bonds = mol.get_JsonMolecule().bonds.len();
        // println!("{:?} jsonMolBonds", jsonMolBonds);
        let json_mol_atoms = mol.get_JsonMolecule().atoms.len();
        // println!("{:?} jsonMolAtoms", jsonMolAtoms);
        // // let coord = (mol.get_coords().iter().map(|c| c.clone()).collect::<Vec<f32>>());
        // let coord = flatten(mol.get_coords()).clone();
        // println!("{:?} coord", coord);
        // // let clog = json.get("CrippenClogP").unwrap().clone();
        // let mut clog = json.as_object().unwrap()["CrippenClogP"].clone();
        //
        // println!("{:?} clog ", clog);
        // for (name, obj) in json.as_object().unwrap().iter() {
        //     println!("{} is {:?}", name, obj);
        // }

        let chi = json.as_object().unwrap()["chi0n"].clone();
        let k1 = json.as_object().unwrap()["kappa1"].clone();
        let k2 = json.as_object().unwrap()["kappa2"].clone();
        let k3 = json.as_object().unwrap()["kappa2"].clone();
        let lab = json.as_object().unwrap()["labuteASA"].clone();

        use serde_json::json;

        let doc = doc!(
                smile => mol.get_smiles(""),
                descriptors => json,
                fragment_parent => fr_parent,
                canonical_tautomer => canonical,
                smarts => mol.get_smarts(""),
                cxsmiles => mol.get_cxsmiles(""),
                numatoms => json_mol_atoms.to_string().as_str(),
                numbonds => json_mol_bonds.to_string().as_str(),
                inchikey => mol.get_inchikey(""),
                chi0n => json!(chi).as_f64().unwrap(),
                kappa1 => json!(k1).as_f64().unwrap(),
                kappa2 => json!(k2).as_f64().unwrap(),
                kappa3 => json!(k3).as_f64().unwrap(),
                labute_asa => json!(lab).as_f64().unwrap(),
        );
        // println!("{:?} doc ", doc);

        index_writer.add_document(doc)?;
    }

    index_writer.commit()?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bad_mol1() {
        let molblock = "THIOL_12\n     RDKit          3D\n\n 25 25  0  0  0  0  0  0  0  0999 V2000\n   -2.2510   -2.6650   -2.0550 S   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.3040   -2.7120   -2.1100 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.7910   -1.5140   -0.7240 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.1270   -2.0270    0.1920 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.4730   -0.6640   -0.8710 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.2780   -0.7500   -0.3280 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.2420   -1.8480   -0.5140 O   0  0  0  0  0  0  0  0  0  0  0  0\n    0.4860    0.3560   -0.2740 N   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0540    1.2670   -0.1190 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.9050    0.2450   -0.7390 C   0  0  1  0  0  0  0  0  0  0  0  0\n    1.9190   -0.1360   -1.7740 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.4830    1.6820   -0.6980 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.4150    2.0880    0.3240 H   0  0  0  0  0  0  0  0  0  0  0  0\n    3.5420    1.6740   -0.9990 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.7270    2.5420   -1.6810 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.7070    2.3180   -2.9770 N   0  0  0  0  0  0  0  0  0  0  0  0\n    2.2060    1.5550   -3.4590 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.9600    3.1990   -3.5600 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.7540    3.2500   -4.6280 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.5040    3.9880   -2.6970 N   0  0  0  0  0  0  0  0  0  0  0  0\n    0.9760    3.6220   -1.3870 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.7740    4.0890   -0.4240 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.7730   -0.7080    0.0810 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.9440   -0.8660   -0.2210 O   0  0  0  0  0  0  0  0  0  0  0  0\n    2.2450   -1.3190    1.1390 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  1  3  1  0\n  3  4  1  0\n  3  5  1  0\n  3  6  1  0\n  6  7  2  0\n  6  8  1  0\n  8  9  1  0\n  8 10  1  0\n 10 11  1  6\n 10 12  1  0\n 10 23  1  0\n 12 13  1  0\n 12 14  1  0\n 12 15  1  0\n 15 16  1  0\n 15 21  2  0\n 16 17  1  0\n 16 18  1  0\n 18 19  1  0\n 18 20  2  0\n 20 21  1  0\n 21 22  1  0\n 23 24  2  0\n 23 25  1  0\nM  CHG  1  25  -1\nM  END\n";
        let mut pkl_mol = Molecule::new(molblock, "").unwrap();
        println!("1: {:?}", pkl_mol);
        pkl_mol.canonical_tautomer("");
        pkl_mol.cleanup(""); // this is needed to avoid exception...
                             //pkl_mol2.remove_all_hs();
        println!("2: {:?}", pkl_mol);
        println!("2.1: fragment {:?}", pkl_mol.fragment_parent(""));
        println!("3.0: {:?}", pkl_mol.get_smiles(""));
        pkl_mol.canonical_tautomer("");
        println!("3: {:?}", pkl_mol.get_smiles(""));
        // println!("4: {:?}", pkl_mol.canonical_tautomer(""));
        assert_eq!(pkl_mol.get_smiles(""), "O=C(CS)NC(Cc1c[nH]cn1)C(=O)[O-]");
    }

    #[test]
    fn find_substructure() {
        let mol = Molecule::new("Cl[C@H](F)C[C@H](F)Cl", "").unwrap();
        let query_mol = Molecule::get_qmol("Cl[C@@H](F)C", "").unwrap();
        let res = mol.get_substruct_match(&query_mol, "");
        assert_eq!(res, "{\"atoms\":[0,1,2,3],\"bonds\":[0,1,2]}");
        let res = mol.get_substruct_matches(&query_mol, "");
        assert_eq!(
            res,
            "[{\"atoms\":[0,1,2,3],\"bonds\":[0,1,2]},{\"atoms\":[6,4,5,3],\"bonds\":[5,4,3]}]"
        );
    }

    #[test]
    fn normalize() {
        let orig_smiles = "CN=N#N";
        let mut pkl_mol = Molecule::new(orig_smiles, "").unwrap();
        pkl_mol.normalize("");
        let smiles = pkl_mol.get_smiles("");
        assert_eq!(smiles, "CN=[N+]=[N-]");
    }

    #[test]
    fn get_smthg() {
        println!("{}", env::consts::OS);
        let orig_smiles = "Oc1c2ccccc2nc2ccncc12";
        let mut pkl_mol = Molecule::new(orig_smiles, "").unwrap();
        println!("original {:?} ", pkl_mol);
        pkl_mol.neutralize("");
        println!("neutro {:?} ", pkl_mol);
        pkl_mol.cleanup("");
        println!("cleanend {:?} ", pkl_mol);
        pkl_mol.reionize("");
        println!("reionize {:?} ", pkl_mol);
        pkl_mol.normalize("");
        println!("normalized {:?} ", pkl_mol);
        pkl_mol.fragment_parent("");
        println!("parent {:?} ", pkl_mol);
        pkl_mol.canonical_tautomer("");
        println!("canonical {:?} ", pkl_mol);
        println!("get_pattern_fp  {:?} ", pkl_mol.get_pattern_fp_as_bytes(""));
        println!("get_rdkit_fp  {:?} ", pkl_mol.get_rdkit_fp_as_bytes(""));
        println!("get_morgan_fp  {:?} ", pkl_mol.get_morgan_fp_as_bytes(""));
        let smile = pkl_mol.get_smiles("");
        // original Oc1c2ccccc2nc2ccncc12
        // neutro Oc1c2ccccc2nc2ccncc12
        // cleanend Oc1c2ccccc2nc2ccncc12
        // reionize Oc1c2ccccc2nc2ccncc12
        // normalized Oc1c2ccccc2nc2ccncc12
        // parent Oc1c2ccccc2nc2ccncc12
        // canonical O=c1c2ccccc2[nH]c2ccncc12
        assert_eq!(&smile, "O=c1c2ccccc2[nH]c2ccncc12");
    }

    // #[test]
    // fn test_fragment_parent() {
    //     let orig_smiles = "CC";
    //     let json_info = CString::new(json_info).unwrap();
    //     let mut pkl_mol = Molecule::new(orig_smiles, "").unwrap();
    //     let result = pkl_mol.fragment_parent("");
    //     let result1 = pkl_mol.canonical_tautomer("");
    //     println!("{:?} can", result1);
    //     println!(" hhhh");
    //     panic!("{} => {:?}", result, pkl_mol);
    // }
    #[test]
    fn hi() {
        println!("jj");
    }
    #[test]
    fn smiles_from_smiles_via_pkl() {
        let orig_smiles = "OCCC#CO";
        // let orig_smiles = "[H]C([H])([H])Sc1nc(N(C([H])([H])C([H])([H])C([H])([H])C([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])[H])c2sc3nc(C([H])([H])C([H])([H])C([H])([H])C([H])([H])[H])c4c(c3c2n1)C([H])([H])C([H])([H])C4([H])[H]";
        let pkl_mol = Molecule::new(orig_smiles, "").unwrap();
        println!(
            "SMILES: {} Canonical SMILES: {}",
            orig_smiles,
            pkl_mol.get_smiles("")
        );
        assert_eq!(pkl_mol.get_smiles(""), "OC#CCCO");
    }
}
