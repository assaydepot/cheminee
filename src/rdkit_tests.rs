#[cfg(test)]
mod tests {
    use rdkit_sys::molecule::*;
    #[test]
    fn test_fragment_parent() {
        let orig_smiles = "CCCC.CC";
        let mut pkl_mol = Molecule::new(orig_smiles, "").unwrap();
        pkl_mol.fragment_parent("");
        let parent_smiles = pkl_mol.get_smiles("");
        println!("{:?}", parent_smiles);
        assert_eq!(pkl_mol.get_smiles(""), "CCCC");
    }

    #[test]
    fn test_canonical_tautomer() {
        let orig_smiles = "Oc1c(cccc3)c3nc2ccncc12";
        let mut pkl_mol = Molecule::new(orig_smiles, "").unwrap();
        pkl_mol.canonical_tautomer("");
        let can_taut_smiles = pkl_mol.get_smiles("");
        println!("{:?}", can_taut_smiles);
        assert_eq!(can_taut_smiles, "O=c1c2ccccc2[nH]c2ccncc12");
    }

    #[test]
    fn test_rdkit_fingerprint() {
        let orig_smiles = "Oc1c(cccc3)c3nc2ccncc12";
        let mut pkl_mol = Molecule::new(orig_smiles, "").unwrap();
        let rdkit_fp = pkl_mol.get_rdkit_fp("");
        println!("{:?}", rdkit_fp);
        assert_eq!(rdkit_fp, "10000110000100000000000001000000001000000000100000000001010001000110000000100001000101000000000010000101100000000000001000001000001100000000001001000000010000000100001000110100000000000001000100110010001000011010000000010000000000010100000000100100000000000000100001001010000000100100010000001011001000000100000000000000000000000111000100000000000001000000000000100000100000000010000000100000000000000000001000000100000001000001001000000011000001000010000000000100001010011011000000000000000010000010011000000100000100000010000000000010000000000001010000000001000100000000010000000000000000000000000000010000101100000001000000100000110000001100000010000000000000000000100100100100000111000001000001000000000000001000100100000000111100000010001000000001000000000000001100000000000010010001000111000001000100000000010000000000000010010000011001100100010001000011000000100010001010010001000000000000000000001000100000000010010100010100011000001000000000001010100110100100000000000100001000100000000000000000100000000000100000000011000000000100001000000001001000001000000000000101000000010000100000100100001000100001000001000100000001100100000000000101000011011000000000100001000010000000100000000001000001001000100000000000000010000000000001100000010000000001000000000100000000000001000000011000000000000000011100000100000100000000000000001000011001000000000100000001000011110001010000000000001010000000100000100001010010101000000000101000000100010000000100001000000000000000000001010011010100001110000011000000000100010000101000000000010000000000000000000000000000000000100000000011000010010000000001000000101000001010010000000101000000000100001000000010000010000000000010000000100011100100001000100000000000000000001000000000000100001010000000000110000000000000001010000000000000100000100010001000000000010000000010010100100100000000100000000000001000000000000101000001000000000000000000000110100011000000000100000000000000000100001100000011010001000000100000100010000000001011100011100000000100010000000100000000000000001000010000000000000001000000")
    }

    #[test]
    fn test_substruct_match() {
        let mol = Molecule::new("CCCCCCC", "").unwrap();
        let query_mol = Molecule::get_qmol("CCCC", "").unwrap();
        let res = mol.get_substruct_match(&query_mol, "");
        println!("{:?}", res);
        assert!(res.len() > 2);
    }

    #[test]
    fn test_tanimoto_similarity() {
        let smiles1 = "c1ccccc1CCCCCCCC";
        let mol1 = Molecule::new(smiles1, "").unwrap();
        let fp1 = mol1.get_rdkit_fp("");
        let fp1 = fp1.as_bytes();

        let smiles2 = "c1ccccc1CCCCCC";
        let mol2 = Molecule::new(smiles2, "").unwrap();
        let fp2 = mol2.get_rdkit_fp("");
        let fp2 = fp2.as_bytes();

        let mut num_matches = 0;
        let mut fp1_ones = 0;
        let mut fp2_ones = 0;
        for idx in 0..fp1.len() {
            if fp1[idx] == 49 { fp1_ones += 1; }
            if fp2[idx] == 49 { fp2_ones += 1; }
            if fp1[idx] == 49 && fp2[idx] == 49 { num_matches += 1; }
        }

        let num_matches = num_matches as f32;
        let fp1_ones = fp1_ones as f32;
        let fp2_ones = fp2_ones as f32;

        let tanimoto = num_matches / (fp1_ones + fp2_ones - num_matches);

        println!("{:?}", tanimoto);
        assert!(tanimoto > 0.94);
    }
}