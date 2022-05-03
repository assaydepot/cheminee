#![allow(soft_unstable)]
#![feature(test)]

extern crate test;
use rdkit_sys::molecule::Molecule;
use test::Bencher;

#[bench]
fn bench_xor_1000_ints(b: &mut Bencher) {
    let smiles1 = "c1ccccc1CCCCCCCC";
    let mol1 = Molecule::new(smiles1, "").unwrap();

    let smiles2 = "c1ccccc1CCCCCC";
    let mol2 = Molecule::new(smiles2, "").unwrap();

    b.iter(|| mol1.get_tanimoto_similarity(&mol2, "", ""));
}

// running 1 test
// test bench_xor_1000_ints ... bench:     132,690 ns/iter (+/- 2,681)

#[bench]
fn bench_molecules(bencher: &mut test::bench::Bencher) {
    let smiles1 = "c1ccccc1CCCCCCCC";

    bencher.iter(|| {
        Molecule::new(smiles1, "").unwrap();
    })
}
