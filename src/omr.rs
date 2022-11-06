use crate::{
    bfv::BfvCipherText,
    modulus::Modulus,
    ntt::NttOperator,
    pvw::{self, PvwCiphertext, PvwParams},
    rq::Poly,
};

fn inner_product(cts: &Vec<BfvCipherText>, polys: &Vec<Poly>) -> BfvCipherText {
    assert!(cts.len() == polys.len());

    let mut product = &cts[0] * &polys[0];
    for i in 1..cts.len() {
        let p = &cts[i] * &polys[i];
        product = &product + &p;
    }

    product
}

fn process() {
    // collection of clues
    let pvw_params = PvwParams::default();
    let clues = Vec::<PvwCiphertext>::new();

    // Process clues for SIMD ops.
    // For now imagine the N < D.
    let D = 64;

    let q_mod = Modulus::new(pvw_params.q);
    let ntt_operator = NttOperator::new(&q_mod, D);

    let mut b_polys = vec![];
    let mut a_polys = vec![];
    for l_index in 0..pvw_params.ell {
        let mut poly = vec![0u64; D];
        for i in 0..clues.len() {
            poly[i] = clues[i].b[l_index];
        }

        // perform SIMD encoding
        ntt_operator.backward(&mut poly);
        b_polys.push(poly);
    }

    for n_index in 0..pvw_params.n {
        let mut poly = vec![0u64; D];
        for i in 0..clues.len() {
            poly[i] = clues[i].a[n_index];
        }

        // perform SIMD encoding
        ntt_operator.backward(&mut poly);
        a_polys.push(poly);
    }
}

// #[tests]
// mod tests {}
