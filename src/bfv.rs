use std::sync::Arc;

use crate::poly::{Context, Modulus};

use super::poly::Poly;
use super::utils::sample_uniform_vec;

struct BfvParameters {
    t: u64,
    q: u64,
    n: usize,
    st_dev: f64,
    poly_ctx: Arc<Context>,
}

struct BfvSecretKey {
    pub params: Arc<BfvParameters>,
    pub poly: Poly,
}

impl BfvSecretKey {
    pub fn new(params: &Arc<BfvParameters>) -> Self {
        let mut sk = Poly::random_gaussian(
            &Arc::new(Context::new(Modulus { q: 2 }, params.n)),
            params.st_dev,
        );
        sk.switch_context(&params.poly_ctx);
        Self {
            params: params.clone(),
            poly: sk,
        }
    }

    pub fn encrypt(&self, plaintext: &BfvPlaintext) {}
}

struct BfvPublicKey {
    pub params: Arc<BfvParameters>,
    pub ciphertext: BfvCipherText,
}

impl BfvPublicKey {
    pub fn new(sk: &BfvSecretKey) -> Self {
        let a = Poly::random_uniform(&sk.params.poly_ctx);
        let e = Poly::random_gaussian(&sk.params.poly_ctx, sk.params.st_dev);

        Self {
            params: sk.params.clone(),
            ciphertext: BfvCipherText::new(&sk.params, vec![-(&a * &sk.poly + e), a]),
        }
    }
}

struct BfvPlaintext {
    pub params: Arc<BfvParameters>,
    pub values: Box<[u64]>,
}

struct BfvCipherText {
    pub params: Arc<BfvParameters>,
    pub c: Vec<Poly>,
}

impl BfvCipherText {
    pub fn new(params: &Arc<BfvParameters>, c: Vec<Poly>) -> Self {
        Self {
            params: params.clone(),
            c,
        }
    }
}
