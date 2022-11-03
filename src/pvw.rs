use crate::modulus::Modulus;
use ndarray::{Array2, Axis};
use rand::thread_rng;

#[derive(Debug)]
struct PvwParams {
    n: usize,
    ell: usize,
    m: usize,
    q: u64,
    variance: f64,
}

impl Default for PvwParams {
    fn default() -> Self {
        Self {
            n: 450,
            ell: 4,
            m: 16000,
            q: 65537,
            variance: 1.3,
        }
    }
}

impl PvwParams {
    pub fn new() -> PvwParams {
        todo!()
    }
}

struct PvwSk(Array2<u64>);

impl PvwSk {
    pub fn gen_sk(params: &PvwParams) -> PvwSk {
        let mut rng = thread_rng();

        // Sample from uniform distribution that is mod `q`
        // and return 2d matrix
        let mut sk: Array2<u64> = Array2::zeros((params.ell, params.n));
        let q_mod = Modulus::new(params.q);
        sk.outer_iter_mut().for_each(|mut v| {
            v.as_slice_mut()
                .unwrap()
                .copy_from_slice(q_mod.random_vec(params.n, &mut rng).as_slice())
        });

        PvwSk(sk)
    }

    /// Encrypts message vector `m`
    /// and returns Ciphertext
    pub fn encrypt(&self, m: Vec<u64>) {}

    /// Generates public key from sk
    pub fn gen_pk(&self) {}
}

struct PvwPk(Vec<PvwCiphertext>);

struct PvwCiphertext {}

mod tests {
    use super::*;

    #[test]
    fn trial() {
        let params = PvwParams::default();

        PvwSk::gen_sk(&params);
    }
}
