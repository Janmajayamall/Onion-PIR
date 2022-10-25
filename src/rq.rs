use crate::{
    ntt::NttOperator,
    poly::{self, Modulus},
    rns::{RnsContext, RnsScaler, ScalingFactor},
    utils::sample_vec_cbd,
};

use itertools::{izip, Itertools};
use ndarray::{s, Array2, ArrayView2, Axis};
use num_bigint::BigUint;
use num_bigint_dig::BigUint as BigUintDig;
use num_traits::ToPrimitive;
use rand::{thread_rng, CryptoRng, RngCore, SeedableRng};
use rand_chacha::ChaCha8Rng;
use sha2::{Digest, Sha256};
use std::{
    clone,
    error::Error,
    fmt::Debug,
    io::Repeat,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
    sync::Arc,
};

/// Substitution exponent
/// substitute x -> x^k, where
/// k is the exponent
///
/// TODO: add support for substitution in NTT form
#[derive(Clone, Debug)]
pub struct Substitution {
    exponent: usize,
}

impl Substitution {
    pub fn new(exponent: usize) -> Self {
        // TODO: Why is this necessary ? https://github.com/tlepoint/fhe.rs/blob/8aafe4396d0b771e6aa25257c7daa61c109eb367/crates/fhe-math/src/rq/mod.rs#L54

        Substitution { exponent }
    }
}

/// Polt scaler
#[derive(Debug, Clone, PartialEq)]
pub struct RqScaler {
    from: Arc<RqContext>,
    to: Arc<RqContext>,
    rns_scaler: RnsScaler,
    scaling_factor: ScalingFactor,
    number_common_moduli: usize,
}

impl RqScaler {
    pub fn new(from: &Arc<RqContext>, to: &Arc<RqContext>, scaling_factor: ScalingFactor) -> Self {
        assert!(from.degree == to.degree);

        let mut number_common_moduli = 0;
        if scaling_factor.is_one() {
            number_common_moduli = izip!(&from.moduli_64, &to.moduli_64)
                .into_iter()
                .fold(0, |acc, (from, to)| if from == to { acc + 1 } else { acc });
        }

        let rns_scaler = RnsScaler::new(&from.rns, &to.rns, scaling_factor.clone());
        Self {
            from: from.clone(),
            to: to.clone(),
            rns_scaler,
            scaling_factor,
            number_common_moduli,
        }
    }

    /// Scale a polynomial
    pub(crate) fn scale(&self, p: &Poly) -> Poly {
        assert!(p.context.as_ref() == self.from.as_ref());
        assert!(p.representation == Representation::PowerBasis);

        let mut representation = p.representation.clone();
        if representation == Representation::NttShoup {
            representation = Representation::Ntt;
        }

        let mut new_coefficients = Array2::<u64>::zeros((self.to.moduli.len(), self.to.degree));

        if self.number_common_moduli > 0 {
            new_coefficients
                .slice_mut(s![..self.number_common_moduli, ..])
                .assign(&p.coefficients.slice(s![..self.number_common_moduli, ..]));
        }

        if self.number_common_moduli < self.to.moduli_64.len() {
            izip!(
                new_coefficients
                    .slice_mut(s![self.number_common_moduli.., ..])
                    .axis_iter_mut(Axis(1)),
                p.coefficients.axis_iter(Axis(1))
            )
            .for_each(|(new_column, column)| {
                self.rns_scaler
                    .scale(column, new_column, self.number_common_moduli)
            });
        }

        Poly {
            context: self.to.clone(),
            coefficients: new_coefficients,
            representation: Representation::PowerBasis,
        }
    }
}

// Context
#[derive(PartialEq, Clone)]
pub struct RqContext {
    pub moduli_64: Vec<u64>,
    pub moduli: Vec<Modulus>,
    pub rns: Arc<RnsContext>,
    pub ntt_ops: Vec<NttOperator>,
    pub degree: usize,
}

impl Debug for RqContext {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Context")
            .field("moduli", &self.moduli_64)
            // .field("q", &self.q)
            // .field("rns", &self.rns)
            // .field("ops", &self.ops)
            // .field("degree", &self.degree)
            // .field("bitrev", &self.bitrev)
            // .field("inv_last_qi_mod_qj", &self.inv_last_qi_mod_qj)
            // .field("inv_last_qi_mod_qj_shoup", &self.inv_last_qi_mod_qj_shoup)
            // .field("next_context", &self.next_context)
            .finish()
    }
}

impl RqContext {
    pub fn new(moduli_64: Vec<u64>, degree: usize) -> Self {
        let rns = Arc::new(RnsContext::new(&moduli_64));
        let moduli = rns.moduli.clone();
        let ntt_ops = moduli.iter().map(|m| NttOperator::new(m, degree)).collect();

        Self {
            moduli_64,
            moduli,
            rns,
            ntt_ops,
            degree,
        }
    }
}

#[derive(Debug, PartialEq, Clone)]
pub enum Representation {
    PowerBasis,
    Ntt,
    NttShoup,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Poly {
    pub context: Arc<RqContext>,
    pub representation: Representation,
    coefficients: Array2<u64>,
}

impl Poly {
    pub fn coefficients(&self) -> ArrayView2<u64> {
        self.coefficients.view()
    }

    pub fn change_representation(&mut self, to: Representation) {
        match self.representation {
            Representation::PowerBasis => match to {
                Representation::Ntt => {
                    izip!(
                        self.coefficients.outer_iter_mut(),
                        self.context.ntt_ops.iter()
                    )
                    .for_each(|(mut coefficients, ntt_ops)| {
                        ntt_ops.forward(coefficients.as_slice_mut().unwrap());
                    });
                    self.representation = Representation::Ntt;
                }
                Representation::PowerBasis => {}
                _ => {
                    todo!()
                }
            },
            Representation::Ntt => match to {
                Representation::PowerBasis => {
                    izip!(
                        self.coefficients.outer_iter_mut(),
                        self.context.ntt_ops.iter()
                    )
                    .for_each(|(mut coefficients, ntt_ops)| {
                        ntt_ops.backward(coefficients.as_slice_mut().unwrap());
                    });
                    self.representation = Representation::PowerBasis;
                }
                Representation::Ntt => {}
                _ => {
                    todo!()
                }
            },
            _ => {
                todo!()
            }
        }
    }

    pub fn substitute(&self, a: &Substitution) -> Poly {
        let mut poly = Poly::zero(&self.context, self.representation.clone());

        match self.representation {
            Representation::PowerBasis => {
                let mut power = 0usize;
                let mask = self.context.degree - 1;
                for j in 0..self.context.degree {
                    izip!(
                        self.context.moduli.iter(),
                        self.coefficients().slice(s![.., j]),
                        poly.coefficients.slice_mut(s![.., power & mask]),
                    )
                    .for_each(|(q, o_p, t_p)| {
                        // Subtract degree coefficient when ceil(power / self.context.degree) is odd,
                        // otherwise add.
                        // Note that: X^N = -1. Therefore, substitution by `k` is
                        // (X^N)^k == -1^k
                        if (power & self.context.degree) != 0 {
                            *t_p = q.sub(*t_p, *o_p);
                        } else {
                            *t_p = q.add(*t_p, *o_p);
                        }
                    });
                    power += a.exponent;
                }
            }
            _ => {
                //FIXME: quick hack
                let mut p = self.clone();
                p.change_representation(Representation::PowerBasis);
                p = p.substitute(a);
                p.change_representation(self.representation.clone());

                poly = p;
            }
        }
        poly
    }

    pub fn zero(ctx: &Arc<RqContext>, representation: Representation) -> Poly {
        Poly {
            context: ctx.clone(),
            representation,
            coefficients: Array2::zeros((ctx.moduli.len(), ctx.degree)),
        }
    }

    pub fn inverse(&self) -> Poly {
        let mut poly = Poly::zero(&self.context, self.representation.clone());
        izip!(
            poly.coefficients.outer_iter_mut(),
            self.coefficients().outer_iter(),
            self.context.moduli.iter()
        )
        .for_each(|(mut inv_p, p, q)| {
            inv_p
                .as_slice_mut()
                .unwrap()
                .copy_from_slice(&q.inv_vec(p.as_slice().unwrap()))
        });
        poly
    }

    pub fn random<R: RngCore + CryptoRng>(
        ctx: &Arc<RqContext>,
        rng: &mut R,
        representation: Representation,
    ) -> Poly {
        let mut poly = Poly::zero(ctx, representation);
        izip!(poly.coefficients.outer_iter_mut(), ctx.moduli.iter()).for_each(
            |(mut coeff_vec, q)| {
                coeff_vec
                    .as_slice_mut()
                    .unwrap()
                    .copy_from_slice(q.random_vec(ctx.degree, rng).as_slice())
            },
        );
        poly
    }

    ///
    /// Ref - https://github.com/Janmajayamall/fhe.rs/blob/8aafe4396d0b771e6aa25257c7daa61c109eb367/crates/fhe-math/src/rq/mod.rs#L243
    pub fn random_from_seed(
        ctx: &Arc<RqContext>,
        representation: Representation,
        seed: <ChaCha8Rng as SeedableRng>::Seed,
    ) -> Poly {
        // hash seed into a ChaCha8Rng seed.
        let mut hasher = Sha256::new();
        hasher.update(seed);

        let mut prng =
            ChaCha8Rng::from_seed(<ChaCha8Rng as SeedableRng>::Seed::from(hasher.finalize()));
        let mut poly = Poly::zero(ctx, representation);
        izip!(poly.coefficients.outer_iter_mut(), ctx.moduli.iter()).for_each(
            |(mut coeff_vec, q)| {
                coeff_vec
                    .as_slice_mut()
                    .unwrap()
                    .copy_from_slice(q.random_vec(ctx.degree, &mut prng).as_slice())
            },
        );
        poly
    }

    pub fn random_small<R: CryptoRng + RngCore>(
        ctx: &Arc<RqContext>,
        representation: Representation,
        variance: usize,
        rng: &mut R,
    ) -> Poly {
        let poly = Poly::zero(ctx, Representation::PowerBasis);

        let values = sample_vec_cbd(ctx.degree, variance, rng);
        let mut poly = Poly::try_from_vec_i64(ctx, values.as_slice());
        poly.change_representation(representation);
        poly
    }

    pub fn try_from_vec_u64(ctx: &Arc<RqContext>, a: &[u64]) -> Poly {
        let mut poly = Poly::zero(ctx, Representation::PowerBasis);
        izip!(poly.coefficients.outer_iter_mut(), ctx.moduli.iter()).for_each(|(mut coeffs, q)| {
            let coeffs = coeffs.as_slice_mut().unwrap();
            coeffs[..a.len()].copy_from_slice(&q.reduce_vec_u64(a))
        });
        poly
    }

    pub fn try_from_vec_i64(ctx: &Arc<RqContext>, a: &[i64]) -> Poly {
        let mut poly = Poly::zero(ctx, Representation::PowerBasis);
        izip!(poly.coefficients.outer_iter_mut(), ctx.moduli.iter()).for_each(|(mut coeffs, q)| {
            let coeffs = coeffs.as_slice_mut().unwrap();
            coeffs[..a.len()].copy_from_slice(&q.reduce_vec_i64(a))
        });
        poly
    }

    pub fn try_from_bigint(ctx: &Arc<RqContext>, a: &[BigUint]) -> Poly {
        let mut poly = Poly::zero(ctx, Representation::PowerBasis);
        izip!(poly.coefficients.outer_iter_mut(), ctx.moduli.iter()).for_each(|(mut coeffs, q)| {
            let coeffs = coeffs.as_slice_mut().unwrap();
            coeffs[..a.len()].copy_from_slice(&q.reduce_vec_biguint(a))
        });
        poly
    }
}

// OPS

impl AddAssign<&Poly> for Poly {
    fn add_assign(&mut self, rhs: &Poly) {
        assert!(self.representation == Representation::Ntt);
        assert!(rhs.representation == Representation::Ntt);

        izip!(
            self.coefficients.outer_iter_mut(),
            rhs.coefficients.outer_iter(),
            self.context.moduli.iter()
        )
        .for_each(|(mut a, b, q)| {
            q.add_vec(a.as_slice_mut().unwrap(), b.as_slice().unwrap());
        })
    }
}

impl Add<&Poly> for &Poly {
    type Output = Poly;
    fn add(self, rhs: &Poly) -> Self::Output {
        let mut tmp = self.clone();
        tmp += rhs;
        tmp
    }
}

impl SubAssign<&Poly> for Poly {
    fn sub_assign(&mut self, rhs: &Poly) {
        assert!(self.representation == Representation::Ntt);
        assert!(rhs.representation == Representation::Ntt);

        izip!(
            self.coefficients.outer_iter_mut(),
            rhs.coefficients.outer_iter(),
            self.context.moduli.iter()
        )
        .for_each(|(mut a, b, q)| {
            q.sub_vec(a.as_slice_mut().unwrap(), b.as_slice().unwrap());
        })
    }
}

impl Sub<&Poly> for &Poly {
    type Output = Poly;
    fn sub(self, rhs: &Poly) -> Self::Output {
        let mut tmp = self.clone();
        tmp -= rhs;
        tmp
    }
}

impl MulAssign<&Poly> for Poly {
    fn mul_assign(&mut self, rhs: &Poly) {
        assert!(self.representation == Representation::Ntt);
        assert!(rhs.representation == Representation::Ntt);

        izip!(
            self.coefficients.outer_iter_mut(),
            rhs.coefficients.outer_iter(),
            self.context.moduli.iter()
        )
        .for_each(|(mut a, b, q)| {
            q.mul_vec(a.as_slice_mut().unwrap(), b.as_slice().unwrap());
        })
    }
}

impl MulAssign<&BigUint> for Poly {
    fn mul_assign(&mut self, rhs: &BigUint) {
        dbg!();
        let mut rhs = Poly::try_from_bigint(&self.context, &[rhs.clone()]);

        rhs.change_representation(Representation::Ntt);
        assert!(self.representation == Representation::Ntt);
        *self *= &rhs;
    }
}

impl Mul<&Poly> for &Poly {
    type Output = Poly;
    fn mul(self, rhs: &Poly) -> Self::Output {
        let mut tmp = self.clone();
        tmp *= rhs;
        tmp
    }
}

impl Neg for &Poly {
    type Output = Poly;
    fn neg(self) -> Self::Output {
        let mut tmp = Poly::zero(&self.context, self.representation.clone());
        izip!(
            self.coefficients.outer_iter(),
            tmp.coefficients.outer_iter_mut(),
            self.context.moduli.iter()
        )
        .for_each(|(a, mut b, q)| {
            b.as_slice_mut()
                .unwrap()
                .copy_from_slice(&q.neg_vec(&a.as_slice().unwrap()));
        });
        tmp
    }
}

impl From<&Poly> for Vec<BigUint> {
    fn from(value: &Poly) -> Self {
        value
            .coefficients()
            .axis_iter(Axis(1))
            .map(|rests| value.context.rns.lift(rests))
            .collect()
    }
}

// TODO: change this to try from
impl From<&Poly> for Vec<u64> {
    fn from(value: &Poly) -> Self {
        assert!(value.context.moduli_64.len() == 1);
        value
            .coefficients()
            .axis_iter(Axis(1))
            .map(|rests| value.context.rns.lift(rests).to_u64().unwrap())
            .collect()
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::bfv::{BfvParameters, SecretKey};

    // Moduli to be used in tests.
    const MODULI: &[u64; 1] = &[
        1153,
        // 4611686018326724609,
        // 4611686018309947393,
        // 4611686018232352769,
        // 4611686018171535361,
    ];

    #[test]
    fn trial() {
        let mut rng = thread_rng();
        let values = sample_vec_cbd(8, 10, &mut rng);

        let rq = Arc::new(RqContext::new(vec![1153u64], 8));
        let v = Poly::try_from_vec_i64(&rq, &values);
        dbg!(v);
    }

    #[test]
    fn change_representation() {
        for _ in 0..100 {
            let params = Arc::new(BfvParameters::default(1, 8));
            let rq = params.rq_context.clone();
            let mut rng = thread_rng();

            let mut p = Poly::random(&rq, &mut rng, Representation::PowerBasis);
            let mut q = p.clone();
            p.change_representation(Representation::Ntt);
            p.change_representation(Representation::PowerBasis);

            assert_eq!(p, q)
        }
    }

    #[test]
    fn substitute() {
        let mut rng = thread_rng();
        let ctx = Arc::new(RqContext::new(MODULI.to_vec(), 8));

        let mut p = Poly::random(&ctx, &mut rng, Representation::PowerBasis);
        let q = p.clone();
        p.change_representation(Representation::Ntt);

        let mut p3 = p.substitute(&Substitution { exponent: 3 });
        p3 = p3.substitute(&Substitution { exponent: 11 });
        p3.change_representation(Representation::PowerBasis);
        dbg!(p3);
        dbg!(q);
    }
}
