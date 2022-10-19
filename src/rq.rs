use crate::{
    ntt::NttOperator,
    poly::{Context, Modulus},
    rns::{RnsContext, RnsScaler, ScalingFactor},
    utils::sample_vec_cbd,
};
use itertools::izip;

use ndarray::{s, Array2, ArrayView2, Axis};
use num_bigint::BigUint;
use rand::{CryptoRng, RngCore, SeedableRng};
use rand_chacha::ChaCha8Rng;
use sha2::{Digest, Sha256};
use std::{
    clone,
    ops::{AddAssign, MulAssign, SubAssign},
    sync::Arc,
};

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

        let number_common_moduli = izip!(&from.moduli_64, &to.moduli_64)
            .into_iter()
            .fold(0, |acc, (from, to)| if from == to { acc + 1 } else { acc });

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
                dbg!(&column);
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
#[derive(Debug, PartialEq, Clone)]
pub struct RqContext {
    pub moduli_64: Vec<u64>,
    pub moduli: Vec<Modulus>,
    pub rns: Arc<RnsContext>,
    pub ntt_ops: Vec<NttOperator>,
    pub degree: usize,
}

impl RqContext {
    pub fn new(moduli_64: Vec<u64>, degree: usize) -> Self {
        let rns = Arc::new(RnsContext::new(moduli_64.clone()));
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
    pub fn new() {}

    pub fn coefficients(&self) -> ArrayView2<u64> {
        self.coefficients.view()
    }

    pub fn change_representation(&mut self, to: Representation) {
        match self.representation {
            Representation::PowerBasis => match to {
                Representation::Ntt => izip!(
                    self.coefficients.outer_iter_mut(),
                    self.context.ntt_ops.iter()
                )
                .for_each(|(mut coefficients, ntt_ops)| {
                    ntt_ops.forward(coefficients.as_slice_mut().unwrap());
                }),
                Representation::PowerBasis => {}
                _ => {
                    todo!()
                }
            },
            Representation::Ntt => match to {
                Representation::PowerBasis => izip!(
                    self.coefficients.outer_iter_mut(),
                    self.context.ntt_ops.iter()
                )
                .for_each(|(mut coefficients, ntt_ops)| {
                    ntt_ops.backward(coefficients.as_slice_mut().unwrap());
                }),
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

    pub fn zero(ctx: &Arc<RqContext>, representation: Representation) -> Poly {
        Poly {
            context: ctx.clone(),
            representation,
            coefficients: Array2::zeros((ctx.moduli.len(), ctx.degree)),
        }
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

    pub fn random_small(
        ctx: &Arc<RqContext>,
        representation: Representation,
        variance: isize,
    ) -> Poly {
        let poly = Poly::zero(&ctx, Representation::PowerBasis);
        let values = sample_vec_cbd(ctx.degree, variance);
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
