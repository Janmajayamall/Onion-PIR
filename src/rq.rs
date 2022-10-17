use crate::{ntt::NttOperator, poly::Modulus, rns::RnsContext};
use itertools::izip;
use ndarray::Array2;
use rand::{CryptoRng, RngCore};
use std::{
    ops::{AddAssign, MulAssign, SubAssign},
    sync::Arc,
};

// Context
pub struct RqContext {
    pub moduli_64: Vec<u64>,
    pub moduli: Vec<Modulus>,
    pub rns: RnsContext,
    pub ntt_ops: Vec<NttOperator>,
    pub degree: usize,
}

impl RqContext {
    pub fn new(moduli_64: Vec<u64>, degree: usize) -> Self {
        let rns = RnsContext::new(moduli_64.clone());
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

#[derive(Debug, PartialEq)]
enum Representation {
    PowerBasis,
    Ntt,
    NttShoup,
}

struct Poly {
    context: Arc<RqContext>,
    representation: Representation,
    coefficients: Array2<u64>,
}

impl Poly {
    pub fn new() {}

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
