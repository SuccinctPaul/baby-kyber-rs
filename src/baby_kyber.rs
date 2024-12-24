use crate::baby_kyber::keygen::{PrivateKey, PublickKey};
use crate::poly::Polynomial;
use crate::poly::ring_poly::RingPolynomial;
use crate::ring::Ring;
use rand::RngCore;
use std::marker::PhantomData;

mod constant;
pub mod decrypto;
pub mod encrypto;
pub mod keygen;

pub struct BabyKyber<P: Polynomial> {
    // dimension of matrices
    pub dimension: usize,
    // Degree of polynomial (256 in kyber)
    pub degree: usize,
    phantom_data: PhantomData<P>,
}

impl<P: Polynomial> BabyKyber<P> {
    pub fn init(dimension: usize, degree: usize) -> BabyKyber<P> {
        Self {
            dimension,
            degree,
            phantom_data: Default::default(),
        }
    }

    pub fn keygen(&self, rng: &mut impl RngCore) -> (PrivateKey<P>, PublickKey<P>) {
        let sk = PrivateKey::<P>::new(rng, self.dimension, self.degree);
        let pk = PublickKey::from_private(rng, self.dimension, self.degree, &sk);
        (sk, pk)
    }

    pub fn compute_scalar() -> P::Coefficient {
        let target = P::Coefficient::MODULUS / 2;

        let value = if target % 2 == 1 { target } else { target + 1 };
        P::Coefficient::from(value)
    }
}

pub struct Ciphtertexts<P: Polynomial> {
    pub u: Vec<P>,
    pub v: P,
}
