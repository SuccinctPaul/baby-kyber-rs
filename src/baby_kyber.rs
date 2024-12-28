use crate::baby_kyber::keygen::{PrivateKey, PublickKey};
use crate::matrix::poly_ring_matrix::PolyRingMatrix;
use crate::ring::{PolynomialRingTrait, Ring};
use rand::RngCore;
use std::marker::PhantomData;

pub mod decrypto;
pub mod encrypto;
pub mod keygen;
pub mod utils;

pub struct BabyKyber<P: PolynomialRingTrait, const MATRIX_DIMENSION: usize> {
    // dimension of matrices
    // pub dimension: usize,
    // Degree of polynomial (256 in kyber)
    // pub degree: usize,
    phantom_data: PhantomData<P>,
}

impl<P: PolynomialRingTrait, const MATRIX_DIMENSION: usize> BabyKyber<P, MATRIX_DIMENSION> {
    pub fn init() -> BabyKyber<P, MATRIX_DIMENSION> {
        Self {
            phantom_data: Default::default(),
        }
    }

    pub fn keygen(&self, rng: &mut impl RngCore) -> (PrivateKey<P>, PublickKey<P>) {
        let sk = PrivateKey::<P>::new(rng, MATRIX_DIMENSION);
        let pk = PublickKey::from_private(rng, MATRIX_DIMENSION, &sk);

        (sk, pk)
    }

    pub fn compute_scalar() -> P::PolyCoeff {
        let target = P::PolyCoeff::MODULUS / 2;

        let value = if target % 2 == 1 { target } else { target + 1 };
        P::PolyCoeff::from(value)
    }

    pub fn compute_half_scalar() -> P::PolyCoeff {
        let target = P::PolyCoeff::MODULUS / 2;

        let value = if target % 2 == 1 { target } else { target + 1 };

        let harf_value = if value % 2 == 1 {
            value / 2 + 1
        } else {
            value / 2
        };
        println!("harf_value = {}", harf_value);
        P::PolyCoeff::from(harf_value)
    }
}

pub struct Ciphtertexts<P: PolynomialRingTrait> {
    pub u: PolyRingMatrix<P>,
    pub v: P,
}
