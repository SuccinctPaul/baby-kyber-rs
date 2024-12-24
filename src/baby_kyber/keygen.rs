use crate::baby_kyber::constant::ERROR_POLY_DEGREE;
use crate::matrix::poly_matrix::PolyMatrix;
use crate::matrix::ring_matrix::RingMatrix;
use crate::matrix::vector_arithmatic::VectorArithmatic;
use crate::poly::ring_poly::RingPolynomial;
use crate::poly::{Polynomial, random_poly_vector};
use crate::ring::Ring;
use rand::RngCore;
use std::ops::Mul;

// The private key of a Kyber key pair consists of polynomials with small coefficients
#[derive(Debug, Clone, Eq, PartialEq, Ord, PartialOrd)]
pub struct PrivateKey<P: Polynomial> {
    pub s: Vec<P>,
}

impl<P: Polynomial> PrivateKey<P> {
    pub fn new(rng: &mut impl rand::RngCore, dimension: usize, degree: usize) -> Self {
        let mut s = random_poly_vector(rng, dimension, degree);
        Self { s }
    }
}

// A Kyber public key consists of two elements.
// 1. matrix of random polynomials `A`
// 2. vector of polynomials `t`
#[derive(Debug, Clone, Eq, PartialEq, Ord, PartialOrd)]
pub struct PublickKey<P: Polynomial> {
    // It's a square matrix.
    pub A: PolyMatrix<P>,
    pub t: Vec<P>,
}
impl<P: Polynomial> PublickKey<P> {
    pub fn from_private(
        rng: &mut impl RngCore,
        dimension: usize,
        degree: usize,
        private_key: &PrivateKey<P>,
    ) -> Self {
        let A = Self::random_A(rng, dimension, degree);
        let e = random_poly_vector(rng, dimension, degree);

        // t= A*s + e
        let t = PolyMatrix::vec_add(&A.mul_vector(&private_key.s), &e);
        Self { A, t }
    }

    // A=[
    //      6x^3 + 16x^2 + 16x + 11, 9x^3 + 4x^2 + 6x + 3,
    //      5x^3 + 3x^2 + 10x + 1,   6x^3 +  x^2 + 9x + 15,
    //   ]
    pub fn random_A(rng: &mut impl RngCore, dimension: usize, degree: usize) -> PolyMatrix<P> {
        let mut values = vec![];
        for _ in 0..dimension {
            values.push(random_poly_vector(rng, dimension, degree));
        }

        PolyMatrix {
            rows: dimension,
            cols: dimension,
            values,
        }
    }

    // // To calculate `t`, needs an additional error vector `e`.
    // // e = [x^2, x^2 - x]
    // fn random_e(rng: &mut impl RngCore, dimension: usize, degree: usize) -> Vec<P> {
    //     let mut s = vec![];
    //     for _ in 0..dimension {
    //         s.push(P::rand(rng, degree));
    //     }
    //     s
    // }
}

#[cfg(test)]
mod tests {
    use crate::baby_kyber::keygen::{PrivateKey, PublickKey};

    // #[test]
    // fn test_keygen() {
    //     let dimentsion = 2;
    //     let degree = 4;
    //     let rng = &mut rand::thread_rng();
    //     let private_key = PrivateKey::new(rng, dimentsion, degree);
    //     let public_key = PublickKey::from_private(rng, dimentsion, degree, &private_key);
    // }
}
