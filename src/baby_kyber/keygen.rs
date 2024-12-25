use crate::baby_kyber::constant::ERROR_POLY_DEGREE;
use crate::baby_kyber::utils::small_poly_vector;
use crate::matrix::poly_matrix::PolyMatrix;
use crate::matrix::ring_matrix::RingMatrix;
use crate::matrix::vector_arithmatic::VectorArithmatic;
use crate::poly::Polynomial;
use crate::poly::uni_poly::UniPolynomial;
use crate::ring::Ring;
use rand::RngCore;

// The private key of a Kyber key pair consists of polynomials with small coefficients
#[derive(Debug, Clone, Eq, PartialEq, Ord, PartialOrd)]
pub struct PrivateKey<P: Polynomial> {
    pub s: PolyMatrix<P>,
}

impl<P: Polynomial> PrivateKey<P> {
    pub fn new(rng: &mut impl rand::RngCore, dimension: usize, degree: usize) -> Self {
        let vec = small_poly_vector(rng, dimension, degree);
        let s = PolyMatrix::from_vector_as_col(vec);
        assert_eq!(s.cols, 1, "S.cols != 1");
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
    pub t: PolyMatrix<P>,
}
impl<P: Polynomial> PublickKey<P> {
    pub fn from_private(
        rng: &mut impl RngCore,
        dimension: usize,
        degree: usize,
        private_key: &PrivateKey<P>,
    ) -> Self {
        let A = Self::random_A(rng, dimension, degree);
        let e = {
            let vector = small_poly_vector(rng, dimension, degree);
            PolyMatrix::from_vector_as_col(vector)
        };

        // t= A*s + e
        let A_s = A.mul_matrix(&private_key.s);

        let t = A_s + e;
        assert_eq!(t.cols, 1, "t.cols != 1");
        Self { A, t }
    }

    pub fn random_A(rng: &mut impl RngCore, dimension: usize, degree: usize) -> PolyMatrix<P> {
        let mut values = vec![];
        for _ in 0..dimension {
            values.push(small_poly_vector(rng, dimension, degree));
        }

        PolyMatrix {
            rows: dimension,
            cols: dimension,
            values,
        }
    }
}
