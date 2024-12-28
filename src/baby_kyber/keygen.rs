use crate::baby_kyber::utils::small_poly_vector;
use crate::matrix::poly_ring_matrix::PolyRingMatrix;
use crate::ring::PolynomialRingTrait;
use rand::RngCore;

// The private key of a Kyber key pair consists of polynomials with small coefficients
#[derive(Debug, Clone, Eq, PartialEq, Ord, PartialOrd)]
pub struct PrivateKey<P: PolynomialRingTrait> {
    // It's a col_vector.
    pub s: PolyRingMatrix<P>,
}

impl<P: PolynomialRingTrait> PrivateKey<P> {
    pub fn new(rng: &mut impl rand::RngCore, dimension: usize) -> Self {
        let vec = small_poly_vector(rng, dimension);
        let s = PolyRingMatrix::from_col_vector(vec);
        assert_eq!(s.cols, 1, "S.cols != 1, it should be a column vector");
        Self { s }
    }
}

// A Kyber public key consists of two elements.
// 1. matrix of random polynomials `A`
// 2. vector of polynomials `t`
#[derive(Debug, Clone, Eq, PartialEq, Ord, PartialOrd)]
pub struct PublickKey<P: PolynomialRingTrait> {
    // It's a square matrix.
    pub A: PolyRingMatrix<P>,
    // It's a col_vector.
    pub t: PolyRingMatrix<P>,
}
impl<P: PolynomialRingTrait> PublickKey<P> {
    pub fn from_private(
        rng: &mut impl RngCore,
        dimension: usize,
        private_key: &PrivateKey<P>,
    ) -> Self {
        // A=n*n
        let A = Self::random_A(rng, dimension);
        // e=1*n
        let e = {
            let vector = small_poly_vector(rng, dimension);
            PolyRingMatrix::from_col_vector(vector)
        };

        // t= A*s + e
        let A_s = A.mul_matrix(&private_key.s);

        let t = A_s + e;
        assert_eq!(t.cols, 1, "cols != 1, it should be a column vector");
        Self { A, t }
    }

    pub fn random_A(rng: &mut impl RngCore, dimension: usize) -> PolyRingMatrix<P> {
        let mut values = vec![];
        for _ in 0..dimension {
            values.push(small_poly_vector(rng, dimension));
        }

        PolyRingMatrix {
            rows: dimension,
            cols: dimension,
            values,
        }
    }
}
