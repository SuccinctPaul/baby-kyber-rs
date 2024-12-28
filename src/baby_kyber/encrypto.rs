use crate::baby_kyber::keygen::PublickKey;
use crate::baby_kyber::utils::{random_poly_vector, small_poly_vector};
use crate::baby_kyber::{BabyKyber, Ciphtertexts};
use crate::debug::debug_poly_matrix;
use crate::matrix::Matrix;
use crate::matrix::poly_ring_matrix::PolyRingMatrix;
use crate::poly::Polynomial;
use crate::ring::{PolynomialRingTrait, Ring};
use crate::utils::{bits_normalize, bytes_to_le_bits};
use rand::RngCore;
use std::ops::Mul;

impl<P: PolynomialRingTrait, const MATRIX_DIMENSION: usize> BabyKyber<P, MATRIX_DIMENSION> {
    // Kyber ciphtertexts consist of those two values: (u,v)
    pub fn encrypto(
        &self,
        rng: &mut impl RngCore,
        msg: u8,
        public_key: &PublickKey<P>,
    ) -> Ciphtertexts<P> {
        // 1. generate msg poly
        let poly_m = {
            // 1.1 encode msg into bits
            let msg_bits = bytes_to_le_bits(&vec![msg]);
            let msg_bits = bits_normalize(&msg_bits);
            // 1.2. generate msg poly according the msg_bits
            let msg_poly_coeffs = msg_bits
                .into_iter()
                .map(|bit| {
                    if bit {
                        P::PolyCoeff::one()
                    } else {
                        P::PolyCoeff::zero()
                    }
                })
                .collect::<Vec<_>>();
            let msg_poly = P::from_coefficients(msg_poly_coeffs);
            // 1.3 scale polynomial by the (q/2)
            let scalar = Self::compute_scalar();
            let scalared_coeffs = msg_poly
                .coefficients()
                .into_iter()
                .map(|coeff| coeff * scalar)
                .collect::<Vec<_>>();
            P::from_coefficients(scalared_coeffs)
        };
        println!("encrypto.poly_m: {:?}", poly_m.to_string());

        // 2. generate random poly and error poly
        //      These polynomial vectors are freshly generated for every encryption.
        let (e1, e2, r) = {
            let poly_e1 = small_poly_vector(rng, MATRIX_DIMENSION);
            let degree = 4;
            let poly_e2 = P::from_coefficients(vec![P::PolyCoeff::one(); degree]);
            let poly_r = small_poly_vector(rng, MATRIX_DIMENSION);

            let r = PolyRingMatrix::from_col_vector(poly_r);
            let e1 = PolyRingMatrix::from_col_vector(poly_e1);

            (e1, poly_e2, r)
        };

        // 3. compute u,v
        let (u, v) = {
            let A_transpose = public_key.A.transpose();
            let t = public_key.t.clone();

            // 3.1 comptue u=A^T * r + e1
            let u = {
                let A_transpose_mul_r = A_transpose * r.clone();
                A_transpose_mul_r + e1.clone()
            };

            // 3.2 compute  v=t^T * r + e2 + m
            let v = {
                let t_tranpose = t.transpose();
                debug_poly_matrix("t_tranpose", &t_tranpose);
                let t_tranpose_mul_r = t_tranpose * r.clone();
                debug_poly_matrix("t_tranpose_mul_r", &t_tranpose_mul_r);
                assert!(
                    t_tranpose_mul_r.rows == t_tranpose_mul_r.cols && t_tranpose_mul_r.rows == 1,
                    "t_tranpose_mul_r is the result of inner product of two vector"
                );
                let t_mul_r_plus_e2 = t_tranpose_mul_r.values[0][0].clone() + e2;

                t_mul_r_plus_e2 + poly_m
            };
            (u, v)
        };

        Ciphtertexts { u, v }
    }
}
