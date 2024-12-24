use crate::baby_kyber::constant::ERROR_POLY_DEGREE;
use crate::baby_kyber::keygen::PublickKey;
use crate::baby_kyber::{BabyKyber, Ciphtertexts};
use crate::matrix::Matrix;
use crate::matrix::poly_matrix::PolyMatrix;
use crate::poly::{Polynomial, small_poly_vector};
use crate::ring::Ring;
use crate::utils::{bits_normalize, bytes_to_le_bits};
use rand::RngCore;
use std::ops::Mul;

impl<P: Polynomial> BabyKyber<P> {
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
            println!("msg_bits {:?}", msg_bits);
            // 1.2. generate msg poly according the msg_bits
            let msg_poly_coeffs = msg_bits
                .into_iter()
                .map(|bit| {
                    if bit {
                        P::Coefficient::one()
                    } else {
                        P::Coefficient::zero()
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
        println!("poly_m: {:?}", poly_m.to_string());

        // 2. generate random poly and error poly
        //      These polynomial vectors are freshly generated for every encryption.
        let (e1, e2, r) = {
            let poly_e1 = small_poly_vector(rng, self.dimension, self.degree);
            // 1 + x
            let poly_e2 = P::from_coefficients(vec![P::Coefficient::one(); 2]);
            // let poly_e2 =
            //     P::from_coefficients(vec![P::Coefficient::from(15), P::Coefficient::from(2)]);
            // let poly_e2 = P::rand(rng, self.degree);
            let poly_r = small_poly_vector(rng, self.dimension, self.degree);
            // println!("poly_e1:{:?}", poly_e1.to_string());
            // println!("poly_e2:{:?}", poly_e2.to_string());
            // println!("poly_r:{:?}", poly_r.to_string());

            let r = PolyMatrix::from_vector_as_col(poly_r);
            let e1 = PolyMatrix::from_vector_as_col(poly_e1);

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
            // assert_eq!(u.rows, self.dimension, "u.len != dimension");
            // assert_eq!(u.cols, self.dimension, "u.len != dimension");

            // 3.2 compute  v=t^T * r + e2 + m
            let v = {
                let t_tranpose = t.transpose();
                let t_tranpose_mul_r = t_tranpose * r.clone();
                // assert_eq!(t_dot_r.len(), 1, "t_dot_r.len() != 1");
                // t_dot_r.u: "4x^3 + 14x^2 + 16x + 6"
                // println!("");
                // println!("t_tranpose_mul_r[0]:{:?}", t_tranpose_mul_r.to_string());
                // t_tranpose_mul_r[0]:"2x^3 + 7x^2 + 8x + 3"
                let t_mul_r_plus_e2 = t_tranpose_mul_r.values[0][0].clone() + e2;
                // println!("t_mul_r_plus_e2: {:?}", t_mul_r_plus_e2.to_string());
                // t_mul_r_plus_e2: "2x^3 + 7x^2 + 10x + 8"
                t_mul_r_plus_e2 + poly_m
            };
            (u, v)
        };

        Ciphtertexts { u, v }
    }
}
