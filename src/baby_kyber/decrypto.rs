use crate::baby_kyber::keygen::PrivateKey;
use crate::baby_kyber::{BabyKyber, Ciphtertexts};
use crate::matrix::Matrix;
use crate::poly::Polynomial;
use crate::ring::{PolynomialRingTrait, Ring};
use crate::utils::{bits_normalize, bytes_to_le_bits};
use rand::RngCore;
use std::ops::{Div, Sub};

impl<P: PolynomialRingTrait, const MATRIX_DIMENSION: usize> BabyKyber<P, MATRIX_DIMENSION> {
    // Kyber ciphtertexts consist of those two values: (u,v)
    pub fn decrypto(msg: u8, ciphter_text: Ciphtertexts<P>, private_key: &PrivateKey<P>) -> bool {
        // 1. First, we compute a noisy result m_n
        //      m_n = v - s^T * u
        let m_n = {
            let s_tranpose = {
                let s = private_key.s.clone();

                let s_transpose = s.transpose();

                s_transpose
            };
            let s_tranpose_dot_u = s_tranpose * ciphter_text.u;
            assert!(
                s_tranpose_dot_u.rows == s_tranpose_dot_u.cols && s_tranpose_dot_u.rows == 1,
                "s_tranpose_dot_u is the result of inner product of two vector"
            );
            ciphter_text.v - s_tranpose_dot_u.values[0][0].clone()
        };
        // decrypto.m_n: "9x^3 + 10x + 10"
        println!("decrypto.m_n: {:?}", m_n.to_string());

        // 2. round and scalar m_n
        let scalaed_m_n = {
            println!("decrypto.m_n.coefficients(): {:?}", m_n.coefficients());
            // 2.1 round m_n and
            let rounded_coeffs = Self::round(m_n.coefficients());
            println!("decrypto.m_n.rounded_coeffs(): {:?}", rounded_coeffs);
            // 2.2 mul 1/scalar to recover the msg_bin
            let scalar = Self::compute_scalar();
            let scalared_coeffs = rounded_coeffs
                .into_iter()
                .map(|m| m.div(scalar.clone()))
                .collect::<Vec<_>>();
            println!("decrypto.m_n.scalared_coeffs(): {:?}", scalared_coeffs);

            P::from_coefficients(scalared_coeffs)
        };
        println!("decrypto.scalaed_m_n: {:?}", scalaed_m_n.to_string());

        // 4. check msg_bin,
        let expect_msg_bin = {
            // 4.1 compute msg_bin
            // 1.1 encode msg into bits
            let msg_bits = bytes_to_le_bits(&vec![msg]);
            let msg_bits = bits_normalize(&msg_bits);
            // 1.2. generate msg poly according the msg_bits
            let coeffs = msg_bits
                .into_iter()
                .map(|bit| {
                    if bit {
                        P::PolyCoeff::one()
                    } else {
                        P::PolyCoeff::zero()
                    }
                })
                .collect::<Vec<_>>();
            P::from_coefficients(coeffs)
        };
        println!("expect_msg_bin: {:?}", expect_msg_bin.to_string());
        // 5. final check
        scalaed_m_n == expect_msg_bin
    }

    // Round the value
    //    closer  module/2 -> scalar
    //    closer  0 or p -> 0
    pub fn round(coeffs: Vec<P::PolyCoeff>) -> Vec<P::PolyCoeff> {
        let scalar = Self::compute_scalar();
        let half_scalar = Self::compute_half_scalar();
        // println!("scalar: {:?}", scalar.to_string());
        // println!("half_scalar: {:?}", half_scalar.to_string());
        coeffs
            .into_iter()
            .map(|c| {
                if c == scalar {
                    return scalar;
                }

                let scalar_diff = if c > scalar {
                    c.clone() - scalar.clone()
                } else {
                    scalar.clone() - c.clone()
                };
                if scalar_diff > half_scalar {
                    P::PolyCoeff::zero()
                } else {
                    scalar.clone()
                }
            })
            .collect::<Vec<_>>()
    }
}
