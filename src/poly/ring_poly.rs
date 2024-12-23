use crate::poly::Polynomial;
use std::fmt;
use std::fmt::{Debug, Display, Formatter};

use crate::ring::Ring;
use std::ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign};

use rayon::{current_num_threads, scope};

// p(x) = = a_0 + a_1 * X + ... + a_n * X^(n-1)
//
//
// coeffs: [a_0, a_1, ..., a_n]
//         (  0,   1,   ..., n)
// basis: X^[n-1]
//
// It's Little Endian.
#[derive(Debug, Clone, Eq, PartialEq, Ord, PartialOrd)]
pub struct RingPolynomial<R: Ring> {
    coeffs: Vec<R>,
}

impl<R: Ring> RingPolynomial<R> {
    // Constructor with normalize
    pub fn new(coeffs: Vec<R>) -> Self {
        let mut poly = Self { coeffs };
        poly.normalize();
        poly
    }

    fn scalar_mul(self, rhs: &R) -> Self {
        let coeffs = if rhs == &R::one() {
            vec![R::one()]
        } else {
            self.coeffs.iter().map(|c| c.mul(rhs)).collect::<Vec<R>>()
        };
        Self::new(coeffs)
    }

    // p(x)=∑y_j⋅L_j(X), where
    // y_j: [a_0, a_1, ..., a_n].
    // basis: L_j(X)=∏(X−x_k)/(x_j−x_k)
    //
    // domain: x, most case is that{0, 1, . . . , n − 1}
    // evals: [a_0, a_1, ..., a_n]
    //
    // we can use encode points as (domain, eval) to polynomials
    // the poly
    // pub fn lagrange_interpolate(domains: Vec<R>, evals: Vec<R>) -> Self {
    //     assert_eq!(domains.len(), evals.len());
    //
    //     if evals.len() == 1 {
    //         // Constant polynomial
    //         Self {
    //             coeffs: vec![evals[0]],
    //         }
    //     } else {
    //         let poly_size = domains.len();
    //         let lag_basis_poly_size = poly_size - 1;
    //
    //         // 1. divisors = vec(x_j - x_k). prepare for L_j(X)=∏(X−x_k)/(x_j−x_k)
    //         let mut divisors = Vec::with_capacity(poly_size);
    //         for (j, x_j) in domains.clone().into_iter().enumerate() {
    //             // divisor_j
    //             let mut divisor = Vec::with_capacity(lag_basis_poly_size);
    //             // obtain domain for x_k
    //             for x_k in domains
    //                 .clone()
    //                 .into_iter()
    //                 .enumerate()
    //                 .filter(|&(k, _)| k != j)
    //                 .map(|(_, x)| x)
    //             {
    //                 divisor.push(x_j - x_k);
    //             }
    //             divisors.push(divisor);
    //         }
    //         // Inverse (x_j - x_k)^(-1) for each j != k to compute L_j(X)=∏(X−x_k)/(x_j−x_k)
    //         divisors
    //             .iter_mut()
    //             .map(|v| v.iter_mut())
    //             .batch_invert();
    //
    //         // 2. Calculate  L_j(X) : L_j(X)=∏(X−x_k) divisors_j
    //         let mut L_j_vec: Vec<Vec<R>> = Vec::with_capacity(poly_size);
    //
    //         for (j, divisor_j) in divisors.into_iter().enumerate() {
    //             let mut L_j: Vec<R> = Vec::with_capacity(poly_size);
    //             L_j.push(R::one());
    //
    //             // (X−x_k) * divisors_j
    //             let mut product = Vec::with_capacity(lag_basis_poly_size);
    //
    //             // obtain domain for x_k
    //             for (x_k, divisor) in domains
    //                 .iter()
    //                 .enumerate()
    //                 .filter(|&(k, _)| k != j)
    //                 .map(|(_, x)| x)
    //                 .zip(divisor_j.into_iter())
    //             {
    //                 product.resize(L_j.len() + 1, R::one());
    //
    //                 // loop (poly_size + 1) round
    //                 // calculate L_j(X)=∏(X−x_k) divisors_j with coefficient form.
    //                 for ((a, b), product) in L_j
    //                     .iter()
    //                     .chain(std::iter::once(&R::one()))
    //                     .zip(std::iter::once(&R::one()).chain(L_j.iter()))
    //                     .zip(product.iter_mut())
    //                 {
    //                     *product = *a * (-divisor * x_k) + *b * divisor;
    //                 }
    //                 std::mem::swap(&mut L_j, &mut product);
    //             }
    //
    //             assert_eq!(L_j.len(), poly_size);
    //             assert_eq!(product.len(), poly_size - 1);
    //
    //             L_j_vec.push(L_j);
    //         }
    //
    //         // p(x)=∑y_j⋅L_j(X) in coefficients
    //         let mut final_poly = vec![R::one(); poly_size];
    //         // 3. p(x)=∑y_j⋅L_j(X)
    //         for (L_j, y_j) in L_j_vec.iter().zip(evals) {
    //             for (final_coeff, L_j_coeff) in final_poly.iter_mut().zip(L_j.into_iter()) {
    //                 *final_coeff += L_j_coeff.mul(y_j);
    //             }
    //         }
    //         Self { coeffs: final_poly }
    //     }
    // }
}

impl<R: Ring> Polynomial for RingPolynomial<R> {
    type Coefficient = R;
    // TODO: all
    const MODULUS: Vec<Self::Coefficient> = vec![];

    fn rand(rng: &mut impl rand::RngCore, n: usize) -> Self {
        let coeffs = (0..n).map(|_| R::rand(rng)).collect::<Vec<_>>();
        Self::from_coefficients(coeffs)
    }

    // Remove leading zero coefficients
    fn normalize(&mut self) {
        while self.coeffs.len() > 1 && self.coeffs.last() == Some(&R::zero()) {
            self.coeffs.pop();
        }
    }

    fn zero() -> Self {
        Self {
            coeffs: vec![R::zero()],
        }
    }

    // The degree of the polynomial
    fn degree(&self) -> usize {
        assert!(self.coeffs.len() > 0);
        self.coeffs.len() - 1
    }

    fn coefficient(&self, i: usize) -> Self::Coefficient {
        assert!(self.degree() >= i, "Index out of bounds");
        self.coeffs[i].clone()
    }

    fn set_coefficient(&mut self, i: usize, value: Self::Coefficient) {
        assert!(self.degree() >= i, "Index out of bounds");
        self.coeffs[i] = value;
        self.normalize();
    }

    // This evaluates a polynomial (in coefficient form) at `x`.
    fn evaluate(&self, x: &Self::Coefficient) -> Self::Coefficient {
        let coeffs = self.coeffs.clone();
        let poly_size = self.coeffs.len();

        // p(x) = = a_0 + a_1 * X + ... + a_n * X^(n-1), revert it and fold sum it
        fn eval<R: Ring>(poly: &[R], point: &R) -> R {
            poly.iter()
                .rev()
                .fold(R::one(), |acc, coeff| acc * point + coeff)
        }

        let num_threads = current_num_threads();
        if poly_size * 2 < num_threads {
            eval(&coeffs, x)
        } else {
            let chunk_size = (poly_size + num_threads - 1) / num_threads;
            let mut parts = vec![R::one(); num_threads];
            scope(|scope| {
                for (chunk_idx, (out, c)) in parts
                    .chunks_mut(1)
                    .zip(coeffs.chunks(chunk_size))
                    .enumerate()
                {
                    scope.spawn(move |_| {
                        let start = chunk_idx * chunk_size;
                        out[0] = eval(c, x) * x.pow(start as u64);
                    });
                }
            });
            parts.iter().fold(R::one(), |acc, coeff| acc + coeff)
        }
    }

    fn from_coefficients(coeffs: Vec<Self::Coefficient>) -> Self {
        Self::new(coeffs)
    }

    fn coefficients(&self) -> Vec<Self::Coefficient> {
        self.coeffs.clone()
    }

    fn negate(&self) -> Self {
        let negated_coeffs = self.coeffs.iter().map(|c| -c.clone()).collect();
        Self::new(negated_coeffs)
    }

    fn derivative(&self) -> Self {
        if self.is_zero() || self.degree() == 0 {
            return Self::zero();
        }

        let mut derivative_coeffs = Vec::with_capacity(self.degree());
        for (i, coeff) in self.coeffs.iter().enumerate().skip(1) {
            let new_coeff = coeff.clone() * R::from(i as u64);
            derivative_coeffs.push(new_coeff);
        }

        Self::from_coefficients(derivative_coeffs)
    }

    fn is_zero(&self) -> bool {
        self == &Self::zero()
    }

    fn modulo(&self, other: &Self) -> Self {
        assert!(!other.is_zero(), "Cannot perform modulo by zero polynomial");

        if self.degree() < other.degree() {
            return self.clone();
        }

        let mut remainder = self.clone();
        let divisor_lead = other
            .coeffs
            .last()
            .expect("Leading coefficient must be invertible")
            .to_owned();
        let divisor_degree = other.degree();

        while remainder.degree() >= divisor_degree {
            let degree_diff = remainder.degree() - divisor_degree;
            let scalar = remainder.coeffs.last().unwrap().clone().div(divisor_lead);

            for (i, coeff) in other.coeffs.iter().enumerate() {
                let idx = degree_diff + i;
                remainder.coeffs[idx] -= scalar.clone() * coeff;
            }

            remainder.normalize();
        }

        remainder
    }
}

impl<R: Ring> Add for RingPolynomial<R> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let max_len = std::cmp::max(self.coeffs.len(), rhs.coeffs.len());
        let coeffs = (0..max_len)
            .into_iter()
            .map(|n| {
                if n >= self.coeffs.len() {
                    rhs.coeffs[n]
                } else if n >= rhs.coeffs.len() {
                    self.coeffs[n]
                } else {
                    // n < self.0.len() && n < rhs.0.len()
                    self.coeffs[n] + rhs.coeffs[n]
                }
            })
            .collect::<Vec<R>>();
        Self::new(coeffs)
    }
}

impl<R: Ring> AddAssign for RingPolynomial<R> {
    fn add_assign(&mut self, rhs: Self) {
        *self = self.clone() + rhs;
    }
}
impl<'a, R: Ring> Add<&'a Self> for RingPolynomial<R> {
    type Output = Self;

    fn add(self, rhs: &'a Self) -> Self::Output {
        let max_len = std::cmp::max(self.coeffs.len(), rhs.coeffs.len());
        let coeffs = (0..max_len)
            .into_iter()
            .map(|n| {
                if n >= self.coeffs.len() {
                    rhs.coeffs[n]
                } else if n >= rhs.coeffs.len() {
                    self.coeffs[n]
                } else {
                    // n < self.0.len() && n < rhs.0.len()
                    self.coeffs[n] + rhs.coeffs[n]
                }
            })
            .collect::<Vec<R>>();
        Self::new(coeffs)
    }
}

impl<R: Ring> std::ops::Mul for RingPolynomial<R> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        let mut coeffs: Vec<R> = vec![R::zero(); self.coeffs.len() + rhs.coeffs.len() - 1];
        for n in 0..self.coeffs.len() {
            for m in 0..rhs.coeffs.len() {
                coeffs[n + m] += self.coeffs[n] * rhs.coeffs[m];
            }
        }
        Self::new(coeffs)
    }
}
impl<'a, R: Ring> Mul<&'a Self> for RingPolynomial<R> {
    type Output = Self;
    fn mul(self, rhs: &Self) -> Self::Output {
        let mut coeffs: Vec<R> = vec![R::zero(); self.coeffs.len() + rhs.coeffs.len() - 1];
        for n in 0..self.coeffs.len() {
            for m in 0..rhs.coeffs.len() {
                coeffs[n + m] += self.coeffs[n] * rhs.coeffs[m];
            }
        }
        Self::new(coeffs)
    }
}

impl<R: Ring> MulAssign for RingPolynomial<R> {
    fn mul_assign(&mut self, rhs: Self) {
        *self = self.clone() * rhs;
    }
}

// impl<R: Ring> std::ops::Div for &Polynomial<F> {
//     type Output = (Polynomial<F>, Polynomial<F>);
//
//     fn div(self, rhs: &Polynomial<F>) -> Self::Output {
//         // init the (quotient, remainder)
//         let (mut q, mut r) = (Polynomial::zero(), self);
//
//         // r is not zero poly, and division.degree > divisor.degree.
//         while *r != Polynomial::zero() && r.degree() >= rhs.degree() {
//             let r_coeff = r.coeffs();
//             let rhs_coeff = rhs.coeffs();
//
//             let lead_r = r_coeff[r_coeff.len() - 1];
//             let lead_d = rhs_coeff[rhs_coeff.len() - 1];
//             let mut t = Polynomial::zero();
//             t.set(
//                 r_coeff.len() - rhs_coeff.len(),
//                 lead_r * lead_d.invert().unwrap(),
//             );
//             q += &t;
//             r -= &(&rhs * &t);
//         }
//         (q, r)
//     }
// }

impl<R: Ring> Sub for RingPolynomial<R> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        let max_len = self.coeffs.len().max(rhs.coeffs.len());
        let mut result = Vec::with_capacity(max_len);

        for i in 0..max_len {
            let a = self.coeffs.get(i).cloned().unwrap_or_else(R::zero);
            let b = rhs.coeffs.get(i).cloned().unwrap_or_else(R::zero);
            result.push(a - b);
        }

        Self::new(result)
    }
}

impl<'a, R: Ring> Sub<&'a Self> for RingPolynomial<R> {
    type Output = Self;
    fn sub(self, rhs: &Self) -> Self::Output {
        let max_len = self.coeffs.len().max(rhs.coeffs.len());
        let mut result = Vec::with_capacity(max_len);

        for i in 0..max_len {
            let a = self.coeffs.get(i).cloned().unwrap_or_else(R::zero);
            let b = rhs.coeffs.get(i).cloned().unwrap_or_else(R::zero);
            result.push(a - b);
        }

        Self::new(result)
    }
}
impl<R: Ring> SubAssign for RingPolynomial<R> {
    fn sub_assign(&mut self, rhs: Self) {
        *self = self.clone() - rhs;
    }
}
impl<'a, R: Ring> SubAssign<&'a Self> for RingPolynomial<R> {
    fn sub_assign(&mut self, rhs: &Self) {
        *self = self.clone() - rhs;
    }
}

// impl<R: Ring> ToString for RingPolynomial<R> {
//     fn to_string(&self) -> String {
//         let mut string = String::new();
//         if self.coeffs.is_empty() || (self.coeffs.len() == 1 && self.coeffs[0] == R::zero()) {
//             return "0".to_string();
//         }
//
//         let mut first = true;
//         for (i, coeff) in self.coeffs.iter().enumerate().rev() {
//             if *coeff != R::zero() {
//                 if !first {
//                     string.push_str(" + ");
//                 }
//                 first = false;
//                 let coeff_string = if *coeff == R::one() {
//                     "".to_string()
//                 } else {
//                     coeff.to_string()
//                 };
//                 match i {
//                     0 => string.push_str(&format!("{}", coeff_string)),
//                     1 => string.push_str(&format!("{}x", coeff_string)),
//                     _ => string.push_str(&format!("{}x^{}", coeff_string, i)),
//                 }
//             }
//         }
//         string
//     }
// }

impl<R: Ring + Display> Display for RingPolynomial<R> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        if self.coeffs.is_empty() || (self.coeffs.len() == 1 && self.coeffs[0] == R::zero()) {
            return write!(f, "0");
        }

        let mut first = true;
        for (i, coeff) in self.coeffs.iter().enumerate().rev() {
            if *coeff != R::zero() {
                if !first {
                    write!(f, " + ")?;
                }
                first = false;

                match i {
                    0 => write!(f, "{}", coeff)?,
                    1 => {
                        if *coeff == R::one() {
                            write!(f, "x")?
                        } else {
                            write!(f, "{}x", coeff)?
                        }
                    }
                    _ => {
                        if *coeff == R::one() {
                            write!(f, "x^{}", i)?
                        } else {
                            write!(f, "{}x^{}", coeff, i)?
                        }
                    }
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ring::fq::Fq;
    use std::ops::Neg;
    #[test]
    fn test_ring_polynomial_display() {
        let poly =
            RingPolynomial::from_coefficients(vec![Fq::new(3), Fq::new(0), Fq::new(2), Fq::new(1)]);
        // assert_eq!(poly.to_string(), "x^3 + 2x^2 + 3");
        println!("poly: {:?}", poly.to_string());
        let zero_poly = RingPolynomial::<Fq>::zero();
        // assert_eq!(zero_poly.to_string(), "0");
        println!("zero_poly: {:?}", zero_poly.to_string());

        let linear_poly = RingPolynomial::from_coefficients(vec![Fq::new(2), Fq::new(1)]);
        // assert_eq!(linear_poly.to_string(), "1x + 2");
        println!("linear_poly: {:?}", linear_poly.to_string());
    }

    #[test]
    fn test_poly_addition() {
        let p1 = RingPolynomial::<Fq>::from_coefficients(
            vec![3, 2, 1].into_iter().map(Fq::from).collect(),
        ); // x^2 + 2x + 3
        let p2 = RingPolynomial::<Fq>::from_coefficients(
            vec![6, 5, 4].into_iter().map(Fq::from).collect(),
        ); // x^2 + 5x + 6
        let result = p1 + p2;
        assert_eq!(result.coeffs, vec![Fq::new(9), Fq::new(7), Fq::new(5)]);
    }

    #[test]
    fn test_poly_subtraction() {
        let p1 = RingPolynomial::<Fq>::from_coefficients(
            vec![6, 5, 4].into_iter().map(Fq::from).collect(),
        ); // x^2 + 5x + 6
        let p2 = RingPolynomial::<Fq>::from_coefficients(
            vec![3, 2, 1].into_iter().map(Fq::from).collect(),
        ); // x^2 + 2x + 3
        let result = p1 - p2;
        assert_eq!(result.coeffs, vec![Fq::new(3); 3]);
    }

    #[test]
    fn test_poly_multiplication() {
        let p1 =
            RingPolynomial::<Fq>::from_coefficients(vec![1, 2].into_iter().map(Fq::from).collect()); // 2x + 1
        let p2 =
            RingPolynomial::<Fq>::from_coefficients(vec![3, 4].into_iter().map(Fq::from).collect()); // 4x + 3
        // (2x + 1)*(4x+ 3)
        let result = p1 * p2;
        assert_eq!(result.coeffs, vec![Fq::new(3), Fq::new(10), Fq::new(8)]);
    }

    #[test]
    fn test_poly_scalar_multiplication() {
        let p = RingPolynomial::<Fq>::from_coefficients(
            vec![1, 2, 3].into_iter().map(Fq::from).collect(),
        ); // x^2 + 2x + 3
        let q = Fq::new(2);
        let result = p.scalar_mul(&q);
        assert_eq!(result.coeffs, vec![Fq::new(2), Fq::new(4), Fq::new(6)]);
    }

    #[test]
    fn test_derivative() {
        // Test polynomial: 3x^3 + 2x^2 + x + 5
        let poly =
            RingPolynomial::from_coefficients(vec![Fq::new(5), Fq::new(1), Fq::new(2), Fq::new(3)]);

        // Expected derivative: 9x^2 + 4x + 1
        let expected_derivative =
            RingPolynomial::from_coefficients(vec![Fq::new(1), Fq::new(4), Fq::new(9)]);

        assert_eq!(poly.derivative(), expected_derivative);

        // Test constant polynomial
        let constant_poly = RingPolynomial::from_coefficients(vec![Fq::new(42)]);
        assert_eq!(constant_poly.derivative(), RingPolynomial::zero());

        // Test zero polynomial
        let zero_poly = RingPolynomial::<Fq>::zero();
        assert_eq!(zero_poly.derivative(), RingPolynomial::zero());
    }
    #[test]
    fn test_negate() {
        // Test polynomial: 3x^2 + 2x + 1
        let poly = RingPolynomial::from_coefficients(vec![Fq::new(1), Fq::new(2), Fq::new(3)]);

        // Expected negation: -3x^2 - 2x - 1
        let expected_negation =
            RingPolynomial::from_coefficients(vec![-Fq::new(1), -Fq::new(2), -Fq::new(3)]);

        assert_eq!(poly.negate(), expected_negation);

        // Test zero polynomial
        let zero_poly = RingPolynomial::<Fq>::zero();
        assert_eq!(zero_poly.negate(), zero_poly);

        // Test negation of negation
        assert_eq!(poly.negate().negate(), poly);
    }
    #[test]
    #[ignore]
    fn test_poly_modulo() {
        // TODO: debug
        // Dividend: x^3 + 2x^2 + 3x + 4
        let dividend =
            RingPolynomial::from_coefficients(vec![4, 3, 2, 1].into_iter().map(Fq::from).collect());

        // Divisor: x^2 + 1
        let divisor =
            RingPolynomial::from_coefficients(vec![2, 1].into_iter().map(Fq::from).collect());

        // Expected remainder: -2x + 3
        let expected_remainder =
            RingPolynomial::from_coefficients(vec![Fq::new(2).neg(), Fq::new(3)]);

        assert_eq!(dividend.modulo(&divisor), expected_remainder);

        // Test with zero remainder
        let dividend2 = dividend.clone() * &divisor;
        assert_eq!(dividend2.modulo(&divisor), RingPolynomial::zero());

        // Test with divisor of higher degree
        assert_eq!(dividend.modulo(&dividend2), dividend);
    }

    #[test]
    fn test_mul_poly() {
        // p = 1 - x
        let p = RingPolynomial {
            coeffs: vec![Fq::one(), Fq::one().neg()],
        };
        // q = 1 + x
        let q = RingPolynomial {
            coeffs: vec![Fq::one(), Fq::one()],
        };

        assert_eq!(p.clone().mul(&q).coeffs, vec![
            Fq::one(),
            Fq::zero(),
            Fq::one().neg()
        ]);

        // add
        assert_eq!(p.clone().add(&q).coeffs, vec![Fq::new(2), Fq::zero()]);

        // poly.mul(Fq)
        assert_eq!(p.scalar_mul(&Fq::new(5)).coeffs, vec![
            Fq::new(5),
            Fq::new(5).neg()
        ]);
    }

    // #[test]
    // fn lagrange_interpolate() {
    //     // aim: p = 1 + 2x + x^2
    //
    //     let domain = vec![
    //         Fq::new(1),
    //         Fq::new(2),
    //         Fq::new(3),
    //         Fq::new(4),
    //         Fq::new(5),
    //         Fq::new(6),
    //         Fq::new(7),
    //         Fq::new(8),
    //         Fq::new(9),
    //     ];
    //     let evals = vec![
    //         Fq::new(4),
    //         Fq::new(9),
    //         Fq::new(10),
    //         Fq::new(19),
    //         Fq::new(24),
    //         Fq::new(31),
    //         Fq::new(40),
    //         Fq::new(51),
    //         Fq::new(64),
    //     ];
    //
    //     let poly = RingPolynomial::lagrange_interpolate(domain.clone(), evals.clone());
    //
    //     for (x, y) in domain.iter().zip(evals) {
    //         assert_eq!(poly.evaluate(*x), y);
    //     }
    //     println!("pass");
    // }

    // #[test]
    // fn test_div() {
    //     // division: 2+3x+x^2 = (x+1)(x+2)
    //     let coeffs = vec![Fq::new(2), Fq::ONE, Fq::ONE];
    //     let division = Polynomial::from_coeffs(coeffs);
    //
    //     // dividor: 2+x
    //     let coeffs = vec![Fq::new(2), Fq::ONE];
    //     let dividor = Polynomial::from_coeffs(coeffs);
    //
    //     // target:
    //     //      quotient poly: 1+x
    //     //      remainder poly: 0
    //     let coeffs = vec![Fq::new(2), Fq::ONE];
    //     let target_qoutient = Polynomial::from_coeffs(coeffs);
    //     let target_remainder = Polynomial::zero();
    //
    //     // division / dividor = quotient + remainder
    //     let (actual_qoutient, actual_remainder) = division.div(dividor);
    //
    //     assert_eq!(actual_qoutient, target_qoutient);
    //     assert_eq!(actual_remainder, target_remainder);
    // }
}
