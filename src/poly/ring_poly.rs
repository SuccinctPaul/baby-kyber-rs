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

    fn scalar_mul(&self, rhs: &R) -> Self {
        let coeffs = if rhs == &R::one() {
            vec![R::one()]
        } else {
            self.coeffs.iter().map(|c| c.mul(rhs)).collect::<Vec<R>>()
        };
        Self::new(coeffs)
    }

    #[ignore]
    pub fn div_rem(&self, other: &Self) -> (Option<Self>, Self) {
        if other.is_zero() {
            panic!("Division by zero polynomial");
        }

        if self.degree() < other.degree() {
            return (None, self.clone());
        }

        let mut quotient = Self::zero();
        let mut remainder = self.clone();

        while remainder.degree() >= other.degree() {
            let lead_r = remainder.coeffs.last().unwrap().clone();
            let lead_d = other.coeffs.last().unwrap().clone();
            let degree_diff = remainder.degree() - other.degree();

            let mut term = Self::zero();
            term.coeffs.resize(degree_diff + 1, R::zero());
            term.coeffs[degree_diff] = lead_r.div(lead_d);

            quotient += term.clone();
            remainder -= (other.clone() * term);
        }

        quotient.normalize();
        remainder.normalize();

        if quotient.is_zero() {
            (None, remainder)
        } else {
            (Some(quotient), remainder)
        }
    }
}

impl<R: Ring> Polynomial for RingPolynomial<R> {
    type Coefficient = R;

    fn rand(rng: &mut impl rand::RngCore, degree: usize) -> Self {
        let coeffs = (0..degree + 1).map(|_| R::rand(rng)).collect::<Vec<_>>();
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

impl<R: Ring> Display for RingPolynomial<R> {
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

        //     s_tranpose_dot_u: "2x^3 + 7x^2 + 8x + 3"
        // ciphter_text.v:
        // "11x^3 + 7x^2 + x + 13"
        let p3 = RingPolynomial::<Fq>::from_coefficients(
            vec![13, 1, 7, 11].into_iter().map(Fq::from).collect(),
        );

        // 2x^3 + 7x^2 + 8x + 3
        let p4 = RingPolynomial::<Fq>::from_coefficients(
            vec![3, 8, 7, 2].into_iter().map(Fq::from).collect(),
        );

        let result = p3 - p4;
        println!("result: {:?}", result.to_string());
        assert_eq!(result.coeffs, vec![
            Fq::new(10),
            Fq::new(7).neg(),
            Fq::zero(),
            Fq::new(9)
        ]);
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
    fn test_div_rem() {
        // Define polynomials
        let p1 = RingPolynomial::from_coefficients(vec![Fq::new(1), Fq::new(2), Fq::new(1)]); // x^2 + 2x + 1
        let p2 = RingPolynomial::from_coefficients(vec![Fq::new(1), Fq::new(1)]); // x + 1

        // Perform division
        let (quotient_opt, remainder) = p1.div_rem(&p2);

        // Check quotient
        assert!(quotient_opt.is_none());
        let quotient = quotient_opt.unwrap();
        assert_eq!(quotient.coefficients(), vec![Fq::new(1), Fq::new(1)]); // x + 1

        // Check remainder
        assert_eq!(remainder.coefficients(), vec![Fq::new(0)]); // 0

        // Test division by higher degree polynomial
        let p3 = RingPolynomial::from_coefficients(vec![Fq::new(1), Fq::new(1), Fq::new(1)]); // x^2 + x + 1
        let (quotient_opt, remainder) = p1.div_rem(&p3);

        // Check quotient is None
        assert!(quotient_opt.is_none());

        // Check remainder is the same as the dividend
        assert_eq!(remainder.coefficients(), p1.coefficients());

        // Test division by zero polynomial
        let zero_poly = RingPolynomial::zero();
        let result = std::panic::catch_unwind(|| p1.div_rem(&zero_poly));
        assert!(result.is_err());
    }
    #[test]
    #[ignore]
    fn test_poly_modulo() {
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
        assert_eq!(dividend2.div_rem(&divisor).1, RingPolynomial::zero());

        // Test with divisor of higher degree
        assert_eq!(dividend.div_rem(&dividend2).0.unwrap(), dividend);
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
