use crate::ring::Ring;
use rand::RngCore;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::ops::{Add, Div, Mul, Neg, Sub};
use std::random::{Random, RandomSource};

//
#[derive(Debug, Copy, Clone, PartialEq, Ord, PartialOrd, Eq, Serialize, Deserialize)]
pub struct Fq {
    value: u64,
}

impl Fq {
    // To ensure the value within the [0, MODULUS - 1].
    pub fn new(value: u64) -> Self {
        Self {
            value: value % Self::MODULUS,
        }
    }

    // Helper function to calculate modular multiplicative inverse
    // The Extended Euclidean Algorithm to find the modular multiplicative inverse of a modulo m.
    // Here's a breakdown:
    fn inverse(a: i64, m: i64) -> Option<u64> {
        let mut t = 0i64;
        let mut newt = 1i64;
        let mut r = m;
        let mut newr = a;

        while newr != 0 {
            let quotient = r / newr;
            (t, newt) = (newt, t - quotient * newt);
            (r, newr) = (newr, r - quotient * newr);
        }

        if r > 1 {
            None // a is not invertible
        } else {
            Some(((t + m as i64) % m as i64) as u64)
        }
    }
}

impl Random for Fq {
    fn random(source: &mut (impl RandomSource + ?Sized)) -> Self {
        Self::new(u64::random(source))
    }
}

impl Ring for Fq {
    const MODULUS: u64 = 17;
    fn rand(rng: &mut impl RngCore) -> Self {
        Self::new(rng.next_u64())
    }
    fn zero() -> Self {
        Self::new(0)
    }
    fn one() -> Self {
        Self::new(1)
    }
}

impl Add for Fq {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::new(self.value + rhs.value)
    }
}

impl Mul for Fq {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self::new(self.value * rhs.value)
    }
}

impl Neg for Fq {
    type Output = Self;

    fn neg(self) -> Self::Output {
        let value = if self.value == 0 {
            0
        } else {
            Self::MODULUS - self.value
        };
        Self { value }
    }
}

impl Sub for Fq {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self::new(self.value + Self::MODULUS - rhs.value)
    }
}

// Division
impl Div for Fq {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        if rhs.value == 0 {
            panic!("Division by zero");
        }

        let inv = Self::inverse(rhs.value as i64, Self::MODULUS as i64)
            .expect("Divisor has no modular inverse");

        Self::new(self.value * inv)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::random::DefaultRandomSource;

    #[test]
    fn test_new_fq() {
        assert_eq!(Fq::zero(), Fq::new(Fq::MODULUS));
        assert_eq!(Fq::one(), Fq::new(1));
    }

    #[test]
    fn test_random_fq() {
        let mut rng = DefaultRandomSource::default();
        let lhs = Fq::random(&mut rng);
        let rhs = Fq::random(&mut rng);
        assert_ne!(lhs, rhs);
    }

    #[test]
    fn test_rand_fq() {
        let mut rng = rand::thread_rng();
        let lhs = Fq::rand(&mut rng);
        let rhs = Fq::rand(&mut rng);
        assert_ne!(lhs, rhs);
    }

    #[test]
    fn test_fq_addition() {
        assert_eq!(Fq::new(8), Fq::new(5) + Fq::new(3)); // 13 mod 17 = 13
        assert_eq!(Fq::new(15), Fq::new(7) + Fq::new(8)); // 15 mod 17 = 15
        assert_eq!(Fq::new(1), Fq::new(9) + Fq::new(9)); // 18 mod 17 = 1
    }

    #[test]
    fn test_fq_subtraction() {
        let a = Fq::new(5);
        let b = Fq::new(3);
        assert_eq!(Fq::new(2), a - b);
        assert_eq!(Fq::new(16), Fq::new(3) - Fq::new(4)); // -1 mod 17 = 16
        assert_eq!(Fq::new(0), Fq::new(7) - Fq::new(7));
    }

    #[test]
    fn test_fq_multiplication() {
        let a = Fq::new(5);
        let b = Fq::new(3);
        assert_eq!(Fq::new(15), a * b);
        assert_eq!(Fq::new(1), Fq::new(6) * Fq::new(3)); // 18 mod 17 = 1
        assert_eq!(Fq::new(0), Fq::new(17) * Fq::new(5)); // 85 mod 17 = 0
    }

    #[test]
    fn test_fq_division() {
        let a = Fq::new(8);
        let b = Fq::new(2);
        assert_eq!(Fq::new(4), a / b);
        assert_eq!(Fq::new(1), Fq::new(5) / Fq::new(5));
        assert_eq!(Fq::new(9), Fq::new(1) / Fq::new(2)); // 1 * 9 = 18 ≡ 1 (mod 17)
    }

    #[test]
    #[should_panic(expected = "Division by zero")]
    fn test_fq_division_by_zero() {
        let a = Fq::new(5);
        let b = Fq::new(0);
        let _ = a / b;
    }

    #[test]
    fn test_fq_negation() {
        assert_eq!(Fq::new(0), -Fq::new(0));
        assert_eq!(Fq::new(12), -Fq::new(5)); // -5 mod 17 = 12
        assert_eq!(Fq::new(1), -Fq::new(16));
    }

    #[test]
    fn test_fq_additive_inverse() {
        for i in 0..17 {
            let a = Fq::new(i);
            assert_eq!(Fq::new(0), a + (-a));
        }
    }

    #[test]
    fn test_fq_multiplicative_inverse() {
        for i in 1..17 {
            let a = Fq::new(i);
            assert_eq!(Fq::new(1), a * (Fq::new(1) / a));
        }
    }

    #[test]
    fn test_fq_associativity() {
        let a = Fq::new(5);
        let b = Fq::new(7);
        let c = Fq::new(11);
        assert_eq!((a + b) + c, a + (b + c));
        assert_eq!((a * b) * c, a * (b * c));
    }

    #[test]
    fn test_fq_commutativity() {
        let a = Fq::new(5);
        let b = Fq::new(7);
        assert_eq!(a + b, b + a);
        assert_eq!(a * b, b * a);
    }

    #[test]
    fn test_fq_distributivity() {
        let a = Fq::new(5);
        let b = Fq::new(7);
        let c = Fq::new(11);
        assert_eq!(a * (b + c), (a * b) + (a * c));
    }

    #[test]
    fn test_fq_serde_and_deserde() {
        let mut rng = rand::thread_rng();
        let lhs = Fq::rand(&mut rng);
        let serde = serde_json::to_string(&lhs).unwrap();
        let rhs = serde_json::from_str::<Fq>(&serde).unwrap();
        assert_eq!(lhs, rhs);
    }
}