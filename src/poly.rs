use crate::ring::Ring;
use rand::RngCore;
use std::fmt::{Debug, Display};
use std::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign};

pub mod ring_poly;

pub trait Polynomial:
    Add<Output = Self>
    + AddAssign
    + for<'a> Add<&'a Self, Output = Self>
    + Mul<Output = Self>
    + for<'a> Mul<&'a Self, Output = Self>
    + MulAssign
    + Sub<Output = Self>
    + for<'a> Sub<&'a Self, Output = Self>
    + SubAssign
    + Sized
    + Clone
    + Debug
    + PartialEq
    + Eq
    + Display
{
    type Coefficient: Ring;

    /// The modulus of the polynomial ring
    // const MODULE: Self;

    const MODULUS: Vec<Self::Coefficient>;

    fn rand(rng: &mut impl rand::RngCore, n: usize) -> Self;
    /// Remove leading zero coefficients
    fn normalize(&mut self);

    /// Create a zero polynomial
    fn zero() -> Self;

    /// Create a polynomial representing 1
    // fn one() -> Self;

    /// Get the degree of the polynomial
    fn degree(&self) -> usize;

    /// Get the coefficient of the x^i term
    fn coefficient(&self, i: usize) -> Self::Coefficient;

    /// Set the coefficient of the x^i term
    fn set_coefficient(&mut self, i: usize, value: Self::Coefficient);

    /// Evaluate the polynomial at a given point
    fn evaluate(&self, x: &Self::Coefficient) -> Self::Coefficient;

    /// Create a polynomial from a list of coefficients
    fn from_coefficients(coeffs: Vec<Self::Coefficient>) -> Self;

    /// Get all coefficients of the polynomial
    fn coefficients(&self) -> Vec<Self::Coefficient>;

    /// Negate the polynomial
    fn negate(&self) -> Self;

    /// Compute the derivative of the polynomial
    fn derivative(&self) -> Self;

    /// Check if the polynomial is zero
    fn is_zero(&self) -> bool;

    /// Reduce the polynomial modulo another polynomial
    fn modulo(&self, other: &Self) -> Self;

    // /// Compute the greatest common divisor of two polynomials
    // fn gcd(a: &Self, b: &Self) -> Self;
    //
    // /// Compute the extended Euclidean algorithm for polynomials
    // fn extended_gcd(a: &Self, b: &Self) -> (Self, Self, Self);
}

// a randomizer polynomial vector
// These polynomial vectors are freshly generated for every encryption.
pub fn random_poly_vector<P: Polynomial>(
    rng: &mut impl RngCore,
    dimension: usize,
    degree: usize,
) -> Vec<P> {
    let mut s = vec![];
    for _ in 0..dimension {
        s.push(P::rand(rng, degree));
    }
    s
}
