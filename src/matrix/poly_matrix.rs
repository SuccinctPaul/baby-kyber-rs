use crate::matrix::Matrix;
use crate::matrix::vector_arithmatic::VectorArithmatic;
use crate::poly::Polynomial;
use crate::ring::Ring;
use std::ops::{Add, AddAssign, Div, Mul, Sub};

/// This defined `matrix` (rows * cols) （m × n）
#[derive(Debug, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub struct PolyMatrix<P: Polynomial> {
    pub rows: usize,
    pub cols: usize,
    pub values: Vec<Vec<P>>,
}

impl<P: Polynomial> PolyMatrix<P> {
    pub fn rand(
        rng: &mut impl rand::RngCore,
        rows: usize,
        cols: usize,
        poly_degree: usize,
    ) -> Self {
        let values = (0..rows)
            .map(|_| {
                (0..cols)
                    .map(|_| P::rand(rng, poly_degree))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();

        Self { cols, rows, values }
    }

    pub fn get_columns(&self, column_index: usize) -> Vec<P> {
        assert!(self.cols > column_index, "Column index out of bounds");

        self.values
            .iter()
            .map(|v| v.get(column_index).unwrap().clone())
            .collect::<Vec<_>>()
    }

    pub fn vec_dot_mul(a: &Vec<P>, b: &Vec<P>) -> P {
        assert_eq!(a.len(), b.len(), "Vectors must have the same length");

        let mut poly = a
            .iter()
            .zip(b.iter())
            .map(|(ai, bi)| ai.clone() * bi.clone())
            .fold(P::zero(), |acc, x| acc + x);
        poly.normalize();
        poly
    }
    // Add between two vectors.
    pub fn vec_add(a: &Vec<P>, b: &Vec<P>) -> Vec<P> {
        assert_eq!(a.len(), b.len(), "Vectors must have the same length");

        a.iter()
            .zip(b.iter())
            .map(|(ai, bi)| ai.clone() + bi.clone())
            .collect::<Vec<_>>()
    }
    // Add between two vectors.
    pub fn vec_sub(a: &Vec<P>, b: &Vec<P>) -> Vec<P> {
        assert_eq!(a.len(), b.len(), "Vectors must have the same length");

        a.iter()
            .zip(b.iter())
            .map(|(ai, bi)| ai.clone() - bi.clone())
            .collect::<Vec<_>>()
    }

    // Mul between two vectors.
    // It's a special case of matrix mul.
    pub fn vec_mul(a: &Vec<P>, b: &Vec<P>) -> Vec<P> {
        assert_eq!(a.len(), b.len(), "Vectors must have the same length");

        a.iter()
            .zip(b.iter())
            .map(|(ai, bi)| ai.clone() * bi.clone())
            .collect::<Vec<_>>()
    }

    /// https://en.wikipedia.org/wiki/Dot_product
    /// Suppose A(m * n), x(n) => A * x = y(m)
    pub fn mul_vector(&self, vector: &Vec<P>) -> Vec<P> {
        assert_eq!(
            self.cols,
            vector.len(),
            "Matrix columns must match vector length"
        );

        self.values
            .iter()
            .map(|row| Self::vec_dot_mul(row, vector))
            .collect()
    }

    /// https://en.wikipedia.org/wiki/Dot_product
    /// Suppose A(m * n), B(n, p) => A * B = C(m * p)
    pub fn mul_matrix(&self, m_b: &Self) -> Self {
        assert!(self.cols > 0 && m_b.rows > 0, "Matrices cannot be empty");
        assert_eq!(
            self.cols, m_b.rows,
            "Matrix dimensions must be compatible for multiplication"
        );

        let m = self.rows;
        let p = m_b.cols;

        let m_b_columns: Vec<Vec<P>> = (0..p).map(|j| m_b.get_columns(j)).collect();

        let matrix: Vec<Vec<P>> = (0..m)
            .map(|i| {
                let row_i = &self.values[i];
                (0..p)
                    .map(|j| Self::vec_dot_mul(row_i, &m_b_columns[j]))
                    .collect()
            })
            .collect();

        Self {
            rows: m,
            cols: p,
            values: matrix,
        }
    }
}

impl<P: Polynomial> Matrix<P> for PolyMatrix<P> {
    fn new(rows: usize, cols: usize) -> Self {
        Self {
            rows,
            cols,
            values: vec![vec![P::zero(); cols]; rows],
        }
    }

    fn rows(&self) -> usize {
        self.rows
    }

    fn cols(&self) -> usize {
        self.cols
    }

    fn get(&self, row: usize, col: usize) -> Option<&P> {
        assert!(row < self.rows, "row out of bounds");
        assert!(col < self.cols, "row out of bounds");
        self.values.get(row).and_then(|row| row.get(col))
    }

    fn set(&mut self, row: usize, col: usize, value: P) -> Result<(), &'static str> {
        assert!(row < self.rows, "row out of bounds");
        assert!(col < self.cols, "row out of bounds");
        self.values[row][col] = value;
        Ok(())
    }

    fn transpose(&self) -> Self {
        let mut transposed = Self::new(self.cols, self.rows);
        for i in 0..self.rows {
            for j in 0..self.cols {
                transposed
                    .set(j, i, self.get(i, j).unwrap().clone())
                    .unwrap();
            }
        }
        transposed
    }

    fn identity(size: usize) -> Self {
        todo!()
    }

    fn scalar_mul(&self, scalar: P) -> Self {
        todo!()
    }

    fn determinant(&self) -> Option<P> {
        todo!()
    }

    fn inverse(&self) -> Option<Self> {
        todo!()
    }
}
impl<P: Polynomial> Add for PolyMatrix<P> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        assert_eq!(
            self.cols, rhs.cols,
            "Matrix rows and cols are not the same length"
        );
        assert_eq!(
            self.rows, rhs.rows,
            "Matrix rows and cols are not the same length"
        );

        let values = self
            .values
            .into_iter()
            .zip(rhs.values.into_iter())
            .map(|(l_row, r_row)| {
                l_row
                    .into_iter()
                    .zip(r_row.into_iter())
                    .map(|(l, r)| l + r)
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        Self {
            rows: self.rows,
            cols: self.cols,
            values,
        }
    }
}

impl<P: Polynomial> Sub for PolyMatrix<P> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        assert_eq!(
            self.cols, rhs.cols,
            "Matrix rows and cols are not the same length"
        );
        assert_eq!(
            self.rows, rhs.rows,
            "Matrix rows and cols are not the same length"
        );

        let values = self
            .values
            .into_iter()
            .zip(rhs.values.into_iter())
            .map(|(l_row, r_row)| {
                l_row
                    .into_iter()
                    .zip(r_row.into_iter())
                    .map(|(l, r)| l - r)
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        Self {
            rows: self.rows,
            cols: self.cols,
            values,
        }
    }
}

impl<P: Polynomial> Mul for PolyMatrix<P> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        self.mul_matrix(&rhs)
    }
}

#[cfg(test)]
mod test {
    use crate::matrix::Matrix;
    use crate::matrix::poly_matrix::PolyMatrix;
    use crate::poly::Polynomial;
    use crate::poly::ring_poly::RingPolynomial;
    use crate::ring::Ring;
    use crate::ring::fq::Fq;
    use std::ops::Neg;

    #[test]
    fn test_ring_matrix_new() {
        let matrix = PolyMatrix::<RingPolynomial<Fq>>::new(3, 4);
        println!("{:?}", matrix);
    }

    #[test]
    pub fn test_matrix_vec_dot_mul() {
        let poly_1 = RingPolynomial::from_coefficients(vec![Fq::one(), Fq::one(), Fq::one()]);
        let vec1 = vec![RingPolynomial::zero(), poly_1.clone()];
        let vec2 = vec![poly_1.clone(), RingPolynomial::zero()];
        assert_eq!(
            PolyMatrix::<RingPolynomial<Fq>>::vec_dot_mul(&vec1, &vec2),
            RingPolynomial::zero()
        );

        // [x+1, x-1]
        let vec3 = vec![
            RingPolynomial::from_coefficients(vec![Fq::one(), Fq::one()]),
            RingPolynomial::from_coefficients(vec![Fq::one().neg(), Fq::one()]),
        ];
        // [3x+1, x+3]
        let mut vec4 = vec![
            RingPolynomial::from_coefficients(vec![Fq::new(1), Fq::new(3)]),
            RingPolynomial::from_coefficients(vec![Fq::new(3), Fq::new(1)]),
        ];
        assert_eq!(
            PolyMatrix::<RingPolynomial<Fq>>::vec_dot_mul(&vec3, &vec4),
            // -2 + 6x + 4x^2
            RingPolynomial::from_coefficients(vec![Fq::new(2).neg(), Fq::new(6), Fq::new(4)])
        );
    }

    #[test]
    pub fn test_matrix_transpose() {
        let m = 2;
        let matrix = PolyMatrix::<RingPolynomial<Fq>> {
            rows: m,
            cols: m,
            values: vec![
                vec![
                    RingPolynomial::zero(),
                    RingPolynomial::from_coefficients(vec![Fq::new(1), Fq::new(3)]),
                ],
                vec![
                    RingPolynomial::from_coefficients(vec![Fq::one(), Fq::new(2)]),
                    RingPolynomial::from_coefficients(vec![Fq::new(3), Fq::new(4)]),
                ],
            ],
        };
        let transposed: PolyMatrix<RingPolynomial<Fq>> = matrix.transpose();
        let expect = PolyMatrix::<RingPolynomial<Fq>> {
            rows: m,
            cols: m,
            values: vec![
                vec![
                    RingPolynomial::zero(),
                    RingPolynomial::from_coefficients(vec![Fq::one(), Fq::new(2)]),
                ],
                vec![
                    RingPolynomial::from_coefficients(vec![Fq::new(1), Fq::new(3)]),
                    RingPolynomial::from_coefficients(vec![Fq::new(3), Fq::new(4)]),
                ],
            ],
        };

        assert_eq!(transposed, expect);

        let recovered: PolyMatrix<RingPolynomial<Fq>> = transposed.transpose();
        assert_eq!(recovered, matrix);
    }

    #[test]
    fn test_matrix_add_and_sub() {
        let cols = 3;
        let rows = 4;
        let degree = 4;
        let rng = &mut rand::thread_rng();
        let lhs = PolyMatrix::<RingPolynomial<Fq>>::rand(rng, rows, cols, degree);
        let rhs = PolyMatrix::<RingPolynomial<Fq>>::rand(rng, rows, cols, degree);

        let sum = lhs.clone() + rhs.clone();

        let expect_lhs = sum.clone() - rhs.clone();
        assert_eq!(expect_lhs, lhs);
        let expect_rhs = sum.clone() - lhs.clone();
        assert_eq!(expect_rhs, rhs);
        let expect_sum = expect_lhs.clone() + expect_rhs.clone();
        assert_eq!(expect_sum, sum);
        let expect_zero = sum - rhs - lhs;
        assert_eq!(expect_zero, Matrix::new(rows, cols));
    }

    #[test]
    fn test_mul_vector() {
        let m = 2;
        let rng = &mut rand::thread_rng();
        let lhs = PolyMatrix {
            rows: m,
            cols: m,
            values: vec![
                // [x+1, x-1]
                vec![
                    RingPolynomial::from_coefficients(vec![Fq::one(), Fq::one()]),
                    RingPolynomial::from_coefficients(vec![Fq::one().neg(), Fq::one()]),
                ],
                // [x+1, x-2]
                vec![
                    RingPolynomial::from_coefficients(vec![Fq::one(), Fq::one()]),
                    RingPolynomial::from_coefficients(vec![Fq::new(2).neg(), Fq::one()]),
                ],
            ],
        };
        // [3x+1, x+3]
        let rhs = vec![
            RingPolynomial::from_coefficients(vec![Fq::new(1), Fq::new(3)]),
            RingPolynomial::from_coefficients(vec![Fq::new(3), Fq::new(1)]),
        ];

        let expect = vec![
            // -2 + 6x + 4x^2
            RingPolynomial::from_coefficients(vec![Fq::new(2).neg(), Fq::new(6), Fq::new(4)]),
            // (4x^2 + 5x - 5)
            RingPolynomial::from_coefficients(vec![Fq::new(5).neg(), Fq::new(5), Fq::new(4)]),
        ];

        let res = lhs.mul_vector(&rhs);
        assert_eq!(expect, res);
    }

    #[test]
    fn test_mul_matrix() {
        let m: usize = 2;
        let lhs = PolyMatrix {
            rows: m,
            cols: m,
            values: vec![
                // [x+1, x-1]
                vec![
                    RingPolynomial::from_coefficients(vec![Fq::one(), Fq::one()]),
                    RingPolynomial::from_coefficients(vec![Fq::one().neg(), Fq::one()]),
                ],
                // [x+1, x-2]
                vec![
                    RingPolynomial::from_coefficients(vec![Fq::one(), Fq::one()]),
                    RingPolynomial::from_coefficients(vec![Fq::new(2).neg(), Fq::one()]),
                ],
            ],
        };
        let rhs = PolyMatrix {
            rows: m,
            cols: m,
            values: vec![
                // [3x+1, x+3]
                vec![
                    RingPolynomial::from_coefficients(vec![Fq::new(1), Fq::new(3)]),
                    RingPolynomial::from_coefficients(vec![Fq::new(3), Fq::new(1)]),
                ],
                // [3x+1, x+1]
                vec![
                    RingPolynomial::from_coefficients(vec![Fq::new(1), Fq::new(3)]),
                    RingPolynomial::from_coefficients(vec![Fq::new(1), Fq::new(1)]),
                ],
            ],
        };

        let expect = PolyMatrix {
            rows: m,
            cols: m,
            values: vec![
                // [2x + 6x^2, 2 + 4x + 2x^2]
                vec![
                    RingPolynomial::from_coefficients(vec![Fq::new(0), Fq::new(2), Fq::new(6)]),
                    RingPolynomial::from_coefficients(vec![Fq::new(2), Fq::new(4), Fq::new(2)]),
                ],
                // [-1 - x + 6x^2, 1 + 3x + 2x^2]
                vec![
                    RingPolynomial::from_coefficients(vec![
                        Fq::new(1).neg(),
                        Fq::new(1).neg(),
                        Fq::new(6),
                    ]),
                    RingPolynomial::from_coefficients(vec![Fq::new(1), Fq::new(3), Fq::new(2)]),
                ],
            ],
        };
        let res = lhs.mul_matrix(&rhs);
        assert_eq!(expect, res);
    }

    #[test]
    fn test_matrix_scalar_mul_and_matrix_mul() {
        let m = 2;
        let mut rng = rand::thread_rng();

        let lhs = PolyMatrix {
            rows: m,
            cols: m,
            values: vec![
                // [x+1, x-1]
                vec![
                    RingPolynomial::from_coefficients(vec![Fq::one(), Fq::one()]),
                    RingPolynomial::from_coefficients(vec![Fq::one().neg(), Fq::one()]),
                ],
                // [x+1, x-2]
                vec![
                    RingPolynomial::from_coefficients(vec![Fq::one(), Fq::one()]),
                    RingPolynomial::from_coefficients(vec![Fq::new(2).neg(), Fq::one()]),
                ],
            ],
        };
        let rhs = PolyMatrix {
            rows: m,
            cols: m,
            values: vec![
                // [3x+1, x+3]
                vec![
                    RingPolynomial::from_coefficients(vec![Fq::new(1), Fq::new(3)]),
                    RingPolynomial::from_coefficients(vec![Fq::new(3), Fq::new(1)]),
                ],
                // [3x+1, x+1]
                vec![
                    RingPolynomial::from_coefficients(vec![Fq::new(1), Fq::new(3)]),
                    RingPolynomial::from_coefficients(vec![Fq::new(1), Fq::new(1)]),
                ],
            ],
        };
        // [3x+1, x+3]
        let x = vec![
            RingPolynomial::from_coefficients(vec![Fq::new(1), Fq::new(3)]),
            RingPolynomial::from_coefficients(vec![Fq::new(3), Fq::new(1)]),
        ];

        // A*B*x
        let res1 = PolyMatrix::<RingPolynomial<Fq>>::mul_matrix(&lhs, &rhs).mul_vector(&x);
        // A*(B*x)
        let res2 = lhs.mul_vector(&rhs.mul_vector(&x));
        assert_eq!(res1, res2);
    }
}
