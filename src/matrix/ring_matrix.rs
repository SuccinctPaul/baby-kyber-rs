use crate::matrix::Matrix;
use crate::ring::Ring;
use std::ops::{Add, AddAssign, Div, Mul, Sub};

/// This define `matrix` (rows * cols) （m × n）
#[derive(Debug, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub struct RingMatrix<R: Ring> {
    rows: usize,
    cols: usize,
    values: Vec<Vec<R>>,
}

impl<R: Ring> RingMatrix<R> {
    pub fn get_columns(&self, column_index: usize) -> Vec<R> {
        assert!(self.cols > column_index, "Column index out of bounds");

        self.values
            .iter()
            .map(|v| v.get(column_index).unwrap().clone())
            .collect::<Vec<_>>()
    }

    // Multiple between two vectors.
    pub fn vec_mul(a: &Vec<R>, b: &Vec<R>) -> R {
        assert_eq!(a.len(), b.len(), "Vectors must have the same length");

        a.iter()
            .zip(b.iter())
            .map(|(ai, bi)| ai.clone() * bi.clone())
            .fold(R::zero(), |acc, x| acc + x)
    }

    /// https://en.wikipedia.org/wiki/Dot_product
    /// Suppose A(m * n), x(n) => A * x = y(n)
    pub fn mul_vector(&self, vector: &Vec<R>) -> Vec<R> {
        assert_eq!(
            self.cols,
            vector.len(),
            "Matrix columns must match vector length"
        );
        let n = self.cols;

        self.values
            .iter()
            .map(|row| Self::vec_mul(row, vector))
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

        let m_b_columns: Vec<Vec<R>> = (0..p).map(|j| m_b.get_columns(j)).collect();

        let matrix: Vec<Vec<R>> = (0..m)
            .map(|i| {
                let row_i = &self.values[i];
                (0..p)
                    .map(|j| Self::vec_mul(row_i, &m_b_columns[j]))
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

impl<R: Ring> Matrix<R> for RingMatrix<R> {
    fn new(rows: usize, cols: usize) -> Self {
        Self {
            rows,
            cols,
            values: vec![vec![R::zero(); cols]; rows],
        }
    }
    fn rand(rng: &mut impl rand::RngCore, rows: usize, cols: usize) -> Self {
        let values = (0..rows)
            .map(|_| (0..cols).map(|_| R::rand(rng)).collect::<Vec<_>>())
            .collect::<Vec<_>>();

        Self { cols, rows, values }
    }

    fn rows(&self) -> usize {
        self.rows
    }

    fn cols(&self) -> usize {
        self.cols
    }

    fn get(&self, row: usize, col: usize) -> Option<&R> {
        assert!(row < self.rows, "row out of bounds");
        assert!(col < self.cols, "row out of bounds");
        self.values.get(row).and_then(|row| row.get(col))
    }

    fn set(&mut self, row: usize, col: usize, value: R) -> Result<(), &'static str> {
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
        let identity = R::one();
        let mut matrix = Self::new(size, size);
        for i in 0..size {
            matrix.set(i, i, identity.clone()).unwrap();
        }
        matrix
    }

    fn scalar_mul(&self, scalar: R) -> Self {
        let mut result = Self::new(self.rows, self.cols);
        for i in 0..self.rows {
            for j in 0..self.cols {
                let value = self.get(i, j).unwrap().clone() * scalar;
                result.set(i, j, value).unwrap();
            }
        }
        result
    }

    fn determinant(&self) -> Option<R> {
        todo!()
    }

    fn inverse(&self) -> Option<Self> {
        todo!()
    }
}
impl<R: Ring> Add for RingMatrix<R> {
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

impl<R: Ring> Sub for RingMatrix<R> {
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

impl<R: Ring> Mul for RingMatrix<R> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        self.mul_matrix(&rhs)
    }
}

#[cfg(test)]
mod test {
    use crate::matrix::Matrix;
    use crate::matrix::ring_matrix::RingMatrix;
    use crate::ring::Ring;
    use crate::ring::fq::Fq;

    #[test]
    fn test_ring_matrix_new() {
        let matrix = RingMatrix::<Fq>::new(3, 4);
        println!("{:?}", matrix);
    }

    #[test]
    pub fn test_matrix_mul_vector() {
        let vec1 = vec![Fq::one(), Fq::zero(), Fq::new(3)];
        let vec2 = vec![Fq::zero(), Fq::one(), Fq::new(2)];
        assert_eq!(RingMatrix::<Fq>::vec_mul(&vec1, &vec2), Fq::new(6));

        let vec3 = vec![Fq::one(), Fq::one()];
        let vec4 = vec![Fq::zero(), Fq::zero()];
        assert_eq!(RingMatrix::<Fq>::vec_mul(&vec3, &vec4), Fq::zero());
    }
    #[test]
    pub fn test_matrix_identity() {
        let m = 2;
        let vector = vec![Fq::one(), Fq::zero()];

        let matrix = RingMatrix::<Fq> {
            rows: m,
            cols: m,
            values: vec![vec![Fq::one(), Fq::zero()], vec![Fq::zero(), Fq::one()]],
        };
        assert_eq!(matrix, RingMatrix::identity(m));

        let rng = &mut rand::thread_rng();
        let len = 5;
        let vector = (0..len)
            .into_iter()
            .map(|_| Fq::rand(rng))
            .collect::<Vec<_>>();

        let identity = RingMatrix::<Fq>::identity(len);

        let actual = identity.mul_vector(&vector);
        assert_eq!(vector, actual);
    }

    #[test]
    pub fn test_matrix_transpose() {
        let m = 2;
        let matrix = RingMatrix::<Fq> {
            rows: m,
            cols: m,
            values: vec![vec![Fq::one(), Fq::new(2)], vec![Fq::new(3), Fq::new(4)]],
        };
        let transposed: RingMatrix<Fq> = matrix.transpose();
        let expect = RingMatrix::<Fq> {
            rows: m,
            cols: m,
            values: vec![vec![Fq::one(), Fq::new(3)], vec![Fq::new(2), Fq::new(4)]],
        };

        assert_eq!(transposed, expect);

        let recovered: RingMatrix<Fq> = transposed.transpose();
        assert_eq!(recovered, matrix);
    }

    #[test]
    fn test_matrix_add_and_sub() {
        let cols = 10;
        let rows = 20;
        let rng = &mut rand::thread_rng();
        let lhs = RingMatrix::<Fq>::rand(rng, rows, cols);
        let rhs = RingMatrix::<Fq>::rand(rng, rows, cols);

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
    #[should_panic(expected = "Vectors must have the same length")]
    fn test_vec_mul_unequal_lengths() {
        let vec1 = vec![Fq::one(), Fq::zero()];
        let vec2 = vec![Fq::zero(), Fq::one(), Fq::new(2)];
        RingMatrix::<Fq>::vec_mul(&vec1, &vec2);
    }
    #[test]
    fn test_mul_vector() {
        let m: usize = 2;
        // | 1 0 |
        // | 0 1 |
        let matrix = vec![vec![Fq::one(), Fq::zero()], vec![Fq::zero(), Fq::one()]];

        let vector = vec![Fq::one(), Fq::zero()];

        let a = RingMatrix::<Fq> {
            rows: m,
            cols: m,
            values: matrix,
        };

        let res = a.mul_vector(&vector);
        assert_eq!(vector, res);
    }

    #[test]
    #[should_panic(expected = "Matrix columns must match vector length")]
    fn test_mul_vector_mismatch() {
        let matrix = RingMatrix::<Fq> {
            rows: 2,
            cols: 2,
            values: vec![vec![Fq::one(), Fq::zero()], vec![Fq::zero(), Fq::one()]],
        };
        let vector = vec![Fq::one()];
        matrix.mul_vector(&vector);
    }
    #[test]
    fn test_mul_matrix() {
        let m: usize = 2;
        let a = RingMatrix::<Fq> {
            rows: m,
            cols: m,
            values: vec![vec![Fq::one(), Fq::zero()], vec![Fq::zero(), Fq::one()]],
        };
        let b = a.clone();

        let res = a.mul_matrix(&b);
        assert_eq!(a.values, res.values);
        let lhs = a * b;
        assert_eq!(lhs.values, res.values);
    }

    #[test]
    #[should_panic(expected = "Matrix dimensions must be compatible for multiplication")]
    fn test_mul_matrix_incompatible() {
        let a = RingMatrix::<Fq> {
            rows: 2,
            cols: 3,
            values: vec![vec![Fq::one(), Fq::zero(), Fq::one()], vec![
                Fq::zero(),
                Fq::one(),
                Fq::zero(),
            ]],
        };
        let b = RingMatrix::<Fq> {
            rows: 2,
            cols: 2,
            values: vec![vec![Fq::one(), Fq::zero()], vec![Fq::zero(), Fq::one()]],
        };
        a.mul_matrix(&b);
    }
    #[test]
    fn test_matrix_scalar_mul_and_matrix_mul() {
        let n = 2;
        let mut rng = rand::thread_rng();

        let A = RingMatrix::<Fq>::rand(&mut rng, n, n);
        let B = RingMatrix::<Fq>::rand(&mut rng, n, n);
        let x = vec![Fq::new(3), Fq::new(5)];

        // A*B*x
        let res1 = RingMatrix::<Fq>::mul_matrix(&A, &B).mul_vector(&x);
        // A*(B*x)
        let res2 = A.mul_vector(&B.mul_vector(&x));
        assert_eq!(res1, res2);
    }
}
