pub mod finite_field;

use finite_field::{FiniteField, FiniteFieldError};
use num_bigint::BigUint;
use num_iter::range;
use num_traits::{One, Zero};
use thiserror::Error;

const TWO: BigUint = BigUint::new_const(2);
const THREE: BigUint = BigUint::new_const(3);
const FOUR: BigUint = BigUint::new_const(4);
const TWENTY_SEVEN: BigUint = BigUint::new_const(27);

#[derive(Eq, PartialEq, Clone, Debug)]
pub enum Point {
    Coordinate { x: BigUint, y: BigUint },
    Identity,
}

// E: y^2 = x^3 + ax + b mod p
pub struct EllipticCurve {
    pub a: BigUint,
    pub b: BigUint,
    f: FiniteField,
}

#[derive(Error, Debug, PartialEq)]
pub enum EllipticCurveError {
    #[error("Invalid curve, determinant is zero")]
    InvalidCurve,
    #[error("Cannot multiply by zero scalar")]
    ZeroScalar,
    #[error("Point ({x}, {y}) is not on the curve")]
    PointNotOnCurve { x: BigUint, y: BigUint },
    #[error("Operation {op} failed: {reason}")]
    OperationFailed { op: String, reason: String },
    #[error(transparent)]
    FieldError(#[from] FiniteFieldError),
}

impl EllipticCurve {
    pub fn new(a: BigUint, b: BigUint, p: BigUint) -> Result<Self, EllipticCurveError> {
        let f = FiniteField { p };
        // Ensure the curve is non-singular: 4a^3 + 27b^2 != 0 (mod p)
        let a_cubed = f.pow(&THREE, &a);
        let b_squared = f.pow(&TWO, &b);
        let discriminant = f.add(&f.mul(&FOUR, &a_cubed), &f.mul(&TWENTY_SEVEN, &b_squared));

        if discriminant.is_zero() {
            return Err(EllipticCurveError::InvalidCurve);
        }

        Ok(Self { a, b, f })
    }

    pub fn add(&self, point1: &Point, point2: &Point) -> Result<Point, EllipticCurveError> {
        match (point1, point2) {
            (Point::Identity, Point::Identity) => Ok(Point::Identity),
            (Point::Identity, p) => Ok(p.clone()),
            (p, Point::Identity) => Ok(p.clone()),
            (Point::Coordinate { x: x1, y: y1 }, Point::Coordinate { x: x2, y: y2 }) => {
                // P + (-P) = Identity
                if x1 == x2 && y1 != y2 {
                    return Ok(Point::Identity);
                }

                let s = self.calculate_slope(point1, point2).map_err(|e| {
                    EllipticCurveError::OperationFailed {
                        op: format!("{:?} + {:?}", point1, point2),
                        reason: e.to_string(),
                    }
                })?;
                // x = s^2 - x1 - x2
                let x = self
                    .f
                    .sub(&self.f.sub(&self.f.pow(&TWO, &s), &x1), &x2);
                // y = s (x1 - x) - y1
                let y = self
                    .f
                    .sub(&self.f.mul(&s, &self.f.sub(&x1, &x)), &y1);
                Ok(Point::Coordinate { x, y })
            }
        }
    }

    // TODO: Improve with bitwise operations
    pub fn scalar_mul(&self, n: &BigUint, point: &Point) -> Result<Point, EllipticCurveError> {
        if n.is_zero() {
            return Err(EllipticCurveError::ZeroScalar);
        }

        let mut result = point.clone();
        for _ in range(BigUint::one(), n.clone()) {
            result = self.add(&result, point)?
        }
        Ok(result)
    }

    pub fn is_on_curve(&self, point: &Point) -> bool {
        match point {
            Point::Coordinate { x, y } => {
                let x_cubed = self.f.pow(&THREE, &x);
                let a_x = self.f.mul(&self.a, &x);
                self.f.pow(&TWO, y) == self.f.add(&self.f.add(&x_cubed, &a_x), &self.b)
            }
            Point::Identity => true
        }
    }

    fn calculate_slope(&self, point1: &Point, point2: &Point) -> Result<BigUint, EllipticCurveError> {
        let f = &self.f;
        match (point1, point2) {
            (Point::Coordinate { x: x1, y: y1 }, Point::Coordinate { x: x2, y: y2 }) => {
                let (numerator, denominator) = if point1 == point2 {
                    (
                        f.add(&f.mul(&THREE, &f.pow(&TWO, &x1)), &self.a),
                        f.mul(&TWO, &y1)
                    )
                } else {
                    (f.sub(&y2, &y1), f.sub(&x2, &x1))
                };

                let inv_denom = f.mul_inv(&denominator)?;
                let s = f.mul(&numerator, &inv_denom);
                Ok(s)
            }
            _ => Err(EllipticCurveError::OperationFailed {
                op: "calculating slope".into(),
                reason: "Cannot calculate slope for Identity point".into()
            })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ecc_add() {
        let curve = EllipticCurve {
            a: BigUint::from(2u32),
            b: BigUint::from(2u32),
            f: FiniteField {
                p: BigUint::from(17u32),
            },
        };
        let p1 = Point::Coordinate {
            x: BigUint::from(5u32),
            y: BigUint::from(1u32),
        };
        let p2 = Point::Coordinate {
            x: BigUint::from(3u32),
            y: BigUint::from(1u32),
        };
        let res = curve.add(&p1, &p2).unwrap();
        assert_eq!(
            res,
            Point::Coordinate {
                x: BigUint::from(9u32),
                y: BigUint::from(16u32)
            }
        );
    }
}
