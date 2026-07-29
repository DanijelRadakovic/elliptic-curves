use crate::finite_field::{FiniteField, FiniteFieldError};
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
#[derive(Debug)]
pub struct EllipticCurve {
    a: BigUint,
    b: BigUint,
    f: FiniteField,
}

#[derive(Error, Debug, PartialEq)]
pub enum EllipticCurveError {
    #[error("Invalid curve, determinant is zero")]
    InvalidCurve,
    #[error("Cannot multiply by zero scalar")]
    ZeroScalar,
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
                let x = self.f.sub(&self.f.sub(&self.f.pow(&TWO, &s), &x1), &x2);
                // y = s (x1 - x) - y1
                let y = self.f.sub(&self.f.mul(&s, &self.f.sub(&x1, &x)), &y1);
                Ok(Point::Coordinate { x, y })
            }
        }
    }

    // TODO: Improve with bitwise operations
    pub fn scalar_mul(&self, n: &BigUint, point: &Point) -> Result<Point, EllipticCurveError> {
        // 0 * P = Identity
        if n.is_zero() {
            return Ok(Point::Identity);
        }
        if n.is_zero() {
            return Ok(Point::Identity);
        }

        // Identity * n = Identity
        if *point == Point::Identity {
            return Ok(Point::Identity);
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
            Point::Identity => true,
        }
    }

    pub fn a(&self) -> &BigUint { &self.a }
    pub fn b(&self) -> &BigUint { &self.b }
    pub fn p(&self) -> &BigUint { &self.f.p }

    fn calculate_slope(
        &self,
        point1: &Point,
        point2: &Point,
    ) -> Result<BigUint, EllipticCurveError> {
        let f = &self.f;
        match (point1, point2) {
            (Point::Coordinate { x: x1, y: y1 }, Point::Coordinate { x: x2, y: y2 }) => {
                let (numerator, denominator) = if point1 == point2 {
                    (
                        f.add(&f.mul(&THREE, &f.pow(&TWO, &x1)), &self.a),
                        f.mul(&TWO, &y1),
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
                reason: "Cannot calculate slope for Identity point".into(),
            }),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn setup_curve() -> EllipticCurve {
        // y^2 = x^3 + 2x + 2 (mod 17)
        EllipticCurve::new(
            2u32.into(),
            2u32.into(),
            17u32.into(),
        ).unwrap()
    }

    #[test]
    fn test_point_on_curve() {
        let curve = setup_curve();
        let p = Point::Coordinate { x: 5u32.into(), y: 1u32.into() };
        assert!(curve.is_on_curve(&p));

        let p_bad = Point::Coordinate { x: 5u32.into(), y: 2u32.into() };
        assert!(!curve.is_on_curve(&p_bad));
    }

    #[test]
    fn test_identity_rules() {
        let curve = setup_curve();
        let p = Point::Coordinate { x: 5u32.into(), y: 1u32.into() };

        assert_eq!(curve.add(&p, &Point::Identity).unwrap(), p);
        assert_eq!(curve.add(&Point::Identity, &p).unwrap(), p);
    }

    #[test]
    fn test_point_addition() {
        let curve = setup_curve();
        let p1 = Point::Coordinate { x: 5u32.into(), y: 1u32.into() };
        let p2 = Point::Coordinate { x: 6u32.into(), y: 3u32.into() };

        let res = curve.add(&p1, &p2).unwrap();
        assert_eq!(res, Point::Coordinate { x: 10u32.into(), y: 6u32.into() });
    }

    #[test]
    fn test_point_doubling() {
        let curve = setup_curve();
        let p = Point::Coordinate { x: 5u32.into(), y: 1u32.into() };
        let res = curve.add(&p, &p).unwrap();

        assert!(curve.is_on_curve(&res));
        assert_eq!(res, Point::Coordinate { x: 6u32.into(), y: 3u32.into() });
    }

    #[test]
    fn test_inverse_addition() {
        let curve = setup_curve();
        let p1 = Point::Coordinate { x: 5u32.into(), y: 1u32.into() };
        let p2 = Point::Coordinate { x: 5u32.into(), y: 16u32.into() };

        assert_eq!(curve.add(&p1, &p2).unwrap(), Point::Identity);
    }

    #[test]
    fn test_scalar_multiplication() {
        let curve = setup_curve();
        let p = Point::Coordinate { x: 5u32.into(), y: 1u32.into() };

        let res2 = curve.scalar_mul(&2u32.into(), &p).unwrap();
        let manual2 = curve.add(&p, &p).unwrap();
        assert_eq!(res2, manual2);

        let res0 = curve.scalar_mul(&0u32.into(), &p).unwrap();
        assert_eq!(res0, Point::Identity);

        let res_identity = curve.scalar_mul(&5u32.into(), &Point::Identity).unwrap();
        assert_eq!(res_identity, Point::Identity);
    }

    #[test]
    fn test_invalid_curve_determinant() {
        // Curve where 4a^3 + 27b^2 = 0
        let res = EllipticCurve::new(0u32.into(), 0u32.into(),17u32.into());
        assert!(res.is_err());
        assert_eq!(res.unwrap_err(), EllipticCurveError::InvalidCurve);
    }
}