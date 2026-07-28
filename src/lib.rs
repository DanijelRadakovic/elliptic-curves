use num_bigint::BigUint;
use num_iter::range;
use num_traits::{One, Zero};
use thiserror::Error;

const TWO: BigUint = BigUint::new_const(2);
const THREE: BigUint = BigUint::new_const(3);
const FOUR: BigUint = BigUint::new_const(4);
const TWENTY_SEVEN: BigUint = BigUint::new_const(27);

#[derive(Eq, PartialEq, Clone, Debug)]
pub struct Point {
    pub x: BigUint,
    pub y: BigUint,
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
}

impl EllipticCurve {
    pub fn new(a: BigUint, b: BigUint, p: BigUint) -> Result<Self, EllipticCurveError> {
        // Ensure the curve is non-singular: 4a^3 + 27b^2 != 0 (mod p)
        let a_cubed = a.modpow(&THREE, &p);
        let b_squared = b.modpow(&TWO, &p);

        let discriminant = (&FOUR * a_cubed + &TWENTY_SEVEN * b_squared) % &p;

        if discriminant.is_zero() {
            return Err(EllipticCurveError::InvalidCurve);
        }

        Ok(Self {
            a,
            b,
            f: FiniteField { p },
        })
    }

    pub fn add(&self, point1: &Point, point2: &Point) -> Result<Point, EllipticCurveError> {
        let s = self.calculate_slope(point1, point2).map_err(|e| {
            EllipticCurveError::OperationFailed {
                op: format!("{:?} + {:?}", point1, point2),
                reason: e.to_string(),
            }
        })?;
        // x = s^2 - x1 - x2
        let x = self
            .f
            .sub(&self.f.sub(&self.f.pow(&TWO, &s), &point1.x), &point2.x);
        // y = s (x1 - x) - y1
        let y = self
            .f
            .sub(&self.f.mul(&s, &self.f.sub(&point1.x, &x)), &point1.y);
        Ok(Point { x, y })
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

    fn calculate_slope(&self, point1: &Point, point2: &Point) -> Result<BigUint, FiniteFieldError> {
        let s: BigUint;
        let numerator: BigUint;
        let denominator: BigUint;
        let f = &self.f;
        if point1 == point2 {
            // (3 * x^2 + a) % p
            numerator = f.add(&f.mul(&THREE, &f.pow(&TWO, &point1.x)), &self.a);
            // (2 * y) % p
            denominator = f.mul(&TWO, &point1.y);

            let inv_denom = f.mul_inv(&denominator)?;
            s = f.mul(&numerator, &inv_denom);
            Ok(s)
        } else {
            numerator = &point2.y - &point1.y;
            denominator = &point2.x - &point1.x;
            let inv_denom = f.mul_inv(&denominator)?;
            s = f.mul(&numerator, &inv_denom);
            Ok(s)
        }
    }
}

pub struct FiniteField {
    pub p: BigUint,
}

#[derive(Error, Debug, PartialEq)]
pub enum FiniteFieldError {
    #[error("Cannot divide by zero in a finite field")]
    DivisionByZero,
    #[error("Zero does not have modular inverse")]
    ZeroHasNoInverse,
    #[error("Modular inverse does not exist (modulus might not be prime)")]
    NoInverse,
}

impl FiniteField {
    pub fn add(&self, a: &BigUint, b: &BigUint) -> BigUint {
        (a + b) % &self.p
    }

    pub fn sub(&self, a: &BigUint, b: &BigUint) -> BigUint {
        // Easy implementation: add(a, add_inv(b)).
        // Improved implementation: a -b = ((a + p - (b mod p)) mod p
        (a + &self.p - (b % &self.p)) % &self.p
    }

    pub fn add_inv(&self, a: &BigUint) -> BigUint {
        // a + x = 0 mod p --> x = p - a
        // There are cases that needs to be aware of:
        // 1. a = 0, p - 0 = p, but should be 0.
        // 2. if a > p, substraction will crash.

        // 1. (a % p) ensures the number is within the field range, handles if a > p.
        // 2. (p - (a % p)) finds the difference
        // 3. The final % p ensures that if the result is p, it becomes 0
        (&self.p - (a % &self.p)) % &self.p
    }

    pub fn mul(&self, a: &BigUint, b: &BigUint) -> BigUint {
        (a * b) % &self.p
    }

    pub fn div(&self, a: &BigUint, b: &BigUint) -> Result<BigUint, FiniteFieldError> {
        if b.is_zero() {
            return Err(FiniteFieldError::DivisionByZero);
        }
        Ok((a * self.mul_inv(b)?) % &self.p)
    }

    pub fn pow(&self, n: &BigUint, a: &BigUint) -> BigUint {
        a.modpow(n, &self.p)
    }

    pub fn mul_inv(&self, a: &BigUint) -> Result<BigUint, FiniteFieldError> {
        if a.is_zero() {
            return Err(FiniteFieldError::ZeroHasNoInverse);
        }
        a.modinv(&self.p).ok_or(FiniteFieldError::NoInverse)
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
        let p1 = Point {
            x: BigUint::from(5u32),
            y: BigUint::from(1u32),
        };
        let p2 = Point {
            x: BigUint::from(3u32),
            y: BigUint::from(1u32),
        };
        let res = curve.add(&p1, &p2).unwrap();
        println!("Result: x={}, y={}", res.x, res.y);
    }
}
