use num_bigint::BigUint;
use num_traits::{Num};
use elliptic_curves::ecdsa::ECDSA;
use elliptic_curves::elliptic_curve::{EllipticCurve, EllipticCurveError, Point};
use elliptic_curves::elliptic_curve::Point::{Coordinate, Identity};

fn setup_secp256k1() -> EllipticCurve {
    let p = BigUint::from_str_radix("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16).unwrap();
    let a = BigUint::from(0u32);
    let b = BigUint::from(7u32);
    EllipticCurve::new(a, b, p).unwrap()
}

fn g() -> Point {
    let x = BigUint::from_str_radix("79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16).unwrap();
    let y = BigUint::from_str_radix("483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8", 16).unwrap();
    Coordinate { x, y }
}

fn n() -> BigUint {
    BigUint::from_str_radix("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141", 16).unwrap()
}

#[test]
fn test_multiplication() {
    let ec = setup_secp256k1();

    let g = g();
    let res = ec.scalar_mul(&1u32.into(), &g).unwrap();
    assert_eq!(res, g);

    // 7 * G
    let point = Coordinate {
        x: BigUint::from_str_radix("41948375291644419605210209193538855353224492619856392092318293986323063962044", 10).unwrap(),
        y: BigUint::from_str_radix("48361766907851246668144012348516735800090617714386977531302791340517493990618", 10).unwrap(),
    };
    let res = ec.scalar_mul(&7u32.into(), &g).unwrap();
    assert_eq!(res, point);

    // n * G = Identity
    let res = ec.scalar_mul(&n(), &g).unwrap();
    assert_eq!(res, Identity);

    let point = ec.scalar_mul(&5u32.into(), &Coordinate {x: 1u32.into(), y: 1u32.into()});
    assert!(point.is_err());
    assert_eq!(point.unwrap_err(), EllipticCurveError::PointNotOnCurve {x: 1u32.into() ,y: 1u32.into(), a: ec.a().clone(), b: ec.b().clone(), p: ec.p().clone()});
}

#[test]
fn test_signature() {
    let ec = setup_secp256k1();
    let ecdsa = ECDSA::new(ec, g(), n()).unwrap();

    let (private_key, public_key) = ecdsa.generate_key_pair().unwrap();
    let msg = b"For the Emperor!";
    let digest = BigUint::from_str_radix(&sha256::digest(msg), 16).unwrap();
    let signature = ecdsa.sign(&digest, &private_key).unwrap();
    let result = ecdsa.verify(&digest, &signature, &public_key).unwrap();
    assert_eq!(result, true);
}