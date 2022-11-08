use std::collections::HashMap;
use std::iter::repeat_with;

use crate::number_type::{NumberOps, NumberType};
use crate::{sieve::SmoothNumber, solver::Solution};
use num_bigint::BigUint;
use num_traits::FromPrimitive;
use rand::{thread_rng, Rng};

pub fn is_square(number: usize) -> bool {
    let t = ((number as f64).sqrt() + 0.5).floor() as usize;
    t * t == number
}

pub fn is_big_square(number: &BigUint) -> bool {
    let t = number.sqrt();
    &(t.clone() * t) == number
}

pub fn euclid(x: &BigUint, y: &BigUint) -> BigUint {
    num_integer::gcd(x.clone(), y.clone())
}

/// `a` is expected to be greater
fn test_factorization(n: &BigUint, a: &BigUint, b: &BigUint) -> Option<(BigUint, BigUint)> {
    let gcd = euclid(n, &(a.clone() - b));
    if &gcd != n && gcd.bits() != 1 {
        let other_term = n / gcd.clone();
        return Some((gcd, other_term));
    }
    None
}

pub fn find_factor(
    n: &NumberType,
    smoothies: &[SmoothNumber<NumberType>],
    solution: &Solution,
) -> Option<(BigUint, BigUint)> {
    //first, try lonely numbers - maybe, we will find perfect square without multiplying

    let n = n.to_varsize();

    let two = BigUint::from_i32(2).unwrap();

    for &lonely_factor in &solution.lonely_variables {
        let candidate = &smoothies[lonely_factor];
        let candidate_number = candidate.number.to_varsize();

        let p2 = candidate.number.to_varsize().modpow(&two, &n);
        if is_big_square(&p2) {
            let Some(solution) = test_factorization(&n, &candidate_number, &p2.sqrt()) else {
                continue;
            };

            return Some(solution);
        }
    }

    if solution.free_variables.is_empty() {
        // no nontrivial except for lonely numbers which all produce trivials
        return None;
    }

    // if it fails, it means that all lonely numbers produce trivial factorization, do not include them

    let attempts = 2usize
        .checked_pow(solution.free_variables.len() as u32)
        .map(|v| v.min(10_000))
        .unwrap_or(10_000);

    let mut rng = thread_rng();

    let one = BigUint::from_i32(1).unwrap();

    let limit = solution.free_variables.len();

    println!("limit: {limit}");

    let mut current_limit = limit;

    let pressure_variants = repeat_with(|| {
        let next_value = current_limit;
        current_limit /= 2;
        next_value as u32
    })
    .take_while(|&n| n >= 2);

    let mut free_mapping = vec![false; solution.free_variables.len()];

    for pressure in pressure_variants {
        println!("trying 1/{}", pressure);
        for _attempt in 0..attempts {
            free_mapping
                .iter_mut()
                .for_each(|position| *position = rng.gen_ratio(1, pressure));

            let inclusion = solution.subsitute(&free_mapping, false);

            let (mut a, mut b) = smoothies
                .iter()
                .zip(inclusion.iter())
                .filter_map(|(a, &b)| if b { Some(a) } else { None })
                .fold((one.clone(), one.clone()), |acc, item| {
                    let left = acc.0 * &item.number.to_varsize();
                    let right = acc.1
                        * item
                            .divisors
                            .iter()
                            .map(|&(divisor, power)| BigUint::from(divisor).pow(power as u32))
                            .product::<BigUint>();
                    (left, right)
                });

            if b > a {
                std::mem::swap(&mut a, &mut b);
            }

            let Some(solution) = test_factorization(&n, &a, &b) else {
                continue;
            };

            return Some(solution);
        }
    }

    None
}
