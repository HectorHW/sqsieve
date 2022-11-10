use std::iter::repeat_with;
use std::ops::Rem;

use crate::number_type::NumberOps;
use crate::{sieve::SmoothNumber, solver::Solution};
use itertools::Itertools;
use lazy_static::lazy_static;
use num_bigint::BigUint;
use num_traits::FromPrimitive;
use rand::{thread_rng, Rng};

pub fn is_big_square(number: &BigUint) -> bool {
    let t = number.sqrt();
    &(t.clone() * t) == number
}

pub fn euclid(x: &BigUint, y: &BigUint) -> BigUint {
    num_integer::gcd(x.clone(), y.clone())
}

/// `a` is expected to be greater
fn test_factorization(n: &BigUint, a: &BigUint, b: &BigUint) -> Option<(BigUint, BigUint)> {
    let mut a = a.rem(n);
    let mut b = b.rem(n);

    if b > a {
        std::mem::swap(&mut a, &mut b);
    }

    let gcd = euclid(n, &(a.clone() - b));
    if &gcd != n && gcd.bits() != 1 {
        let other_term = n / gcd.clone();
        return Some((gcd, other_term));
    }
    None
}

fn is_zero(vec: &[bool]) -> bool {
    vec.iter().all(|&item| !item)
}

#[inline]
fn increase(vec: &mut [bool]) {
    let mut carry = true;
    for digit in vec {
        let new_carry = *digit && carry;
        *digit = *digit != carry;
        carry = new_carry;
    }
}

lazy_static! {
    static ref ONE: BigUint = BigUint::from_i32(1).unwrap();
    static ref TWO: BigUint = BigUint::from_i32(2).unwrap();
}

fn build_numbers<NT: NumberOps>(
    n: &BigUint,
    smoothies: &[SmoothNumber<NT>],
    positions: &[usize],
) -> (BigUint, BigUint) {
    let (mut a, mut b) = (ONE.clone(), ONE.clone());

    for &index in positions {
        a *= smoothies[index].number.to_varsize();
        b *= smoothies[index].number.to_varsize().modpow(&TWO, n);
    }

    b = b.sqrt();

    (a, b)
}

#[inline]
fn attempt_factorization<NT: NumberOps>(
    n: &BigUint,
    smoothies: &[SmoothNumber<NT>],
    solution: &Solution,
    substitution_vector: &[bool],
) -> Option<(BigUint, BigUint)> {
    let inclusion = solution
        .subsitute(substitution_vector, false)
        .into_iter()
        .enumerate()
        .filter_map(|(i, p)| if p { Some(i) } else { None })
        .collect_vec();

    let (a, b) = build_numbers(n, smoothies, &inclusion);

    test_factorization(n, &a, &b)
}

fn search_lonelies<NT: NumberOps>(
    n: &NT,
    smoothies: &[SmoothNumber<NT>],
    solution: &Solution,
) -> Option<(BigUint, BigUint)> {
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
    None
}

pub fn find_factor_exhaustive<NT: NumberOps>(
    n: &NT,
    smoothies: &[SmoothNumber<NT>],
    solution: &Solution,
) -> Option<(BigUint, BigUint)> {
    //first, try lonely numbers - maybe, we will find perfect square without multiplying

    println!("using exhaustive search");

    if let Some(answ) = search_lonelies(n, smoothies, solution) {
        return Some(answ);
    }

    let n = n.to_varsize();

    if solution.free_variables.is_empty() {
        // no nontrivial except for lonely numbers which all produce trivials
        return None;
    }

    // if it fails, it means that all lonely numbers produce trivial factorization, do not include them

    let limit = solution.free_variables.len();

    println!("variables for exhaustive search: {limit}");

    let mut free_mapping = vec![false; solution.free_variables.len()];

    increase(&mut free_mapping);

    while !is_zero(&free_mapping) {
        if let Some(sol) = attempt_factorization(&n, smoothies, solution, &free_mapping) {
            return Some(sol);
        }
        increase(&mut free_mapping);
    }

    None
}

pub fn find_factor_simple<NT: NumberOps>(
    n: &NT,
    smoothies: &[SmoothNumber<NT>],
    solution: &Solution,
) -> Option<(BigUint, BigUint)> {
    //first, try lonely numbers - maybe, we will find perfect square without multiplying

    println!("trying base vector search");

    if let Some(answ) = search_lonelies::<NT>(n, smoothies, solution) {
        return Some(answ);
    }

    let n = n.to_varsize();

    if solution.free_variables.is_empty() {
        // no nontrivial except for lonely numbers which all produce trivials
        return None;
    }

    // if it fails, it means that all lonely numbers produce trivial factorization, do not include them

    let limit = solution.free_variables.len();

    println!("variables for base vector search: {limit}");

    let mut free_mapping = vec![false; solution.free_variables.len()];

    for i in 0..free_mapping.len() {
        if i > 0 {
            free_mapping[i - 1] = false;
        }
        free_mapping[i] = true;

        if let Some(sol) = attempt_factorization(&n, smoothies, solution, &free_mapping) {
            return Some(sol);
        }
    }

    None
}

pub fn find_factors_random<NT: NumberOps>(
    n: &NT,
    smoothies: &[SmoothNumber<NT>],
    solution: &Solution,
) -> Option<(BigUint, BigUint)> {
    //first, try lonely numbers - maybe, we will find perfect square without multiplying

    println!("trying random search");

    if let Some(answ) = search_lonelies(n, smoothies, solution) {
        return Some(answ);
    }

    let n = n.to_varsize();

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

            b = b.sqrt();

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

pub fn find_factors_from_pivots<NT: NumberOps>(
    n: &NT,
    smoothies: &[SmoothNumber<NT>],
    solution: &[Vec<usize>],
) -> Option<(BigUint, BigUint)> {
    let n = n.to_varsize();
    for vector in solution {
        let (a, b) = build_numbers(&n, smoothies, vector);

        if let Some(answ) = test_factorization(&n, &a, &b) {
            return Some(answ);
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::{increase, is_zero};

    #[test]
    fn should_turn_0_into_1() {
        let mut v = vec![false, false, false];
        increase(&mut v);
        assert_eq!(v, vec![true, false, false]);
    }

    #[test]
    fn should_apply_carry_bit() {
        let mut v = vec![true, true, true, false];
        increase(&mut v);
        assert_eq!(v, vec![false, false, false, true]);
    }

    #[test]
    fn should_increase() {
        let mut v = vec![true, false, true];
        increase(&mut v);
        assert_eq!(v, vec![false, true, true]);
    }

    #[test]
    fn should_be_zero() {
        let v = vec![false, false, false];
        assert!(is_zero(&v))
    }
}
