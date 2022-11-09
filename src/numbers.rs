use crate::number_type::{NumberOps, NumberType};

use crypto_bigint::Zero;
use itertools::Itertools;
use num_bigint::BigUint;

pub fn small_eratosphenes(upper_limit: usize) -> Vec<usize> {
    assert!(
        upper_limit.checked_add(1).is_some(),
        "upper limit overflows usize"
    );
    // TODO: run this for large numbers when computers will have 2**63 bytes of ram
    let mut candidates = vec![false; upper_limit + 1];
    let mut result = vec![];
    for number in 2..=upper_limit {
        if candidates[number] {
            continue;
        }

        result.push(number);

        if let Some(lower_bound) = number.checked_mul(number) {
            let mut multiple = lower_bound;
            while multiple <= upper_limit {
                candidates[multiple] = true;
                multiple += number;
            }
        }
    }
    result
}

fn legendre(a: BigUint, p: BigUint) -> isize {
    use num_traits::Zero;
    let power = (p.clone() - BigUint::from(1usize)) / BigUint::from(2usize);

    let res = a.modpow(&power, &p);

    if res == (p - BigUint::from(1usize)) {
        return -1;
    }

    if res.is_zero() {
        return 0;
    }
    1
}

pub fn build_factor_base(primes: Vec<usize>, n: &NumberType) -> Vec<usize> {
    primes
        .into_iter()
        .filter(|&prime| legendre(n.to_varsize(), BigUint::from(prime)) == 1)
        .collect_vec()
}

pub fn small_smooth(upper_limit: usize, b_value: usize) -> Vec<usize> {
    // same algorithm as above but instead we divide numbers
    assert!(
        upper_limit.checked_add(1).is_some(),
        "upper limit overflows usize"
    );
    let mut candidates = (0..=upper_limit).into_iter().collect_vec();
    for number in 2..=b_value {
        if candidates[number] == 1 {
            continue;
        }

        candidates[number] = 1;

        if let Some(lower_bound) = number.checked_mul(2) {
            let mut multiple = lower_bound;
            while multiple <= upper_limit {
                while candidates[multiple] % number == 0 {
                    candidates[multiple] /= number;
                }

                multiple += number;
            }
        }
    }
    candidates
        .into_iter()
        .enumerate()
        .skip(1) // do not include 0
        .filter_map(|(i, n)| if n == 1 { Some(i) } else { None })
        .collect_vec()
}

pub type DenseMultiplierMap = Vec<(usize, usize)>;

pub fn factor_smooth(n: &NumberType, prime_table: &[usize]) -> Option<DenseMultiplierMap> {
    let mut result = vec![];

    let mut n = n.clone();

    'outer: for &prime in prime_table {
        loop {
            let (d, r) = n.divmod(prime);
            if r.is_zero().unwrap_u8() == 0 {
                break;
            }
            match result.last_mut() {
                Some((n, d)) if n == &prime => {
                    *d += 1;
                }
                _ => result.push((prime, 1)),
            }

            n = d;

            if &n == NumberType::one() {
                break 'outer;
            }
        }
    }
    if &n == NumberType::one() {
        Some(result)
    } else {
        None
    }
}
