use crate::number_type::NumberOps;

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

pub fn legendre(a: BigUint, p: BigUint) -> isize {
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

pub fn gcd(a: usize, b: usize) -> usize {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

fn upcast_modpow(x: usize, e: usize, m: usize) -> usize {
    BigUint::modpow(&BigUint::from(x), &BigUint::from(e), &BigUint::from(m))
        .to_u64_digits()
        .pop()
        .unwrap_or_default() as usize
}

/// finds x from x^2 === n (mod p)
#[allow(non_snake_case)]
pub fn tonelli_shanks(n: usize, p: usize) -> Option<usize> {
    if gcd(n, p) != 1 || legendre(BigUint::from(n), BigUint::from(p)) == -1 {
        return None;
    }

    let mut q = p - 1;
    let mut s = 0;
    while q % 2 == 0 {
        s += 1;
        q /= 2;
    }

    let mut z_value = 2;
    loop {
        if upcast_modpow(z_value, (p - 1) / 2, p) == p - 1 {
            break;
        }
        z_value += 1;
    }

    let mut M = s;
    let mut c = upcast_modpow(z_value, q, p);
    let mut t = upcast_modpow(n, q, p);

    let mut R = upcast_modpow(n, (q + 1) / 2, p);

    loop {
        if t == 0 {
            return Some(0);
        }

        if t == 1 {
            return Some(R);
        }

        let mut i = 1;
        while i < M {
            if upcast_modpow(t, 2usize.pow(i), p) == 1 {
                break;
            }

            i += 1;
        }

        let b = upcast_modpow(c, 2usize.pow(M - i - 1), p);

        M = i;
        c = upcast_modpow(b, 2, p);
        t = t * b * b % p;
        R = R * b % p;
    }
}

pub fn build_factor_base<NT: NumberOps>(primes: Vec<usize>, n: &NT) -> Vec<usize> {
    primes
        .into_iter()
        .filter(|&prime| prime == 2 || legendre(n.to_varsize(), BigUint::from(prime)) == 1)
        .collect_vec()
}

pub type DenseMultiplierMap = Vec<(usize, usize)>;

pub fn trial_divide<NT: NumberOps>(n: &NT, prime_table: &[usize]) -> Option<DenseMultiplierMap> {
    let mut result = vec![];

    let mut n = *n;

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

            if &n == NT::one() {
                break 'outer;
            }
        }
    }
    if &n == NT::one() {
        Some(result)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use num_bigint::BigUint;

    use super::tonelli_shanks;

    //upcast modpow relies on that
    #[test]
    fn bigint_returns_nothing_in_state_of_zero() {
        assert_eq!(BigUint::from(0usize).to_u64_digits(), vec![])
    }

    #[test]
    fn should_compute_value() {
        assert_eq!(tonelli_shanks(3, 13), Some(9))
    }

    #[test]
    fn should_solve_example() {
        assert_eq!(tonelli_shanks(5, 41), Some(28))
    }
}
