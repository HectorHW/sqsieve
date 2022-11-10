use std::{env::args, error::Error, iter::repeat, str::FromStr};

use crate::factor_building::find_factors_from_pivots;
use crate::number_type::NumberOps;
use crate::sieve::SmoothNumber;
use crypto_bigint::UInt;
use itertools::Itertools;
use num_bigint::BigUint;
use num_integer::Roots;
use num_traits::ToPrimitive;

use crate::{
    factor_building::{find_factor_exhaustive, find_factor_simple, find_factors_random},
    numbers::{build_factor_base, small_eratosphenes},
    sieve::{compute_b_limit, BlockSieve, LogSieve, TestDivisionSieve},
    solver::{produce_solution, CongruenceSystem},
};
fn gaussian_multistage<NT: NumberOps>(
    n: &NT,
    table: Vec<SmoothNumber<NT>>,
    mut system: CongruenceSystem,
) -> Option<(BigUint, BigUint)> {
    system.diagonalize();

    #[cfg(feature = "verbose")]
    println!("{}", system.print_as_dense());

    #[cfg(feature = "verbose")]
    println!("---------");

    let solution = produce_solution(&system);

    #[cfg(feature = "verbose")]
    println!("linear system {:?}", solution);

    println!(
        "number of dependencies in solution: {}",
        solution.dependencies.len()
    );

    println!("built solution dependencies, searching for factors");

    find_factor_simple::<NT>(n, &table, &solution)
        .or_else(|| find_factors_random(n, &table, &solution))
        .or_else(|| find_factor_exhaustive(n, &table, &solution))
}

fn pivot_search<NT: NumberOps>(
    n: &NT,
    table: Vec<SmoothNumber<NT>>,
    mut system: CongruenceSystem,
) -> Option<(BigUint, BigUint)> {
    println!("using fast pivot algorithm");
    let vectors = system.fast_pivot();
    println!("produced {} candidates", vectors.len());
    if let Some(answ) = find_factors_from_pivots(n, &table, &vectors) {
        return Some(answ);
    }

    None
}

fn run_factor<NT: NumberOps>(n: &NT, prime_bound: usize) -> Option<(BigUint, BigUint)> {
    let primes = small_eratosphenes(prime_bound);

    println!("primes until bound: {}", primes.len());

    let factor_base = build_factor_base(primes, n);

    println!("built factor base of size {}", factor_base.len(),);

    if factor_base.len() < 2 {
        println!("this is too small");
        return None;
    }

    #[cfg(feature = "truncate")]
    if factor_base.len() >= 250 {
        println!("factor base is too huge, truncating");
        factor_base.truncate(250);
    }

    #[cfg(feature = "verbose")]
    println!("factor base: {:?}", factor_base);

    let mut ratio = 1.05;

    const NUM_ATTEMPTS: usize = 5;

    let mut table = vec![];

    let mut sieve = LogSieve::new(*n, factor_base.clone());

    for _ in 0..NUM_ATTEMPTS {
        let sieving_limit = usize::max(
            (factor_base.len() as f64 * ratio).ceil() as usize,
            factor_base.len() + 5,
        );

        println!("need about {sieving_limit} numbers");

        let mut additional_table = sieve.run(sieving_limit.saturating_sub(table.len()));

        table.append(&mut additional_table);

        println!("done collecting, building solution");

        #[cfg(feature = "verbose")]
        println!("{:?}", table);

        #[cfg(feature = "verbose")]
        for (i, item) in table.iter().enumerate() {
            println!("{i} -> {}", item.number);
        }

        let row_labels = (0..table.len()).collect_vec();
        let rows = table.iter().map(|row| row.divisors.clone()).collect_vec();

        let system =
            CongruenceSystem::with_labels(&rows, factor_base.clone(), row_labels).transpose();

        #[cfg(feature = "verbose")]
        println!("{}", system.print_as_dense());

        #[cfg(feature = "verbose")]
        println!("---------");

        if cfg!(feature = "multistage") {
            if let Some(answ) = gaussian_multistage(n, table.clone(), system) {
                return Some(answ);
            }
        } else if let Some(answ) = pivot_search(n, table.clone(), system) {
            return Some(answ);
        }

        if factor_base.len() < 50 {
            ratio *= 1.5;
        } else {
            ratio += 0.05;
        }

        println!("increasing number of smoothies to find")
    }
    None
}

fn run_factorization_generic<const LIMBS: usize>(bytes: Vec<u8>) -> Result<Vec<BigUint>, String>
where
    UInt<LIMBS>: NumberOps,
{
    println!(
        "used integer size: {}",
        LIMBS * std::mem::size_of::<usize>()
    );

    let bytes = repeat(0u8)
        .take(LIMBS * std::mem::size_of::<usize>() - bytes.len())
        .chain(bytes.into_iter())
        .collect_vec();

    let n = <UInt<LIMBS>>::from_be_slice(&bytes);

    let limit = compute_b_limit(&n).min(10_000);

    let mut prime_bound = 100;

    #[cfg(feature = "test-small")]
    while prime_bound <= limit {
        println!("trying prime bound of {prime_bound}");
        let factorization = run_factor(&n, prime_bound);

        println!("factorization: {factorization:?}");

        if let Some((a, b)) = factorization {
            let prod = a.clone() * b.clone();
            println!("{} * {} = {}", a, b, prod);
            return Ok(());
        }

        prime_bound *= 2;
    }

    prime_bound = limit;

    println!("trying prime bound of {prime_bound}");
    let factorization = run_factor(&n, prime_bound);

    println!("factorization: {factorization:?}");

    if let Some((a, b)) = factorization {
        let prod = a.clone() * b.clone();
        println!("{} * {} = {}", a, b, prod);
        return Ok(vec![a, b]);
    }

    Err("could not factorize".to_string())
}

pub fn factorize(number_repr: String) -> Result<Vec<BigUint>, String> {
    let n = BigUint::from_str(&number_repr).map_err(|e| e.to_string())?;

    println!("n: {}", n);

    println!("base 10 digits: {}", n.to_string().len());

    println!("bit size: {}", n.bits());

    let bytes = n.to_bytes_be();

    if bytes.len() >= 64 {
        return Err("number is too big".into());
    }

    if bytes.len() < 4 {
        return trial_divide(n.to_u64().unwrap() as usize)
            .map(|ok| ok.into_iter().map(BigUint::from).collect_vec());
    }

    if bytes.len() < 8 {
        return run_factorization_generic::<1>(bytes);
    }

    if bytes.len() < 16 {
        return run_factorization_generic::<2>(bytes);
    }

    if bytes.len() < 32 {
        return run_factorization_generic::<4>(bytes);
    }

    run_factorization_generic::<8>(bytes)
}

fn trial_divide(number: usize) -> Result<Vec<usize>, String> {
    let mut to_factor = number;

    let mut divisors = vec![];

    for i in 2..number {
        while to_factor % i == 0 {
            divisors.push(i);
            to_factor /= i;
        }
    }

    if divisors.is_empty() {
        Err("number is prime (tested all divisors up to n)".to_string())
    } else {
        Ok(divisors)
    }
}
