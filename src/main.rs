use std::{env::args, error::Error, iter::repeat, str::FromStr};

use crypto_bigint::prelude::*;
use itertools::Itertools;
use num_bigint::BigUint;

use crate::{
    factor_building::{find_factor_exhaustive, find_factor_simple, find_factors_random},
    number_type::{NumberOps, NumberType},
    numbers::{build_factor_base, small_eratosphenes},
    sieve::{compute_b_limit, sieve_for_smoothies},
    solver::{produce_solution, CongruenceSystem},
};

mod number_type;
mod numbers;

mod factor_building;
mod sieve;
mod solver;

fn run_factor(n: &NumberType, prime_bound: usize) -> Option<(BigUint, BigUint)> {
    let primes = small_eratosphenes(prime_bound);

    println!("primes until bound: {}", primes.len());

    let factor_base = build_factor_base(primes, n);

    println!("built factor base of size {}", factor_base.len(),);

    if factor_base.len() < 2 {
        println!("this is too small");
        return None;
    }

    let sieving_limit = factor_base.len() + 33;

    println!("need about {sieving_limit} numbers");

    let table = sieve_for_smoothies(n, &factor_base, |_, _| sieving_limit).unwrap();

    println!("done collecting, building solution");

    #[cfg(feature = "verbose")]
    println!("{:?}", table);

    #[cfg(feature = "verbose")]
    for (i, item) in table.iter().enumerate() {
        println!("{i} -> {}", item.number);
    }

    let row_labels = (0..table.len()).collect_vec();
    let rows = table.iter().map(|row| row.divisors.clone()).collect_vec();

    let system = CongruenceSystem::new(&rows, row_labels);

    #[cfg(feature = "verbose")]
    println!("{}", system.print_as_dense());

    #[cfg(feature = "verbose")]
    println!("---------");

    let mut system = system.transpose();
    #[cfg(feature = "verbose")]
    println!("{}", system.print_as_dense());

    #[cfg(feature = "verbose")]
    println!("---------");

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

    find_factor_simple(n, &table, &solution)
        .or_else(|| find_factors_random(n, &table, &solution))
        .or_else(|| find_factor_exhaustive(n, &table, &solution))
}

fn main() -> Result<(), Box<dyn Error>> {
    //let n = BigUint::from_str("1577271624417732056618338337651").unwrap();

    let args = args().collect::<Vec<String>>();

    if args.len() != 2 {
        return Err("please provide number as single argument".into());
    }

    let n = BigUint::from_str(&args[1]).unwrap();

    println!("n: {}", n);

    println!("base 10 digits: {}", n.to_string().len());

    println!("bit size: {}", n.bits());

    let bytes = n.to_bytes_be();

    if bytes.len() > 64 {
        return Err("number is too big".into());
    }

    let bytes = repeat(0u8)
        .take(64 - bytes.len())
        .chain(bytes.into_iter())
        .collect_vec();

    let n = NumberType::from_be_slice(&bytes);

    //let prime_bound = compute_b_limit(&n);
    //let prime_bound = compute_b_limit(&n).min(10usize.pow(3));

    let limit = compute_b_limit(&n);

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
        return Ok(());
    }

    Ok(())
}
