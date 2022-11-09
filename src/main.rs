use std::{env::args, error::Error, iter::repeat, str::FromStr};

use factor_building::find_factors_from_pivots;
use itertools::Itertools;
use num_bigint::BigUint;
use sieve::SmoothNumber;

use crate::{
    factor_building::{find_factor_exhaustive, find_factor_simple, find_factors_random},
    number_type::NumberType,
    numbers::{build_factor_base, small_eratosphenes},
    sieve::{block_division_sieve, compute_b_limit, test_division_sieve},
    solver::{produce_solution, CongruenceSystem},
};

mod number_type;
mod numbers;

mod factor_building;
mod sieve;
mod solver;

fn gaussian_multistage(
    n: &NumberType,
    table: Vec<SmoothNumber<NumberType>>,
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

    find_factor_simple(n, &table, &solution)
        .or_else(|| find_factors_random(n, &table, &solution))
        .or_else(|| find_factor_exhaustive(n, &table, &solution))
}

fn pivot_search(
    n: &NumberType,
    table: Vec<SmoothNumber<NumberType>>,
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

fn run_factor(n: &NumberType, prime_bound: usize) -> Option<(BigUint, BigUint)> {
    let primes = small_eratosphenes(prime_bound);

    println!("primes until bound: {}", primes.len());

    let factor_base = build_factor_base(primes, n);

    println!("built factor base of size {}", factor_base.len(),);

    if factor_base.len() < 2 {
        println!("this is too small");
        return None;
    }

    #[cfg(feature = "verbose")]
    println!("factor base: {:?}", factor_base);

    let mut ratio = 1.05;

    const NUM_ATTEMPTS: usize = 5;

    let mut cont = None;

    let mut table = vec![];

    for _ in 0..NUM_ATTEMPTS {
        let sieving_limit = usize::max(
            (factor_base.len() as f64 * ratio).ceil() as usize,
            factor_base.len() + 5,
        );

        println!("need about {sieving_limit} numbers");

        let (mut additional_table, new_cont) = test_division_sieve(
            n,
            &factor_base,
            |_, _| sieving_limit.saturating_sub(table.len()),
            cont,
        )
        .unwrap();

        table.append(&mut additional_table);

        cont = Some(new_cont);

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
