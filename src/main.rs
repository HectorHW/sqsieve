use std::{env::args, error::Error, iter::repeat, str::FromStr};

use crypto_bigint::prelude::*;
use itertools::Itertools;
use num_bigint::BigUint;

use crate::{
    factor_building::find_factor,
    number_type::{NumberOps, NumberType},
    numbers::small_eratosphenes,
    sieve::{compute_b_limit, sieve_for_smoothies},
    solver::{produce_solution, CongruenceSystem},
};

mod number_type;
mod numbers;

mod factor_building;
mod sieve;
mod solver;

fn main() -> Result<(), Box<dyn Error>> {
    //let n = BigUint::from_str("1577271624417732056618338337651").unwrap();

    let args = args().collect::<Vec<String>>();

    if args.len() != 2 {
        return Err("please provide number as single argument".into());
    }

    let n = BigUint::from_str(&args[1]).unwrap();

    let bytes = n.to_bytes_be();

    if bytes.len() > 64 {
        return Err("number is too big".into());
    }

    let bytes = repeat(0u8)
        .take(64 - bytes.len())
        .chain(bytes.into_iter())
        .collect_vec();

    let n = NumberType::from_be_slice(&bytes);

    println!("n: {}", n.to_varsize());

    //let prime_bound = compute_b_limit(&n);
    let prime_bound = compute_b_limit(&n).min(10usize.pow(3));

    println!("prime bound is {prime_bound}");

    let primes = small_eratosphenes(prime_bound);

    let sieving_limit = compute_b_limit(&n);

    println!("need about {sieving_limit} numbers");

    let table = sieve_for_smoothies(&n, &primes, |_, _| sieving_limit)?;

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

    println!("built solution dependencies, searching for factors");

    let factorization = find_factor(&n, &table, &solution);

    println!("factorization: {factorization:?}");

    if let Some((a, b)) = factorization {
        let prod = a.clone() * b.clone();
        println!("{} * {} = {}", a, b, prod);
    }

    Ok(())
}
