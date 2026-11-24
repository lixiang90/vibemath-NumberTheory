// Independently compute D = psi * chi * conjugate(chi) by divisor sums.
// The Python certificate instead constructs D from local prime powers.
use std::env;

const SCALE: i64 = 100_000_000_000_000;
const PSI: [i32; 5] = [0, 1, -1, -1, 1];

fn coefficients(limit: usize) -> Vec<i32> {
    let mut divisors = vec![0; limit + 1];
    for d in 1..=limit {
        let value = PSI[d % 5];
        if value != 0 {
            for n in (d..=limit).step_by(d) {
                divisors[n] += value;
            }
        }
    }
    let mut result = vec![0i32; limit + 1];
    for (a, &sum) in divisors.iter().enumerate().skip(1) {
        // chi * conjugate(chi) is conjugate(chi(a)) times (1 * psi)(a)
        // on units. At a nonsquare unit the divisor sum vanishes.
        let value = match a % 5 {
            1 => sum,
            4 => -sum,
            2 | 3 => {
                assert_eq!(sum, 0);
                0
            }
            _ => 0,
        };
        if value != 0 {
            for b in 1..=limit / a {
                result[a * b] += value * PSI[b % 5];
            }
        }
    }
    result
}

fn main() {
    let limit = env::args().nth(1).map(|n| n.parse().expect("integer endpoint"))
        .unwrap_or(10_000_000usize);
    assert!((40..=10_000_000).contains(&limit));
    let values = coefficients(limit);
    let mut minima = [(0usize, i64::MAX); 2];
    let mut total = 0i64;
    for (n, &d) in values.iter().enumerate().skip(1) {
        total += (SCALE * i64::from(d)).div_euclid(n as i64);
        let block = usize::from(n >= 40);
        if total < minima[block].1 {
            minima[block] = (n, total);
        }
    }
    println!("checked {limit} terms, minima {minima:?}");
    if limit == 10_000_000 {
        assert_eq!(minima, [(3, 16_666_666_666_666), (143, 30_434_989_347_753)]);
    }
}
