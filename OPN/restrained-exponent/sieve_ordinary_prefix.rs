// A prime s above cutoff can be ordinary only when s^2+s+1 has no prime
// factor from 5 through cutoff. Units below are 10^-18 in the reciprocal sum.
use std::env;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::thread;

const SCALE: u64 = 1_000_000_000_000_000_000;
const SEGMENT: usize = 1 << 21;

fn primes(limit: usize) -> Vec<usize> {
    let mut is_prime = vec![true; limit + 1];
    for p in 2..=limit.isqrt() {
        if is_prime[p] {
            for j in (p * p..=limit).step_by(p) {
                is_prime[j] = false;
            }
        }
    }
    (2..=limit).filter(|&p| is_prime[p]).collect()
}

fn modular_power(mut a: u64, mut n: usize, p: u64) -> usize {
    let mut value = 1;
    while n != 0 {
        if n % 2 == 1 {
            value = value * a % p;
        }
        a = a * a % p;
        n /= 2;
    }
    value as usize
}

struct Sieve {
    first: usize,
    total: usize,
    wheel: Vec<u8>,
    progressions: Vec<(usize, usize)>,
}

impl Sieve {
    fn new(cutoff: usize, upper: usize) -> Self {
        assert!(cutoff >= 19 && cutoff < upper && upper <= cutoff * cutoff);
        let first = (cutoff + 1) | 1;
        let total = (upper + 2 - first) / 2;
        let mut wheel = vec![1; 3 * 5 * 7 * 11 * 13 * 17 * 19];
        let mut progressions = Vec::new();
        for p in primes(cutoff).into_iter().filter(|&p| p != 2) {
            let mut forbid = |residue: usize| {
                // The array index i represents the odd number first + 2*i.
                let start = ((residue + p - first % p) % p) * p.div_ceil(2) % p;
                if p <= 19 {
                    for i in (start..wheel.len()).step_by(p) {
                        wheel[i] = 0;
                    }
                } else {
                    progressions.push((p, start));
                }
            };
            if p <= upper.isqrt() {
                forbid(0);
            }
            if p % 3 == 1 {
                let root = (2..p)
                    .map(|g| modular_power(g as u64, (p - 1) / 3, p as u64))
                    .find(|&root| root != 1)
                    .expect("missing nontrivial cube root");
                forbid(root);
                forbid(p - 1 - root);
            }
        }
        Self { first, total, wheel, progressions }
    }

    fn segment(&self, start: usize) -> (u64, u128) {
        let len = SEGMENT.min(self.total - start);
        let at = start % self.wheel.len();
        let split = len.min(self.wheel.len() - at);
        let mut keep = Vec::with_capacity(len);
        keep.extend_from_slice(&self.wheel[at..at + split]);
        keep.extend_from_slice(&self.wheel[..len - split]);
        for &(p, first) in &self.progressions {
            let offset = (first + p - start % p) % p;
            for i in (offset..len).step_by(p) {
                keep[i] = 0;
            }
        }
        let mut count = 0;
        let mut reciprocal = 0;
        for (i, &yes) in keep.iter().enumerate() {
            if yes != 0 {
                let prime = (self.first + 2 * (start + i)) as u64;
                count += 1;
                reciprocal += SCALE.div_ceil(prime - 1) as u128;
            }
        }
        (count, reciprocal)
    }
}

fn main() {
    let args: Vec<_> = env::args().skip(1).collect();
    assert!(args.len() == 3, "expected cutoff upper threads");
    let cutoff: usize = args[0].parse().unwrap();
    let upper: usize = args[1].parse().unwrap();
    let threads: usize = args[2].parse().unwrap();
    assert!(threads >= 1);
    let sieve = Sieve::new(cutoff, upper);
    let next = AtomicUsize::new(0);
    let totals = thread::scope(|scope| {
        let handles: Vec<_> = (0..threads).map(|_| scope.spawn(|| {
            let mut local = (0u64, 0u128);
            loop {
                let start = next.fetch_add(SEGMENT, Ordering::Relaxed);
                if start >= sieve.total {
                    return local;
                }
                let (count, sum) = sieve.segment(start);
                local.0 += count;
                local.1 += sum;
            }
        })).collect();
        handles.into_iter().map(|h| h.join().unwrap())
            .fold((0u64, 0u128), |(n, s), (k, t)| (n + k, s + t))
    });
    println!("{{\"cutoff\":{cutoff},\"upper\":{upper},\"survivors\":{},\"scaled_reciprocal_upper\":{},\"scale\":{SCALE}}}", totals.0, totals.1);
}

#[cfg(test)]
mod tests {
    use super::*;

    fn trial_prime(n: usize) -> bool {
        n >= 2 && (2..=n.isqrt()).all(|d| n % d != 0)
    }

    fn direct(cutoff: usize, lower: usize, upper: usize) -> (u64, u128) {
        let divisors: Vec<_> = (5..=cutoff).filter(|&p| trial_prime(p)).collect();
        let mut answer = (0u64, 0u128);
        for s in lower..=upper {
            if trial_prime(s) && divisors.iter().all(|p| (s*s+s+1) % p != 0) {
                answer.0 += 1;
                answer.1 += SCALE.div_ceil((s-1) as u64) as u128;
            }
        }
        answer
    }

    #[test]
    fn segmented_sum_matches_direct_trial_division() {
        for (cutoff, upper) in [(19, 20), (20, 31), (31, 961), (97, 7938), (100, 10000)] {
            let sieve = Sieve::new(cutoff, upper);
            assert_eq!(sieve.segment(0), direct(cutoff, cutoff+1, upper));
        }
    }

    #[test]
    fn segment_crossing_a_wheel_boundary() {
        let cutoff = 4001;
        let start = 3*5*7*11*13*17*19 - 64;
        let first = 4003 + 2*start;
        let upper = first + 254;
        let sieve = Sieve::new(cutoff, upper);
        assert_eq!(sieve.segment(start), direct(cutoff, first, upper));
    }
}
