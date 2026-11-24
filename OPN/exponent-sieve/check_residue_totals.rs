// Directly sum rounded reciprocal weights in the seven reduced classes other than 1.
const START: u64 = 100_000_000;
const STOP: u64 = 8_000_000_000;
const PRECISION: u64 = 1_000_000_000_000_000_000;

fn main() {
    let limit = STOP.isqrt() as usize;
    let mut prime = vec![true; limit + 1];
    for p in 2..=limit.isqrt() {
        if prime[p] {
            for multiple in (p * p..=limit).step_by(p) {
                prime[multiple] = false;
            }
        }
    }
    let mut total = 0;
    let mut weight = 0u64;
    for residue in [7, 11, 13, 17, 19, 23, 29] {
        // For p > 5, determine where 30k+residue is zero modulo p.
        let mut progressions = Vec::new();
        for p in (7..=limit).filter(|&p| prime[p]) {
            let p = p as u64;
            let k = (0..30).find(|&j| (j * p + 1) % 30 == 0).unwrap();
            let inverse = (k * p + 1) / 30;
            let offset = (p - (residue * inverse) % p) % p;
            progressions.push((p, offset));
        }
        let mut begin = (START - residue) / 30 + 1;
        let end = (STOP - residue) / 30 + 1;
        while begin < end {
            let next = (begin + 32768).min(end);
            let mut allowed = vec![true; (next - begin) as usize];
            for &(p, root) in &progressions {
                if p * p > 30 * (next - 1) + residue {
                    break;
                }
                let first = ((root + p - begin % p) % p) as usize;
                for j in (first..allowed.len()).step_by(p as usize) {
                    allowed[j] = false;
                }
            }
            for (j, keep) in allowed.into_iter().enumerate() {
                if keep {
                    let denominator = 30 * (begin + j as u64) + residue - 1;
                    weight += PRECISION.div_ceil(denominator);
                    total += 1;
                }
            }
            begin = next;
        }
    }
    println!("admissible primes {total}, reciprocal upper numerator {weight}, scale {PRECISION}");
}
