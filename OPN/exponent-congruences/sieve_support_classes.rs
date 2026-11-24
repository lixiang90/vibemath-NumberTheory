// Independent residue-class sieve for the four finite prefixes in the paper.
// Enumerate only q = 3,000,000 + residue + 6*i. Obtain the roots of Phi
// modulo p from a square root of -3, using Tonelli-Shanks.
use std::{env, thread};

const P: u64 = 3_000_000;
const SCALE: u64 = 1_000_000_000_000_000_000;
const SEGMENT: u64 = 1 << 20;
const SMALL: [u64; 6] = [5, 7, 11, 13, 17, 19];

fn power(mut base: u64, mut exponent: u64, p: u64) -> u64 {
    let mut result = 1;
    while exponent > 0 {
        if exponent % 2 == 1 {
            result = result * base % p;
        }
        base = base * base % p;
        exponent /= 2;
    }
    result
}

fn square_root(a: u64, p: u64) -> u64 {
    if p % 4 == 3 {
        return power(a, (p + 1) / 4, p);
    }
    let mut m = (p - 1).trailing_zeros();
    let odd_part = (p - 1) >> m;
    let nonresidue = (2..p).find(|&g| power(g, (p - 1) / 2, p) == p - 1).unwrap();
    let mut c = power(nonresidue, odd_part, p);
    let mut t = power(a, odd_part, p);
    let mut root = power(a, (odd_part + 1) / 2, p);
    while t != 1 {
        let mut v = t;
        let mut i = 0;
        while v != 1 {
            v = v * v % p;
            i += 1;
        }
        assert!(i < m);
        let b = power(c, 1 << (m - i - 1), p);
        c = b * b % p;
        root = root * b % p;
        t = t * c % p;
        m = i;
    }
    root
}

fn progressions(first: u64, upper: u64) -> Vec<(u64, u64)> {
    let mut composite = vec![false; P as usize + 1];
    let mut steps = Vec::new();
    for p in 2..=P {
        if composite[p as usize] {
            continue;
        }
        for k in (p*p..=P).step_by(p as usize) {
            composite[k as usize] = true;
        }
        if p <= 19 {
            continue;
        }
        let inverse_six = power(6, p-2, p);
        let offset = |r: u64| (p + r - first % p) % p * inverse_six % p;
        if p*p <= upper {
            steps.push((p, offset(0)));
        }
        if p % 3 == 1 {
            let root = square_root(p-3, p);
            for sign in [root, p-root] {
                let r = (p-1 + sign) * ((p+1)/2) % p;
                assert_eq!((r*r + r + 1) % p, 0);
                steps.push((p, offset(r)));
            }
        }
    }
    steps
}

fn main() {
    let args: Vec<u64> = env::args().skip(1).map(|s| s.parse().unwrap()).collect();
    assert!(args.len() == 3, "expected upper residue_mod_6 threads");
    let (upper, residue, threads) = (args[0], args[1], args[2]);
    assert!(matches!(residue, 1 | 5) && threads > 0 && threads <= usize::MAX as u64);
    let first = P + residue;
    assert!(first <= upper && upper <= P*P);
    let total = (upper - first) / 6 + 1;
    let period: u64 = SMALL.iter().product();
    let wheel: Vec<u8> = (0..period).map(|i| {
        let q = first + 6*i;
        SMALL.iter().all(|&p| {
            let r = q % p;
            r != 0 && (r*r+r+1) % p != 0
        }) as u8
    }).collect();
    let steps = progressions(first, upper);

    let (count, weight) = thread::scope(|scope| {
        let workers: Vec<_> = (0..threads).map(|id| {
            let steps = &steps;
            let wheel = &wheel;
            scope.spawn(move || {
                let mut count = 0u64;
                let mut weight = 0u128;
                for block in (id..total.div_ceil(SEGMENT)).step_by(threads as usize) {
                    let start = block * SEGMENT;
                    let size = SEGMENT.min(total - start) as usize;
                    let at = (start % period) as usize;
                    let split = size.min(wheel.len() - at);
                    let mut keep = Vec::with_capacity(size);
                    keep.extend_from_slice(&wheel[at..at+split]);
                    keep.extend_from_slice(&wheel[..size-split]);
                    for &(p, offset) in steps {
                        let pos = ((offset + p - start % p) % p) as usize;
                        for i in (pos..size).step_by(p as usize) {
                            keep[i] = 0;
                        }
                    }
                    for (i, flag) in keep.into_iter().enumerate() {
                        if flag != 0 {
                            let q = first + 6*(start+i as u64);
                            count += 1;
                            weight += (1+(SCALE-1)/(q-1)) as u128;
                        }
                    }
                }
                (count, weight)
            })
        }).collect();
        workers.into_iter().map(|w| w.join().unwrap())
            .fold((0, 0), |(a, b), (c, d)| (a+c, b+d))
    });
    println!("{upper} {residue} {count} {weight}");
}
