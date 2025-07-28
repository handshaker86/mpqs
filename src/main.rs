use num_bigint::{BigInt, Sign};
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};
use std::collections::HashMap;
use std::str::FromStr;

/// Legendre symbol
fn legendre_symbol(a: &BigInt, p: &BigInt) -> i8 {
    let ls = a.modpow(&((p - BigInt::one()) / 2), p);
    if ls == p - BigInt::one() {
        -1
    } else {
        1
    }
}

/// calculate the square root modulo a prime p
fn mod_sqrt(n: &BigInt, p: &BigInt) -> Vec<BigInt> {
    if legendre_symbol(n, p) != 1 {
        return vec![];
    }

    if p.mod_floor(&BigInt::from(4)) == BigInt::from(3) {
        let x = n.modpow(&((p + BigInt::one()) / 4), p);
        return vec![x.clone(), p - x];
    }

    // For other cases, use a simple search
    let mut x = BigInt::from(2);
    while &x < p {
        if x.modpow(&BigInt::from(2), p) == *n {
            return vec![x.clone(), p - x];
        }
        x += BigInt::one();
    }
    vec![]
}

/// calculate the square root modulo p^k
fn mod_sqrt_prime_power(n: &BigInt, p: &BigInt, k: u32) -> Vec<BigInt> {
    if k == 1 {
        return mod_sqrt(n, p);
    }

    let prev_roots = mod_sqrt_prime_power(n, p, k - 1);
    if prev_roots.is_empty() {
        return vec![];
    }

    let p_k_minus_1 = p.pow(k - 1);
    let mut roots = Vec::new();
    let two = BigInt::from(2);

    for r in &prev_roots {
        // Solve (r + x*p^(k-1))^2 ≡ n (mod p^k)
        let val = (n - r.modpow(&two, &(p_k_minus_1.clone() * p))) / &p_k_minus_1;

        if let Some(inv_2r) = (&two * r).modinv(p) {
            let x0 = (val * inv_2r).mod_floor(p);
            roots.push(r + x0 * &p_k_minus_1);
        }
    }
    roots
}

fn find_factor_base(n: &BigInt, b_limit: u64) -> Vec<BigInt> {
    let mut factor_base = vec![BigInt::from(-1)];
    for p_u64 in 2..=b_limit {
        if primal::is_prime(p_u64) {
            let p = BigInt::from(p_u64);
            if legendre_symbol(n, &p) == 1 {
                factor_base.push(p);
            }
        }
    }
    factor_base
}

fn solve_poly_mod_p(a: &BigInt, b: &BigInt, _c: &BigInt, p: &BigInt, n: &BigInt) -> Vec<BigInt> {
    if a.is_multiple_of(p) {
        return vec![];
    }
    if let Some(inv_a) = a.modinv(p) {
        let sqrt_n_mod_p = mod_sqrt(&n.mod_floor(p), p);
        let mut solutions = Vec::new();
        for r in sqrt_n_mod_p {
            let sol = ((r - b) * inv_a.clone()).mod_floor(p);
            solutions.push(sol);
        }
        return solutions;
    }
    vec![]
}

fn trial_division(num: &BigInt, factor_base: &[BigInt]) -> Option<Vec<u8>> {
    let mut num_abs = num.abs();
    let mut factors_exp = vec![0; factor_base.len()];

    if num.sign() == Sign::Minus {
        factors_exp[0] = 1;
    }

    for (i, p) in factor_base.iter().enumerate() {
        if p.is_negative() {
            continue;
        }
        while num_abs.is_multiple_of(p) {
            factors_exp[i] = (factors_exp[i] + 1) % 2;
            num_abs /= p;
        }
    }

    if num_abs.is_one() {
        Some(factors_exp)
    } else {
        None
    }
}

fn build_and_solve_matrix(matrix: &[Vec<u8>]) -> Vec<Vec<u8>> {
    if matrix.is_empty() {
        return vec![];
    }
    let num_relations = matrix.len();
    let num_primes = matrix[0].len();

    let mut m: Vec<Vec<u8>> = (0..num_primes)
        .map(|i| (0..num_relations).map(|j| matrix[j][i]).collect())
        .collect();

    let mut pivot_cols = Vec::new();
    let mut lead = 0;
    for r in 0..num_primes {
        if lead >= num_relations {
            break;
        }
        let mut i = r;
        while m[i][lead] == 0 {
            i += 1;
            if i == num_primes {
                i = r;
                lead += 1;
                if lead == num_relations {
                    break;
                }
            }
        }
        if lead >= num_relations {
            break;
        }
        m.swap(i, r);
        for i in 0..num_primes {
            if i != r && m[i][lead] == 1 {
                for j in lead..num_relations {
                    m[i][j] ^= m[r][j];
                }
            }
        }
        pivot_cols.push(lead);
        lead += 1;
    }

    let mut dependencies = Vec::new();
    let free_cols: Vec<usize> = (0..num_relations)
        .filter(|c| !pivot_cols.contains(c))
        .collect();

    for free_col in free_cols {
        let mut dependency = vec![0; num_relations];
        dependency[free_col] = 1;
        for (i, &pivot_col) in pivot_cols.iter().enumerate() {
            if m[i][free_col] == 1 {
                dependency[pivot_col] = 1;
            }
        }
        dependencies.push(dependency);
    }
    dependencies
}

fn mpqs(n: &BigInt, m_val: i64, fb_limit: u64) -> Option<(BigInt, BigInt)> {
    println!("Starting MPQS factorization for: {}", n);

    let factor_base = find_factor_base(n, fb_limit);
    let fb_size = factor_base.len();
    println!("Sieving interval: [{}, {}]", -m_val, m_val);
    println!("Factor base (size {}): {:?}\n", fb_size, factor_base);

    let mut smooth_relations = Vec::new();
    let mut smooth_x_vals = Vec::new();
    let mut poly_params_list = Vec::new();
    let needed_relations = fb_size + 5;

    let q_candidate_float =
        (n * BigInt::from(2)).to_f64().unwrap().powf(0.25) / (m_val as f64).powf(0.5);
    let mut q_candidate = q_candidate_float.ceil() as u64;

    let a_prime = loop {
        if primal::is_prime(q_candidate) && legendre_symbol(n, &BigInt::from(q_candidate)) == 1 {
            break BigInt::from(q_candidate);
        }
        q_candidate += 1;
    };
    let a = &a_prime * &a_prime;
    println!("Selected polynomial coefficient a = {}^2 = {}", a_prime, a);

    let b_sol = mod_sqrt_prime_power(&n.mod_floor(&a), &a_prime, 2);
    let b = if let Some(b_val) = b_sol.get(0) {
        b_val.clone()
    } else {
        println!("Could not find b for selected a={}", a);
        return None;
    };
    let c = (b.pow(2) - n) / &a;
    println!("Polynomial Q(x) = {}x^2 + 2*{}x + {}\n", a, b, c);

    let ln_2 = std::f64::consts::LN_2;
    let mut sieving_array: HashMap<i64, f64> = HashMap::new();
    for x_i64 in -m_val..=m_val {
        let x = BigInt::from(x_i64);
        let q_x = &a * &x * &x + BigInt::from(2) * &b * &x + &c;
        if !q_x.is_zero() {
            sieving_array.insert(x_i64, q_x.abs().bits() as f64 * ln_2);
        }
    }

    let threshold = if fb_size > 1 {
        (factor_base.last().unwrap().bits() as f64 * ln_2) * 2.0
    } else {
        2.0
    };

    for p in &factor_base {
        if p.is_negative() {
            continue;
        }
        if let Some(log_p) = p.to_f64() {
            if log_p > 0.0 {
                // Ensure log_p > 0 to avoid issues
                let solutions = solve_poly_mod_p(&a, &b, &c, p, n);
                for s_big in solutions {
                    if let Some(s) = s_big.to_i64() {
                        if let Some(p_i64) = p.to_i64() {
                            if p_i64 == 0 {
                                continue;
                            }
                            // Sieve positive side
                            for x_val in (s..=m_val).step_by(p_i64 as usize) {
                                if let Some(val) = sieving_array.get_mut(&x_val) {
                                    *val -= log_p.ln();
                                }
                            }
                            // Sieve negative side
                            let mut x_val = s - p_i64;
                            while x_val >= -m_val {
                                if let Some(val) = sieving_array.get_mut(&x_val) {
                                    *val -= log_p.ln();
                                }
                                x_val -= p_i64;
                            }
                        }
                    }
                }
            }
        }
    }

    println!("Sieving complete, collecting smooth relations...");
    for x_i64 in -m_val..=m_val {
        if smooth_relations.len() >= needed_relations {
            break;
        }
        if let Some(rem_log) = sieving_array.get(&x_i64) {
            if *rem_log < threshold {
                let x = BigInt::from(x_i64);
                let q_x = &a * x.pow(2) + BigInt::from(2) * &b * &x + &c;
                if q_x.is_zero() {
                    continue;
                }
                if let Some(exp_vector) = trial_division(&q_x, &factor_base) {
                    smooth_relations.push(exp_vector);
                    smooth_x_vals.push(x);
                    poly_params_list.push((a.clone(), b.clone()));
                }
            }
        }
    }

    if smooth_relations.len() < fb_size + 1 {
        println!("Failed to find enough smooth relations.");
        return None;
    }
    println!("Found {} smooth relations.\n", smooth_relations.len());

    let dependencies = build_and_solve_matrix(&smooth_relations);
    if dependencies.is_empty() {
        println!("No linear dependencies found.");
        return None;
    }
    println!("Found {} linear dependencies.\n", dependencies.len());

    for (i, dep) in dependencies.iter().enumerate() {
        println!("--- Trying dependency {} ---", i + 1);
        let mut cap_x = BigInt::one();
        let mut y_sq_factors: HashMap<BigInt, u32> = HashMap::new();
        let mut k: u32 = 0;

        for (rel_idx, &bit) in dep.iter().enumerate() {
            if bit == 1 {
                k += 1;
                let x_k = &smooth_x_vals[rel_idx];
                let (a_k, b_k) = &poly_params_list[rel_idx];
                cap_x = (cap_x * (a_k * x_k + b_k)).mod_floor(n);

                let q_x_k = a_k * x_k.pow(2) + BigInt::from(2) * b_k * x_k + (b_k.pow(2) - n) / a_k;

                let mut q_x_abs = q_x_k.abs();
                if q_x_k.sign() == Sign::Minus {
                    *y_sq_factors.entry(BigInt::from(-1)).or_insert(0) += 1;
                }
                for p_fb in &factor_base {
                    if p_fb.is_negative() {
                        continue;
                    }
                    while q_x_abs.is_multiple_of(p_fb) {
                        *y_sq_factors.entry(p_fb.clone()).or_insert(0) += 1;
                        q_x_abs /= p_fb;
                    }
                }
            }
        }

        let mut y_prime = BigInt::one();
        for (p_fb, exp) in &y_sq_factors {
            y_prime = (y_prime * p_fb.modpow(&(BigInt::from(*exp / 2)), n)).mod_floor(n);
        }

        let sqrt_a_k = a_prime.modpow(&BigInt::from(k), n);
        let cap_y = (y_prime * sqrt_a_k).mod_floor(n);

        let factor = (cap_x.clone() - &cap_y).abs().gcd(n);
        if !factor.is_one() && &factor != n {
            println!("Factor found successfully!");
            return Some((factor.clone(), n / factor));
        }

        let factor = (cap_x + cap_y).abs().gcd(n);
        if !factor.is_one() && &factor != n {
            println!("Factor found successfully!");
            return Some((factor.clone(), n / factor));
        }
        println!("This dependency only yielded a trivial factor, trying next...");
    }

    println!("\nAll dependencies resulted in trivial factors. Algorithm failed.");
    None
}

fn main() {
    let number_to_factor = BigInt::from_str("10666351").unwrap(); // 1237 * 8623

    let factors = mpqs(&number_to_factor, 30000, 130);

    if let Some((p, q)) = factors {
        println!(
            "\nFactorization successful: {} = {} * {}",
            number_to_factor, p, q
        );
    }
}
