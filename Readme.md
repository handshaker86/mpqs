# MPQS Factoring Algorithm in Rust

## Overview

This project is a Rust implementation of the Multiple Polynomial Quadratic Sieve (MPQS) algorithm, an integer factorization method for large composite numbers. The implementation is based on the principles outlined in the paper "The Quadratic Sieve Factoring Algorithm" . [Here is the paper](http://www.cs.virginia.edu/crab/QFS_Simple.pdf)

The algorithm's goal is to find a non-trivial factor of a number `n` by finding integers `x` and `y` that satisfy the congruence $x^2 \equiv y^2 \pmod{n}$.

## Implementation Details


### 1. Setup and Polynomial Selection
* **Function:** `mpqs()`, `find_factor_base()`, `mod_sqrt_prime_power()`
* **Description:** The main `mpqs()` function begins by setting parameters. It calls `find_factor_base()` to generate a set of small prime numbers `p` (the factor base) where `n` is a quadratic residue modulo `p`. It then selects a polynomial `Q(x) = ax² + 2bx + c` by choosing `a` as a perfect square and solving for `b` such that `b² ≡ n (mod a)`. This modular square root for a prime-power modulus is handled by `mod_sqrt_prime_power()`.

### 2. Sieving for Smooth Relations
* **Function:** `mpqs()` , `solve_poly_mod_p()`, `trial_division()`
* **Description:** The algorithm sieves over a given interval `[-M, M]` to find "smooth" numbers, which are numbers that factor completely over the established factor base. For each prime `p` in the base, `solve_poly_mod_p()` finds the roots of `Q(x) ≡ 0 (mod p)` to identify arithmetic progressions where `Q(x)` is divisible by `p`. After an efficient logarithmic sieve identifies likely candidates, `trial_division()` confirms if they are truly smooth.

### 3. Linear Algebra
* **Function:** `build_and_solve_matrix()`
* **Description:** For each smooth relation found, the exponents of its prime factors are collected into a binary vector (modulo 2). These vectors form a matrix where `build_and_solve_matrix()` performs Gaussian elimination. The goal is to find linear dependencies, which is a set of rows that sum to the zero vector and corresponds to products that are perfect squares.

### 4. Factor Extraction
* **Function(s):** `mpqs()` 
* **Description:** Each dependency we have found ensures that the product of the corresponding `Q(xᵢ)` values is a perfect square, and let's call it `Y²`. This leads to the final congruence `X² ≡ Y² (mod n)`. The program then computes `gcd(X - Y, n)` to find a non-trivial factor. Since any dependency has at least a 50% chance of success, the program iterates through all found dependencies until a non-triival factor is discovered. 

### Example Usage
In the main function, we provide a simple example to factorize 10666351. The final answer should be `10333351 = 2521 * 4231`. Feel free to try other big composite numbers! A kind note is that you can try larger `m_val` and `fb_limit` when you fail to factorize some very large composite numbers. `m_val` represents the boundary value of sieving interval `[-M,M]` and `fb_limit` represents the upper limit for the prime base.