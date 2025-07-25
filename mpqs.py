import math
from collections import defaultdict


# --- 1. 数学辅助函数 (包含新的mod_sqrt_prime_power) ---
def is_prime(num):
    if num < 2:
        return False
    for i in range(2, int(math.sqrt(num)) + 1):
        if num % i == 0:
            return False
    return True


def gcd(a, b):
    while b:
        a, b = b, a % b
    return a


def legendre_symbol(a, p):
    ls = pow(a, (p - 1) // 2, p)
    return -1 if ls == p - 1 else ls


def mod_sqrt(n, p):
    """计算模素数p下的平方根 (Tonelli-Shanks简化版)"""
    if legendre_symbol(n, p) != 1:
        return []
    if p % 4 == 3:
        x = pow(n, (p + 1) // 4, p)
        return [x, p - x]
    for x in range(2, p):
        if (x * x) % p == n:
            return [x, p - x]
    return []


def mod_sqrt_prime_power(n, p, k):
    """使用亨泽尔引理计算模 p^k 下的平方根"""
    if k == 1:
        return mod_sqrt(n, p)

    # 递归计算模 p^(k-1) 下的解
    prev_roots = mod_sqrt_prime_power(n, p, k - 1)
    if not prev_roots:
        return []

    p_k_minus_1 = pow(p, k - 1)
    p_k = pow(p, k)

    roots = []
    for r in prev_roots:
        # 求解 (r + x*p^(k-1))^2 ≡ n (mod p^k)
        # 2*r*x*p^(k-1) ≡ n - r^2 (mod p^k)
        # 2*r*x ≡ (n - r^2)/p^(k-1) (mod p)
        val = (n - r * r) // p_k_minus_1

        # 求解 2*r*k ≡ val (mod p)
        inv_2r = pow(2 * r, -1, p)
        x0 = (val * inv_2r) % p

        # 新的解是 r + x0 * p^(k-1)
        roots.append(r + x0 * p_k_minus_1)

    return roots


# --- 2. MPQS 核心实现 (使用最终修正逻辑) ---
def find_factor_base(n, b_limit):
    factor_base = [-1]
    for p in range(2, b_limit + 1):
        if is_prime(p) and legendre_symbol(n, p) == 1:
            factor_base.append(p)
    return factor_base


def solve_poly_mod_p(a, b, c, p, n):
    if a % p == 0:
        if (2 * b) % p != 0:
            inv_2b = pow(2 * b, -1, p)
            return [(-c * inv_2b) % p]
        return []
    inv_a = pow(a, -1, p)
    sqrt_n_mod_p = mod_sqrt(n, p)
    solutions = []
    for r in sqrt_n_mod_p:
        sol = ((r - b) * inv_a) % p
        solutions.append(sol)
    return solutions


def trial_division(num, factor_base):
    num_abs = abs(num)
    factors_exp = [0] * len(factor_base)
    if num < 0:
        factors_exp[0] = 1
    for i, p in enumerate(factor_base):
        if p == -1:
            continue
        while num_abs % p == 0:
            factors_exp[i] = (factors_exp[i] + 1) % 2
            num_abs //= p
    return factors_exp if num_abs == 1 else None


def build_and_solve_matrix(matrix):
    if not matrix:
        return []
    num_relations = len(matrix)
    num_primes = len(matrix[0])
    M = [[matrix[j][i] for j in range(num_relations)] for i in range(num_primes)]
    pivot_cols, lead = [], 0
    for r in range(num_primes):
        if lead >= num_relations:
            break
        i = r
        while M[i][lead] == 0:
            i += 1
            if i == num_primes:
                i = r
                lead += 1
                if lead == num_relations:
                    break
        if lead >= num_relations:
            break
        M[i], M[r] = M[r], M[i]
        for i in range(num_primes):
            if i != r and M[i][lead] == 1:
                for j in range(lead, num_relations):
                    M[i][j] ^= M[r][j]
        pivot_cols.append(lead)
        lead += 1
    dependencies = []
    free_cols = [i for i in range(num_relations) if i not in pivot_cols]
    for free_col in free_cols:
        dependency = [0] * num_relations
        dependency[free_col] = 1
        for i, pivot_col in enumerate(pivot_cols):
            if M[i][free_col] == 1:
                dependency[pivot_col] = 1
        dependencies.append(dependency)
    return dependencies


def mpqs(n, m_val=10000, fb_limit=100):
    print(f"开始使用MPQS分解: {n}")
    M = m_val
    factor_base = find_factor_base(n, fb_limit)
    fb_size = len(factor_base)
    print(f"筛选区间: [{-M}, {M}]")
    print(f"因子基 (大小 {fb_size}): {factor_base}\n")

    smooth_relations = []
    smooth_x_vals = []
    poly_params = []
    needed_relations = fb_size + 5

    # --- ✨ 修正点 1: 选择 a 为一个平方数 ---
    q_candidate = int((n * 2) ** 0.25 / M**0.5)
    a_prime = 0
    while True:
        q_candidate += 1
        if is_prime(q_candidate) and legendre_symbol(n, q_candidate) == 1:
            a_prime = q_candidate
            break

    a = a_prime * a_prime
    print(f"选择多项式系数 a = {a_prime}^2 = {a}")

    # --- ✨ 修正点 2: 使用新的求解器计算 b ---
    b_sol = mod_sqrt_prime_power(n % a, a_prime, 2)
    if not b_sol:
        print(f"无法为选定的a={a}找到b")
        return None
    b = b_sol[0]

    c = (b * b - n) // a
    print(f"多项式 Q(x) = {a}x^2 + 2*{b}x + {c}\n")

    sieving_array = {
        x: math.log(abs(a * x**2 + 2 * b * x + c) or 1) for x in range(-M, M + 1)
    }
    threshold = math.log(factor_base[-1] if len(factor_base) > 1 else 2)

    for p in factor_base:
        if p == -1:
            continue
        log_p = math.log(p)
        solutions = solve_poly_mod_p(a, b, c, p, n)
        for s in solutions:
            for x_val in range(s, M + 1, p):
                if x_val in sieving_array:
                    sieving_array[x_val] -= log_p
            for x_val in range(s - p, -M - 1, -p):
                if x_val in sieving_array:
                    sieving_array[x_val] -= log_p

    print("筛选完成，开始收集平滑关系...")
    for x, rem_log in sieving_array.items():
        if len(smooth_relations) >= needed_relations:
            break
        if rem_log < threshold * 2:
            q_x = a * x**2 + 2 * b * x + c
            if q_x == 0:
                continue
            exp_vector = trial_division(q_x, factor_base)
            if exp_vector is not None:
                smooth_relations.append(exp_vector)
                smooth_x_vals.append(x)
                poly_params.append({"a": a, "b": b})

    if len(smooth_relations) < fb_size + 1:
        print(f"未能找到足够的平滑关系。")
        return None
    print(f"已找到 {len(smooth_relations)} 个平滑关系。\n")

    dependencies = build_and_solve_matrix(smooth_relations)
    if not dependencies:
        print("未找到线性相关性。")
        return None
    print(f"找到 {len(dependencies)} 个线性相关性。\n")

    for i, dep in enumerate(dependencies):
        print(f"--- 尝试第 {i+1} 个依赖关系 ---")
        X = 1
        Y_sq_factors = defaultdict(int)
        k = 0  # 关系计数器

        for rel_idx, bit in enumerate(dep):
            if bit == 1:
                k += 1
                x_k = smooth_x_vals[rel_idx]
                params = poly_params[rel_idx]
                a_k, b_k = params["a"], params["b"]
                X = (X * (a_k * x_k + b_k)) % n
                q_x_k = a_k * x_k**2 + 2 * b_k * x_k + (b_k**2 - n) // a_k
                if q_x_k < 0:
                    Y_sq_factors[-1] += 1
                q_x_abs = abs(q_x_k)
                for j, p_fb in enumerate(factor_base):
                    if p_fb == -1:
                        continue
                    while q_x_abs % p_fb == 0:
                        Y_sq_factors[p_fb] += 1
                        q_x_abs //= p_fb
        Y_prime = 1
        for p_fb, exp in Y_sq_factors.items():
            Y_prime = (Y_prime * pow(p_fb, exp // 2, n)) % n

        # --- ✨ 修正点 3: 将 a 的平方根乘回 Y ---
        # 因为 a 是平方数 a = a_prime^2, 所以 sqrt(a^k) = a_prime^k
        sqrt_a_k = pow(a_prime, k, n)
        Y = (Y_prime * sqrt_a_k) % n

        factor = gcd(abs(X - Y), n)
        if factor != 1 and factor != n:
            print("成功找到因子！")
            return factor, n // factor

        factor = gcd(X + Y, n)
        if factor != 1 and factor != n:
            print("成功找到因子！")
            return factor, n // factor

        print("这个依赖关系只得到平凡因子，继续尝试下一个...")

    print("\n所有依赖关系都只得到平凡因子。算法失败。")
    return None


if __name__ == "__main__":
    # 再次尝试 n = 10666351
    number_to_factor = 10403
    factors = mpqs(number_to_factor, m_val=30000, fb_limit=1300)

    if factors:
        p, q = factors
        print(f"\n分解成功: {number_to_factor} = {p} * {q}")
        print(f"验证: {p * q == number_to_factor}")
