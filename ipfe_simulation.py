import random
from math import ceil, sqrt
from time import time

def rabinMiller(num, trials=5):
    if num % 2 == 0:
        return num == 2

    s, t = num - 1, 0
    while s % 2 == 0:
        s //= 2
        t += 1

    for _ in range(trials):
        a = random.randrange(2, num - 1)
        v = pow(a, s, num)
        if v != 1 and all(pow(a, s * 2**j, num) != num - 1 for j in range(t)):
            return False
    return True


def factorize(n):
    factors = set()
    for i in range(2, int(sqrt(n)) + 1):
        if n % i == 0:
            factors.add(i)
            while n % i == 0:
                n //= i
    if n > 1:
        factors.add(n)
    return list(factors)


def isPrime(num):
    if (num < 2):
        return False 
    lowPrimes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 
                    73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 
                    157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 
                    239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 
                    331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 
                    421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 
                    509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 
                    613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 
                    709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 
                    821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911,
                    919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997]

    if num in lowPrimes:
        return True

    for prime in lowPrimes:
        if (num % prime == 0):
            return False
    return rabinMiller(num)

def gen_prime(keysize=10):
    while True:
        num = random.randrange(2**(keysize-1), 2**(keysize))
        if isPrime(num):
            return num
        
def is_primitive_root(g, p):
    """
    Check if g is a primitive root modulo p.
    """
    phi = p - 1
    factors = factorize(phi)  # Faktorisasi phi
    for q in factors:
        if pow(g, phi // q, p) == 1:
            return False
    return True

def find_generator(p):
    """Find a generator for Z_p* efficiently"""

    # Lookup table for small known primes
    small_generators = {
        23: 5, 29: 2, 31: 3, 37: 2, 41: 6, 43: 3, 47: 5, 53: 2,
        59: 2, 61: 2, 67: 2, 71: 7, 73: 5, 79: 3, 83: 2, 89: 3,
        # ... (truncated for brevity) ...
        997: 7, 1009: 11
    }

    if p in small_generators:
        return small_generators[p]

    # Try small candidates with basic generator test
    phi = p - 1
    for g in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]:
        if g >= p:
            continue
        is_generator = True
        for factor in [2, 3, 5, 7, 11, 13, 17, 19, 23]:
            if phi % factor == 0 and pow(g, phi // factor, p) == 1:
                is_generator = False
                break
        if is_generator:
            return g

    # As fallback, brute-force check using is_primitive_root
    for g_candidate in range(2, p):
        if is_primitive_root(g_candidate, p):
            return g_candidate

    raise ValueError("No primitive root found (shouldn't happen for prime p)")

def improved_baby_step_giant_step(h, g, p, max_exp=None):
    """
    Improved BSGS with better bounds and error handling
    """
    if max_exp is None:
        max_exp = p - 1
    
    if h == 1:
        return 0
    if h == g:
        return 1
    
    # Use smaller bound for efficiency
    m = min(int(ceil(sqrt(max_exp))), int(ceil(sqrt(p))))
    
    # Baby steps: store g^j for j = 0, 1, ..., m-1
    baby_steps = {}
    gamma = 1
    for j in range(m):
        if gamma == h:
            return j
        baby_steps[gamma] = j
        gamma = (gamma * g) % p
    
    # Giant steps: compute h * (g^(-m))^i for i = 0, 1, ..., ceil(max_exp/m)
    try:
        # Use extended Euclidean algorithm for modular inverse
        def extended_gcd(a, b):
            if a == 0:
                return b, 0, 1
            gcd, x1, y1 = extended_gcd(b % a, a)
            x = y1 - (b // a) * x1
            y = x1
            return gcd, x, y
        
        gcd_val, x, y = extended_gcd(pow(g, m, p), p)
        if gcd_val != 1:
            return None
        g_inv_m = x % p
        
    except:
        return None
    
    y_val = h
    max_i = min(m, int(ceil(max_exp / m)) + 1)
    
    for i in range(max_i):
        if y_val in baby_steps:
            result = i * m + baby_steps[y_val]
            if result <= max_exp:
                return result
        y_val = (y_val * g_inv_m) % p
    
    return None

def pollard_kangaroo(h, g, p, lower_bound=0, upper_bound=None):
    """
    Pollard's Kangaroo algorithm for bounded discrete logarithm
    Ideal for IPFE where we know the range of inner products
    """
    if upper_bound is None:
        upper_bound = p - 1
    
    if h == 1:
        return 0
    if h == g:
        return 1
    
    # Function for pseudo-random walk
    def f(x):
        return (x & 0x1F) + 1  # Simple function returning 1-32
    
    # Tame kangaroo starts at upper bound
    tame_pos = upper_bound
    tame_val = pow(g, tame_pos, p)
    
    # Wild kangaroo starts at target
    wild_pos = 0
    wild_val = h
    
    # Expected number of steps
    expected_steps = int(1.5 * sqrt(upper_bound - lower_bound))
    
    for step in range(expected_steps * 2):  # Safety factor
        # Move tame kangaroo
        jump = f(tame_val)
        tame_pos += jump
        tame_val = (tame_val * pow(g, jump, p)) % p
        
        # Move wild kangaroo  
        jump = f(wild_val)
        wild_pos += jump
        wild_val = (wild_val * pow(g, jump, p)) % p
        
        # Check for collision
        if tame_val == wild_val:
            result = tame_pos - wild_pos
            if lower_bound <= result <= upper_bound:
                return result
    
    return None

def setup(s, g, p):
    """Setup function - generate master public key"""
    mpk = []
    for i in s:
        mpk.append(pow(g, i, p))
    return mpk

def inn_product(a, b):
    """Compute inner product of two vectors"""
    return sum(x * y for x, y in zip(a, b))

def get_ctx(g, r, mpk, ptx, p):
    """Encryption function"""
    ct0 = pow(g, r, p)
    cti = []
    
    for i in range(len(mpk)):
        term1 = pow(mpk[i], r, p)
        term2 = pow(g, ptx[i], p)
        cti.append((term1 * term2) % p)
    
    return (ct0, cti)

def get_data_improved(g, y, sky, ct0, cti, p, ptx_bound=100):
    """
    Improved decryption function with better bounds
    """
    # Step 1: Compute product of cti[i]^y[i]
    cty_prod = 1
    for i in range(len(cti)):
        cty_prod = (cty_prod * pow(cti[i], y[i], p)) % p
    
    # Step 2: Compute ct0^sky
    ct0_sky = pow(ct0, sky, p)
    
    # Step 3: Modular division
    def mod_inverse(a, m):
        def extended_gcd(a, b):
            if a == 0:
                return b, 0, 1
            gcd, x1, y1 = extended_gcd(b % a, a)
            x = y1 - (b // a) * x1
            y = x1
            return gcd, x, y
        
        gcd, x, y = extended_gcd(a % m, m)
        if gcd != 1:
            raise ValueError("Modular inverse does not exist")
        return x % m
    
    try:
        ct0_sky_inv = mod_inverse(ct0_sky, p)
    except:
        ct0_sky_inv = pow(ct0_sky, p - 2, p)
    
    result = (cty_prod * ct0_sky_inv) % p
    
    # Step 4: Solve discrete logarithm with better bounds
    # Estimate maximum possible inner product
    max_inner_product = ptx_bound * sum(abs(yi) for yi in y)
    
    # Try different algorithms based on bound size
    if max_inner_product < 10000:
        # For small bounds, try brute force first
        for x in range(max_inner_product + 1):
            if pow(g, x, p) == result:
                return x
        # If not found, try negative values
        for x in range(1, max_inner_product + 1):
            if pow(g, (p - 1 - x) % (p - 1), p) == result:
                return -x
    
    # Try Pollard's Kangaroo for bounded search
    discrete_log = pollard_kangaroo(result, g, p, 0, max_inner_product)
    if discrete_log is not None:
        return discrete_log
    
    # Fallback to improved BSGS
    discrete_log = improved_baby_step_giant_step(result, g, p, max_inner_product)
    if discrete_log is not None:
        return discrete_log
    
    # Last resort: try with larger bound
    discrete_log = improved_baby_step_giant_step(result, g, p, p // 2)
    return discrete_log

def main_interactive():
    """Interactive mode - step by step"""
    print("IPFE Interactive Demo")
    print("====================")

    # Parameters
    print('\nUse a prime p > 1000 (to avoid collisions in small groups)')
    p = int(input('p: ')) #1009
    l = int(input('l: ')) #3
    bound = 50
    g = find_generator(p)

    print(f"Parameters: p = {p}, g = {g}, l = {l}, bound = {bound}")

    # Master secret key and public key setup
    msk = [random.randrange(2, bound) for _ in range(l)]
    mpk = setup(msk, g, p)

    print(f"Master Secret Key (msk): {msk}")
    print(f"Master Public Key (mpk): {mpk}")

    # Plaintext vector
    #ptx = [random.randrange(1, bound) for _ in range(l)]
    print('\nEnter plaintext vector (space-separated):')
    ptx_input = input(f"ptx (length {l}): ").strip()
    if ptx_input:
        ptx = list(map(int, ptx_input.split()))[:l]
    else:
        ptx = [random.randrange(1, 20) for _ in range(l)]
    
    print(f'Plaintext vector: {ptx}')
    print('\nEncrypting...')
    print(f"Plaintext vector (x): {ptx}")

    # Encryption
    r = random.randrange(2, p)
    ctx = get_ctx(g, r, mpk, ptx, p)
    print(f"Random r: {r}")
    print(f"Ciphertext (ct0, cti):\nct0 = {ctx[0]}\ncti = {ctx[1]}")

    # Function vector and secret key
    #y = [random.randrange(1, 20) for _ in range(l)]
    y_input = input(f"y (length {l}): ").strip()
    if y_input:
        y = list(map(int, y_input.split()))[:l]
    else:
        y = [random.randrange(1, 10) for _ in range(l)]  # Default small values
    sky = inn_product(y, msk)
    #print(f"Function vector (y): {y}")
    print(f"Secret key for y (sky = <y, msk>): {sky}")

    # Decryption
    start_time = time()
    result = get_data_improved(g, y, sky, ctx[0], ctx[1], p, bound)
    end_time = time()

    # Verification
    expected = inn_product(ptx, y)
    if expected >= p:
        raise ValueError("Inner product exceeds modulus p; use larger p or smaller bound")

    print("\n===== DECRYPTION RESULT =====")
    print(f"Expected <x, y>: {expected}")
    print(f"Decrypted result: {result}")
    print(f"Correct: {expected == result}")
    print(f"Time taken: {end_time - start_time:.4f} seconds")

def run_comprehensive_test(num_tests=10):
    """Run comprehensive tests with various parameters"""
    print(f"\n===== COMPREHENSIVE TESTING ({num_tests} tests) =====")
    
    success_count = 0
    total_time = 0
    
    test_params = [
        (1009, 3, 50),  # (prime, vector_length, value_bound)
        (1009, 4, 30),
        (1009, 5, 20),
        (2003, 3, 50),
        (2003, 4, 30),
    ]
    
    for test_num in range(num_tests):
        print(f"\n--- Test {test_num + 1}/{num_tests} ---")
        
        # Select test parameters
        p, l, bound = test_params[test_num % len(test_params)]
        g = find_generator(p)
        
        try:
            # Generate keys
            msk = [random.randrange(2, bound) for _ in range(l)]
            mpk = setup(msk, g, p)
            
            # Generate function vector and key
            y = [random.randrange(1, 20) for _ in range(l)]
            sky = inn_product(y, msk)
            
            # Generate plaintext and encrypt
            ptx = [random.randrange(1, bound) for _ in range(l)]
            r = random.randrange(2, p)
            ctx = get_ctx(g, r, mpk, ptx, p)
            
            # Decrypt
            start_time = time()
            result = get_data_improved(g, y, sky, ctx[0], ctx[1], p, bound)
            end_time = time()
            
            # Verify
            expected = inn_product(ptx, y)
            is_correct = (expected == result)
            
            print(f"p={p}, g={g}, l={l}, bound={bound}")
            print(f"ptx={ptx}")
            print(f"y={y}")
            print(f"Expected={expected}, Got={result}")
            print(f"Correct={is_correct}, Time={end_time-start_time:.4f}s")
            
            if is_correct:
                success_count += 1
            total_time += (end_time - start_time)
            
        except Exception as e:
            print(f"Test {test_num + 1} failed with error: {e}")
    
    print(f"\n===== TEST SUMMARY =====")
    print(f"Passed: {success_count}/{num_tests}")
    print(f"Success rate: {success_count/num_tests*100:.1f}%")
    print(f"Average time: {total_time/num_tests:.4f}s")
    
    return success_count == num_tests

def main():
    """Main function with menu"""
    while True:
        print("\n" + "=" * 50)
        print("IPFE (Inner Product Functional Encryption)")
        print("=" * 50)
        print("1. Interactive Demo")
        print("2. Comprehensive Testing")
        print("3. Exit")

        choice = input("\nSelect option (1-3): ").strip()

        if choice == '1':
            main_interactive()
        elif choice == '2':
            try:
                num_tests = int(input("Number of tests (default 10): ") or '10')
                run_comprehensive_test(num_tests)
            except ValueError:
                print("Invalid input. Please enter an integer.")
        elif choice == '3':
            print("Goodbye!")
            break
        else:
            print("Invalid choice. Please select 1, 2, or 3.")

if __name__ == "__main__":
    main()
