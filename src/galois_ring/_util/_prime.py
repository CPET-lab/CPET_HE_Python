import random

def _miller_rabin(n, k=40):
    """
    Miller-Rabin prime test algorithm
    - n: test prime
    - k: number of test
    """
    if n <= 1: return False
    if n <= 3: return True
    if n % 2 == 0: return False

    r, d = 0, n - 1
    while d % 2 == 0:
        r += 1
        d //= 2

    for _ in range(k):
        a = random.randint(2, n - 2)
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True

def _generate_prime(bit_length : int, n_degree):
    """
    - bit_length: bit length of prime
    - n_degree: degree of modulus polynomial (1024, 2048, 4096...)
    """

    m = 2 * n_degree
    
    lower_bound = (1 << (bit_length - 1))
    upper_bound = (1 << bit_length) - 1
    
    while True:
        candidate = random.randint(lower_bound, upper_bound)
        candidate = (candidate // m) * m + 1
        while candidate >= lower_bound:
            if _miller_rabin(candidate):
                return candidate
            candidate -= m

def _generate_rns_bases(bit_lengths : list, n_degree):
    """
    - bit_lengths: list of bit length
    - n_degree: degree of modulus polynomial (1024, 2048, 4096...)
    """
    m = 2 * n_degree
    generated_primes = set()
    primes_list = []

    for bit_len in bit_lengths:
        lower_bound = (1 << (bit_len - 1))
        upper_bound = (1 << bit_len) - 1
        
        found = False
        while not found:
            candidate = random.randint(lower_bound, upper_bound)
            candidate = (candidate // m) * m + 1
            
            while candidate >= lower_bound:
                if candidate not in generated_primes:
                    if _miller_rabin(candidate):
                        generated_primes.add(candidate)
                        primes_list.append(candidate)
                        found = True
                        break
                candidate -= m
            
            if not found:
                continue
                
    return primes_list

if __name__ == "__main__":
    n = 4096
    bits = 60
    prime = _generate_prime(bits, n)
    
    print(f"prime: {prime}")
    print(f"bit length: {prime.bit_length()} bits")
    print(f"NTT condition check (p % {2*n} == 1): {prime % (2*n) == 1}")

    bit_lengths = [60, 60, 60]
    n_degree = 2048
    rns_bases = _generate_rns_bases(bit_lengths, n_degree)
    print(rns_bases)