from typing import Self
from he.galois_ring.poly import Poly
from he.galois_ring.rns_poly import RNS_Poly
from he.galois_ring._util import _random
from he.galois_ring._util import _modulus

def _gcd(a, b):
    while b:
        a, b = b, a % b
    return a

# return prime factor with n in 1 ~ n
def _get_primte_factors(num : int) -> int:
    factors = set()
    d = 2
    temp = num
    while d * d <= temp:
        if temp % d == 0:
            factors.add(d)
            while temp % d == 0:
                temp //= d
        d += 1
    if temp > 1:
        factors.add(temp)
    return factors

# find primitive root of F_q*
def _find_primitive_root(prime : int) -> int:
    phi = prime - 1
    factors = _get_primte_factors(phi)
    for i in range(2, prime):
        if _gcd(i, prime) != 1:
            continue
        is_primitive = True
        for q in factors:
            if pow(i, phi // q, prime) == 1:
                is_primitive = False
                break
        if is_primitive:
            return i
    return None

# return monic irreducible polynomial of degree degree
def _sample_irreducible_poly(self, prime : int, degree : int) -> Poly:
    if degree <= 0:
        raise Exception("degree must be positive integer")
    # Randoma coins
    rho1 = [ _random._random_int(prime) for _ in range(degree)]
    rho2 = [ _random._random_int(prime) for _ in range(degree)]
    generator = _find_primitive_root(prime)
    generator_pow_prime = _modulus._modulus(generator ** prime, prime)
    generator_list = [ generator ]
    for i in range(1, degree):
        generator_list.append(_modulus._modulus(generator_list[-1] * generator_pow_prime, prime))