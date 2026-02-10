def _centered_modulus(n : int, modular : int) -> int:
    temp = n % modular
    if temp > modular / 2:
        temp -= modular
    return temp

def _modulus(n : int, modulus : int) -> int:
    return n % modulus

def extended_gcd(a, b):
    if a == 0:
        return b, 0, 1
    else:
        gcd, x, y = extended_gcd(b % a, a)
        return gcd, y - (b // a) * x, x

def _mod_inverse(a, m):
    gcd, x, y = extended_gcd(a, m)
    if gcd != 1:
        raise Exception('Modular inverse does not exist')
    return (x % m + m) % m