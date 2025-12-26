from galois_ring._util import _modulus

import random

def _random_mod_int(mod : int) -> int:
    return random.randint(0, mod)