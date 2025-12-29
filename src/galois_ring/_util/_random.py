import random

def _random_mod_int(mod : int) -> int:
    return random.randint(0, mod)

def _random_centered_mod_int(mod : int) -> int:
    t = random.randint(0, mod)
    sign = random.randint(0, 1)
    if sign == 1:
        return t
    else:
        return t * -1