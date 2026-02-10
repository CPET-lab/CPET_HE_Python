import random

# sample integer on Z_q
def _random_int(q : int) -> int:
    return random.randint(0, q - 1)

def _random_bound_int(bound : int) -> int:
    t = random.randint(0, bound)
    sign = random.randint(0, 1)
    if sign == 1:
        return t
    else:
        return t * -1
    
def _random_centered_mod_int(mod : int) -> int:
    temp = _random_int(mod)
    if temp > (mod // 2):
        return temp - mod
    return temp