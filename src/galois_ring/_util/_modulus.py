def _modulus(n : int, modular : int):
    if n < 0:
        if 0 <= n + modular:
            return n + modular
        else:
            _modulus(n + modular, modular)
    else:
        return n % modular