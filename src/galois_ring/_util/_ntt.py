class _NTT_Engine:
    # n: poly_modulus, q: coeff_modulus
    def __init__(self, n, q):
        self._n = n
        self._q = q
        g = self._find_primitive_root(q)
        self._psi = pow(g, (q - 1) // (2 * n), q)
        self._psi_inv = pow(self._psi, q - 2, q)
        self._n_inv = pow(n, q - 2, q)
        self._tables = [0] * n
        self._inv_tables = [0] * n
        self._precompute_tables()

    # find generator of finite field Z_q
    def _find_primitive_root(self, q):
        phi = q - 1
        factors = [2]
        d = 3
        temp = phi // 2
        while d * d <= temp:
            if temp % d == 0:
                factors.append(d)
                while temp % d == 0: temp //= d
            d += 2
        if temp > 1: factors.append(temp)
        for res in range(2, q):
            if all(pow(res, phi // f, q) != 1 for f in factors):
                return res
        return None
    
    def _precompute_tables(self):
        for i in range(self._n):
            rev_i = self._bit_rev_val(i, self._n.bit_length() - 1)
            self._tables[i] = pow(self._psi, rev_i, self._q)
            self._inv_tables[i] = pow(self._psi_inv, rev_i, self._q)
    
    def _bit_rev_val(self, i, bits):
        res = 0
        for _ in range(bits):
            res = (res << 1) | (i & 1)
            i >>= 1
        return res
    
    def _transform_to_ntt_form(self, a):
        """Cooley-Tukey NTT: Input(Normal) -> Output(Bit-reversed)"""
        n = self._n
        t = n
        m = 1
        while m < n:
            t >>= 1
            for i in range(m):
                w = self._tables[m + i]
                for j in range(2 * i * t, (2 * i + 1) * t):
                    u = a[j]
                    v = (a[j + t] * w) % self._q
                    a[j] = (u + v) % self._q
                    a[j + t] = (u - v + self._q) % self._q
            m <<= 1
        return a
    
    def _transform_from_ntt_form(self, a):
        """Gentleman-Sande INTT: Input(Bit-reversed) -> Output(Normal)"""
        n = self._n
        t = 1
        m = n >> 1
        while m > 0:
            for i in range(m):
                w = self._inv_tables[m + i]
                for j in range(2 * i * t, (2 * i + 1) * t):
                    u = a[j]
                    v = a[j + t]
                    a[j] = (u + v) % self._q
                    a[j + t] = ((u - v + self._q) * w) % self._q
            t <<= 1
            m >>= 1
        for i in range(n):
            a[i] = (a[i] * self._n_inv) % self._q
        return a

if __name__ == "__main__":
    n = 8
    q = 12289
    engine = _NTT_Engine(n, q)

    poly = [1, 2, 3, 4, 0, 0, 0, 0]
    print(f"poly: {poly}")

    # Forward (Normal -> Bit-reversed)
    curr = engine._transform_to_ntt_form(poly[:])
    print(f"NTT (Bit-reversed order): {curr}")

    # Inverse (Bit-reversed -> Normal)
    recovered = engine._transform_from_ntt_form(curr)
    print(f"INTT: {recovered}")