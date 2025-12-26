from _util import _modulus as mod
import _poly as poly

class _RNS_Poly:

    def __init__(self, rns_base : list):
        self._rns_base = rns_base
        self._rns_poly = { e: poly._Poly(e, n) for e in rns_base }

    def __add__(self, other):
        if self.keys() != other.keys():
            raise Exception(f"Except RNS base {self.keys()} but {other.keys()}")
        ret = _RNS_Poly(self._rns_base)
        for rns_base in ret:
            ret._rns_poly[rns_base] = self._rns_poly[rns_base] + other._rns_poly[rns_base]
        return ret
    
    def __neg__(self):
        ret = _RNS_Poly(self._rns_base)
        for rns_base in ret:
            ret[rns_base] = -ret[rns_base]
        return ret
    
    def __sub__(self, other):
        if self.keys() != other.keys():
            raise Exception(f"Except RNS base {self.keys()} but {other.keys()}")
        ret = _RNS_Poly(self._rns_base)
        for rns_base in ret:
            ret._rns_poly[rns_base] = self._rns_poly[rns_base] - other._rns_poly[rns_base]
        return ret
    
    def __mul__(self, other):
        if self.keys() != other.keys():
            raise Exception(f"Except RNS base {self.keys()} but {other.keys()}")
        ret = _RNS_Poly(self._rns_base)
        for rns_base in ret:
            ret._rns_poly[rns_base] = self._rns_poly[rns_base] * other._rns_poly[rns_base]
        return ret

    def _set_poly(self, rns_poly : poly._Poly):
        if rns_poly._coeff_modulus not in self._rns_poly:
            raise Exception(f"no RNS modulus {rns_poly._coeff_modulus}")
        self._rns_poly[rns_poly._coeff_modulus] = rns_poly
        return self
    
    def toString(self, print_zero = True) -> str:
        ret = ""
        for rns_base, rns_poly in self._rns_poly.items():
            ret += f"{rns_base}: {rns_poly.toString(print_zero)}\n"
        return ret

if __name__ == "__main__":
    n = 8
    q = 12289
    rp1 = _RNS_Poly([7, 11, 12289])
    p1 = poly._Poly(q, n, [2, 8, 12, 0, 7])
    rp1._set_poly(p1)
    print(rp1.toString())