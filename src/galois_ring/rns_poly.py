from typing import Self
from galois_ring.poly import Poly
from galois_ring._util._ntt import _NTT_Engine

class RNS_Poly:

    def __init__(self, rns_base : list):
        self._rns_base = rns_base
        self._rns_poly = { e: Poly(e, n) for e in rns_base }
        self._is_ntt_form = False

    def __add__(self, other : Self):
        if self._rns_poly.keys() != other._rns_poly.keys():
            raise Exception(f"Except RNS base {self.keys()} but {other.keys()}")
        ret = RNS_Poly(self._rns_base)
        for rns_base in ret:
            ret._rns_poly[rns_base] = self._rns_poly[rns_base] + other._rns_poly[rns_base]
        return ret
    
    def __neg__(self):
        ret = RNS_Poly(self._rns_base)
        for rns_base in ret:
            ret[rns_base] = -ret[rns_base]
        return ret
    
    def __sub__(self, other : Self):
        if self.keys() != other.keys():
            raise Exception(f"Except RNS base {self.keys()} but {other.keys()}")
        ret = RNS_Poly(self._rns_base)
        for rns_base in ret:
            ret._rns_poly[rns_base] = self._rns_poly[rns_base] - other._rns_poly[rns_base]
        return ret
    
    def __mul__(self, other : Self):
        if self.keys() != other.keys():
            raise Exception(f"Except RNS base {self.keys()} but {other.keys()}")
        ret = RNS_Poly(self._rns_base)
        for rns_base in ret:
            ret._rns_poly[rns_base] = self._rns_poly[rns_base] * other._rns_poly[rns_base]
        return ret

    def _set_poly(self, rns_poly : Poly):
        if rns_poly._coeff_modulus not in self._rns_poly:
            raise Exception(f"no RNS modulus {rns_poly._coeff_modulus}")
        if rns_poly.is_ntt_form() != self._is_ntt_form:
            raise Exception(f"Expect is_ntt_form {self._is_ntt_form} but {rns_poly.is_ntt_form()}")
        self._rns_poly[rns_poly._coeff_modulus] = rns_poly
        return self
    
    def _set_ntt_engines(self, prime : int, ntt_engine : _NTT_Engine):
        if prime not in self._rns_base:
            raise Exception(f"Invalid RNS modulus")
        self._rns_poly[prime]._set_ntt_engine(ntt_engine)
        return self
    
    def transform_to_ntt_form(self):
        if self._is_ntt_form == True:
            raise Exception(f"polynomial is already ntt form")
        for _, poly in self._rns_poly:
            poly.transform_to_ntt_form()
        return self
    
    def transform_from_ntt_form(self):
        if self._is_ntt_form == False:
            raise Exception(f"polynomial is already basic form")
        for _, poly in self._rns_poly.items():
            poly.transform_from_ntt_form()
        return self
    
    def copy(self) -> Self:
        ret = self.__init__(self._rns_base)
        for _, poly in self._rns_poly:
            ret._set_poly(poly.copy())
        return ret
    
    def add_inplace(self, other : Self):
        if self._rns_base != other._rns_base:
            raise Exception("Invalid parameter")
        for base, poly in self._rns_poly.items():
            poly.add_inplace(other._rns_poly[base])
        return self
    
    def sub_inplace(self, other : Self):
        if self._rns_base != other._rns_base:
            raise Exception("Invalid parameter")
        for base, poly in self._rns_poly.items():
            poly.sub_inplace(other._rns_poly[base])
        return self
    
    def mul_inplace(self, other : Self):
        if self._rns_base != other._rns_base:
            raise Exception("Invalid parameter")
        for base, poly in self._rns_poly.items():
            poly.mul_inplace(other._rns_poly[base])
        return self
    
    def toString(self, print_zero = True) -> str:
        ret = ""
        for rns_base, rns_poly in self._rns_poly.items():
            ret += f"{rns_base}: {rns_poly.toString(print_zero)}\n"
        return ret

if __name__ == "__main__":
    n = 8
    q = 12289
    rp1 = RNS_Poly([7, 11, 12289])
    p1 = Poly(q, n, [2, 8, 12, 0, 7])
    rp1._set_poly(p1)
    print(rp1.toString())