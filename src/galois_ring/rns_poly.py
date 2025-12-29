from typing import Self
from galois_ring.poly import Poly
from galois_ring._util import _modulus
from galois_ring._util._ntt import _NTT_Engine

class RNS_Poly:

    def __init__(self, rns_base : list[int], poly_modulus : int, is_ntt_form=False):
        self._rns_base = rns_base
        self._poly_modulus = poly_modulus
        self._rns_poly = { e: Poly(e, poly_modulus) for e in rns_base }
        self._is_ntt_form = is_ntt_form

    def __add__(self, other : Self):
        if self._rns_poly.keys() != other._rns_poly.keys():
            raise Exception(f"Expect RNS base {self._rns_poly.keys()} but {other._rns_poly.keys()}")
        if self.is_ntt_form() != other.is_ntt_form():
            raise Exception(f"form is different")
        ret = RNS_Poly(self._rns_base, self._poly_modulus, self._is_ntt_form)
        for rns_base in self._rns_base:
            ret._rns_poly[rns_base] = self._rns_poly[rns_base] + other._rns_poly[rns_base]
        return ret
    
    def __neg__(self):
        ret = RNS_Poly(self._rns_base)
        for rns_base in self._rns_base:
            ret[rns_base] = -ret[rns_base]
        return ret
    
    def __sub__(self, other : Self):
        if self._rns_poly.keys() != other._rns_poly.keys():
            raise Exception(f"Expect RNS base {self._rns_poly.keys()} but {other._rns_poly.keys()}")
        if self.is_ntt_form() != other.is_ntt_form():
            raise Exception(f"form is different")
        ret = RNS_Poly(self._rns_base, self._poly_modulus, self._is_ntt_form)
        for rns_base in self._rns_base:
            ret._rns_poly[rns_base] = self._rns_poly[rns_base] - other._rns_poly[rns_base]
        return ret
    
    def __mul__(self, other : Self):
        if self._rns_poly.keys() != other._rns_poly.keys():
            raise Exception(f"Expect RNS base {self._rns_poly.keys()} but {other._rns_poly.keys()}")
        if self.is_ntt_form() != other.is_ntt_form():
            raise Exception(f"form is different")
        ret = RNS_Poly(self._rns_base, self._poly_modulus, self._is_ntt_form)
        for rns_base in self._rns_base:
            ret._rns_poly[rns_base] = self._rns_poly[rns_base] * other._rns_poly[rns_base]
        return ret

    def _set_poly(self, rns_poly : Poly):
        if rns_poly._coeff_modulus not in self._rns_poly:
            raise Exception(f"no RNS modulus {rns_poly._coeff_modulus}")
        if rns_poly.is_ntt_form() != self.is_ntt_form():
            raise Exception(f"Expect is_ntt_form {self._is_ntt_form} but {rns_poly.is_ntt_form()}")
        self._rns_poly[rns_poly._coeff_modulus] = rns_poly
        return self
    
    def _eval_rns(self, poly : Poly):
        self._is_ntt_form = poly.is_ntt_form()
        for rns_base in self._rns_base:
            self._rns_poly[rns_base] = Poly(rns_base, poly._poly_modulus, [], poly.is_ntt_form())
            for idx, coeff in enumerate(poly._data):
                if idx >= len(self._rns_poly[rns_base]._data):
                    self._rns_poly[rns_base]._data.append(_modulus._centered_modulus(coeff, rns_base))
                else:
                    self._rns_poly[rns_base]._data[idx] = _modulus._centered_modulus(coeff, rns_base)
        return self

    def _set_ntt_engine(self, prime : int, ntt_engines : _NTT_Engine):
        if prime not in self._rns_base:
            raise Exception(f"Invalid RNS modulus")
        self._rns_poly[prime]._set_ntt_engine(ntt_engines)
        return self
    
    def _set_ntt_engines(self, ntt_engine : dict[int, _NTT_Engine]):
        for base in self._rns_base:
            if base not in dict.keys():
                raise Exception("Invalid ntt engine dictionary")
            self._rns_poly[base]._set_ntt_engine(ntt_engine[base])
        return self
    
    def transform_to_ntt_form(self):
        if self._is_ntt_form == True:
            raise Exception(f"polynomial is already ntt form")
        for _, poly in self._rns_poly.items():
            poly.transform_to_ntt_form()
        self._is_ntt_form = True
        return self
    
    def transform_from_ntt_form(self):
        if self._is_ntt_form == False:
            raise Exception(f"polynomial is already basic form")
        for _, poly in self._rns_poly.items():
            poly.transform_from_ntt_form()
        self._is_ntt_form = False
        return self
    
    def copy(self) -> Self:
        ret = RNS_Poly(self._rns_base, self._poly_modulus)
        ret._is_ntt_form = self._is_ntt_form
        for _, poly in self._rns_poly.items():
            ret._set_poly(poly.copy())
        return ret
    
    def add_inplace(self, other : Self):
        if self._rns_poly.keys() != other._rns_poly.keys():
            raise Exception(f"Expect RNS base {self._rns_poly.keys()} but {other._rns_poly.keys()}")
        if self.is_ntt_form() != other.is_ntt_form():
            raise Exception(f"form is different")
        for base, poly in self._rns_poly.items():
            poly.add_inplace(other._rns_poly[base])
        return self
    
    def add_poly_inplace(self, other : Poly):
        if self.is_ntt_form() != other.is_ntt_form():
            raise Exception("form is different")
        other_rns = RNS_Poly(self._rns_base, self._poly_modulus, self.is_ntt_form())
        other_rns._eval_rns(other)
        self.add_inplace(other_rns)
    
    def sub_inplace(self, other : Self):
        if self._rns_poly.keys() != other._rns_poly.keys():
            raise Exception(f"Expect RNS base {self._rns_poly.keys()} but {other._rns_poly.keys()}")
        if self.is_ntt_form() != other.is_ntt_form():
            raise Exception(f"form is different")
        for base, poly in self._rns_poly.items():
            poly.sub_inplace(other._rns_poly[base])
        return self
    
    def sub_poly_inplace(self, other : Poly):
        if self.is_ntt_form() != other.is_ntt_form():
            raise Exception("form is different")
        other_rns = RNS_Poly(self._rns_base, self._poly_modulus, self.is_ntt_form())
        other_rns._eval_rns(other)
        self.sub_inplace(other_rns)
    
    def neg_inplace(self):
        for _, poly in self._rns_poly.items():
            poly.neg_inplace()
        return self
    
    def mul_inplace(self, other : Self):
        if self._rns_poly.keys() != other._rns_poly.keys():
            raise Exception(f"Expect RNS base {self._rns_poly.keys()} but {other._rns_poly.keys()}")
        if self.is_ntt_form() != other.is_ntt_form():
            raise Exception(f"form is different")
        for base, poly in self._rns_poly.items():
            poly.mul_inplace(other._rns_poly[base])
        return self
    
    def mul_poly_inplace(self, other : Poly):
        if self.is_ntt_form() != other.is_ntt_form():
            raise Exception("form is different")
        other_rns = RNS_Poly(self._rns_base, self._poly_modulus, self.is_ntt_form())
        other_rns._eval_rns(other)
        self.mul_inplace(other_rns)
    
    def is_ntt_form(self) -> bool:
        return self._is_ntt_form
    
    def toString(self, length=-1, print_zero = True) -> str:
        ret = ""
        for rns_base, rns_poly in self._rns_poly.items():
            ret += f"{rns_base}: {rns_poly.toString(length, print_zero)}\n"
        return ret

if __name__ == "__main__":
    n = 8
    q = 12289
    rp1 = RNS_Poly([7, 11, 12289])
    p1 = Poly(q, n, [2, 8, 12, 0, 7])
    rp1._set_poly(p1)
    print(rp1.toString())