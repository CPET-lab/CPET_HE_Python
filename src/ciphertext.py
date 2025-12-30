from typing import Self
from galois_ring.poly import Poly
from galois_ring.rns_poly import RNS_Poly
from he_parameter import HE_Parameter

class Ciphertext:
    # data: list of RNS_Poly
    def __init__(self, param : HE_Parameter, data : list[RNS_Poly], error_bound : int, is_ntt_form=False):
        self._param = param
        self._coeff_modulus = param.coeff_modulus
        self._poly_modulus = param.poly_modulus
        self._debug = param._debug
        self._data = data
        self._is_ntt_form = is_ntt_form
        self._error_bound = error_bound
        if self._debug:
            for d in data:
                if d._rns_base != self._param.coeff_modulus:
                    raise Exception(f"data rns base isn't match")
                if d.is_ntt_form() != is_ntt_form:
                    raise Exception(f"data form isn't match")
    
    def __add__(self, other : Self):
        if self.size() != other.size():
            raise Exception(f"ciphertext size is different")
        if self.is_ntt_form() != other.is_ntt_form():
            raise Exception(f"ciphertext form is defferent")
        ret = Ciphertext(self._param, [], self._error_bound + other._error_bound, self._is_ntt_form)
        for idx in range(self.size()):
            ret._data.append(self._data[idx] + other._data[idx])
        return ret
    
    def __sub__(self, other : Self):
        if self.size() != other.size():
            raise Exception(f"ciphertext size is different")
        if self.is_ntt_form() != other.is_ntt_form():
            raise Exception(f"ciphertext form is defferent")
        ret = Ciphertext(self._param, [], self._error_bound + other._error_bound, self._is_ntt_form)
        for idx in range(self.size()):
            ret._data.append(self._data[idx] - other._data[idx])
        return ret
    
    def __neg__(self):
        ret = self.copy()
        for rns_poly in ret._data:
            rns_poly.neg_inplace()
        return ret
    
    def __mul__(self, other : Self):
        if self.is_ntt_form() != other.is_ntt_form():
            raise Exception(f"ciphertext form is different")
        temp_data = []
        for _ in range(self.size() + other.size() - 1):
            temp_poly = RNS_Poly(self._param.coeff_modulus, self._param.poly_modulus, self._is_ntt_form)
            temp_poly._set_ntt_engines(self._param.ntt_engines)
            temp_data.append(temp_poly)
        for i in range(self.size()):
            for j in range(other.size()):
                temp_data[i + j].add_inplace(self._data[i] * other._data[j])
        temp_error_bound = self._error_bound * other._error_bound\
             + (self._error_bound + other._error_bound) * self._param.plain_modulus
        return Ciphertext(self._param, temp_data, temp_error_bound, self._is_ntt_form)
    
    def add_inplace(self, other : Self):
        if self.size() != other.size():
            raise Exception(f"ciphertext size is different")
        for idx, rns_poly in enumerate(self._data):
            rns_poly.add_inplace(other._data[idx])
        return self
    
    def sub_inplace(self, other : Self):
        if self.size() != other.size():
            raise Exception(f"ciphertext size is different")
        for idx, rns_poly in enumerate(self._data):
            rns_poly.sub_inplace(other._data[idx])
        return self
    
    def neg_inplace(self):
        for rns_poly in self._data:
            rns_poly.neg_inplace()
        return self
    
    def mul_inplace(self, other : Self):
        if self.is_ntt_form() != other.is_ntt_form():
            raise Exception(f"ciphertext form is defferent")
        temp_data = []
        for _ in range(self.size() + other.size):
            temp_poly = RNS_Poly(self._param.coeff_modulus, self._param.poly_modulus, self._is_ntt_form)
            temp_poly._set_ntt_engines(self._param.ntt_engines)
            temp_data.append(temp_poly)
        for i in range(self.size()):
            for j in range(other.size()):
                temp_data[i + j].add_inplace(self._data[i] * other._data[j])
        self._data = temp_data
        self._error_bound = self._error_bound * other._error_bound\
             + (self._error_bound + other._error_bound) * self._param.plain_modulus
        return self
    
    def add_plain(self, other : Poly) -> Self:
        if other.is_ntt_form():
            raise Exception("plain must be coeff form")
        ret = self.copy()
        ret._data[0].add_poly_inplace(other)
        return ret
    
    def sub_plain(self, other : Poly) -> Self:
        if other.is_ntt_form():
            raise Exception("plain must be coeff form")
        ret = self.copy()
        ret._data[0].sub_poly_inplace(other)
        return ret
    
    def mul_plain(self, other : Poly) -> Self:
        if other.is_ntt_form():
            raise Exception("plain must be coeff form")
        ret = self.copy()
        for rns_poly in ret._data:
            rns_poly.mul_poly_inplace(other)
        ret._error_bound = self._error_bound * other.norm()
        return ret
    
    def add_plain_inplace(self, other : Poly) -> Self:
        if other.is_ntt_form():
            raise Exception("plain must be coeff form")
        self._data[-1].add_poly_inplace(other)
        return self
    
    def sub_plain_inplace(self, other : Poly) -> Self:
        if other.is_ntt_form():
            raise Exception("plain must be coeff form")
        self._data[-1].sub_poly_inplace(other)
        return self
    
    def mul_plain_inplace(self, other : Poly) -> Self:
        if other.is_ntt_form():
            raise Exception("plain must be coeff form")
        for rns_poly in self._data:
            rns_poly.mul_poly_inplace(other)
        self._error_bound = self._error_bound * other.norm()
        return self

    def size(self) -> int:
        return len(self._data)
    
    def copy(self) -> Self:
        temp_data = []
        for d in self._data:
            temp_data.append(d.copy())
        return Ciphertext(self._param, temp_data, self._error_bound, self._is_ntt_form)
    
    def is_ntt_form(self) -> bool:
        return self._is_ntt_form
    
    def transform_to_ntt_form(self):
        if self.is_ntt_form():
            raise Exception("alread NTT form")
        for poly in self._data:
            poly.transform_to_ntt_form()
        self._is_ntt_form = True
        return self
    
    def transform_from_ntt_form(self):
        if not self.is_ntt_form():
            raise Exception("alread basic form")
        for poly in self._data:
            poly.transform_from_ntt_form()
        self._is_ntt_form = False
        return self

    def toString(self, length=-1, print_zero = True) -> str:
        ret = ""
        for idx, poly in enumerate(self._data):
            ret += f"Poly {idx}\n"
            ret += poly.toString(length, print_zero)
        return ret