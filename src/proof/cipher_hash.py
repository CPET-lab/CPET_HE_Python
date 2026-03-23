from typing import Self
from _util import _random
from he.galois_ring.poly import Poly
from he.galois_ring.rns_poly import RNS_Poly
from he.he_parameter import HE_Parameter
from he.ciphertext import Ciphertext
from he.encoder import Encoder
from he.evaluator import Evaluator
from he.key_generator import Key_Generator
from he.encryptor import Encryptor
from he.decryptor import Decryptor
from proof.circuit import Field

class HomHash_Manager:
    def __init__(self, param : HE_Parameter):
        self.param = param
        self.r_poly, self.min_coeff = self.sampling_r()

    def sampling_r(self) -> tuple[Poly, int]:
        param = self.param
        min_coeff = min(param.coeff_modulus)
        ranint = _random._random_int(min_coeff - 1) + 1
        r_data = [ ranint for _ in range(param.poly_modulus)]
        r_poly = Poly(param.plain_modulus, param.poly_modulus, r_data, True)
        r_poly._set_ntt_engine(param.ntt_engines[param.plain_modulus])
        return r_poly, min_coeff

    def cipher_hash(self, cipher : Ciphertext) -> RNS_Poly:
        r_poly = self.r_poly
        ret = cipher._data[-1].copy()
        r_pow = r_poly._data[0]
        if not cipher.is_ntt_form():
            cipher.transform_to_ntt_form()
        # if not r_pow.is_ntt_form():
        #     r_pow.transform_to_ntt_form()
        for rns_poly in reversed(cipher._data[:-1]):
            rns_poly_copy = rns_poly.copy()
            rns_poly_copy.mul_scalar_inplace(r_pow)
            ret += rns_poly_copy
            r_pow *= r_pow
        return ret
    
    def to_field(self, val : int) -> Field:
        return Field(val, self.min_coeff)
    
    def sampling_field(self) -> Field:
        return Field(_random._random_int(self.min_coeff), self.min_coeff)
    
    def sampling_field_range(self, start : int) -> Field:
        return Field(_random._random_int(self.min_coeff - start) + start, self.min_coeff)