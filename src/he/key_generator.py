from he.galois_ring._util import _random
from he.galois_ring._util import _modulus
from he.galois_ring.poly import Poly
from he.galois_ring.rns_poly import RNS_Poly
from he.ciphertext import Ciphertext
from he.he_parameter import HE_Parameter

class Key_Generator:
    def __init__(self, parameter : HE_Parameter):
        if parameter.generate_context == False:
            raise Exception("parameter setting is not complete")
        self.param = parameter

    def _generate_random_poly(self, bound : int) -> list[int]:
        if bound < 0:
            raise Exception("bound must be positive integer")
        ret = []
        for _ in range(self.param.poly_modulus):
            ret.append(_random._random_bound_int(bound))
        return ret

    def generate_bound_poly(self, modulus : int, bound : int) -> Poly:
        if bound < 0:
            raise Exception("bound must be positive integer")
        return Poly(modulus, self.param.poly_modulus, self._generate_random_poly(bound))\
            ._set_ntt_engine(self.param.ntt_engines[modulus])

    # for public key, secret key, error
    def generate_bound_rns_poly(self, bound : int) -> RNS_Poly:
        if bound < 0:
            raise Exception("bound must be positive integer")
        ret = RNS_Poly(self.param.coeff_modulus, self.param.poly_modulus)
        poly = self.generate_bound_poly(self.param.plain_modulus, bound)
        ret._eval_rns(poly)
        for coeff, poly in ret._rns_poly.items():
            poly._set_ntt_engine(self.param.ntt_engines[coeff])
        return ret
    
    def generate_secret_key(self) -> RNS_Poly:
        return self.generate_bound_rns_poly(self.param._secret_key_bound).transform_to_ntt_form()
    
    def generate_public_key(self, secret_key : RNS_Poly) -> Ciphertext:
        if not secret_key.is_ntt_form():
            raise Exception("Secret Key must be NTT form")
        c0 = self.generate_bound_rns_poly(0)
        c0.transform_to_ntt_form()
        for base in self.param.coeff_modulus:
            temp_poly = self.generate_bound_poly(base, 0)
            temp_poly._is_ntt_form = True
            c0._set_poly(temp_poly.copy())
        c1 = c0.copy()
        c1.mul_inplace(secret_key)
        loc_error = [ self.param.plain_modulus * \
            _random._random_bound_int(self.param._first_error_bound)\
            for _ in range(self.param.poly_modulus) ]
        loc_error_plain = Poly(self.param.plain_modulus, self.param.poly_modulus, loc_error)
        loc_error_plain._set_ntt_engine(self.param.ntt_engines[self.param.plain_modulus])
        c1.add_poly_inplace(loc_error_plain)
        ret = Ciphertext(self.param, [c0, c1.neg_inplace()], self.param._first_error_bound, True)
        return ret