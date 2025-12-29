from galois_ring._util import _random
from galois_ring.poly import Poly
from galois_ring.rns_poly import RNS_Poly
from ciphertext import Ciphertext
from he_parameter import HE_Parameter

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
            ret.append(_random._random_centered_mod_int(bound))
        return ret

    def generate_bound_poly(self, modulus : int, bound : int) -> Poly:
        if bound < 0:
            raise Exception("bound must be positive integer")
        return Poly(modulus, self.param.poly_modulus, self._generate_random_poly(bound))\
            ._set_ntt_engine(self.param.ntt_engines[self.param.plain_modulus])

    # for public key, secret key, error
    def generate_bound_rns_poly(self, bound : int) -> RNS_Poly:
        if bound < 0:
            raise Exception("bound must be positive integer")
        ret = RNS_Poly(self.param.coeff_modulus, self.param.poly_modulus)
        if bound == 0:
            return ret
        poly = self.generate_bound_poly(bound, bound)
        ret._eval_rns(poly)
        for coeff, poly in ret._rns_poly.items():
            poly._set_ntt_engine(self.param.ntt_engines[coeff])
        return ret
    
    def generate_secret_key(self) -> RNS_Poly:
        return self.generate_bound_rns_poly(self.param._secret_key_bound)
    
    def generate_public_key(self, secret_key : RNS_Poly) -> Ciphertext:
        if not secret_key.is_ntt_form():
            raise Exception("Secret Key must be NTT form")
        c1 = self.generate_bound_rns_poly(0)
        c1._is_ntt_form = True
        for base in self.param.coeff_modulus:
            temp_poly = self.generate_bound_poly(base, base - 1)
            if secret_key.is_ntt_form():
                temp_poly._is_ntt_form = True
            c1._set_poly(temp_poly.copy())
        c0 = c1.copy()
        c0.mul_inplace(secret_key)
        error_poly = self.generate_bound_rns_poly(self.param._first_error_bound)
        error_poly.transform_to_ntt_form()
        loc_error = [ self.param.plain_modulus for _ in range(self.param.poly_modulus) ]
        loc_error_plain = Poly(self.param.plain_modulus, self.param.poly_modulus, loc_error)
        loc_error_plain._set_ntt_engine(self.param.ntt_engines[self.param.plain_modulus])
        loc_error_plain.transform_to_ntt_form()
        # c0.add_inplace(error_poly)
        ret = Ciphertext(self.param, [c0, c1.neg_inplace()], self.param._first_error_bound, True)
        return ret