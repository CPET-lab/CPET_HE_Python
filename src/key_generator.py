from galois_ring._util import _random
from galois_ring.poly import Poly
from galois_ring.rns_poly import RNS_Poly
from he_parameter import HE_Parameter

class Key_Generator:
    def __init__(self, parameter : HE_Parameter):
        if parameter.generate_context == False:
            raise Exception("parameter setting is not complete")
        self.param = parameter

    def generate_zero_poly(self, bound : int) -> Poly:
        if bound < 0:
            raise Exception("bound must be positive integer")
        data = []
        for _ in range(self.param.plain_modulus):
            data.append(_random._random_mod_int(bound))
        return Poly(self.param.plain_modulus, self.param.poly_modulus, data)\
            ._set_ntt_engine(self.param.ntt_engines[self.param.plain_modulus])

    # for public key, secret key, error
    def generate_zero_rns_poly(self, bound : int) -> RNS_Poly:
        if bound < 0:
            raise Exception("bound must be positive integer")