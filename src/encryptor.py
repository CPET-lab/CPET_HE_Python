from galois_ring.poly import Poly
from galois_ring.rns_poly import RNS_Poly
from he_parameter import HE_Parameter

class Encryptor:
    def __init__(self, parameter : HE_Parameter):
        if parameter.generate_context == False:
            raise Exception("parameter setting is not complete")
        self.param = parameter

    def encrypt(self, poly : Poly) -> RNS_Poly:
        ret = RNS_Poly(self.param.coeff_modulus)