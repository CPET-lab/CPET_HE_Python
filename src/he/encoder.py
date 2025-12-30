from he.galois_ring.poly import Poly
from he.he_parameter import HE_Parameter

class Encoder:
    def __init__(self, parameter : HE_Parameter):
        self.param = parameter
        self.poly_modulus = parameter.poly_modulus
        self.plain_modulus = parameter.plain_modulus

    def slot_encode(self, data : list[int]):
        if len(data) > self.poly_modulus:
            raise Exception("Too many datas")
        ret = Poly(self.plain_modulus, self.poly_modulus, data, True)
        ret._set_ntt_engine(self.param.ntt_engines[self.plain_modulus])
        return ret
    
    def coeff_encode(self, data : list[int]):
        if len(data) > self.poly_modulus:
            raise Exception("Too many datas")
        ret = Poly(self.plain_modulus, self.poly_modulus, data, False)
        ret._set_ntt_engine(self.param.ntt_engines[self.plain_modulus])
        return ret