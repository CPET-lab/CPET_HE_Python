from galois_ring.poly import Poly
from he_parameter import HE_Parameter

class Evaluator:
    def __init__(self, parameter : HE_Parameter):
        if parameter.generate_context == False:
            raise Exception("parameter setting is not complete")
        self.param = parameter