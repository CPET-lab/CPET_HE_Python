from galois_ring.poly import Poly
from galois_ring.rns_poly import RNS_Poly
from he_parameter import HE_Parameter
from ciphertext import Ciphertext

class Encryptor:
    def __init__(self, parameter : HE_Parameter, public_key : Ciphertext):
        if public_key.size() != 2:
            raise Exception("Invalid Public Key")
        if parameter._setup_complete == False:
            raise Exception("parameter setting is not complete")
        if not public_key.is_ntt_form():
            raise Exception("Public Key must be NTT form")
        self._param = parameter
        self.public_key = public_key.copy()

    def encrypt(self, poly : Poly) -> Ciphertext:
        if poly.is_ntt_form():
            raise Exception("plaintext must be basic form")
        encrypted = self.public_key.copy()
        encrypted.add_plain_inplace(poly)
        return encrypted