from galois_ring._util import _modulus
from galois_ring.poly import Poly
from galois_ring.rns_poly import RNS_Poly
from he_parameter import HE_Parameter
from ciphertext import Ciphertext

class Decryptor:
    def __init__(self, parameter : HE_Parameter, secret_key : RNS_Poly):
        if parameter._setup_complete == False:
            raise Exception("parameter setting is not complete")
        self._param = parameter
        secret_key_copy = secret_key.copy()
        if not secret_key_copy.is_ntt_form():
            secret_key_copy.transform_to_ntt_form()
        self._secret_keys = [secret_key_copy]
        for idx in range(1, 40):
            self._secret_keys.append(self._secret_keys[idx - 1] * secret_key_copy)

    def _recover_rns(self, rns_poly : RNS_Poly) -> Poly:
        ret = Poly(self._param.plain_modulus, self._param.poly_modulus)\
            ._set_ntt_engine(self._param.ntt_engines[self._param.plain_modulus])
        ret._is_ntt_form = rns_poly.is_ntt_form()
        ret._data = [ 0 for _ in range(self._param.poly_modulus) ]
        for deg in range(self._param.poly_modulus):
            acc = 0
            for coeff_modulus in self._param.coeff_modulus:
                if deg >= len(rns_poly._rns_poly[coeff_modulus]._data):
                    break
                residue = rns_poly._rns_poly[coeff_modulus]._data[deg]
                # print(f"acc, residue, basis: {acc}, {residue}, {self._param._basis[coeff_modulus]}")
                acc += residue * self._param._basis[coeff_modulus]
            ret._data[deg] = _modulus._centered_modulus(ret._data[deg] + acc, self._param._total_modulus)
        return ret

    def decrypt(self, encrypted : Ciphertext) -> Poly:
        encrypted_copy = encrypted.copy()
        if not encrypted_copy.is_ntt_form():
            encrypted_copy.transform_to_ntt_form()
        enc_data = encrypted_copy._data # list of rns_poly
        poly_sum = enc_data[0] # rns_poly
        # print("poly_sum\n" + poly_sum.toString(10, False))
        for i in range(1, len(enc_data)):
            poly_sum.add_inplace(enc_data[i] * self._secret_keys[i - 1])
        poly_sum.transform_from_ntt_form()
        # print("poly_sum\n" + poly_sum.toString(10, False))
        return self._recover_rns(poly_sum) # poly