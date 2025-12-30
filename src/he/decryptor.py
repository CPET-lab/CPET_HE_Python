from he.galois_ring._util import _modulus
from he.galois_ring.poly import Poly
from he.galois_ring.rns_poly import RNS_Poly
from he.he_parameter import HE_Parameter
from he.ciphertext import Ciphertext

class Decryptor:
    def __init__(self, parameter : HE_Parameter, secret_key : RNS_Poly):
        if parameter._setup_complete == False:
            raise Exception("parameter setting is not complete")
        self._param = parameter
        secret_key_copy = secret_key.copy()
        if not secret_key_copy.is_ntt_form():
            raise Exception("Secret key must be NTT form")
        self._secret_keys = [secret_key_copy]
        for idx in range(1, 40):
            self._secret_keys.append(self._secret_keys[idx - 1] * secret_key_copy)
            # print("secret key\n" + self._secret_keys[-1].transform_from_ntt_form().toString(self._param.poly_modulus, False))
            # self._secret_keys[-1].transform_to_ntt_form()

    def _recover_rns(self, rns_poly : RNS_Poly) -> Poly:
        ret = Poly(self._param.plain_modulus, self._param.poly_modulus)\
            ._set_ntt_engine(self._param.ntt_engines[self._param.plain_modulus])
        ret._is_ntt_form = False
        ret._data = [ 0 for _ in range(self._param.poly_modulus) ]
        for deg in range(self._param.poly_modulus):
            acc = 0
            for coeff_modulus in self._param.coeff_modulus:
                if deg >= len(rns_poly._rns_poly[coeff_modulus]._data):
                    break
                residue = rns_poly._rns_poly[coeff_modulus]._data[deg]
                acc += residue * self._param._basis[coeff_modulus]
            ret._data[deg] = _modulus._centered_modulus(_modulus._centered_modulus(\
                acc, self._param._total_modulus), self._param.plain_modulus)
        return ret

    def decrypt(self, encrypted : Ciphertext) -> Poly:
        encrypted_copy = encrypted.copy()
        if not encrypted_copy.is_ntt_form():
            encrypted_copy.transform_to_ntt_form()
        enc_data = encrypted_copy._data # list of rns_poly
        poly_sum = enc_data[-1] # rns_poly
        for i in range(len(enc_data) - 1):
            poly_sum.add_inplace(enc_data[i] * self._secret_keys[len(enc_data) - i - 2])
            print("i: " + str(i) + "\n" + (enc_data[i] * self._secret_keys[len(enc_data) - i - 2]).toString(20, False))
        poly_sum.transform_from_ntt_form()
        print(f"poly_sum\n{poly_sum.toString(10, False)}")
        return self._recover_rns(poly_sum) # poly