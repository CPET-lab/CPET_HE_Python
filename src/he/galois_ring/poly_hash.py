# Homomorphic Hash Function for Ciphertext
# Ciphertext(list of RNS_Poly) -> RNS_Poly

from he.galois_ring._util._ntt import _NTT_Engine
from he.galois_ring._util import _random
from he.galois_ring.poly import Poly
from he.galois_ring.rns_poly import RNS_Poly
from he.ciphertext import Ciphertext

def _sampling_poly_r(poly_modulus : int, coeff_modulus : int, ntt_engine : _NTT_Engine) -> Poly:
    data = [ _random._random_centered_mod_int(coeff_modulus) for _ in range(poly_modulus) ]
    return Poly(coeff_modulus, poly_modulus, data)._set_ntt_engine(ntt_engine)

# def cipher_hash(cipher : Ciphertext, num=-1) -> set[RNS_Poly]:
#     ret = set()