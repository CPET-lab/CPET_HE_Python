from he.galois_ring.poly import Poly
from he.he_parameter import HE_Parameter
from he.encoder import Encoder
from he.evaluator import Evaluator
from he.key_generator import Key_Generator
from he.encryptor import Encryptor
from he.decryptor import Decryptor

if __name__ == "__main__":
    parms = HE_Parameter("bv").set_poly_modulus(13).set_coeff_modulus([30, 30, 40]).set_plain_modulus(18).set_bound(1, 2)
    parms.generate_context()
    print(parms.toString())

    encoder = Encoder(parms)
    keygen = Key_Generator(parms)
    secret_key = keygen.generate_secret_key()

    secret_key.transform_from_ntt_form()
    temp_plain = encoder.coeff_encode([1, 1, 0])
    secret_key._data = secret_key._eval_rns(temp_plain)
    secret_key.transform_to_ntt_form()

    public_key = keygen.generate_public_key(secret_key)

    encryptor = Encryptor(parms, public_key)
    decryptor = Decryptor(parms, secret_key)

    plain1 = encoder.coeff_encode([1, 2])
    plain2 = encoder.coeff_encode([1, 1, 0, 0, 0, 0, 0, 0, 0, 0])
    
    c1 = encryptor.encrypt(plain1)
    c2 = encryptor.encrypt(plain2)
    pk_copy = public_key.copy()
    pk_copy._data[0].add_poly_inplace(plain2)
    c3 = c1 * c2 * c2

    d3 = decryptor.decrypt(c3)
    print("d3\n" + d3.toString(10, False))