from he.galois_ring.poly import Poly
from he.he_parameter import HE_Parameter
from he.encoder import Encoder
from he.evaluator import Evaluator
from he.key_generator import Key_Generator
from he.encryptor import Encryptor
from he.decryptor import Decryptor
from proof.cipher_hash import *
from proof.circuit import *

if __name__ == "__main__":

    parms = HE_Parameter("bv").set_poly_modulus(13).set_coeff_modulus([30, 30, 40]).set_plain_modulus(18).set_bound(1, 2)
    parms.generate_context()
    print(parms.toString())

    encoder = Encoder(parms)
    keygen = Key_Generator(parms)
    secret_key = keygen.generate_secret_key()
    public_key = keygen.generate_public_key(secret_key)

    encryptor = Encryptor(parms, public_key)
    decryptor = Decryptor(parms, secret_key)

    #############################################
    #           Homomorphic Hash Test           #
    #############################################
    # r_poly, min_coeff = sampling_r(parms)
    # plain1 = encoder.coeff_encode([1, 2])
    # plain2 = encoder.coeff_encode([1, 1, 0, 0, 0, 0, 0, 0, 0, 0])
    
    # c1 = encryptor.encrypt(plain1)
    # c2 = encryptor.encrypt(plain2)
    # pk_copy = public_key.copy()
    # pk_copy._data[0].add_poly_inplace(plain2)
    # c3 = c1 * c2 * c2

    # # d3 = decryptor.decrypt(c3)
    # # print("d3\n" + d3.toString(10, False))
    # # c3.transform_to_ntt_form()
    # c3_hash = cipher_hash(c3, r_poly)
    # print("c3_hash\n" + c3_hash.toString(10, False))

    # # c1.transform_to_ntt_form()
    # # c2.transform_to_ntt_form()
    # c1_hash = cipher_hash(c1, r_poly)
    # c2_hash = cipher_hash(c2, r_poly)
    # c3_hash2 = c1_hash * c2_hash * c2_hash
    # print("c3_hash2\n" + c3_hash2.toString(10, False))

    # print(f"equal : {c3_hash.equal(c3_hash2)}")

    #############################################
    #                Circuit Test               #
    #############################################
    plain1 = encoder.coeff_encode([1, 1])
    plain2 = encoder.coeff_encode([1, 2, 0, 0, 0, 0, 0, 0, 0, 0])
    
    c1 = encryptor.encrypt(plain1)
    c2 = encryptor.encrypt(plain2)
    c3 = c1 * c2 + c2

    d3 = decryptor.decrypt(c3)
    print("d3\n" + d3.toString(10, False))
    print()

    c = parse_circuit("0*1+2")
    print(c.toString())

    ret = c.compute_poly([c1, c2, c2])[0]
    d4 = decryptor.decrypt(ret)
    print("d4\n" + d4.toString(10, False))