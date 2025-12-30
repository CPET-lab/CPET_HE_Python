from galois_ring.poly import Poly
from he_parameter import HE_Parameter
from encoder import Encoder
from evaluator import Evaluator
from key_generator import Key_Generator
from encryptor import Encryptor
from decryptor import Decryptor

if __name__ == "__main__":
    parms = HE_Parameter("bv").set_poly_modulus(13).set_coeff_modulus([30, 30, 40]).set_plain_modulus(18).set_bound(1, 2)
    parms.generate_context()
    print(parms.toString())

    encoder = Encoder(parms)
    plain1 = encoder.coeff_encode([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    plain2 = encoder.slot_encode([11, 12, 13, 14, 15, 16, 17, 18, 19, 20])
    plain3 = plain1.copy()
    print("plain 1 : " + plain1.toString(10, False))
    print("plain 2 : " + plain2.toString(10, False))
    print("plain 3 : " + plain3.toString(10, False))

    # plain1.transform_from_ntt_form()
    plain2.transform_from_ntt_form()
    plain4 = plain1 + plain2
    plain1.transform_to_ntt_form()
    plain2.transform_to_ntt_form()
    plain4.transform_to_ntt_form()
    
    print("plain 1 : " + plain1.toString(10, False))
    print("plain 2 : " + plain2.toString(10, False))
    print("plain 4 : " + plain4.toString(10, False))


    keygen = Key_Generator(parms)
    rns_plain = keygen.generate_bound_rns_poly(2)
    print("\nrns_plain\n" + rns_plain.toString(10, False))

    print("HE")
    secret_key = keygen.generate_secret_key()
    temp_plain = encoder.coeff_encode([1, 0, 0])
    secret_key._data = secret_key._eval_rns(temp_plain)
    print("secret key\n" + secret_key.toString(10, False))
    secret_key.transform_to_ntt_form()
    public_key = keygen.generate_public_key(secret_key)
    # public_key._data[0] = secret_key.copy()
    public_key._data[1] = public_key._data[0].copy().neg_inplace() * secret_key
    print("public key\n" + public_key.toString(10, False))

    encryptor = Encryptor(parms, public_key)
    decryptor = Decryptor(parms, secret_key)

    plain5 = encoder.slot_encode([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    print("\nplain5: " + plain5.toString(10, False))
    # plain5.transform_to_ntt_form()
    cipher1 = encryptor.encrypt(plain5)
    print("cipher1\n" + cipher1.toString(10, False))
    print((cipher1._data[0] * secret_key).toString(10, False))
    decrypt_plain5 = decryptor.decrypt(cipher1)
    decrypt_plain5.transform_to_ntt_form()
    print("decrypt plain: " + decrypt_plain5.toString(10, False), end="\n\n")

    decrypt_public_key = decryptor.decrypt(public_key)
    print("decrypt public key: " + decrypt_public_key.toString(10, False), end="\n\n")

    # plain2 = encoder.slot_encode([11, 12, 13, 14, 15, 16, 17, 18, 19, 20])
    plain1 = encoder.coeff_encode([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    plain2 = encoder.coeff_encode([1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    plain1.transform_to_ntt_form()
    plain2.transform_to_ntt_form()
    c1 = encryptor.encrypt(plain1)
    c2 = encryptor.encrypt(plain2)
    pk_copy = public_key.copy()
    pk_copy._data[0].add_poly_inplace(plain2)
    c3 = c1 * c2
    d3 = decryptor.decrypt(c3)
    print("c3\n" + c3.toString(10, False))
    print("d3\n" + d3.toString(10, False))