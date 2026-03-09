from he.galois_ring.poly import Poly
from he.he_parameter import HE_Parameter
from he.encoder import Encoder
from he.evaluator import Evaluator
from he.key_generator import Key_Generator
from he.encryptor import Encryptor
from he.decryptor import Decryptor
from proof.cipher_hash import *
from proof.circuit import *
from proof import giraffe
import time

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
    # hasher = HomHash_Manager(parms)
    # plain1 = encoder.coeff_encode([1, 2])
    # plain2 = encoder.coeff_encode([1, 1, 0, 0, 0, 0, 0, 0, 0, 0])
    
    # c1 = encryptor.encrypt(plain1)
    # c2 = encryptor.encrypt(plain2)
    # pk_copy = public_key.copy()
    # pk_copy._data[0].add_poly_inplace(plain2)
    # c3 = c1 * c2 + c2

    # # d3 = decryptor.decrypt(c3)
    # # print("d3\n" + d3.toString(10, False))
    # # c3.transform_to_ntt_form()
    # c3_hash = hasher.cipher_hash(c3)
    # print("c3_hash\n" + c3_hash.toString(10, False))

    # # c1.transform_to_ntt_form()
    # # c2.transform_to_ntt_form()
    # c1_hash = hasher.cipher_hash(c1)
    # c2_hash = hasher.cipher_hash(c2)
    # c3_hash2 = c1_hash * c2_hash + c2_hash
    # print("c3_hash2\n" + c3_hash2.toString(10, False))

    # print(f"equal : {c3_hash.equal(c3_hash2)}")


    #############################################
    #                Circuit Test               #
    #############################################
    # plain1 = encoder.coeff_encode([1, 1])
    # plain2 = encoder.coeff_encode([1, 2, 0, 0, 0, 0, 0, 0, 0, 0])
    
    # c1 = encryptor.encrypt(plain1)
    # c2 = encryptor.encrypt(plain2)
    # c = parse_circuit("((0*1)+(2+3))*((4*5)+(6*7))")

    # print(c.toString())

    # ret = c.compute_poly([c1, c1, c1, c1, c2, c2, c2, c2])
    # d4 = decryptor.decrypt(ret[-1][0])
    # print("d4\n" + d4.toString(10, False))

    # for i, l in enumerate(ret):
    #     print(f"\nlayer {i}:")
    #     for j, p in enumerate(l):
    #         d = decryptor.decrypt(p)
    #         print(f"{j} poly : {d.toString(10, False)}")


    #############################################
    #                Giraffe Test               #
    #############################################
    demo = giraffe.Demo()

    # field list giraffe test
    # fielder = Fielder(1009)
    # data = [ [ fielder.to_field(i) for _ in range(8)] for i in range(8) ]
    # demo.giraffe_basic(c, data, fielder, True)

    # ciphertext giraffe test
    # plain_data = [encoder.slot_encode([i for _ in range(10)]) for i in range(8)]
    # start = time.perf_counter()
    # cipher_data = [encryptor.encrypt(_plain.transform_from_ntt_form()) for _plain in plain_data]
    # end = time.perf_counter()
    # print(f"encode + encrypt time: {end - start}")
    # demo.giraffe_cipher(c, cipher_data, parms, False)

    # matrix multiplication giraffe test
    # data1 = [[i for j in range(4)] for i in range(3)]
    # data2 = [[j for j in range(2)] for i in range(4)]
    # demo.matrix_mult(data1, data2, parms, False)

    demo.poly_func([1, 3, 5, 7, 9], [2, 1], parms, False)