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
from _util import _modulus
import time
import random

if __name__ == "__main__":

    # test 15, [50, 50, 50, 60], 29, 6
    parms = HE_Parameter("bv").set_poly_modulus(15).set_coeff_modulus([50, 50, 55, 60]).set_plain_modulus(29).set_bound(1, 2)
    parms.generate_context()
    print(parms.toString())

    encoder = Encoder(parms)
    keygen = Key_Generator(parms)
    secret_key = keygen.generate_secret_key()
    # secret_key._rns_poly[parms.coeff_modulus[0]]._data = [0 for _ in range(parms.poly_modulus)]
    public_key = keygen.generate_public_key(secret_key)

    encryptor = Encryptor(parms, public_key)
    decryptor = Decryptor(parms, secret_key)

    demo = giraffe.Demo()

    mod = min(parms.coeff_modulus)
    fielder = Fielder(mod)
    # data = [random.randint(2, 15) for _ in range(parms.poly_modulus)]
    data = [i % 14 + 2 for i in range(parms.poly_modulus)]
    data_backup = copy.deepcopy(data)
    # coeff = [2, 1, 0, 6, 2, 1]
    # coeff = [15, 14, 13, 12, 11]
    # coeff = [6, 5, 4, 3, 2]
    coeff = [6, 14, 2, 2, 15, 11, 15]
    circuit = build_poly_circuit(coeff)
    # print(circuit.toString())

    input_plain = [
        encoder.slot_encode(data).transform_from_ntt_form(),
        *[encoder.coeff_encode([_coeff]) for _coeff in coeff],
        encoder.coeff_encode([0])
    ]

    input_enc = [encryptor.encrypt(p) for p in input_plain]
    _result_enc = circuit.compute_poly(input_enc)
    # print("---\n")
    # for rrr in _result_enc:
    #     for _rrr in rrr:
    #         print(_rrr.toString(10, False))
    #     print()
    # for rrr in _result_enc:
    #     for _rrr in rrr:
    #         print(decryptor.decrypt(_rrr).transform_to_ntt_form().toString(10, False))
    #     print()
    # print("---\n")
    result_enc = _result_enc[-1]
    print(len(result_enc[0]._data))
    ddd = decryptor.decrypt(result_enc[0]).transform_to_ntt_form()
    print(ddd.toString(10, False))
    # result_plain = parms.plain_modulus
    # for i, c in enumerate(coeff[1:]):
    #     result_plain += (data_backup[0] ** (i + 1)) * c
    #     # result_plain %= parms.plain_modulus
    #     # print(c, data_backup[0] ** (i + 1), result_plain)
    # print(result_plain)

    rt = []
    for i in range(parms.poly_modulus):
        r = coeff[0]
        for j, c in enumerate(coeff[1:]):
            r += (data_backup[i] ** (j + 1)) * c
            # r %= parms.plain_modulus
        # if ddd._data[i] != r:
        #     print(f"{i} false")
        rt.append(r)
    print(rt[:10])

    # print(coeff)
    # print(coeff[::-1])
    rrt = [giraffe._horner(coeff[::-1], d) for d in data_backup]
    print(rrt[:10])

    same = all(rrt[i] == ddd._data[i] for i in range(parms.poly_modulus))
    if not same:
        print("plain - cipher false")

    # demo = giraffe.Demo()
    # encrypt_time, decrypt_time, prover_time, verifier_time, cipher_eval_time, hash_time, decrypt_result\
    #     = demo.poly_func(copy.deepcopy(data_backup), copy.deepcopy(coeff), parms, False)
    # same = all(rt[i] == decrypt_result._data[i] for i in range(parms.poly_modulus))
    # if not same:
    #     print("plain - cipher false 22")


    ppp = [[Field(d, parms.plain_modulus) for d in data_backup],
           *[[Field(_coeff, parms.plain_modulus) for _ in range(parms.poly_modulus)] for _coeff in coeff], 
           [Field(0, parms.plain_modulus) for _ in range(parms.poly_modulus)]]
    # ppp = [
    #     data_backup,
    #     *[[c for _ in range(parms.poly_modulus)] for c in coeff],
    #     [0] * parms.poly_modulus
    # ]
    print(len(ppp))
    # pppp = circuit.compute_list(ppp)
    # print("---")
    # for l in pppp:
    #     for g in l:
    #         for i in range(10):
    #             print(f"{g[i]}, ", end=" ")
    #         print()
    #     print()
    # print("---")

    # for e in input_plain:
    #     e.transform_to_ntt_form()
    # qwer = circuit.compute_poly(input_plain)[0]
    # print(qwer.toString(10, False))

    # print(_result_enc[2][0].toString(10, False))
    # print((_result_enc[2][0]*_result_enc[2][0]).toString(10, False))
    # print(_result_enc[3][1].toString(10, False))

    # print(decryptor.decrypt(_result_enc[2][0]).transform_to_ntt_form().toString(10, False))
    # print(decryptor.decrypt(_result_enc[2][0]*_result_enc[2][0]).transform_to_ntt_form().toString(10, False))
    # print(decryptor.decrypt(_result_enc[3][1]).transform_to_ntt_form().toString(10, False))

    # print(_result_enc[3][1].toString(10, False))
    # print(_result_enc[3][6].toString(10, False))
    # print((_result_enc[3][1]*_result_enc[3][6]).toString(10, False))
    # print(_result_enc[4][1].toString(10, False))

    # print(decryptor.decrypt(_result_enc[3][1]).transform_to_ntt_form().toString(10, False))
    # print(decryptor.decrypt(_result_enc[3][6]).transform_to_ntt_form().toString(10, False))
    # print(decryptor.decrypt(_result_enc[3][1]*_result_enc[3][6]).transform_to_ntt_form().toString(10, False))
    # print(decryptor.decrypt(_result_enc[4][1]).transform_to_ntt_form().toString(10, False))