from homhash.circuit import *
from he.galois_ring.poly import Poly
from he.he_parameter import HE_Parameter
from he.encoder import Encoder
from he.evaluator import Evaluator
from he.key_generator import Key_Generator
from he.encryptor import Encryptor
from he.decryptor import Decryptor
from homhash.sumcheck import *
from homhash.cipher_hash import *

# string with (+, *)
# "+*+" -> [a0 + a1, a2 * a3, a4 + a5]
def string_to_layer(string : str) -> Layer:
    gate_list = []
    for i, c in enumerate(list(string)):
        if c == '+':
            gate_list.append(Gate(OpType.ADD, i * 2, i * 2 + 1))
        elif c == '*':
            gate_list.append(Gate(OpType.MULT, i * 2, i * 2 + 1))
    return Layer(gate_list)

if __name__ == "__main__":

    # HE Setting
    parms = HE_Parameter("bv").set_poly_modulus(10).set_coeff_modulus([30, 30, 40]).set_plain_modulus(18).set_bound(1, 2)
    parms.generate_context()
    encoder = Encoder(parms)
    keygen = Key_Generator(parms)
    secret_key = keygen.generate_secret_key()
    public_key = keygen.generate_public_key(secret_key)
    encryptor = Encryptor(parms, public_key)
    decryptor = Decryptor(parms, secret_key)

    plain_list = []
    plain_list.append(encoder.slot_encode([1, 2]).transform_from_ntt_form())
    plain_list.append(encoder.slot_encode([1, 1]).transform_from_ntt_form())
    plain_list.append(encoder.slot_encode([2, 1]).transform_from_ntt_form())
    plain_list.append(encoder.slot_encode([3, 7]).transform_from_ntt_form())
    plain_list.append(encoder.coeff_encode([1]))
    plain_list.append(encoder.coeff_encode([0]))

    cipher_list = []
    for i in range(len(plain_list)):
        cipher_list.append(encryptor.encrypt(plain_list[i]))


    # CirCuit Setting
    layer_list = []
    layer_list.append(string_to_layer("+")) # output layer
    layer_list.append(string_to_layer("+*")) # input layer, 4 data
    # layer_list.append(string_to_layer("+++*"))
    circuit = Circuit(layer_list, cipher_list[-2])

    res = circuit.compute_circuit(cipher_list)
    decrypted = decryptor.decrypt(res)
    print(f"result: {decrypted.transform_to_ntt_form().toString(10, False)}")


    # Sumcheck Test
    r_poly, min_coeff = sampling_r(parms)
    test_layer = circuit.layers[-1]
    prev_data = []
    for cipher in cipher_list:
        print("t")
        prev_data.append(cipher_hash(cipher, r_poly))
    print(1)
    prover = Sumcheck_Prover(test_layer, r_poly, prev_data)
    print(2)
    initial_claim = prover.initial_claim()
    print(3)
    verifier = Sumcheck_Verifier(test_layer, r_poly, initial_claim, min_coeff)

    for i in range(3):
        print(f"round {i}")
        g_value = prover.round_claim()
        r = verifier.round_verify(g_value)
        prover.fold(r)
    if verifier.final_verify(prover.add[0], prover.mult[0], prover.left[0], prover.right[0]):
        print("Success")
    else:
        print("Fail")