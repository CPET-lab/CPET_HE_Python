from he_parameter import HE_Parameter
from encoder import Encoder
from evaluator import Evaluator
from key_generator import Key_Generator

if __name__ == "__main__":
    parms = HE_Parameter("bv").set_poly_modulus(13).set_coeff_modulus([30, 30, 40]).set_plain_modulus(18)
    parms.generate_context()
    print(parms.toString())

    encoder = Encoder(parms)
    plain1 = encoder.slot_encode([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    plain2 = encoder.slot_encode([11, 12, 13, 14, 15, 16, 17, 18, 19, 20])
    plain3 = plain1.copy()
    print("plain 1 : " + plain1.toString(10, False))
    print("plain 2 : " + plain2.toString(10, False))
    print("plain 3 : " + plain3.toString(10, False))

    plain1.transform_from_ntt_form()
    plain2.transform_from_ntt_form()
    plain4 = plain1 + plain2
    plain1.transform_to_ntt_form()
    plain2.transform_to_ntt_form()
    plain4.transform_to_ntt_form()
    
    print("plain 1 : " + plain1.toString(10, False))
    print("plain 2 : " + plain2.toString(10, False))
    print("plain 4 : " + plain4.toString(10, False))

    keygen = Key_Generator(parms)
    plain5 = keygen.generate_zero_poly(2)
    print("plain 5 : " + plain5.toString(10, True))