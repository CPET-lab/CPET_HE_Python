from galois_ring._util import _prime as prime
from galois_ring._util._ntt import _NTT_Engine

HE_SCHEME = {"bv", "bgv", "bfv"}

class HE_Parameter:
    def __init__(
        self, 
        scheme="bv", 
        poly_modulus_deg=0,
        coeff_modulus_bits=[0],
        plain_modulus_bit=0
    ):
        self.set_scheme(scheme)
        if poly_modulus_deg != 0:
            self.set_poly_modulus(poly_modulus_deg)
        else:
            self.poly_modulus = 0
        if coeff_modulus_bits != [0]:
            self.set_coeff_modulus(coeff_modulus_bits)
        if plain_modulus_bit != 0:
            self.set_plain_modulus(plain_modulus_bit)
        self.setup_complete = False
        
    def set_scheme(self, scheme : str):
        if scheme not in HE_SCHEME:
            raise Exception(f"scheme \"{scheme}\" is not exist")
        self.scheme = scheme
        return self
    
    def set_poly_modulus(self, poly_modulus_deg : int):
        if poly_modulus_deg <= 0:
            raise Exception("Invalid Parameter: poly_modulus_deg must be positive integer")
        self.poly_modulus = 2 ** poly_modulus_deg
        return self
    
    def set_coeff_modulus(self, coeff_modulus_bits : list):
        if len(coeff_modulus_bits) <= 0:
            raise Exception("Invalid Parameter: coeff_modulus_bit is empty")
        if self.poly_modulus == 0:
            raise Exception("poly modulus must be set before coeff modulus")
        self.coeff_modulus_bits = coeff_modulus_bits
        self.coeff_modulus = prime._generate_rns_bases(coeff_modulus_bits, self.poly_modulus)
        return self
    
    def set_plain_modulus(self, plain_modulus_bit : int):
        if plain_modulus_bit <= 0:
            raise Exception("Invalid Parameter: plain_modulus_bit must be positive integer")
        if self.poly_modulus == 0:
            raise Exception("poly modulus must be set before coeff modulus")
        self.plain_modulus_bit = plain_modulus_bit
        self.plain_modulus = prime._generate_prime(plain_modulus_bit, self.poly_modulus)
        return self
    
    def generate_context(self):
        # generate ntt tables
        if self.poly_modulus == 0 or self.coeff_modulus == None or self.plain_modulus == None:
            raise Exception("set parameters before generate context")
        self.ntt_engines = dict()
        for base in self.coeff_modulus:
            self.ntt_engines[base] = _NTT_Engine(self.poly_modulus, base)
        self.ntt_engines[self.plain_modulus] = _NTT_Engine(self.poly_modulus, self.plain_modulus)
        self.setup_complete = True
        return self
    
    def toString(self):
        ret = ""
        ret += f"scheme: {self.scheme}\n"
        if self.poly_modulus == 0:
            return ret + f"setup complete: {self.setup_complete}"
        ret += f"poly modulus: {self.poly_modulus}\n"
        if self.coeff_modulus != None:
            ret += f"coeff modulus bits: {self.coeff_modulus_bits}\n"
        if self.plain_modulus != None:
            ret += f"plain modulus: {self.plain_modulus}\n"
        ret += f"setup complete: {self.setup_complete}\n"
        return ret