import copy
from typing import Self
from homhash.circuit import *
from homhash.cipher_hash import cipher_hash
from he.galois_ring.rns_poly import RNS_Poly
from he.galois_ring.poly import Poly
from he.galois_ring._util import _modulus
from he.galois_ring._util import _random
from he.ciphertext import Ciphertext

class Sumcheck_Prover:
    def __init__(self, layer : Layer, r_poly : Poly, prev_data : list[RNS_Poly]):
        self.layer = layer
        self.r_poly = r_poly
        self.output_hash = []
        for output_cipher in layer.output:
            self.output_hash.append(cipher_hash(output_cipher, r_poly))
        self.zero_poly = self.output_hash[0] - self.output_hash[0]
        self.claim_sum = self.zero_poly.copy()
        for rns_poly in self.output_hash:
            self.claim_sum.add_inplace(rns_poly)
        self.gates = layer.gates
        self.degree = layer.max_idx
        for _ in range(len(self.gates), self.degree):
            self.gates.append(Gate(OpType.NONE, -1, -1))
        for _ in range(len(self.output_hash), self.degree):
            self.output_hash.append(self.zero_poly.copy())
        self.add = [ 0 for _ in range(self.degree) ]
        self.mult = [ 0 for _ in range(self.degree) ]
        self.left = [ _ for _ in range(self.degree) ]
        self.right = [ _ for _ in range(self.degree) ]
        for idx, gate in enumerate(layer.gates):
            if gate.op == OpType.NONE:
                continue
            if gate.op == OpType.ADD:
                self.add[idx] = 1
            elif gate.op == OpType.MULT:
                self.mult[idx] = 1
            self.left[idx] = prev_data[gate.left]
            self.right[idx] = prev_data[gate.right]
    
    def initial_claim(self) -> RNS_Poly:
        return self.claim_sum
    
    def round_claim(self) -> list[RNS_Poly]:
        g_value = [ self.zero_poly.copy() for _ in range(4) ]
        for i in range(self.degree // 2):
            add_zero = self.add[i * 2]
            mult_zero = self.mult[i * 2]
            left_zero = self.left[i * 2]
            right_zero = self.right[i * 2]
            add_one = self.add[i * 2 + 1]
            mult_one = self.mult[i * 2 + 1]
            left_one = self.left[i * 2 + 1]
            right_one = self.right[i * 2 + 1]
            g_value[0] += (left_zero + right_zero).mul_scalar(add_zero)
            g_value[0] += (left_zero * right_zero).mul_scalar(mult_zero)
            g_value[1] += (left_one + right_one).mul_scalar(add_one)
            g_value[1] += (left_one * right_one).mul_scalar(mult_one)
            add_two = add_one * 2 - add_zero
            add_three = add_one * 3 - add_zero * 2
            mult_two = mult_one * 2 - mult_zero
            mult_three = mult_one * 3 - mult_zero * 2
            left_two = left_one.mul_scalar(2) - left_zero
            left_three = left_one.mul_scalar(3) - left_zero.mul_scalar(2)
            right_two = right_one.mul_scalar(2) - right_zero
            right_three = right_one.mul_scalar(3) - right_zero.mul_scalar(2)
            g_value[2] += (left_two + right_two).mul_scalar(add_two)
            g_value[2] += (left_two * right_two).mul_scalar(mult_two)
            g_value[3] += (left_three + right_three).mul_scalar(add_three)
            g_value[3] += (left_three * right_three).mul_scalar(mult_three)
        return g_value
    
    def fold(self, r : int):
        self.degree = self.degree // 2
        op_list = [ self.add, self.mult ]
        value_list = [ self.left, self.right ]
        new_list = []
        for j in range(2):
            temp_op = []
            temp_value = []
            for i in range(self.degree):
                even = op_list[j][2 * i]
                odd = op_list[j][2 * i + 1]
                temp_op.append(even * (1 - r) + odd * r)
                even = value_list[j][2 * i]
                odd = value_list[j][2 * i + 1]
                temp_value.append(even.mul_scalar(1 - r) + odd.mul_scalar(r))
            new_list.append(copy.deepcopy(temp_op))
            new_list.append(copy.deepcopy(temp_value))
        self.add = new_list[0]
        self.left = new_list[1]
        self.mult = new_list[2]
        self.right = new_list[3]
        return self


class Sumcheck_Verifier:
    # mod: minimum RNS base
    def __init__(self, layer : Layer, r_poly : Poly, claim_sum : RNS_Poly, mod : int):
        self.layer = layer
        self.r_poly = r_poly
        self.claim_sum = claim_sum.copy()
        self.degree = layer.max_idx
        self.mod = mod
        self.r_history = []
    
    def _interpolate_deg3(self, g_value : list[RNS_Poly], r : int):
        y0, y1, y2, y3 = g_value
        t0 = y0.mul_scalar((r - 1) * (r - 2) * (r - 3))
        t0.mul_scalar_inplace(_modulus._mod_inverse(-6, self.mod))
        t1 = y1.mul_scalar(r * (r - 2) * (r - 3))
        t1.mul_scalar_inplace(_modulus._mod_inverse(2, self.mod))
        t2 = y2.mul_scalar(r * (r - 1) * (r - 3))
        t2.mul_scalar_inplace(_modulus._mod_inverse(-2, self.mod))
        t3 = y3.mul_scalar(r * (r - 1) * (r - 2))
        t3.mul_scalar_inplace(_modulus._mod_inverse(6, self.mod))
        return t0 + t1 + t2 + t3

    def round_verify(self, g_value : list[RNS_Poly]) -> int:
        new_claim_sum = g_value[0] + g_value[1]
        if not self.claim_sum.equal(new_claim_sum):
            raise Exception("sum is different")
        r = _modulus._centered_modulus(_random._random_int(self.mod - 4) + 4, self.mod)
        self.r_history.append(r)
        self.claim_sum = self._interpolate_deg3(g_value, r)
        return r
    
    def final_verify(self, add : int, mult : int, left : RNS_Poly, right : RNS_Poly) -> bool:
        final_claim = (left + right).mul_scalar(add) + (left * right).mul_scalar(mult)
        if self.claim_sum.equal(final_claim):
            return True
        else:
            return False