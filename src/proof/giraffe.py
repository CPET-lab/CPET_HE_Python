# Giraffe for BV

import copy
import math
from enum import Enum
from typing import Self
from _util import _modulus
from _util import _random
from proof.circuit import *
from proof.cipher_hash import *

class Giraffe:

    def __init__(self, circuit : Circuit, input : list[list[int]], hasher : HomHash_Manager):
        self.circuit = circuit
        self.input = []
        for l in input:
            temp = []
            for e in l:
                temp.append(hasher.to_field(e))
            self.input.append(temp)
        self.hasher = hasher
        self.subAC_num = len(input[0])
        self.demo()

    # eval circuit + generate witness
    def eval_circuit(self) -> list[Field]:
        self.witness = self.circuit.compute_list(self.input)
        return self.witness[-1][0]

    # BV Giraffe Example
    def demo(self):
        # I. Prover/Verifier Setup
        # Both : circuit C, list of input ciphertext x
        # Prover : witness (all layer's output), initial random sub-AC index
        # Verifier : claimed output y
        prover = self.Prover(self.circuit, self.input, self.hasher)
        claimed_output = prover.eval_circuit()
        verifier = self.Verifier(self.circuit, self.input, self.hasher, claimed_output)
        prover.init_proof(verifier.init_proof())

        # II. Sum-Check
        # II-1. layer loop
        for ridx, layer in enumerate(reversed(self.circuit.layers)):
            layer_idx = len(self.circuit.layers) - ridx - 1
            prover.init_sumcheck(layer_idx)

            # II-2. phase 1
            for round in range(math.ceil(math.log2(self.subAC_num))):
                eval_list, coeff_list = prover.phase1(layer_idx)
                rande = verifier.phase1(eval_list, coeff_list)
                print(f"rande: ", rande.val)
                prover.collapse(rande)










    

    class Prover:
        def __init__(self, circuit : Circuit, input : list[list[Field]], hasher : HomHash_Manager):
            self.circuit = circuit
            self.input = input
            self.hasher = hasher
            self.subAC_num = len(input[0])
        
        # eval circuit + generate witness
        def eval_circuit(self) -> list[Field]:
            self.witness = self.circuit.compute_list(self.input)
            return self.witness[-1][0]
        
        def _collapse(self, vec : list[Field], len : int, r : Field):
            for sigma in range(len // 2):
                vec[sigma] = (-r).add(1) * vec[2 * sigma] + r * vec[2 * sigma + 1]
        
        # len(r) = log |vec|
        def _multi_collapse(self, vec : list[Field], r : list[Field]):
            loglen = math.ceil(math.log2(len(vec)))
            for j in range(1, loglen + 1):
                self._collapse(vec, 2 ** (loglen + 1 - j), r[j])
            return vec[0]
        
        # def _dotp_mult_collapse(self, vec : list[Field], r : list[Field]):
        #     loglen = math.ceil(math.log2(len(vec)))

        def init_proof(self, init_subAC_idx : list[Field]):
            self.prev_gate = [ self.hasher.to_field(0) ] * 2
            self.prev_subAC = init_subAC_idx
        
        # def init_sumcheck(self, layer_idx : int):
        #     self.W = copy.deepcopy(self.witness[layer_idx])
            
        #     # init termP1
        #     self.A1 = [ self.hasher.to_field(0) for _ in range(2 ** len(self.prev_subAC)) ]
        #     self.A1[0] = (-self.prev_subAC[0]).add(1)
        #     self.A1[1] = self.prev_subAC[0]
        #     for l in range(1, len(self.prev_subAC)):
        #         for k in range(2 ** (len(self.prev_subAC) - l - 1)):
        #             self.A1[2 * k] = ((-self.prev_subAC[l]).add(1)) * self.A1[k]
        #             self.A1[2 * k + 1] = self.prev_subAC[l] * self.A1[k]
        #     print("A1")
        #     for a in self.A1:
        #         print(a)

        #     # init termP2
        #     self.A2 = [ self.hasher.to_field(0) for _ in range(2 ** len(self.prev_gate)) ]
        #     if len(self.prev_gate) == 0:
        #         self.A2[0] = self.hasher.to_field(1)
        #     else:
        #         self.A2[0] = (-self.prev_gate[-1]).add(1)
        #         self.A2[1] = self.prev_gate[-1]
        #         for l in range(len(self.prev_gate) - 2, -1, -1):
        #             for k in range(2 ** (len(self.prev_gate) - l - 1) - 1, -1, -1):
        #                 self.A2[2 * k] = ((-self.prev_gate[l]).add(1)) * self.A2[k]
        #                 self.A2[2 * k + 1] = self.prev_gate[l] * self.A2[k]
        #     print("A2")
        #     for a in self.A2:
        #         print(a)
        #     return

        def init_sumcheck(self, layer_idx : int):
            self.W = copy.deepcopy(self.witness[layer_idx])

            # init termP1
            self.A1 = [ self.hasher.to_field(0) for _ in range(2 ** len(self.prev_subAC)) ]
            self.A1[0] = (-self.prev_subAC[-1]).add(1)
            self.A1[1] = self.prev_subAC[-1]
            for l in range(len(self.prev_subAC) - 2, -1, -1):
                for k in range(2 ** (len(self.prev_subAC) - l - 1) - 1, -1, -1):
                    self.A1[2 * k] = ((-self.prev_subAC[l]).add(1)) * self.A1[k]
                    self.A1[2 * k + 1] = self.prev_subAC[l] * self.A1[k]

            # init termP2
            self.A2 = [ self.hasher.to_field(0) for _ in range(2 ** len(self.prev_gate)) ]
            self.A2[0] = (-self.prev_gate[-1]).add(1)
            self.A2[1] = self.prev_gate[-1]
            for l in range(len(self.prev_gate) - 2, -1, -1):
                for k in range(2 ** (len(self.prev_gate) - l - 1) - 1, -1, -1):
                    self.A2[2 * k] = ((-self.prev_gate[l]).add(1)) * self.A2[k]
                    self.A2[2 * k + 1] = self.prev_gate[l] * self.A2[k]
        
        def get_coefficients_d3(self, evals):
            """
            evals: [F(0), F(1), F(2), F(-1)] 순서로 담긴 Field 객체들의 리스트
            반환값: 3차 다항식의 계수 [a, b, c, d] (Field 객체들의 리스트)
            """
            # 0, 1, 2, -1 순서의 평가값 언패킹
            y_0, y_1, y_2, y_m1 = evals
            mod = y_0.mod
            
            # 6과 2의 유한체 내 곱셈 역원 계산 (mod가 소수임을 가정)
            inv6 = pow(6, mod - 2, mod)
            inv2 = pow(2, mod - 2, mod)
            
            # 1. d = y_0
            d = y_0.copy()
            
            # 2. a = (y_2 - 3*y_1 + 3*y_0 - y_minus1) * inv6
            # Field 클래스의 mul(int) 메서드와 -, + 연산자 활용
            a_numerator = y_2 - y_1.mul(3) + y_0.mul(3) - y_m1
            a = a_numerator.mul(inv6)
            
            # 3. b = (y_1 + y_minus1 - 2*y_0) * inv2
            b_numerator = y_1 + y_m1 - y_0.mul(2)
            b = b_numerator.mul(inv2)
            
            # 4. c = (y_1 - y_minus1) * inv2 - a
            c_numerator = y_1 - y_m1
            c = c_numerator.mul(inv2) - a
            
            return [d, c, b, a]

        # return (eval value, coeff)
        # [F(0), F(1), F(2), F(-1)], [X^0, X^1, X^2, X^3]
        def phase1(self, layer_idx : int) -> tuple[list[Field], list[Field]]:
            ret = [ self.hasher.to_field(0) for _ in range(4) ]
            layer = self.circuit.layer(layer_idx)
            for gate_i, gate in enumerate(layer.gates):
                for i in range(len(self.W[0]) // 2):
                    # termL, R
                    termL = []
                    termR = []
                    termL.append(self.W[gate.left][2 * i])
                    termL.append(self.W[gate.left][2 * i + 1])
                    termL.append(self.W[gate.left][2 * i + 1].mul(2) - self.W[gate.left][2 * i])
                    termL.append(self.W[gate.left][2 * i].mul(2) - self.W[gate.left][2 * i + 1])
                    termR.append(self.W[gate.right][2 * i])
                    termR.append(self.W[gate.right][2 * i + 1])
                    termR.append(self.W[gate.right][2 * i + 1].mul(2) - self.W[gate.right][2 * i])
                    termR.append(self.W[gate.right][2 * i].mul(2) - self.W[gate.right][2 * i + 1])

                    # termP1
                    termP1 = []
                    termP1.append(self.A1[2 * i])
                    termP1.append(self.A1[2 * i + 1])
                    termP1.append(self.A1[2 * i + 1].mul(2) - self.A1[2 * i])
                    termP1.append(self.A1[2 * i].mul(2) - self.A1[2 * i + 1])

                    # sum
                    if gate.op == OpType.ADD:
                        for j in range(4):
                            ret[j] += termP1[j] * self.A2[gate_i] * (termL[j] + termR[j])
                    elif gate.op == OpType.MULT:
                        for j in range(4):
                            ret[j] += termP1[j] * self.A2[gate_i] * (termL[j] * termR[j])
            coeff = self.get_coefficients_d3(ret)
            return ret, coeff
        
        def collapse(self, r : Field):
            # termL, R
            for j, sub_val in enumerate(self.W):
                for i in range(len(sub_val) // 2):
                    t0, t1 = sub_val[2 * i], sub_val[2 * i + 1]
                    sub_val[i] = (-r).add(1) * t0 + r * t1
                self.W[j] = self.W[j][:len(sub_val) // 2]

            # term P1
            for i in range(len(self.A1) // 2):
                t0, t1 = self.A1[2 * i], self.A1[2 * i + 1]
                self.A1[i] = (-r).add(1) * t0 + r * t1
            self.A1 = self.A1[:len(self.A1) // 2]
        

    class Verifier:
        def __init__(self, circuit : Circuit, input : list[list[Field]], \
                     hasher : HomHash_Manager, claimed_output : list[Field]):
            self.circuit = circuit
            self.input = input
            self.hasher = hasher
            self.subAC_num = len(input[0])
            self.claimed_output = claimed_output

        # return q' (sub-AC index)
        def init_proof(self) -> list[Field]:
            self.qp = [ self.hasher.sampling_field() for _ in range(math.ceil(math.log2(self.subAC_num))) ]
            self.qp = [ self.hasher.to_field(1) for _ in range(math.ceil(math.log2(self.subAC_num))) ]
            temp = copy.deepcopy(self.claimed_output)
            for e in temp:
                print(e)
            for j in range(len(self.qp)):
                length = 2 ** (len(self.qp) - j)
                for i in range(length // 2):
                    temp[i] = (-self.qp[j]).add(1) * temp[2 * i] + self.qp[j] * temp[2 * i + 1]
            for e in temp:
                print(e)
            self.prev_claim = temp[0]
            return self.qp
        
        def phase1(self, eval : list[Field], coeff : list[Field]) -> Field:
            # integrity check
            print(f"phase1: ", self.prev_claim, eval[0] + eval[1])
            if self.prev_claim != eval[0] + eval[1]:
                raise Exception("not valid")
            for i in range(-1, 3):
                if eval[i] != coeff[0] + coeff[1].mul(i) + coeff[2].mul(i ** 2) + coeff[3].mul(i ** 3):
                    raise Exception("not valid")
            rande = self.hasher.sampling_field()
            while rande.val < 2 or rande.val == self.hasher.min_coeff - 1:
                rande = self.hasher.sampling_field()
            self.prev_claim = coeff[0] + coeff[1] * rande + coeff[2] * rande * rande \
                              + coeff[3] * rande * rande * rande
            return rande