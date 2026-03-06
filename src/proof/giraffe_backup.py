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

        # II. Layer loop
        for ridx, layer in enumerate(reversed(self.circuit.layers)):
            layer_idx = len(self.circuit.layers) - ridx - 1
            print(f"\n\n** layer {layer_idx} **")
            verifier.init_sumcheck(layer_idx)

            # II-1. Sumcheck phase 1
            prover.init_phase1(layer_idx)
            for i in range(math.ceil(math.log2(self.subAC_num))):
                eval_list, coeff_list = prover.phase1(layer_idx)
                rande = verifier.phase1(eval_list, coeff_list, i)
                prover.collapse_phase1(rande)

            # II-2. Sumcheck phase 2
            prover.init_phase2(layer_idx)
            if layer_idx == 0:
                next_layer_length = layer.length() * 2
            else:
                next_layer_length = self.circuit.layer(layer_idx - 1).length()
            print("next", next_layer_length)
            for i in range(2 * math.ceil(math.log2(next_layer_length))):
                _, coeff_list = prover.phase2(layer_idx, i)
                rande = verifier.phase2(layer_idx, coeff_list, i)
                prover.collapse_phase2(layer_idx, i, rande)
            
            # II-3. End Sumcheck
            tau = verifier.end_sumcheck(layer_idx, prover.end_sumcheck(layer_idx))
            prover.init_sumcheck(tau)
        
        print("Accept")
            









    

    class Prover:
        def __init__(self, circuit : Circuit, input : list[list[Field]], hasher : HomHash_Manager):
            self.circuit = circuit
            self.input = input
            self.hasher = hasher
            self.subAC_num = len(input[0])
        
        # eval circuit + generate witness
        def eval_circuit(self) -> list[Field]:
            self.witness = self.circuit.compute_list(self.input)  # [layer][gate][sub-AC]
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

        def init_proof(self, init_subAC_idx : list[Field]):
            self.prev_gate = [ self.hasher.to_field(0) ] * 2
            self.prev_subAC = init_subAC_idx

        # EvalTermLR
        def init_phase1(self, layer_idx : int):
            self.W = copy.deepcopy(self.witness[layer_idx])
            self.r_prime = []

            # # init termP1
            self.A1 = [ self.hasher.to_field(0) for _ in range(2 ** len(self.prev_subAC)) ]
            self.A1[0] = (-self.prev_subAC[-1]).add(1)
            self.A1[1] = self.prev_subAC[-1]

            for l in range(len(self.prev_subAC) - 2, -1, -1):
                for k in range(2 ** (len(self.prev_subAC) - l - 1) - 1, -1, -1):
                    temp = self.A1[k]
                    self.A1[2 * k] = ((-self.prev_subAC[l]).add(1)) * temp
                    self.A1[2 * k + 1] = self.prev_subAC[l] * temp

            # # init termP2
            self.A2 = [ self.hasher.to_field(0) for _ in range(2 ** len(self.prev_gate)) ]
            self.A2[0] = (-self.prev_gate[-1]).add(1)
            self.A2[1] = self.prev_gate[-1]

            for l in range(len(self.prev_gate) - 2, -1, -1):
                for k in range(2 ** (len(self.prev_gate) - l - 1) - 1, -1, -1):
                    temp = self.A2[k]
                    self.A2[2 * k] = ((-self.prev_gate[l]).add(1)) * temp
                    self.A2[2 * k + 1] = self.prev_gate[l] * temp
        
        def init_phase2(self, layer_idx : int):
            layer = self.circuit.layer(layer_idx)
            eq_val = self.hasher.to_field(1)
            for l in range(len(self.prev_subAC)):
                q_bit = self.prev_subAC[l]
                r_bit = self.r_prime[l]
                term1 = q_bit * r_bit
                term2 = (self.hasher.to_field(1) - q_bit) * (self.hasher.to_field(1) - r_bit)
                eq_val *= (term1 + term2)
            
            self.U = [eq_val * self.A2[g] for g in range(layer.length())]
            self.current_W_L = [self.W[i][0] for i in range(len(self.W))]
            self.current_W_R = copy.deepcopy(self.current_W_L)
            self.V_r_prime = copy.deepcopy(self.current_W_L)
            self.r0 = []
            self.r1 = []
        
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
        
        def get_coefficients_d2(self, evals):
            """
            evals: [F(0), F(1), F(-1)] 순서로 담긴 Field 객체들의 리스트
            반환값: 2차 다항식의 계수 [c, b, a] (Field 객체들의 리스트, F(x) = ax^2 + bx + c)
            """
            # 0, 1, -1 순서의 평가값 언패킹
            y_0, y_1, y_m1 = evals
            mod = y_0.mod
            
            # 2의 유한체 내 곱셈 역원 계산 (mod가 소수임을 가정)
            inv2 = pow(2, mod - 2, mod)
            
            # 1. c = F(0) = y_0 (상수항)
            c = y_0.copy()
            
            # 2. a = (F(1) + F(-1) - 2*F(0)) / 2 (2차항 계수)
            # F(1) + F(-1) = (a + b + c) + (a - b + c) = 2a + 2c
            a_numerator = y_1 + y_m1 - y_0.mul(2)
            a = a_numerator.mul(inv2)
            
            # 3. b = (F(1) - F(-1)) / 2 (1차항 계수)
            # F(1) - F(-1) = (a + b + c) - (a - b + c) = 2b
            b_numerator = y_1 - y_m1
            b = b_numerator.mul(inv2)
            
            return [c, b, a]

        def get_coefficients_general(self, evals):
            """
            evals: [F(0), F(1), ..., F(n-1)] 순서로 담긴 Field 객체들의 리스트
            반환값: n-1차 다항식의 계수 [c_0, c_1, ..., c_{n-1}] 
                (Field 객체들의 리스트, F(x) = c_0 + c_1*x + ... + c_{n-1}*x^{n-1})
            """
            n = len(evals)
            if n == 0:
                return []
                
            mod = evals[0].mod
            
            # 1. 평가점(xs) 생성: 0부터 n-1까지
            xs = list(range(n))
                    
            # 2. V: Vandermonde 행렬 생성 (n x n, 파이썬 정수로 계산)
            # V[i][j] = (xs[i] ** j) % mod
            V = [[pow(x, j, mod) if j > 0 else 1 for j in range(n)] for x in xs]
            
            # 3. Y: evals의 복사본 (원본 데이터 보호 및 연산 결과 저장용)
            Y = [y.copy() for y in evals]
            
            # 4. 가우스 소거법 (Forward Elimination)
            for i in range(n):
                # Pivot 찾기
                pivot_row = i
                while pivot_row < n and V[pivot_row][i] == 0:
                    pivot_row += 1
                    
                if pivot_row == n:
                    raise ValueError("평가점이 중복되었거나 행렬식이 0이어서 계수를 구할 수 없습니다.")
                    
                # Pivot 교환
                if pivot_row != i:
                    V[i], V[pivot_row] = V[pivot_row], V[i]
                    Y[i], Y[pivot_row] = Y[pivot_row], Y[i]
                    
                # Pivot을 1로 만들기 (모듈로 역원 곱셈)
                pivot_val = V[i][i]
                inv_pivot = pow(pivot_val, mod - 2, mod) 
                
                for j in range(i, n):
                    V[i][j] = (V[i][j] * inv_pivot) % mod
                Y[i] = Y[i].mul(inv_pivot) 
                
                # 아래 행 소거
                for k in range(i + 1, n):
                    factor = V[k][i]
                    if factor != 0:
                        for j in range(i, n):
                            V[k][j] = (V[k][j] - factor * V[i][j]) % mod
                        Y[k] = Y[k] - Y[i].mul(factor)
                        
            # 5. 후진 대입법 (Back Substitution)
            for i in range(n - 1, -1, -1):
                for k in range(i - 1, -1, -1):
                    factor = V[k][i]
                    if factor != 0:
                        V[k][i] = 0
                        Y[k] = Y[k] - Y[i].mul(factor)
                        
            return Y

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
        
        def phase2(self, layer_idx : int, round : int) -> tuple[list[Field], list[Field]]:
            ret = [ self.hasher.to_field(0) for _ in range(3) ]
            layer = self.circuit.layer(layer_idx)
            if layer_idx == 0:
                next_layer_length = layer.length() * 2
            else:
                next_layer_length = self.circuit.layer(layer_idx - 1).length()
            b_G = math.ceil(math.log2(next_layer_length))
            # b_G = int(math.log2(layer.length()))
            for gate_i, gate in enumerate(layer.gates):
                U_val = self.U[gate_i]
                if round < b_G:
                    round_j = round
                    # termL
                    rem_L = gate.left >> (round_j + 1)
                    termL_0 = self.current_W_L[2 * rem_L]
                    termL_1 = self.current_W_L[2 * rem_L + 1]
                    termL_m1 = termL_0.mul(2) - termL_1
                    # termR
                    termR_val = self.V_r_prime[gate.right]
                    termR_0 = termR_val
                    termR_1 = termR_val
                    termR_m1 = termR_val
                    # termP
                    b_L = (gate.left >> round_j) & 1
                    if b_L == 0:
                        termP_0 = U_val
                        termP_1 = self.hasher.to_field(0)
                        termP_m1 = U_val.mul(2)
                    else:
                        termP_0 = self.hasher.to_field(0)
                        termP_1 = U_val
                        termP_m1 = self.hasher.to_field(0) - U_val
                else:
                    round_j = round - b_G
                    # termL
                    termL_val = self.final_W_L[0]
                    termL_0 = termL_val
                    termL_1 = termL_val
                    termL_m1 = termL_val
                    # termR
                    rem_R = gate.right >> (round_j + 1)
                    termR_0 = self.current_W_R[2 * rem_R]
                    termR_1 = self.current_W_R[2 * rem_R + 1]
                    termR_m1 = termR_0.mul(2) - termR_1
                    # termP
                    b_R = (gate.right >> round_j) & 1
                    if b_R == 0:
                        termP_0 = U_val
                        termP_1 = self.hasher.to_field(0)
                        termP_m1 = U_val.mul(2)
                    else:
                        termP_0 = self.hasher.to_field(0)
                        termP_1 = U_val
                        termP_m1 = self.hasher.to_field(0) - U_val
                termP = [termP_0, termP_1, termP_m1]
                termL = [termL_0, termL_1, termL_m1]
                termR = [termR_0, termR_1, termR_m1]
                if gate.op == OpType.ADD:
                    for k_idx in range(3):
                        ret[k_idx] += termP[k_idx] * (termL[k_idx] + termR[k_idx])
                elif gate.op == OpType.MULT:
                    for k_idx in range(3):
                        ret[k_idx] += termP[k_idx] * (termL[k_idx] * termR[k_idx])
            coeff = self.get_coefficients_d2(ret) 
            return ret, coeff
        
        def collapse_phase1(self, r : Field):
            self.r_prime.append(r)
            # termL, R
            for j, sub_val in enumerate(self.W):
                for i in range(len(sub_val) // 2):
                    t0, t1 = sub_val[2 * i], sub_val[2 * i + 1]
                    sub_val[i] = (-r).add(1) * t0 + r * t1
                self.W[j] = self.W[j][:len(sub_val) // 2]
            # print(len(self.W), len(self.W[0]))

            # term P1
            for i in range(len(self.A1) // 2):
                t0, t1 = self.A1[2 * i], self.A1[2 * i + 1]
                self.A1[i] = (-r).add(1) * t0 + r * t1
            self.A1 = self.A1[:len(self.A1) // 2]
        
        def collapse_phase2(self, layer_idx : int, round : int, r : Field):
            layer = self.circuit.layer(layer_idx)
            # b_G = math.ceil(math.log2(layer.length()))
            # if b_G == 0:
            #     b_G = 1
            if layer_idx == 0:
                next_layer_length = layer.length() * 2
            else:
                next_layer_length = self.circuit.layer(layer_idx - 1).length()
            b_G = math.ceil(math.log2(next_layer_length))
            
            one_minus_r_j = self.hasher.to_field(1) - r
            
            if round < b_G:
                self.r0.append(r)
            else:
                self.r1.append(r)

            # U update
            for gate_i, gate in enumerate(layer.gates):
                if round < b_G:
                    
                    round_j = round
                    b = (gate.left >> round_j) & 1
                else:
                    round_j = round - b_G
                    b = (gate.right >> round_j) & 1
                if b == 0:
                    self.U[gate_i] *= one_minus_r_j
                else:
                    self.U[gate_i] *= r
            
            # collapse
            if round < b_G:
                self.current_W_L = self._collapse_1d(self.current_W_L, r)
                if round == b_G - 1:
                    self.final_W_L = self.current_W_L
            else:
                self.current_W_R = self._collapse_1d(self.current_W_R, r)
            
        def _collapse_1d(self, array: list[Field], r: Field) -> list[Field]:
            half_len = len(array) // 2
            new_array = [self.hasher.to_field(0) for _ in range(half_len)]
            one_minus_r = self.hasher.to_field(1) - r
            for i in range(half_len):
                new_array[i] = one_minus_r * array[2 * i] + r * array[2 * i + 1]
            return new_array
        
        def end_sumcheck(self, layer_idx : int) -> list[Field]:
            layer = self.circuit.layer(layer_idx)
            # b_G = max(math.ceil(math.log2(layer.length())), 1)
            # b_G = math.ceil(math.log2(layer.length()))
            if layer_idx == 0:
                next_layer_length = layer.length() * 2
            else:
                next_layer_length = self.circuit.layer(layer_idx - 1).length()
            b_G = math.ceil(math.log2(next_layer_length))
            eval_H = []
            for l in range(b_G + 1):
                l_field = self.hasher.to_field(l)
                target_r = []
                for k in range(b_G):
                    diff = self.r1[k] - self.r0[k]
                    val = diff * l_field + self.r0[k]
                    target_r.append(val)
                eval_val = copy.deepcopy(self.V_r_prime)
                for r_k in target_r:
                    eval_val = self._collapse_1d(eval_val, r_k)
                eval_H.append(eval_val[0])
            coeff_H = self.get_coefficients_general(eval_H)
            return coeff_H

        def init_sumcheck(self, tau : Field):
            q_i = []
            for k in range(len(self.r0)):
                diff = self.r1[k] - self.r0[k]
                q_i.append(diff * tau + self.r0[k])
            self.prev_gate = q_i
            self.prev_subAC = self.r_prime






    class Verifier:
        def __init__(self, circuit : Circuit, input : list[list[Field]], \
                     hasher : HomHash_Manager, claimed_output : list[Field]):
            self.circuit = circuit
            self.input = input
            self.hasher = hasher
            self.subAC_num = len(input[0])
            self.claimed_output = claimed_output

        def _collapse(self, vec : list[Field], len : int, r : Field):
            for sigma in range(len // 2):
                vec[sigma] = (-r).add(1) * vec[2 * sigma] + r * vec[2 * sigma + 1]
        
        # len(r) = log |vec|
        def _multi_collapse(self, vec : list[Field], r : list[Field]):
            loglen = math.ceil(math.log2(len(vec)))
            for j in range(loglen):
                self._collapse(vec, 2 ** (loglen - j), r[j])
            return vec[0]

        # return q' (sub-AC index)
        def init_proof(self) -> list[Field]:
            self.q = [self.hasher.to_field(0), self.hasher.to_field(0)]
            self.qp = [ self.hasher.sampling_field() for _ in range(math.ceil(math.log2(self.subAC_num))) ]
            temp = copy.deepcopy(self.claimed_output)
            for j in range(len(self.qp)):
                length = 2 ** (len(self.qp) - j)
                for i in range(length // 2):
                    temp[i] = (-self.qp[j]).add(1) * temp[2 * i] + self.qp[j] * temp[2 * i + 1]
            self.prev_claim = temp[0]
            # self.prev_claim = self._multi_collapse(temp, self.qp)
            return self.qp
        
        # degree 3
        def phase1(self, eval : list[Field], coeff : list[Field], round : int) -> Field:
            # integrity check
            if self.prev_claim != eval[0] + eval[1]:
                raise Exception("not valid")
            for i in range(-1, 3):
                if eval[i] != coeff[0] + coeff[1].mul(i) + coeff[2].mul(i ** 2) + coeff[3].mul(i ** 3):
                    raise Exception("not valid")
            rande = self.rp[round]
            self.prev_claim = coeff[0] + coeff[1] * rande + coeff[2] * rande * rande \
                              + coeff[3] * rande * rande * rande
            return rande
        
        def phase2(self, layer_idx : int, coeff : list[Field], round : int) -> Field:
            layer = self.circuit.layers[layer_idx]
            eval0 = coeff[0]
            eval1 = coeff[0] + coeff[1] + coeff[2]
            if self.prev_claim != eval0 + eval1:
                raise Exception("not valid")
            # b_G = math.ceil(math.log2(layer.length()))
            if layer_idx == 0:
                next_layer_length = layer.length() * 2
            else:
                next_layer_length = self.circuit.layer(layer_idx - 1).length()
            b_G = math.ceil(math.log2(next_layer_length))
            if b_G == 0:
                b_G = 1
            if round < b_G:
                rande = self.r0[round]
            else:
                rande = self.r1[round - b_G]
            self.prev_claim = coeff[0] + coeff[1] * rande + coeff[2] * rande * rande
            return rande
        
        # verifier precomputation
        def init_sumcheck(self, layer_idx : int):
            layer = self.circuit.layers[layer_idx]
            # next_layer = self.circuit.layers[layer_idx - 1]
            if layer_idx == 0:
                next_layer_length = layer.length() * 2
            else:
                next_layer_length = self.circuit.layer(layer_idx - 1).length()
            self.rp = [self.hasher.sampling_field_range(2) for _ in range(math.ceil(math.log2(self.subAC_num)))]
            if next_layer_length == 1:
                self.r0 = [self.hasher.sampling_field_range(2)]
                self.r1 = [self.hasher.sampling_field_range(2)]
            else:
                self.r0 = [self.hasher.sampling_field_range(2) for _ in range(math.ceil(math.log2(next_layer_length)))]
                self.r1 = [self.hasher.sampling_field_range(2) for _ in range(math.ceil(math.log2(next_layer_length)))]
            self.Aq = [None for _ in range(layer.length())]
            self.Aq[0] = (-self.q[-1]).add(1)
            if layer.length() > 1:
                self.Aq[1] = self.q[-1]
                for l in range(len(self.q) - 2, -1, -1):
                    for k in range(2 ** (len(self.q) - l - 1) - 1, -1, -1):
                        temp = self.Aq[k]
                        self.Aq[2 * k] = ((-self.q[l]).add(1)) * temp
                        self.Aq[2 * k + 1] = self.q[l] * temp
            self.Ar0 = [None for _ in range(next_layer_length)]
            self.Ar0[0] = (-self.r0[-1]).add(1)
            if next_layer_length > 1:
                self.Ar0[1] = self.r0[-1]
                for l in range(len(self.r0) - 2, -1, -1):
                    for k in range(2 ** (len(self.r0) - l - 1) - 1, -1, -1):
                        temp = self.Ar0[k]
                        self.Ar0[2 * k] = ((-self.r0[l]).add(1)) * temp
                        self.Ar0[2 * k + 1] = self.r0[l] * temp
            self.Ar1 = [None for _ in range(next_layer_length)]
            self.Ar1[0] = (-self.r1[-1]).add(1)
            if next_layer_length > 1:
                self.Ar1[1] = self.r1[-1]
                for l in range(len(self.r1) - 2, -1, -1):
                    for k in range(2 ** (len(self.r1) - l - 1) - 1, -1, -1):
                        temp = self.Ar1[k]
                        self.Ar1[2 * k] = ((-self.r1[l]).add(1)) * temp
                        self.Ar1[2 * k + 1] = self.r1[l] * temp
        
        # return tau
        def end_sumcheck(self, layer_idx : int, coeff_H : list[Field]):
            layer = self.circuit.layer(layer_idx)
            v0 = coeff_H[0]
            v1 = self.hasher.to_field(0)
            for c in coeff_H:
                v1 += c
            
            add_val = self.hasher.to_field(0)
            mult_val = self.hasher.to_field(0)

            for gate_i, gate in enumerate(layer.gates):
                term = self.Aq[gate_i] * self.Ar0[gate.left] * self.Ar1[gate.right]
                if gate.op == OpType.ADD:
                    add_val += term
                elif gate.op == OpType.MULT:
                    mult_val += term
            
            eq_val = self.hasher.to_field(1)
            for l in range(len(self.qp)):
                q_b = self.qp[l]
                r_b = self.rp[l]
                term1 = q_b * r_b
                term2 = (self.hasher.to_field(1) - q_b) * (self.hasher.to_field(1) - r_b)
                eq_val *= (term1 + term2)
            
            expected_e = eq_val * (add_val * (v0 + v1) + mult_val * (v0 * v1))
            print(self.prev_claim, expected_e)
            if self.prev_claim != expected_e:
                raise Exception("Sumcheck Reject")
            
            tau_i = self.hasher.sampling_field()
            a_i = self.hasher.to_field(0)
            for c in reversed(coeff_H):
                a_i = a_i * tau_i + c
            
            q_i = []
            for k in range(len(self.r0)):
                diff = self.r1[k] - self.r0[k]
                q_i.append(diff * tau_i + self.r0[k])
            
            self.prev_claim = a_i
            self.q = q_i
            self.qp = self.rp

            return tau_i