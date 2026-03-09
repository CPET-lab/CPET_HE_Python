import copy
import math
import time
from he.galois_ring.poly import Poly
from he.galois_ring.rns_poly import RNS_Poly
from he.he_parameter import HE_Parameter
from he.ciphertext import Ciphertext
from proof.circuit import *
from proof.cipher_hash import *

def _next_gate_len(circuit : Circuit, layer_idx : int, input) -> int:
    if layer_idx == 0:
        return len(input)
    return circuit.layer(layer_idx - 1).length()

def _calc_bg(circuit : Circuit, layer_idx : int, input) -> int:
    return _bit_len(_next_gate_len(circuit, layer_idx, input))

def _bit_len(i : int) -> int:
    return math.ceil(math.log2(i))

def _init_chi_vec(a : list[Field], label : list[Field]):
    a[0] = (-label[-1]).add(1)
    if len(a) > 1:
        a[1] = label[-1]
        for l in range(len(label) - 2, -1, -1):
            for k in range(2 ** (len(label) - l - 1) -1, -1, -1):
                temp = a[k]
                a[2 * k] = (-label[l]).add(1) * temp
                a[2 * k + 1] = label[l] * temp

# MLE of eq(q, r) = q * r + (1 - q) * (1 - r)
# if q == r == 1 then return 1
def _eq(q : Field, r : Field):
    mod = q.mod
    return q * r + (Field(1, mod) - q) * (Field(1, mod) - r)

# MLE of V_i(q', q)

def _collapse(vec : list[Field], len : int, r : Field) -> list[Field]:
    for sigma in range(len // 2):
        vec[sigma] = (-r).add(1) * vec[2 * sigma] + r * vec[2 * sigma + 1]

# len(r) = _bit_len(len(vec))
def _multi_collapse(vec : list[Field], r : list[Field]) -> list[Field]:
    for i, _r in enumerate(r):
        _collapse(vec, 2 ** (len(r) - i), _r)
    return vec[0]

class Prover:
    def __init__(self, circuit : Circuit, input : list[list[Field]], fielder : Fielder, debug=False):
        self.circuit = circuit
        self.input = input
        self.fielder = fielder
        self.debug = debug
        self.subAC_num = len(input[0])
        self.running_time = 0

    def log(self, *objects, sep=' ', end='\n'):
        if self.debug == True:
            print(*objects, sep, end)

    # input: f(0), f(1) / output f(2)
    def _eval_f2(self, f0 : Field, f1 : Field) -> Field:
        return f1.mul(2) - f0
    
    # input: f(0), f(1) / output f(-1)
    def _eval_fm1(self, f0 : Field, f1 : Field) -> Field:
        return f0.mul(2) - f1
        
    # eval circuit + generate witness
    def eval_circuit(self) -> list[Field]:
        self.witness = self.circuit.compute_list(self.input)  # [layer][gate][sub-AC]
        return self.witness[-1][0]

    def init_proof(self, init_gate : list[Field], init_subAC : list[Field]):
        start = time.perf_counter()
        self.prev_gate = init_gate
        self.prev_subAC = init_subAC
        end = time.perf_counter()
        self.running_time += (end - start)
    
    def init_phase1(self, layer_idx : int):
        start = time.perf_counter()

        self.W = copy.deepcopy(self.witness[layer_idx])
        self.rand_subAC = []

        # init termP1, P2
        self.A1 = [None for _ in range(self.subAC_num)]  # for termP1
        self.A2 = [None for _ in range(self.circuit.layer(layer_idx).length())]  # for termP2
        _init_chi_vec(self.A1, self.prev_subAC)
        _init_chi_vec(self.A2, self.prev_gate)

        end = time.perf_counter()
        self.running_time += (end - start)
    
    def init_phase2(self, layer_idx : int):
        start = time.perf_counter()

        eq_val = self.fielder.to_field(1)
        for i in range(len(self.prev_subAC)):
            eq_val *= _eq(self.prev_subAC[i], self.rand_subAC[i])
        
        # U: eq(q, r) * chi_j(q)
        self.U = [eq_val * self.A2[j] for j in range(self.circuit.layer(layer_idx).length())]
        self.W_backup = [self.W[j][0] for j in range(len(self.W))]
        self.wl, self.wr = copy.deepcopy(self.W_backup), copy.deepcopy(self.W_backup)
        self.r0, self.r1 = [], []

        end = time.perf_counter()
        self.running_time += (end - start)
    
    # return coeff (x^3, x^2, x^1, x^0)
    def phase1(self, layer_idx : int) -> list[Field]:
        start = time.perf_counter()

        eval = [self.fielder.to_field(0) for _ in range(4)]
        layer = self.circuit.layer(layer_idx)
        for gate_idx, gate in enumerate(layer.gates):
            for i in range(len(self.W[0]) // 2):
                termL = [
                    self.W[gate.left][2 * i],
                    self.W[gate.left][2 * i + 1]
                ]
                termL.append(self._eval_f2(termL[0], termL[1]))
                termL.append(self._eval_fm1(termL[0], termL[1]))
                termR = [
                    self.W[gate.right][2 * i],
                    self.W[gate.right][2 * i + 1]
                ]
                termR.append(self._eval_f2(termR[0], termR[1]))
                termR.append(self._eval_fm1(termR[0], termR[1]))

                termP1 = [
                    self.A1[2 * i],
                    self.A1[2 * i + 1]
                ]
                termP1.append(self._eval_f2(termP1[0], termP1[1]))
                termP1.append(self._eval_fm1(termP1[0], termP1[1]))

                for j in range(4):
                    if gate.op == OpType.ADD:
                        eval[j] += termP1[j] * self.A2[gate_idx] * (termL[j] + termR[j])
                    elif gate.op == OpType.MULT:
                        eval[j] += termP1[j] * self.A2[gate_idx] * (termL[j] * termR[j])
        ret =  self.fielder.get_coeff_d3(eval)

        end = time.perf_counter()
        self.running_time += (end - start)
        return ret
    
    def phase1_update(self, rand_req : Field):
        start = time.perf_counter()

        self.rand_subAC.append(rand_req)
        for i, _w in enumerate(self.W):
            _collapse(_w, len(_w), rand_req)
            self.W[i] = self.W[i][:len(_w) // 2]
        _collapse(self.A1, len(self.A1), rand_req)
        self.A1 = self.A1[:len(self.A1) // 2]

        end = time.perf_counter()
        self.running_time += (end - start)

    # return coeff (x^2, x^1, x^0)
    def phase2(self, layer_idx : int, round : int) -> list[Field]:
        start = time.perf_counter()

        eval = [self.fielder.to_field(0) for _ in range(3)]
        layer = self.circuit.layer(layer_idx)
        bg = _calc_bg(self.circuit, layer_idx, self.input)
        for gate_idx, gate in enumerate(layer.gates):
            if round < bg:
                rem = gate.left >> (round + 1)
                termL = [self.wl[2 * rem], self.wl[2 * rem + 1]]
                termL.append(self._eval_fm1(termL[0], termL[1]))
                termR = [self.W_backup[gate.right] for _ in range(3)]
                p_bit = (gate.left >> round) & 1
            else:
                round_p = round - bg
                termL = [self.fixed_left for _ in range(3)]
                rem = gate.right >> (round_p + 1)
                termR = [self.wr[2 * rem], self.wr[2 * rem + 1]]
                termR.append(self._eval_fm1(termR[0], termR[1]))
                p_bit = (gate.right >> round_p) & 1
            u_val = self.U[gate_idx]
            if p_bit == 0:
                termP = [u_val, self.fielder.to_field(0), u_val.mul(2)]
            else:
                termP = [self.fielder.to_field(0), u_val, self.fielder.to_field(0) - u_val]
            for i in range(3):
                if gate.op == OpType.ADD:
                    eval[i] += termP[i] * (termL[i] + termR[i])
                elif gate.op == OpType.MULT:
                    eval[i] += termP[i] * (termL[i] * termR[i])
        ret = self.fielder.get_coeff_d2(eval)

        end = time.perf_counter()
        self.running_time += (end - start)
        return ret
    
    # V_i = (1 - r) * V_i(0) + r * V_i(1)
    def phase2_update(self, layer_idx : int, round : int, rand_req : Field):
        start = time.perf_counter()

        bg = _calc_bg(self.circuit, layer_idx, self.input)
        if round < bg:
            self.r0.append(rand_req)
        else:
            self.r1.append(rand_req)
        
        # update U
        layer = self.circuit.layer(layer_idx)
        for gate_idx, gate in enumerate(layer.gates):
            if round < bg:
                p_bit = (gate.left >> round) & 1
            else:
                p_bit = (gate.right >> (round - bg)) & 1
            if p_bit == 0:
                self.U[gate_idx] *= self.fielder.to_field(1) - rand_req
            else:
                self.U[gate_idx] *= rand_req
        
        # collapse
        if round < bg:
            _collapse(self.wl, len(self.wl), rand_req)
            self.wl = self.wl[:len(self.wl) // 2]
            if round == bg - 1:
                self.fixed_left = self.wl[0]
        else:
            _collapse(self.wr, len(self.wr), rand_req)
            self.wr = self.wr[:len(self.wr) // 2]
        
        end = time.perf_counter()
        self.running_time += (end - start)
    
    # return layer function coeff
    def end_sumcheck(self, layer_idx : int) -> list[Field]:
        start = time.perf_counter()

        bg = _calc_bg(self.circuit, layer_idx, self.input)
        eval = []
        for l in range(bg + 1):
            new_rand = []
            for k in range(bg):
                new_rand.append((self.r1[k] - self.r0[k]) * self.fielder.to_field(l) + self.r0[k])
            eval_val = copy.deepcopy(self.W_backup)
            for r in new_rand:
                _collapse(eval_val, len(eval_val), r)
                eval_val = eval_val[:len(eval_val) // 2]
            eval.append(eval_val[0])
        ret = self.fielder.get_coefficients_general(eval)

        end = time.perf_counter()
        self.running_time += (end - start)
        return ret
    
    def init_sumcheck(self, tau : Field):
        start = time.perf_counter()

        self.prev_gate = []
        for k in range(len(self.r0)):
            self.prev_gate.append((self.r1[k] - self.r0[k]) * tau + self.r0[k])
        self.prev_subAC = self.rand_subAC

        end = time.perf_counter()
        self.running_time += (end - start)


class Verifier:
    def __init__(self, circuit : Circuit, input : list[list[Field]], fielder : Fielder, debug=False):
        self.circuit = circuit
        self.input = input
        self.fielder = fielder
        self.debug = debug
        self.subAC_num = len(input[0])
        self.running_time = 0

    def log(self, *objects, sep=' ', end='\n'):
        if self.debug == True:
            print(*objects, sep, end)
    
    # return F(val)
    def horner(self, coeff : list[Field], val : Field) -> Field:
        ret = self.fielder.to_field(0)
        for c in coeff:
            ret = ret * val + c
        return ret
    
    # return initial random gate, sub-AC index
    def init_proof(self, claimed_output : list[Field]) -> tuple[list[Field], list[Field]]:
        start = time.perf_counter()

        self.prev_gate = [self.fielder.to_field(0)]
        self.prev_subAC = [self.fielder.sampling_field() for _ in range(_bit_len(self.subAC_num))]
        self.prev_claim = _multi_collapse(claimed_output, self.prev_subAC)

        end = time.perf_counter()
        self.running_time += (end - start)
        return self.prev_gate, self.prev_subAC
    
    def init_sumcheck(self, layer_idx : int):
        start = time.perf_counter()

        gate_num = self.circuit.layer(layer_idx).length()
        next_gate_len = _next_gate_len(self.circuit, layer_idx, self.input)
        bg = _calc_bg(self.circuit, layer_idx, self.input)

        # rp: random sub-AC index
        # r0, r1: random left / right input index
        self.rp = [self.fielder.sampling_field() for _ in range(_bit_len(self.subAC_num))]
        self.r0 = [self.fielder.sampling_field() for _ in range(bg)]
        self.r1 = [self.fielder.sampling_field() for _ in range(bg)]

        # Aq, Ar0, Ar1: for eval gate location
        self.Aq = [None for _ in range(gate_num)]
        self.Ar0 = [None for _ in range(next_gate_len)]
        self.Ar1 = [None for _ in range(next_gate_len)]
        _init_chi_vec(self.Aq, self.prev_gate)
        _init_chi_vec(self.Ar0, self.r0)
        _init_chi_vec(self.Ar1, self.r1)

        end = time.perf_counter()
        self.running_time += (end - start)
    
    def phase1(self, coeff : list[Field], round : int) -> Field:
        start = time.perf_counter()

        if self.prev_claim != self.horner(coeff, self.fielder.to_field(0)) + self.horner(coeff, self.fielder.to_field(1)):
            raise Exception("not valid")
        rand_req = self.rp[round]
        self.prev_claim = self.horner(coeff, rand_req)

        end = time.perf_counter()
        self.running_time += (end - start)
        return rand_req
    
    def phase2(self, layer_idx : int, coeff : list[Field], round : int) -> Field:
        start = time.perf_counter()

        if self.prev_claim != self.horner(coeff, self.fielder.to_field(0)) + self.horner(coeff, self.fielder.to_field(1)):
            raise Exception("not valid")
        bg = _calc_bg(self.circuit, layer_idx, self.input)
        if round < bg:
            rand_req = self.r0[round]
        else:
            rand_req = self.r1[round - bg]
        self.prev_claim = self.horner(coeff, rand_req)

        end = time.perf_counter()
        self.running_time += (end - start)
        return rand_req

    # return tau
    def end_sumcheck(self, layer_idx : int, coeff : list[Field]) -> Field:
        start = time.perf_counter()

        layer = self.circuit.layer(layer_idx)
        add_val, mult_val = self.fielder.to_field(0), self.fielder.to_field(0)
        for gate_idx, gate in enumerate(layer.gates):
            term = self.Aq[gate_idx] * self.Ar0[gate.left] * self.Ar1[gate.right]
            if gate.op == OpType.ADD: add_val += term
            elif gate.op == OpType.MULT: mult_val += term
        
        eq_val = self.fielder.to_field(1)
        for l in range(len(self.prev_subAC)):
            term1 = self.prev_subAC[l] * self.rp[l]
            term2 = (self.fielder.to_field(1) - self.prev_subAC[l]) * (self.fielder.to_field(1) - self.rp[l])
            eq_val *= (term1 + term2)
        v0, v1 = coeff[-1], self.horner(coeff, self.fielder.to_field(1))
        calc_claim = eq_val * (add_val * (v0 + v1) + mult_val * (v0 * v1))
        if self.prev_claim != calc_claim:
            raise Exception("Sumcheck Reject")
        
        tau = self.fielder.sampling_field()
        self.prev_claim = self.horner(coeff, tau)
        self.prev_subAC = self.rp
        self.prev_gate = []
        for k in range(len(self.r0)):
            self.prev_gate.append((self.r1[k] - self.r0[k]) * tau + self.r0[k])
        
        end = time.perf_counter()
        self.running_time += (end - start)
        return tau


class Demo:

    def log(self, *objects, sep=' ', end='\n'):
        if self.debug == True:
            print(*objects, sep, end)
    
    # return claimed_output
    # if timecheck is true, (claimed_output, prover time, verifier time)
    def giraffe_basic(self, circuit : Circuit, input : list[list[Field]], fielder : Fielder, debug=False, timecheck=False):
        self.debug = debug
        subAC_num = len(input[0])
        self.log("giraffe basic")

        prover = Prover(circuit, input, fielder, debug)
        verifier = Verifier(circuit, input, fielder, debug)
        claimed_output = prover.eval_circuit()
        ret_output = copy.deepcopy(claimed_output)
        init_gate, init_subAC = verifier.init_proof(claimed_output)
        prover.init_proof(init_gate, init_subAC)

        # Layer Loop
        for layer_idx in reversed(range(circuit.length())):
            self.log(f"** layer {layer_idx} **")
            verifier.init_sumcheck(layer_idx)

            # Phase 1
            # F_j(k) = sum(n, g, gL, gR) termP1(j,n,k) * termP2(g) * OP(g, termL(j,n,gL,k), termR(j,n,gR,k))
            self.log("\nphase 1")
            prover.init_phase1(layer_idx)
            for round in range(_bit_len(subAC_num)):
                self.log(f"    Round 1-{round}")
                coeff = prover.phase1(layer_idx)
                rand_req = verifier.phase1(coeff, round)
                self.log(f"        random request: {rand_req}")
                prover.phase1_update(rand_req)
            
            # Phase 2
            # F_j(k) = sum(g, gL, gR)(termP(j, g, k) * OP_g(termL(j, gL, k), termR(j, gR, k)))
            self.log("\nphase 2")
            prover.init_phase2(layer_idx)
            bg = _calc_bg(circuit, layer_idx, input)
            for round in range(2 * bg):
                self.log(f"    Round 2-{round}")
                coeff = prover.phase2(layer_idx, round)
                rand_req = verifier.phase2(layer_idx, coeff, round)
                self.log(f"        random request: {rand_req}")
                prover.phase2_update(layer_idx, round, rand_req)
            
            # End Sumcheck
            coeff = prover.end_sumcheck(layer_idx)
            tau = verifier.end_sumcheck(layer_idx, coeff)
            self.log(f"tau: {tau}")
            prover.init_sumcheck(tau)
        
        self.log("Accept")

        if timecheck == False:
            return ret_output
        else:
            return ret_output, prover.running_time, verifier.running_time
        
    def giraffe_cipher(self, circuit : Circuit, cipher_input : list[Ciphertext], param : HE_Parameter, debug=False):
        self.debug = debug

        # ciphertext -> Field list
        hasher = HomHash_Manager(param)
        fielder_dict = {key: Fielder(key) for key in param.coeff_modulus}
        input_dict = {key: [] for key in param.coeff_modulus}
        start = time.perf_counter()
        for _cipher in cipher_input:
            hashed_cipher = hasher.cipher_hash(_cipher)
            for _bais in param.coeff_modulus:
                input_dict[_bais].append(fielder_dict[_bais].to_field_list(hashed_cipher._rns_poly[_bais]._data))
        end = time.perf_counter()
        print(f"input hashing time: {end - start}")

        prover_time = 0
        verifier_time = 0
        output_dict = {}
        for _bais in param.coeff_modulus:
            to, tp, tv = self.giraffe_basic(circuit, input_dict[_bais], fielder_dict[_bais], debug, True)
            prover_time += tp
            verifier_time += tv
            output_dict[_bais] = to
        
        print(f"prover time: {prover_time}")
        print(f"verifier time: {verifier_time}")
        
        # check result
        start = time.perf_counter()
        eval_result = circuit.compute_poly(cipher_input)
        end = time.perf_counter()
        print(f"circuit evaluation time: {end - start}")
        start = time.perf_counter()
        hashed_result = hasher.cipher_hash(eval_result)
        end = time.perf_counter()
        print(f"result hashing time: {end - start}")
        is_same = True
        for _bais in param.coeff_modulus:
            hashed_result_field = fielder_dict[_bais].to_field_list(hashed_result._rns_poly[_bais]._data)
            same = all(hashed_result_field[i] == output_dict[_bais][i] for i in range(len(output_dict[_bais])))
            if not same: is_same = False
        self.log(is_same)