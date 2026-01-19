from enum import Enum
from typing import Self
from he.ciphertext import Ciphertext

class OpType(Enum):
    ADD = 0
    MUL = 1

class Gate:
    def __init__(self, op : OpType, left : int, right : int):
        self.op = op
        self.left = left
        self.right = right
    
    def toString(self):
        ret = ""
        if self.op == OpType.ADD:
            ret += "ADD"
        else:
            ret += "MUL"
        ret += f" {self.left}, {self.right}"
        return ret

class Layer:
    # gates[-1] = 0, gates[-2] = 1
    # gates[-1] : op = MUL, left = 0, right = 1 (0)
    # gates[-2] : op = ADD, left = 0, right = 1 (1)
    def __init__(self, gates : list[Gate]):
        self.gates = gates + [Gate(OpType.ADD, -1, -2), Gate(OpType.MUL, -1, -2)]
        self.max_idx = 1
        while self.max_idx < len(gates):
            self.max_idx *= 2
    
    def compute_layer(self, prev_data : list[Ciphertext], zero_cipher : Ciphertext) -> list[Ciphertext]:
        ret = []
        for gate in self.gates:
            if gate.op == OpType.ADD:
                # print("add")
                ret.append((prev_data[gate.left] + prev_data[gate.right]) * zero_cipher)
            elif gate.op == OpType.MUL:
                # print("mul")
                ret.append(prev_data[gate.left] * prev_data[gate.right])
        return ret
    
    def toString(self):
        ret = ""
        for gate in self.gates:
            ret += gate.toString() + "\n"
        return ret

class Circuit:
    # layers[0] : output layer
    def __init__(self, layers : list[Layer], enc1 : Ciphertext):
        self.layers = layers
        self.one_ciphers = [enc1.copy()]
        for _ in range(len(layers) - 1):
            self.one_ciphers.append(self.one_ciphers[-1] * self.one_ciphers[-1])
    
    # input[-1] = Enc(0), input[-2] = Enc(1)
    def compute_circuit(self, input : list[Ciphertext]) -> Ciphertext:
        layer_output = input
        for i, layer in enumerate(reversed(self.layers)):
            layer_output = layer.compute_layer(layer_output, self.one_ciphers[i])
        return layer_output[0]