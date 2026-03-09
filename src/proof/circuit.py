import copy
import heapq
import itertools
from enum import Enum
from typing import Self
from he.galois_ring.rns_poly import RNS_Poly
from he.ciphertext import Ciphertext
from _util import _modulus
from _util import _random
from _util._field import *
from collections import deque

class OpType(Enum):
    ADD = 0
    MULT = 1

class Gate:
    def __init__(self, op : OpType, left : int, right : int):
        self.op = op
        if left == -1:
            self.left = right
            self.right = -1
        elif right == -1:
            self.left = left
            self.right = -1
        else:
            self.left = min(left, right)
            self.right = max(left, right)
    
    def toString(self):
        ret = ""
        if self.op == OpType.ADD:
            ret += "ADD"
        else:
            ret += "MULT"
        ret += f" [{self.left}], [{self.right}]"
        return ret
    
    # poly, rns_poly, ciphertext...
    def compute_poly(self, data : list):
        if len(data) < self.right:
            raise Exception("data list is to small")
        if self.op == OpType.ADD:
            return data[self.left] + data[self.right]
        elif self.op == OpType.MULT:
            return data[self.left] * data[self.right]
        
    def compute_int(self, mod : int, data : list[int]) -> int:
        if len(data) < self.right:
            raise Exception("data list is to small")
        if self.op == OpType.ADD:
            return _modulus._centered_modulus(data[self.left] + data[self.right], mod)
        elif self.op == OpType.MULT:
            return _modulus._centered_modulus(data[self.left] * data[self.right], mod)
    
    def compute_field(self, data : list[Field]) -> Field:
        if len(data) < self.right:
            raise Exception("data list is to small")
        if self.op == OpType.ADD:
            return data[self.left] + data[self.right]
        elif self.op == OpType.MULT:
            return data[self.left] * data[self.right]
        
    def compute_list(self, data : list[list[Field]]) -> list[Field]:
        if len(data) < self.right:
            raise Exception("data list is to small")
        ret = []
        if self.op == OpType.ADD:
            for i in range(len(data[0])):
                ret.append(data[self.left][i] + data[self.right][i])
        elif self.op == OpType.MULT:
            for i in range(len(data[0])):
                ret.append(data[self.left][i] * data[self.right][i])
        return ret
            
class Layer:
    def __init__(self, gates : list[Gate]):
        self.gates = gates
        self.max_idx = 1
        while self.max_idx < len(gates):
            self.max_idx *= 2
        while len(self.gates) < self.max_idx:
            self.gates.append(Gate(OpType.ADD, -1, -1))

    # def __init__(self, gates : list[Gate]):
    #     self.gates = gates
    #     self.max_idx = 1
    #     while self.max_idx < len(gates) + 1:
    #         self.max_idx *= 2
    #     while len(self.gates) < self.max_idx:
    #         self.gates.append(Gate(OpType.ADD, -1, -1))
        
    def toString(self):
        ret = ""
        for idx, gate in enumerate(self.gates):
            ret += f"{idx} gate : {gate.toString()}\n"
        return ret
    
    def compute_poly(self, data : list) -> list:
        data.append(data[0] - data[0])  # zero poly
        ret = []
        for gate in self.gates:
            ret.append(gate.compute_poly(data))
        return ret
    
    def compute_int(self, mod : int, data : list[int]) -> list[int]:
        data.append(0)
        ret = []
        for gate in self.gates:
            ret.append(gate.compute_int(mod, data))
        return ret

    def compute_field(self, data : list[Field]) -> list[Field]:
        data.append(data[0] - data[0])
        ret = []
        for gate in self.gates:
            ret.append(gate.compute_field(data))
        return ret
    
    def compute_list(self, data : list[list[Field]]) -> list[list[Field]]:
        ret = []
        for gate in self.gates:
            ret.append(gate.compute_list(data))
        return ret

    def length(self) -> int:
        return len(self.gates)

class Circuit:
    # layer 0 - input layer
    def __init__(self, layers : list[Layer]):
        self.layers = layers
    
    def toString(self):
        ret = ""
        for idx, layer in enumerate(self.layers):
            ret += f"  ** Layer {idx} **\n"
            ret += layer.toString()
            ret += "\n"
        return ret
    
    # input : list of polynomial (plain, cipher...)
    # output : list of layer's output (output[layer][gates]), each elemenent is poly
    def compute_poly(self, data : list):
        self._modify_label_m1(data)
        data_copy = copy.deepcopy(data)
        ret = [ copy.deepcopy(data_copy) ]
        for layer in self.layers:
            data_copy = layer.compute_poly(data_copy)
            ret.append(copy.deepcopy(data_copy))
        return ret[-1]
        
    def compute_int(self, mod : int, data : list[int]) -> list[int]:
        self._modify_label_m1(data)
        data1 = copy.deepcopy(data)
        for layer in self.layers:
            data1 = layer.compute_int(mod, data1)
        return data1
    
    def compute_field(self, data : list[Field]) -> list[Field]:
        self._modify_label_m1(data)
        data1 = copy.deepcopy(data)
        for layer in self.layers:
            data1 = layer.compute_field(data1)
        return data1
    
    def compute_list(self, data : list[list[Field]]) -> list[list[list[Field]]]:
        self._modify_label_m1(data)
        data_copy = copy.deepcopy(data)
        ret = [copy.deepcopy(data)]
        for layer in self.layers:
            data_copy = layer.compute_list(data_copy)
            ret.append(copy.deepcopy(data_copy))
        return ret
        
    def layer(self, layer_idx : int) -> Layer:
        return self.layers[layer_idx]
    
    def length(self) -> int:
        return len(self.layers)

    def _modify_label_m1(self, input):
        input_length = len(input)
        for gate in self.layer(0).gates:
            if gate.left == -1:
                gate.left = input_length - 1
            if gate.right == -1:
                gate.right = input_length - 1
        for layer_idx in range(1, self.length()):
            for gate in self.layer(layer_idx).gates:
                if gate.left == -1:
                    gate.left = self.layer(layer_idx - 1).length() - 1
                if gate.right == -1:
                    gate.right = self.layer(layer_idx - 1).length() - 1


# # str : '(', ')', '+', '*', input index (0 ~ i)
# # ex) "0*1+(2+3*4)+5*6"
# def parse_circuit(expression: str) -> Circuit:
    
#     # 내부용 DAG 노드 클래스
#     class DagNode:
#         def __init__(self, node_type, value=None, left=None, right=None):
#             self.type = node_type  # 'INPUT', 'OP'
#             self.value = value     # 'ADD', 'MULT' or Input Index (int)
#             self.left = left       # Left Child Node
#             self.right = right     # Right Child Node
#             self.depth = 0         # Computation Depth
#             self.id = id(self)     # Unique ID for tracking
            
#         def __repr__(self):
#             return f"Node({self.type}, {self.value}, d={self.depth})"

#     # 1. Tokenize (제공해주신 로직 그대로 사용)
#     def tokenize_expression(text):
#         operators = "()+*"
#         for op in operators:
#             text = text.replace(op, f" {op} ")
#         tokens = text.split()
#         result = []
#         for token in tokens:
#             if token.isdigit():
#                 result.append(int(token))
#             else:
#                 result.append(token)
#         return result

#     tokens = tokenize_expression(expression)

#     # 2. Shunting-yard Algorithm (Infix -> Postfix)
#     precedence = {'+': 1, '*': 2, '(': 0}
#     output_queue = deque()
#     operator_stack = []

#     for token in tokens:
#         if isinstance(token, int):  # Operand (Input Index)
#             output_queue.append(DagNode('INPUT', value=token))
#         elif token == '(':
#             operator_stack.append(token)
#         elif token == ')':
#             while operator_stack and operator_stack[-1] != '(':
#                 output_queue.append(operator_stack.pop())
#             operator_stack.pop()  # Pop '('
#         else:  # Operators + or *
#             while (operator_stack and operator_stack[-1] != '(' and
#                    precedence.get(operator_stack[-1], 0) >= precedence[token]):
#                 output_queue.append(operator_stack.pop())
#             operator_stack.append(token)
    
#     while operator_stack:
#         output_queue.append(operator_stack.pop())

#     # 3. Build DAG (Calculate Depth)
#     evaluation_stack = []
#     all_nodes = [] # To keep track of all nodes for usage analysis

#     # 입력 문자열에 등장하는 최대 인덱스를 찾아 초기 입력 크기 결정
#     max_input_idx = 0
#     for token in tokens:
#         if isinstance(token, int):
#             max_input_idx = max(max_input_idx, token)
    
#     # 초기 입력 노드들 미리 생성 (입력 인덱스 0 ~ max_input_idx)
#     input_nodes = [DagNode('INPUT', value=i) for i in range(max_input_idx + 1)]
#     all_nodes.extend(input_nodes)

#     for item in output_queue:
#         if isinstance(item, DagNode): # It's an input operand from queue
#             # 큐에 있는건 껍데기일 수 있으므로 실제 관리하는 input_nodes에서 가져옴
#             evaluation_stack.append(input_nodes[item.value])
#         else: # Operator string
#             right_node = evaluation_stack.pop()
#             left_node = evaluation_stack.pop()
            
#             new_node = DagNode('OP', value=(OpType.ADD if item == '+' else OpType.MULT), left=left_node, right=right_node)
#             new_node.depth = max(left_node.depth, right_node.depth) + 1
            
#             evaluation_stack.append(new_node)
#             all_nodes.append(new_node)

#     final_output_node = evaluation_stack.pop()
#     max_circuit_depth = final_output_node.depth

#     # 4. DAG -> Leveled Circuit Construction
    
#     # 어떤 노드가 몇 깊이 이상의 노드에서 사용되는지(부모가 누군지) 파악
#     # (Passthrough가 필요한지 판단하기 위함)
#     node_usage_max_depth = {node.id: 0 for node in all_nodes}
    
#     def mark_usage(node, user_depth):
#         if node is None: return
#         node_usage_max_depth[node.id] = max(node_usage_max_depth[node.id], user_depth)
#         # 재귀적으로 내려갈 필요 없음, 상향식 구성 시 부모가 자식을 참조하므로 
#         # 자식 입장에서 자신의 부모 중 가장 깊은 depth를 알면 됨.
#         # DAG 구성 시 부모->자식 링크만 있으므로, 
#         # 전체 노드를 순회하며 자신의 자식들에게 내 depth를 알려주는 방식이 효율적.

#     for node in all_nodes:
#         if node.type == 'OP':
#             mark_usage(node.left, node.depth)
#             mark_usage(node.right, node.depth)
            
#     # 최종 결과 노드는 회로의 끝까지 살아남아야 함
#     mark_usage(final_output_node, max_circuit_depth)

#     layers = []
    
#     # 현재 레이어(이전 레이어의 출력)에 존재하는 노드들과 그 인덱스 매핑
#     # 초기 상태: 입력 데이터들
#     current_nodes = input_nodes 
#     # current_nodes 리스트의 i번째 원소가 실제 데이터 리스트의 i번째 값에 해당함

#     for d in range(1, max_circuit_depth + 1):
#         gates = []
#         next_nodes = []
        
#         # 현재 처리해야 할 노드(Depth == d)와 전달해야 할 노드(Depth < d, but used later) 식별
        
#         # 1. 계산 게이트 생성 (Depth가 현재 d인 노드들)
#         # 이 노드들은 반드시 current_nodes(이전 레이어 출력)에 자식들이 존재함
#         nodes_at_depth = [n for n in all_nodes if n.depth == d]
        
#         # 2. 패스스루 게이트 생성 (Depth < d 이지만, 미래(depth > d)에 사용되는 노드들)
#         # 이 노드들은 현재 current_nodes에 존재해야 함
#         pass_nodes = []
#         for n in current_nodes:
#             # 내 depth는 이미 d보다 작음. 
#             # 내가 사용되는 최대 깊이가 d보다 크다면 다음 레이어로 넘겨야 함
#             if node_usage_max_depth[n.id] > d:
#                 pass_nodes.append(n)
#             # 만약 내가 최종 결과 노드이고 현재 depth에 도달하지 못했다면 넘겨야 함 (예: 입력이 바로 출력인 경우 등)
#             elif n == final_output_node and d < max_circuit_depth:
#                  pass_nodes.append(n)

#         # 게이트 구성 및 다음 레이어 노드 리스트 작성
        
#         # 1) 연산 게이트 추가
#         for node in nodes_at_depth:
#             # current_nodes에서 자식 노드의 인덱스를 찾음
#             try:
#                 l_idx = current_nodes.index(node.left)
#                 r_idx = current_nodes.index(node.right)
#                 gates.append(Gate(node.value, l_idx, r_idx))
#                 next_nodes.append(node)
#             except ValueError:
#                 raise Exception(f"Circuit construction failed: Dependencies for node {node} not found in layer {d-1}")

#         # 2) 패스스루(Identity) 게이트 추가
#         # Gate(ADD, idx, -1) -> data[idx] + 0 (값 복사)
#         for node in pass_nodes:
#             idx = current_nodes.index(node)
#             gates.append(Gate(OpType.ADD, idx, -1))
#             next_nodes.append(node)

#         layers.append(Layer(gates))
#         current_nodes = next_nodes

#     return Circuit(layers)

import heapq
import itertools
from collections import deque

# str : '(', ')', '+', '*', input index (0 ~ i)
# ex) "0*1+(2+3*4)+5*6"
def parse_circuit(expression: str) -> Circuit:
    
    class DagNode:
        def __init__(self, node_type, value=None, left=None, right=None):
            self.type = node_type  # 'INPUT', 'OP'
            self.value = value     # '+', '*' or Input Index (int)
            self.left = left       
            self.right = right     
            self.depth = 0         
            self.id = id(self)     
            
        def __repr__(self):
            return f"Node({self.type}, {self.value}, d={self.depth})"

    # 1. Tokenize
    def tokenize_expression(text):
        operators = "()+*"
        for op in operators:
            text = text.replace(op, f" {op} ")
        tokens = text.split()
        result = []
        for token in tokens:
            if token.isdigit():
                result.append(int(token))
            else:
                result.append(token)
        return result

    tokens = tokenize_expression(expression)

    # 2. Shunting-yard Algorithm (Infix -> Postfix)
    precedence = {'+': 1, '*': 2, '(': 0}
    output_queue = deque()
    operator_stack = []

    max_input_idx = 0
    for token in tokens:
        if isinstance(token, int): 
            output_queue.append(DagNode('INPUT', value=token))
            max_input_idx = max(max_input_idx, token)
        elif token == '(':
            operator_stack.append(token)
        elif token == ')':
            while operator_stack and operator_stack[-1] != '(':
                output_queue.append(operator_stack.pop())
            operator_stack.pop() 
        else: 
            while (operator_stack and operator_stack[-1] != '(' and
                   precedence.get(operator_stack[-1], 0) >= precedence[token]):
                output_queue.append(operator_stack.pop())
            operator_stack.append(token)
    
    while operator_stack:
        output_queue.append(operator_stack.pop())

    # 3. Build Initial DAG (Raw Tree)
    evaluation_stack = []
    input_nodes = [DagNode('INPUT', value=i) for i in range(max_input_idx + 1)]

    for item in output_queue:
        if isinstance(item, DagNode):
            evaluation_stack.append(input_nodes[item.value])
        else: 
            right_node = evaluation_stack.pop()
            left_node = evaluation_stack.pop()
            
            # 임시로 문자열 '+', '*'를 value로 저장 (최적화 편의를 위해)
            new_node = DagNode('OP', value=item, left=left_node, right=right_node)
            new_node.depth = max(left_node.depth, right_node.depth) + 1
            evaluation_stack.append(new_node)

    raw_root = evaluation_stack.pop()

    # =========================================================================
    # 🌟 NEW: 4. AST Balancing (Tree Optimization) 🌟
    # 연속된 동일 연산자(+, *)를 평탄화하고, 최소 Depth를 가지도록 트리를 재조립
    # =========================================================================
    def optimize_tree(node):
        if node.type == 'INPUT':
            return node
        
        # Bottom-up 방식으로 자식 노드부터 최적화
        node.left = optimize_tree(node.left)
        node.right = optimize_tree(node.right)
        
        # 같은 연산자로 묶인 모든 자식 노드들을 평탄화(Flatten)하여 수집
        def gather_same_ops(n, op_val):
            if n.type == 'OP' and n.value == op_val:
                return gather_same_ops(n.left, op_val) + gather_same_ops(n.right, op_val)
            else:
                return [n]
        
        operands = gather_same_ops(node, node.value)
        
        # 피연산자가 2개 이하면 재조립할 필요 없이 깊이만 갱신
        if len(operands) <= 2:
            node.depth = max(node.left.depth, node.right.depth) + 1
            return node
            
        # 피연산자가 3개 이상이면 허프만 코딩(Min-Heap) 방식으로 최소 깊이 트리 구성
        counter = itertools.count() # id 충돌 방지용 카운터
        pq = []
        for op_node in operands:
            heapq.heappush(pq, (op_node.depth, next(counter), op_node))
            
        while len(pq) > 1:
            depth1, _, left_n = heapq.heappop(pq)
            depth2, _, right_n = heapq.heappop(pq)
            
            # 깊이가 가장 얕은 두 노드를 합쳐서 새로운 부모 노드 생성
            new_node = DagNode('OP', value=node.value, left=left_n, right=right_n)
            new_node.depth = max(left_n.depth, right_n.depth) + 1
            
            heapq.heappush(pq, (new_node.depth, next(counter), new_node))
            
        _, _, optimized_root = heapq.heappop(pq)
        return optimized_root

    final_output_node = optimize_tree(raw_root)
    max_circuit_depth = final_output_node.depth

    # 최적화된 트리를 순회하며 all_nodes 리스트 재구성 (위상 정렬 순서 보장)
    all_nodes = list(input_nodes)
    visited_ids = set(n.id for n in all_nodes)
    
    def collect_nodes(n):
        if n is None or n.id in visited_ids:
            return
        collect_nodes(n.left)
        collect_nodes(n.right)
        visited_ids.add(n.id)
        all_nodes.append(n)
        
    collect_nodes(final_output_node)

    # =========================================================================
    # 5. DAG -> Leveled Circuit Construction (기존 로직과 동일)
    # =========================================================================
    node_usage_max_depth = {node.id: 0 for node in all_nodes}
    
    def mark_usage(node, user_depth):
        if node is None: return
        node_usage_max_depth[node.id] = max(node_usage_max_depth[node.id], user_depth)

    for node in all_nodes:
        if node.type == 'OP':
            mark_usage(node.left, node.depth)
            mark_usage(node.right, node.depth)
            
    mark_usage(final_output_node, max_circuit_depth)

    layers = []
    current_nodes = list(input_nodes) 

    for d in range(1, max_circuit_depth + 1):
        gates = []
        next_nodes = []
        
        nodes_at_depth = [n for n in all_nodes if n.depth == d]
        
        pass_nodes = []
        for n in current_nodes:
            if node_usage_max_depth[n.id] > d:
                pass_nodes.append(n)
            elif n == final_output_node and d < max_circuit_depth:
                 pass_nodes.append(n)

        # 1) 연산 게이트 추가
        for node in nodes_at_depth:
            try:
                l_idx = current_nodes.index(node.left)
                r_idx = current_nodes.index(node.right)
                # 이 단계에서 문자열 연산자를 OpType Enum으로 치환
                op_type = OpType.ADD if node.value == '+' else OpType.MULT
                gates.append(Gate(op_type, l_idx, r_idx))
                next_nodes.append(node)
            except ValueError:
                raise Exception(f"Circuit construction failed: Dependencies for node {node} not found in layer {d-1}")

        # 2) 패스스루(Identity) 게이트 추가
        for node in pass_nodes:
            idx = current_nodes.index(node)
            gates.append(Gate(OpType.ADD, idx, -1))
            next_nodes.append(node)

        # =====================================================================
        # 🌟 NEW: 3) 2의 거듭제곱으로 패딩(Padding) 추가 🌟
        # =====================================================================
        current_gate_count = len(gates)
        if current_gate_count > 0:
            # 현재 개수보다 크거나 같은 가장 작은 2의 거듭제곱 계산 (예: 3 -> 4, 5 -> 8)
            target_gate_count = 1 << (current_gate_count - 1).bit_length()
            padding_needed = target_gate_count - current_gate_count
            
            for _ in range(padding_needed):
                # 아무나 참조해도 상관없는 0번 인덱스에 0을 더하는 무의미한 연산 (Pass-through)
                gates.append(Gate(OpType.ADD, 0, -1))
                
                # 다음 레이어에서 인덱스가 꼬이지 않도록 자리만 차지하는 더미 노드 생성
                dummy_node = DagNode('DUMMY', value='PAD')
                # 이 노드는 미래에 절대 사용되지 않으므로 max_depth를 0으로 설정
                node_usage_max_depth[dummy_node.id] = 0 
                next_nodes.append(dummy_node)
        # =====================================================================

        layers.append(Layer(gates))
        current_nodes = next_nodes

    return Circuit(layers)

def build_poly_circuit(coeffs: list) -> Circuit:
    """
    coeffs: 오름차순 정렬된 다항식 계수 리스트 [c_0, c_1, ..., c_n]
    Circuit Input: [x, c_0, c_1, ..., c_n, 0] 
                   (길이 = n + 3)
    """
    class DagNode:
        def __init__(self, node_type, value=None, left=None, right=None):
            self.type = node_type  # 'INPUT', 'OP', 'DUMMY'
            self.value = value     
            self.left = left       
            self.right = right     
            self.depth = 0         
            self.id = id(self)     
            
        def __repr__(self):
            return f"Node({self.type}, {self.value}, d={self.depth})"

    n_coeffs = len(coeffs)
    zero_index = n_coeffs + 1  # 맨 마지막 '0'의 인덱스
    
    # =========================================================================
    # 1. 초기 입력 노드 매핑
    # =========================================================================
    input_nodes = [DagNode('INPUT', value=i) for i in range(zero_index + 1)]
    
    x_node = input_nodes[0]                 # 인덱스 0: 입력 데이터 x
    c_nodes = input_nodes[1:n_coeffs + 1]   # 인덱스 1 ~ n: 계수 c_0 ~ c_n
    zero_node = input_nodes[zero_index]     # 인덱스 n+1: 상수 0

    all_nodes = list(input_nodes)

    # =========================================================================
    # 2. x의 거듭제곱 트리 (x^2, x^3, ... x^n)
    # =========================================================================
    powers = {1: x_node}
    for i in range(2, n_coeffs):
        a = i // 2
        b = i - a
        node = DagNode('OP', value='*', left=powers[a], right=powers[b])
        node.depth = max(powers[a].depth, powers[b].depth) + 1
        powers[i] = node
        all_nodes.append(node)

    # =========================================================================
    # 3. 각 항 계산 (c_i * x^i)
    # =========================================================================
    terms = []
    for i in range(n_coeffs):
        if coeffs[i] == 0: 
            continue # 계수가 0인 항은 회로 연산에서 아예 제외
        
        if i == 0:
            # c_0 항은 곱셈 없이 계수 입력 노드 자체를 항으로 취급 (Depth = 0)
            terms.append(c_nodes[0])
        else:
            # c_i * x^i (두 입력 노드를 곱함)
            term_node = DagNode('OP', value='*', left=powers[i], right=c_nodes[i])
            term_node.depth = max(powers[i].depth, c_nodes[i].depth) + 1
            terms.append(term_node)
            all_nodes.append(term_node)

    # 모든 계수가 0일 경우 (결과가 무조건 0)
    if not terms:
        return Circuit([Layer([Gate(OpType.ADD, zero_index, -1)])])

    # =========================================================================
    # 4. 항들을 모두 더하기 (Min-Heap 덧셈 트리 밸런싱)
    # =========================================================================
    counter = itertools.count()
    pq = []
    for t in terms:
        heapq.heappush(pq, (t.depth, next(counter), t))
        
    while len(pq) > 1:
        d1, _, left_n = heapq.heappop(pq)
        d2, _, right_n = heapq.heappop(pq)
        
        new_node = DagNode('OP', value='+', left=left_n, right=right_n)
        new_node.depth = max(left_n.depth, right_n.depth) + 1
        all_nodes.append(new_node)
        heapq.heappush(pq, (new_node.depth, next(counter), new_node))
        
    _, _, final_output_node = heapq.heappop(pq)
    max_circuit_depth = final_output_node.depth

    # =========================================================================
    # 5. DAG -> Leveled Circuit 변환 및 패딩(Padding)
    # =========================================================================
    node_usage_max_depth = {n.id: 0 for n in all_nodes}
    
    def mark_usage(node, user_depth):
        if node is None: return
        node_usage_max_depth[node.id] = max(node_usage_max_depth[node.id], user_depth)

    for node in all_nodes:
        if node.type == 'OP':
            mark_usage(node.left, node.depth)
            mark_usage(node.right, node.depth)
            
    mark_usage(final_output_node, max_circuit_depth)
    
    # ★ 핵심 로직: 맨 끝의 '0' 노드는 매 레이어마다 패딩에 쓰여야 하므로 
    # 수명을 회로의 끝(max_circuit_depth)까지 강제로 연장합니다.
    node_usage_max_depth[zero_node.id] = max_circuit_depth

    layers = []
    current_nodes = list(input_nodes) 

    for d in range(1, max_circuit_depth + 1):
        gates = []
        next_nodes = []
        
        nodes_at_depth = [n for n in all_nodes if n.type == 'OP' and n.depth == d]
        
        pass_nodes = []
        for n in current_nodes:
            if node_usage_max_depth[n.id] > d:
                pass_nodes.append(n)
            elif n == final_output_node and d < max_circuit_depth:
                pass_nodes.append(n)

        # 1) 연산 게이트 추가
        for node in nodes_at_depth:
            l_idx = current_nodes.index(node.left)
            r_idx = current_nodes.index(node.right)
            op_type = OpType.ADD if node.value == '+' else OpType.MULT
            gates.append(Gate(op_type, l_idx, r_idx))
            next_nodes.append(node)

        # 2) 패스스루 게이트 추가
        for node in pass_nodes:
            idx = current_nodes.index(node)
            gates.append(Gate(OpType.ADD, idx, -1))
            next_nodes.append(node)

        # 3) 2의 거듭제곱으로 패딩 (안전한 0 노드 참조)
        current_gate_count = len(gates)
        if current_gate_count > 0:
            target_gate_count = 1 << (current_gate_count - 1).bit_length()
            padding_needed = target_gate_count - current_gate_count
            
            # 현재 레이어에서 zero_node가 위치한 인덱스를 찾음
            # (수명을 강제 연장했으므로 항상 존재함이 보장됨!)
            zero_idx = current_nodes.index(zero_node)
            
            for _ in range(padding_needed):
                gates.append(Gate(OpType.ADD, zero_idx, -1))
                
                dummy_node = DagNode('DUMMY', value='PAD')
                node_usage_max_depth[dummy_node.id] = 0 
                next_nodes.append(dummy_node)

        layers.append(Layer(gates))
        current_nodes = next_nodes

    return Circuit(layers)