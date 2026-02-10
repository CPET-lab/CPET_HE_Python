import copy
from enum import Enum
from typing import Self
from he.galois_ring.rns_poly import RNS_Poly
from he.ciphertext import Ciphertext
from _util import _modulus
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
        
class Layer:
    def __init__(self, gates : list[Gate]):
        self.gates = gates
        self.max_idx = 1
        while self.max_idx < len(gates) + 1:
            self.max_idx *= 2
        while len(self.gates) < self.max_idx:
            self.gates.append(Gate(OpType.ADD, -1, -1))
        
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

class Circuit:
    def __init__(self, layers : list[Layer]):
        self.layers = layers
    
    def toString(self):
        ret = ""
        for idx, layer in enumerate(self.layers):
            ret += f"  ** Layer {idx} **\n"
            ret += layer.toString()
            ret += "\n"
        return ret
    
    def compute_poly(self, data : list) -> list:
        data1 = copy.deepcopy(data)
        data2 = []
        for layer in self.layers:
            data2 = copy.deepcopy(data1)
            data1 = layer.compute_poly(data2)
        return data1
    
    def compute_int(self, mod : int, data : list[int]) -> list[int]:
        data1 = copy.deepcopy(data)
        data2 = []
        for layer in self.layers:
            data2 = copy.deepcopy(data1)
            data1 = layer.compute_int(mod, data2)
        return data1
            

# str : '(', ')', '+', '*', input index (0 ~ i)
# ex) "0*1+(2+3*4)+5*6"
def parse_circuit(expression: str) -> Circuit:
    
    # 내부용 DAG 노드 클래스
    class DagNode:
        def __init__(self, node_type, value=None, left=None, right=None):
            self.type = node_type  # 'INPUT', 'OP'
            self.value = value     # 'ADD', 'MULT' or Input Index (int)
            self.left = left       # Left Child Node
            self.right = right     # Right Child Node
            self.depth = 0         # Computation Depth
            self.id = id(self)     # Unique ID for tracking
            
        def __repr__(self):
            return f"Node({self.type}, {self.value}, d={self.depth})"

    # 1. Tokenize (제공해주신 로직 그대로 사용)
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

    for token in tokens:
        if isinstance(token, int):  # Operand (Input Index)
            output_queue.append(DagNode('INPUT', value=token))
        elif token == '(':
            operator_stack.append(token)
        elif token == ')':
            while operator_stack and operator_stack[-1] != '(':
                output_queue.append(operator_stack.pop())
            operator_stack.pop()  # Pop '('
        else:  # Operators + or *
            while (operator_stack and operator_stack[-1] != '(' and
                   precedence.get(operator_stack[-1], 0) >= precedence[token]):
                output_queue.append(operator_stack.pop())
            operator_stack.append(token)
    
    while operator_stack:
        output_queue.append(operator_stack.pop())

    # 3. Build DAG (Calculate Depth)
    evaluation_stack = []
    all_nodes = [] # To keep track of all nodes for usage analysis

    # 입력 문자열에 등장하는 최대 인덱스를 찾아 초기 입력 크기 결정
    max_input_idx = 0
    for token in tokens:
        if isinstance(token, int):
            max_input_idx = max(max_input_idx, token)
    
    # 초기 입력 노드들 미리 생성 (입력 인덱스 0 ~ max_input_idx)
    input_nodes = [DagNode('INPUT', value=i) for i in range(max_input_idx + 1)]
    all_nodes.extend(input_nodes)

    for item in output_queue:
        if isinstance(item, DagNode): # It's an input operand from queue
            # 큐에 있는건 껍데기일 수 있으므로 실제 관리하는 input_nodes에서 가져옴
            evaluation_stack.append(input_nodes[item.value])
        else: # Operator string
            right_node = evaluation_stack.pop()
            left_node = evaluation_stack.pop()
            
            new_node = DagNode('OP', value=(OpType.ADD if item == '+' else OpType.MULT), left=left_node, right=right_node)
            new_node.depth = max(left_node.depth, right_node.depth) + 1
            
            evaluation_stack.append(new_node)
            all_nodes.append(new_node)

    final_output_node = evaluation_stack.pop()
    max_circuit_depth = final_output_node.depth

    # 4. DAG -> Leveled Circuit Construction
    
    # 어떤 노드가 몇 깊이 이상의 노드에서 사용되는지(부모가 누군지) 파악
    # (Passthrough가 필요한지 판단하기 위함)
    node_usage_max_depth = {node.id: 0 for node in all_nodes}
    
    def mark_usage(node, user_depth):
        if node is None: return
        node_usage_max_depth[node.id] = max(node_usage_max_depth[node.id], user_depth)
        # 재귀적으로 내려갈 필요 없음, 상향식 구성 시 부모가 자식을 참조하므로 
        # 자식 입장에서 자신의 부모 중 가장 깊은 depth를 알면 됨.
        # DAG 구성 시 부모->자식 링크만 있으므로, 
        # 전체 노드를 순회하며 자신의 자식들에게 내 depth를 알려주는 방식이 효율적.

    for node in all_nodes:
        if node.type == 'OP':
            mark_usage(node.left, node.depth)
            mark_usage(node.right, node.depth)
            
    # 최종 결과 노드는 회로의 끝까지 살아남아야 함
    mark_usage(final_output_node, max_circuit_depth)

    layers = []
    
    # 현재 레이어(이전 레이어의 출력)에 존재하는 노드들과 그 인덱스 매핑
    # 초기 상태: 입력 데이터들
    current_nodes = input_nodes 
    # current_nodes 리스트의 i번째 원소가 실제 데이터 리스트의 i번째 값에 해당함

    for d in range(1, max_circuit_depth + 1):
        gates = []
        next_nodes = []
        
        # 현재 처리해야 할 노드(Depth == d)와 전달해야 할 노드(Depth < d, but used later) 식별
        
        # 1. 계산 게이트 생성 (Depth가 현재 d인 노드들)
        # 이 노드들은 반드시 current_nodes(이전 레이어 출력)에 자식들이 존재함
        nodes_at_depth = [n for n in all_nodes if n.depth == d]
        
        # 2. 패스스루 게이트 생성 (Depth < d 이지만, 미래(depth > d)에 사용되는 노드들)
        # 이 노드들은 현재 current_nodes에 존재해야 함
        pass_nodes = []
        for n in current_nodes:
            # 내 depth는 이미 d보다 작음. 
            # 내가 사용되는 최대 깊이가 d보다 크다면 다음 레이어로 넘겨야 함
            if node_usage_max_depth[n.id] > d:
                pass_nodes.append(n)
            # 만약 내가 최종 결과 노드이고 현재 depth에 도달하지 못했다면 넘겨야 함 (예: 입력이 바로 출력인 경우 등)
            elif n == final_output_node and d < max_circuit_depth:
                 pass_nodes.append(n)

        # 게이트 구성 및 다음 레이어 노드 리스트 작성
        
        # 1) 연산 게이트 추가
        for node in nodes_at_depth:
            # current_nodes에서 자식 노드의 인덱스를 찾음
            try:
                l_idx = current_nodes.index(node.left)
                r_idx = current_nodes.index(node.right)
                gates.append(Gate(node.value, l_idx, r_idx))
                next_nodes.append(node)
            except ValueError:
                raise Exception(f"Circuit construction failed: Dependencies for node {node} not found in layer {d-1}")

        # 2) 패스스루(Identity) 게이트 추가
        # Gate(ADD, idx, -1) -> data[idx] + 0 (값 복사)
        for node in pass_nodes:
            idx = current_nodes.index(node)
            gates.append(Gate(OpType.ADD, idx, -1))
            next_nodes.append(node)

        layers.append(Layer(gates))
        current_nodes = next_nodes

    return Circuit(layers)