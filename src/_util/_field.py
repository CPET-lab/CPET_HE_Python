from _util import _modulus
from _util import _random

class Field:
    def __init__(self, val, mod):
        self.mod = mod
        self.val = _modulus._centered_modulus(val, mod)
    
    def __add__(self, other):
        return Field(self.val + other.val, self.mod)
    
    def __sub__(self, other):
        return Field(self.val - other.val, self.mod)
    
    def __mul__(self, other):
        return Field(self.val * other.val, self.mod)
    
    def __neg__(self):
        return Field(-1 * self.val, self.mod)
    
    def __str__(self):
        return f"{self.val}"
    
    def __eq__(self, other):
        if self.mod != other.mod:
            raise Exception("asdf")
        if self.val == other.val:
            return True
        return False
    
    def __ne__(self, other):
        if self.mod != other.mod:
            raise Exception("asdf")
        if self.val != other.val:
            return True
        return False
    
    def add(self, other : int):
        return Field(self.val + other, self.mod)
    
    def sub(self, other : int):
        return Field(self.val - other, self.mod)
    
    def mul(self, other : int):
        return Field(self.val * other, self.mod)
    
    def eq(self, other : int):
        if self.val == other:
            return True
        return False
    
    def copy(self):
        return Field(self.val, self.mod)
    
class Fielder:
    def __init__(self, mod : int):
        self.mod = mod
    
    def to_field(self, val : int) -> Field:
        return Field(val, self.mod)
    
    def to_field_list(self, data : list[int]) -> list[Field]:
        return [self.to_field(d) for d in data]
    
    def sampling_field(self) -> Field:
        return Field(_random._random_int(self.mod), self.mod)
    
    def sampling_field_range(self, start : int, end = None) -> Field:
        if end == None:
            end = self.mod
        return Field(_random._random_int(end - start) + start, self.mod)
    
    # evals: [F(0), F(1), F(2), F(-1)]
    # return coeff of F (x^3, x^2, x^1, x^0)
    def get_coeff_d3(self, evals : list[Field]) -> list[Field]:
        y_0, y_1, y_2, y_m1 = evals
        inv6 = pow(6, self.mod - 2, self.mod)
        inv2 = pow(2, self.mod - 2, self.mod)
        d = y_0.copy()
        a_numerator = y_2 - y_1.mul(3) + y_0.mul(3) - y_m1
        a = a_numerator.mul(inv6)
        b_numerator = y_1 + y_m1 - y_0.mul(2)
        b = b_numerator.mul(inv2)
        c_numerator = y_1 - y_m1
        c = c_numerator.mul(inv2) - a
        return [a, b, c, d]
    
    # evals: [F(0), F(1), F(-1)]
    # return coeff of F (x^2, x^1, x^0)
    def get_coeff_d2(self, evals : list[Field]) -> list[Field]:
        y_0, y_1, y_m1 = evals
        inv2 = pow(2, self.mod - 2, self.mod)
        c = y_0.copy()
        a_numerator = y_1 + y_m1 - y_0.mul(2)
        a = a_numerator.mul(inv2)
        b_numerator = y_1 - y_m1
        b = b_numerator.mul(inv2)
        return [a, b, c]
    
    def get_coefficients_general(self, evals):
        """
        evals: [F(0), F(1), ..., F(n-1)] 순서로 담긴 Field 객체들의 리스트
        반환값: n-1차 다항식의 계수 [c_0, c_1, ..., c_{n-1}] 
            (Field 객체들의 리스트, F(x) = c_0 + c_1*x + ... + c_{n-1}*x^{n-1})
        """
        n = len(evals)
        if n == 0:
            return []
                    
        # 1. 평가점(xs) 생성: 0부터 n-1까지
        xs = list(range(n))
                
        # 2. V: Vandermonde 행렬 생성 (n x n, 파이썬 정수로 계산)
        # V[i][j] = (xs[i] ** j) % mod
        V = [[pow(x, j, self.mod) if j > 0 else 1 for j in range(n)] for x in xs]
        
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
            inv_pivot = pow(pivot_val, self.mod - 2, self.mod) 
            
            for j in range(i, n):
                V[i][j] = (V[i][j] * inv_pivot) % self.mod
            Y[i] = Y[i].mul(inv_pivot) 
            
            # 아래 행 소거
            for k in range(i + 1, n):
                factor = V[k][i]
                if factor != 0:
                    for j in range(i, n):
                        V[k][j] = (V[k][j] - factor * V[i][j]) % self.mod
                    Y[k] = Y[k] - Y[i].mul(factor)
                    
        # 5. 후진 대입법 (Back Substitution)
        for i in range(n - 1, -1, -1):
            for k in range(i - 1, -1, -1):
                factor = V[k][i]
                if factor != 0:
                    V[k][i] = 0
                    Y[k] = Y[k] - Y[i].mul(factor)
                    
        return Y[::-1]