from typing import Self
from homhash.circuit import *

class Prover:
    def __init__(self, eval_table):
        self.state = eval_table

class Verifier:
    def __init__(self, claim_sum):
        self.claim_sum = claim_sum
        self.r_list = []