from sympy import Symbol, init_printing, adjoint, expand
# from .rules import __RULES
import copy

class operator:
    def __init__(self, label, hilbertspace, hermitian=False):
        self.basis = hilbertspace

        if isinstance(label, str):
            self.symbol = Symbol(label, commutative=False)
        else:
            self.symbol = label

        self.ishermitian = hermitian

    def dagger(self):
        if self.ishermitian:
            return self
        else:
            return operator(adjoint(self.symbol), self.basis)

    def hermitian(self):
        self.ishermitian = True

    def __add__(self, other):
        if isinstance(other, operator):
            basis = self.basis + other.basis
            return operator(self.symbol + other.symbol, basis,
                self.ishermitian and other.ishermitian)
        else:
            return operator(self.symbol + other, self.basis, self.ishermitian)

    def __radd__(self, other):
        return self + other

    def __neg__(self):
        return operator(-self.symbol, self.basis, self.ishermitian)

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return other + (-self)

    def __mul__(self, other):
        if isinstance(other, operator):
            basis = self.basis*other.basis
            return operator(self.symbol*other.symbol, basis)
        else:
            return operator(self.symbol*other, self.basis, self.ishermitian)

    def __rmul__(self, other):
        if isinstance(other, operator):
            basis = other.basis*self.basis
            return operator(other.symbol*self.symbol, basis)
        else:
            return operator(other*self.symbol, self.basis, self.ishermitian)

    def __pow__(self, n):
        return operator(self.symbol**n, self.basis, self.ishermitian)

    def __eq__(self, other):
        if isinstance(other, operator):
            return self.symbol == other.symbol and samebasis(self, other) and (self.ishermitian == other.ishermitian)
        else:
            return False

    def __neq__(self, other):
        return not self == other

    def __str__(self):
        string = str(self.symbol)
        return string.replace("adjoint", "dagger")

    def __repr__(self):
        return self.__str__()


def dagger(x):
    return x.dagger()

def commutator(a, b):
    return a*b - b*a

def samebasis(a, b):
    return a.basis == b.basis

def samebasis_nofactors(a, b):
    """
    Check if all consituents of a basis expression are the same.
    Useful for trivial commutation relations, i.e. if none of
    them are the same (check=False) the operators commute.
    """
    check = True
    for a_base in a.basis.free_symbols:
        for b_base in b.basis.free_symbols:
            check = a_base == b_base
    return check

def hermitian(x):
    return x.hermitian()

def copy_operator(x, new_basis):
    y = copy.copy(x)
    if isinstance(new_basis, str):
        y.basis = Symbol(new_basis, commutative=False)
    else:
        y.basis = new_basis
    return y

def pretty(x):
    init_printing()
    return x.symbol

def apply_rules(op, max_iter=1, iteration=1):
    sym = expand(op.symbol).subs(__RULES)
    basis = expand(op.basis).subs(__RULES_BASIS)
    op_out = operator(sym, basis)
    if iteration < max_iter:
        return apply_rules(op_out, max_iter, iteration+1)
    else:
        return op_out

def basis(label):
    assert isinstance(label, str)
    return Symbol(label, commutative=False)


# Rules
__RULES = {}
__RULES_BASIS = {}

def add_rule(op, res, add_conjugate=True):
    if isinstance(res, operator):
        __RULES[op.symbol] = res.symbol
        __RULES_BASIS[op.basis] = res.basis
        if add_conjugate and not op.ishermitian:
            op_dag = dagger(op)
            res_dag = dagger(res)
            __RULES[op_dag.symbol] = res_dag.symbol
            __RULES_BASIS[op_dag.basis] = res_dag.basis
    else:
        __RULES[op.symbol] = res
        __RULES_BASIS[op.basis] = res
        if add_conjugate and not op.ishermitian:
            op_dag = dagger(op)
            __RULES[op_dag.symbol] = res
            __RULES_BASIS[op_dag.basis] = res

def commutator_set(operators):
    for i in range(len(operators) - 1):
        op1 = operators[i]
        for j in range(i+1, len(operators)):
            op2 = operators[j]
            if not samebasis_nofactors(op1, op2):
                commutator_rule(op1, op2, 0, False)
                if not op1.ishermitian:
                    commutator_rule(op1.dagger(), op2, 0, False)
                    if not op2.ishermitian:
                        commutator_rule(op1.dagger(), op2.dagger(), 0, False)
                        commutator_rule(op1, op2.dagger(), 0, False)

def commutator_rule(a, b, res, add_conjugate=True):
    add_rule(a*b, res + b*a, add_conjugate)

def show_rules():
    print(__RULES)

def clear_rules():
    __RULES = {}
    __RULES_BASIS = {}
