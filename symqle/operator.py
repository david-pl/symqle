from sympy import Symbol, init_printing, adjoint, expand, conjugate, prod, factorial, binomial
import sympy
# from sympy.core import add.Add, mul.Mul
import copy
import itertools

# Global set of rules and operators
__RULES = {}
__RULES_BASIS = {}
OPERATORS_ = []

class operator:
    def __init__(self, label, hilbertspace, hermitian=False):
        self.basis = hilbertspace

        if isinstance(label, str):
            self.symbol = Symbol(label, commutative=False)
        else:
            self.symbol = label

        self.ishermitian = hermitian

        # Add to global operator list
        if isinstance(self.symbol, Symbol):
            global OPERATORS_
            if not self in OPERATORS_:
                OPERATORS_.append(self)
                commutator_set(OPERATORS_)

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
            return operator(other.symbol*self.symbol, basis,
                self.ishermitian and other.ishermitian)
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

def show_operators():
    print(OPERATORS_)

def clear_operators():
    global OPERATORS_
    OPERATORS_ = []

def basis(label):
    assert isinstance(label, str)
    return Symbol(label, commutative=False)

def show_basis(op):
    assert isinstance(op, operator)
    return op.basis.free_symbols

# Average
def average(op_, order=1, divide_cumulants=False, max_iter=1, hard_cutoff=False):
    assert isinstance(op_, operator)
    op = apply_rules(op_, max_iter=max_iter)

    sym = op.symbol
    if isinstance(sym, sympy.mul.Mul):
        return factorize_prod(sym, order, divide_cumulants, hard_cutoff)
    elif isinstance(sym, sympy.add.Add):
        return factorize_sum(sym, order, divide_cumulants, hard_cutoff)
        # number, summands = sym.as_coeff_add()
    elif isinstance(sym, sympy.symbol.Symbol):
        if is_operator(sym):
            return Symbol("<" + str(sym) + ">")
        else:
            return sym
    elif isinstance(sym, adjoint):
        op_str = str(sym)
        op_str = op_str.replace("adjoint(", "")[0:-1]
        return conjugate(Symbol("<" + op_str + ">"))
    elif isinstance(sym, sympy.power.Pow):
        return factorize_pow(sym, order)
    else:
        raise("ERROR: Unsupported math operation!")


def factorize_prod(sym, order, divide_cumulants=False, hard_cutoff=False):
    factor, coeffs = sym.as_coeff_mul()

    # Get which coefficients are operators and parameters (symbols)
    ops = []
    syms = []
    for i in range(len(coeffs)):
        if is_operator(coeffs[i]):
            ops.append(coeffs[i])
        elif isinstance(coeffs[i], sympy.power.Pow):
            if is_operator(coeffs[i].base):
                tmp = [coeffs[i].base for k in range(coeffs[i].exp)]
                ops.extend(tmp)
            else:
                syms.append(coeffs[i])
        else:
            syms.append(coeffs[i])

    if len(syms) == 0:
        syms = [1]

    # Factorize if product has more constituents than order
    if len(ops) > order:
        if hard_cutoff:
            ops_new = 0
        else:
            bases_ = find_bases(ops)
            decision = decide_factor(bases_, order)
            if not decision:
                ops_new = cumulant_expansion(ops, order, divide_cumulants)

    # Return average if product is smaller order
    else:
        if len(ops) > 0:
            op_str = str(prod(ops))
            if len(ops) == 1 and isinstance(ops[0], adjoint):
                op_str = op_str.replace("adjoint(", "")[:-1]
                ops_new = conjugate(Symbol("<" + op_str + ">"))
            else:
                ops_new = Symbol("<" + op_str + ">")
        else:
            ops_new = 1

    return factor*prod(syms)*ops_new

def factorize_sum(sym, order, divide_cumulants=False, hard_cutoff=False):
    number, ops = sym.as_coeff_add()
    expr = []
    # TODO: Clean-up
    for i in range(len(ops)):
        if isinstance(ops[i], sympy.mul.Mul):
            prod_new = factorize_prod(ops[i], order, hard_cutoff=hard_cutoff)
            expr.append(prod_new)
        elif isinstance(ops[i], sympy.add.Add):
            sum_new = factorize_sum(ops[i], order, hard_cutoff=hard_cutoff)
            expr.append(sum_new)
        elif isinstance(ops[i], sympy.symbol.Symbol):
            if is_operator(ops[i]):
                op_ = Symbol("<" + str(ops[i]) + ">")
            else:
                op_ = ops[i]
            expr.append(op_)
        elif isinstance(ops[i], sympy.numbers.ImaginaryUnit):
            expr.append(ops[i])
        elif isinstance(ops[i], adjoint):
            op_str = str(ops[i])
            op_str = op_str.replace("adjoint(", "")[0:-1]
            op_new_ = conjugate(Symbol("<" + op_str + ">"))
            expr.append(op_new_)
        elif isinstance(ops[i], sympy.power.Pow):
            expr.append(factorize_pow, hard_cutoff)
        else:
            print("WARNING:" + str(type(ops[i])))
    return number + sum(expr)

def factorize_pow(sym, order, hard_cutoff=False):
    n, op = sym.exp, sym.base

    if is_operator(op):
        if n > order:
            if hard_cutoff:
                return 0
            else:
                factorize_prod(prod([op for i in range(n)]), order, hard_cutoff=hard_cutoff)
        else:
            op_new = Symbol("<" + str(op) + ">")
            return op_new**n
    elif isinstance(op, adjoint):
        op_str = str(op)
        op_str = op_str.replace("adjoint(", "")[0:-1]
        return conjugate(Symbol("<" + op_str + ">")**n)
    else:
        return sym


def find_bases(ops):
    bases_ = []
    for o in ops:
        if isinstance(o, sympy.symbol.Symbol):
            if is_operator(o):
                bases_.append(find_basis(o))
            else:
                bases_.append(0)
        elif isinstance(o, adjoint):
            bases_.append(find_basis_adjoint(o))
        elif isinstance(o, sympy.power.Pow):
            n, op_ = o.exp, o.base
            if is_operator(op_):
                bases_.append(find_basis(op_))
            elif isinstance(op_, adjoint):
                bases_.append(find_basis_adjoint(op_))
            else:
                bases_.append(0)
        else:
            bases_.append(0)
    return bases_

def find_basis(op):
    for OP_ in OPERATORS_:
        if op == OP_.symbol:
            return OP_.basis
    print("WARNING: Untracked operator!")

def find_basis_adjoint(op):
    for OP_ in OPERATORS_:
        if op == adjoint(OP_.symbol):
            return OP_.basis
    print("WARNING: Untracked operator!")

def is_operator(op):
    for OP_ in OPERATORS_:
        if op == OP_.symbol or op == adjoint(OP_.symbol):
            return True
    return False

def decide_factor(bases, order):
    # TODO: Make more sophisticated choice
    return 0

def cumulant_expansion(ops, order, divide=False):
    # Get all combinations of length order
    combs_order = list(itertools.combinations(ops, order))

    # Get all remainding combinations
    combs_remain = list(itertools.combinations(ops, len(ops) - order))
    combs_remain.reverse()

    ops_new = []
    if order == 1:
        for c in combs_order:
            op_str = str(c[0])
            if "adjoint" in op_str:
                op_str = op_str.replace("adjoint(", "")[0:-1]
                op_str = "<" + op_str + ">"
                ops_new.append(conjugate(Symbol(op_str)))
            else:
                op_str = "<" + op_str + ">"
                ops_new.append(Symbol(op_str))
        return prod(ops_new)

    else:
        for i in range(len(combs_order)):
            ops_tmp = []

            op_str = str(prod(combs_order[i]))
            if len(combs_order[i]) == 1 and "adjoint" in op_str:
                op_str = op_str.replace("adjoint(", "")[:-1]
                ops_tmp.append(conjugate(Symbol("<" + op_str + ">")))
            else:
                ops_tmp.append(Symbol("<" + op_str + ">"))

            op_str2 = str(prod(combs_remain[i]))
            if "adjoint" in op_str2 and len(combs_remain[i]) == 1:
                    op_str2 = op_str2.replace("adjoint(", "")[:-1]
                    ops_tmp.append(conjugate(Symbol("<" + op_str2 + ">")))
                #TODO: replace "adjoint()" by "^\dagger" for longer averages
            else:
                ops_tmp.append(Symbol("<" + op_str2 + ">"))

            ops_new.append(prod(ops_tmp))

        if divide:
            n = binomial(len(ops), order)
        else:
            n = 1
        return sum(ops_new)/n



# Rules
def apply_rules(op, max_iter=1, iteration=1):
    sym = expand(op.symbol).subs(__RULES)
    basis = expand(op.basis).subs(__RULES_BASIS)
    op_out = operator(sym, basis)
    if iteration < max_iter:
        return apply_rules(op_out, max_iter, iteration+1)
    else:
        return op_out

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
    global __RULES, __RULES_BASIS
    __RULES = {}
    __RULES_BASIS = {}

def clear():
    clear_operators()
    clear_rules()
