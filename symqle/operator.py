from sympy import *
import copy


class operator:
    def __init__(self, label, hilbertspace, hermitian=False):

        if isinstance(hilbertspace, str):
            self.basis = Symbol(hilbertspace, commutative=False)
        else:
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
