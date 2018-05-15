from .operator import commutator, apply_rules
from sympy import sqrt
from sympy.matrices.dense import MutableDenseMatrix

def langevin(op, H, L=[], rates=[], noise=[], max_iter=1):
    L_op = []
    if isinstance(rates, list):
        for i in range(0, len(L)):
            L_op.append(-rates[i]*commutator(op, L[i].dagger())*0.5*L[i]+
                rates[i]*0.5*L[i].dagger()*commutator(op, L[i]))

        for i in range(0, len(noise)):
            L_op.append(-sqrt(rates[i])*commutator(op, L[i].dagger())*noise[i]+
                sqrt(rates[i])*noise[i].dagger()*commutator(op, L[i]))

    if isinstance(rates, MutableDenseMatrix):
        assert len(noise) == 0
        dims = rates.shape
        assert dims[0] == dims[1]
        for i in range(dims[0]):
            for j in range(dims[1]):
                L_op.append(-rates[i,j]*commutator(op, L[i].dagger())*0.5*L[j]+
                    rates[j,i]*0.5*L[i].dagger()*commutator(op, L[j]))

    LE = 1j*commutator(H, op) + sum(L_op)
    return apply_rules(LE, max_iter)
