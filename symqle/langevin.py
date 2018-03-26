from .operator import commutator, apply_rules
from sympy import sqrt

def langevin(op, H, L=[], rates=[], noise=[], max_iter=1):
    L_op = []
    for i in range(0, len(L)):
        L_op.append(-rates[i]*commutator(op, L[i].dagger())*0.5*L[i]+
            rates[i]*0.5*L[i].dagger()*commutator(op, L[i]))

    for i in range(0, len(noise)):
        L_op.append(-sqrt(rates[i])*commutator(op, L[i].dagger())*noise[i]+
            sqrt(rates[i])*noise[i].dagger()*commutator(op, L[i]))

    LE = 1j*commutator(H, op) + sum(L_op)
    return apply_rules(LE, max_iter)
