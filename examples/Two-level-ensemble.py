# coding: utf-8
from symqle import *

# Number of atoms
N = 2

# Corresponding bases
bases = []
for i in range(1, N+1):
    bases.append('atom-%i' %i)

# Operators
sm = []
for i in range(1, N+1):
    sm.append(operator('\sigma_{%i}' %i, bases[i-1]))

# Rules
rs = ruleset(sm)
projector_rules(rs, sm)

# Hamiltonian
om = Symbol('ω')
H = om*sum([dagger(s)*s for s in sm])

print(pretty(H))

# Lindblad
gam = Symbol('γ')
rates = [gam for i in range(N)]
L = sm

print(pretty(langevin(sm[1], H, L, rates, ruleset=rs, max_iter=2)))

