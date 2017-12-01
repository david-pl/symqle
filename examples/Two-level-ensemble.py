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
    sm.append(operator('σ%i' %i, bases[i-1]))

# Rules
rs = ruleset(sm)
for i in range(0, N):
    rs.add_rule(sm[i]*sm[i], 0)
    rs.add_rule(sm[i]*dagger(sm[i])*sm[i], sm[i])

print(rs.rules)

# Hamiltonian
om = Symbol('ω')
H = om*sum([dagger(s)*s for s in sm])

print(H)

# Lindblad
gam = Symbol('γ')
rates = [gam for i in range(N)]
L = sm

print(langevin(sm[1], H, L, rates, ruleset=rs, max_iter=3))

