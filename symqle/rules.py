from symqle.operator import *
import copy

class ruleset:
    def __init__(self, operators=[]):
        self.rules = {}
        self.rules_basis = {}
        self.operators = operators
        for i in range(len(operators) - 1):
            op1 = operators[i]
            for j in range(i+1, len(operators)):
                op2 = operators[j]
                if not samebasis_nofactors(op1, op2):
                    self.commutator(op1, op2, 0, False)
                    if not op1.ishermitian:
                        self.commutator(op1.dagger(), op2, 0, False)
                        if not op2.ishermitian:
                            self.commutator(op1.dagger(), op2.dagger(), 0, False)
                            self.commutator(op1, op2.dagger(), 0, False)

    def add_rule(self, op, res, add_conjugate=True):
        if isinstance(res, operator):
            self.rules[op.symbol] = res.symbol
            self.rules_basis[op.basis] = res.basis
            if add_conjugate and not op.ishermitian:
                op_dag = dagger(op)
                res_dag = dagger(res)
                self.rules[op_dag.symbol] = res_dag.symbol
                self.rules_basis[op_dag.basis] = res_dag.basis
        else:
            self.rules[op.symbol] = res
            self.rules_basis[op.basis] = res
            if add_conjugate and not op.ishermitian:
                op_dag = dagger(op)
                self.rules[op_dag.symbol] = res
                self.rules_basis[op_dag.basis] = res

    def commutator(self, a, b, res, add_conjugate=True):
        self.add_rule(a*b, res + b*a, add_conjugate)


def apply_rules(op, ruleset, max_iter=1, iteration=1):
    sym = expand(op.symbol).subs(ruleset.rules)
    basis = expand(op.basis).subs(ruleset.rules_basis)
    op_out = operator(sym, basis)
    if iteration < max_iter:
        return apply_rules(op_out, ruleset, max_iter, iteration+1)
    else:
        return op_out

def copy_rules(old_rules, old_basis, new_basis, new_ops):
    new_rules = copy.copy(old_rules)
    old_ops = old_rules.operators

    old_keys_vals = []
    for key, value in new_rules.rules.items():
        old_keys_vals.append((key, value))

    for old in old_keys_vals:
        key = old[0]
        value = old[1]
        for i in range(len(old_ops)):
            new_key = key.subs({old_ops[i].symbol : new_ops[i].symbol})
            if isinstance(value, Symbol):
                new_value = value.subs({old_ops[i].symbol : new_ops[i].symbol})
            else:
                new_value = value

            new_rules.rules[new_key] = new_value

    new_rules.operators.extend(new_ops)
    return new_rules
