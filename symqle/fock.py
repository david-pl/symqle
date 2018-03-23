from .operator import operator, commutator_rule

def field_operators(b, label="a"):
    assert isinstance(label, str)
    a = operator(label, b)
    ad = a.dagger()
    commutator_rule(a, ad, 1)
    return a, ad
