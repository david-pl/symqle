from .operator import operator, commutator_rule

def field_operators(b, label="a", index=None):
    assert isinstance(label, str)
    assert not index or isinstance(index, int)
    if isinstance(index, int):
        assert index > 0

    if not index:
        a = operator(label, b)
    else:
        a = operator(label + "_{%i}" %index)
    ad = a.dagger()

    commutator_rule(a, ad, 1)
    return a, ad
