from .operator import operator, dagger, add_rule

def two_level(b, label="s", index=None):
    assert isinstance(label, str)
    assert not index or isinstance(index, int)
    if isinstance(index, int):
        assert index > 0

    if not index:
        sm = operator(label, b)
    else:
        sm = operator(label + "_{%i}" %index, b)
    sp = dagger(sm)
    add_rule(sm**2, 0)
    add_rule(sm*sp*sm, sm)
    add_rule(sm*sp, 1 - sp*sm, add_conjugate=False)

    return sm, sp
