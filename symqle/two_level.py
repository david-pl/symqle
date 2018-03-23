from .operator import operator, dagger, add_rule

def two_level(b, label="s"):
    assert isinstance(label, str)

    sm = operator("s", b)
    sp = dagger(sm)
    add_rule(sm**2, 0)
    add_rule(sm*sp*sm, sm)

    return sm, sp
