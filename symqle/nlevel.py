from .operator import operator, dagger, add_rule

def nlevel(b, N, label="s", index=None):
    assert isinstance(N, int)
    assert 1 <= N
    assert isinstance(label, str)
    assert not index or isinstance(index, int)
    if isinstance(index, int):
        assert index > 0

    sigmas = []
    bra = []
    ket = []
    for i in range(0, N):
        for j in range(i+1, N):
            if not index:
                sigmas.append(operator("%s_{%i%i}" %(label, i+1, j+1), b))
            else:
                sigmas.append(operator("%s_{%i%i,(%i)}" %(label, i+1, j+1, index), b))
            bra.append(j+1)
            ket.append(i+1)

    for i in range(0, len(sigmas)):
        add_rule(sigmas[i]**2, 0)
        add_rule(sigmas[i]*sigmas[i].dagger()*sigmas[i], sigmas[i])
        for j in range(i+1, len(sigmas)):
            str1 = str(sigmas[i].symbol)
            str2 = str(sigmas[j].symbol)

            if not index:
                bra1 = str1[-2]
                bra2 = str2[-2]
                ket1 = str1[-3]
                ket2 = str2[-3]

                if bra1 == ket2:
                    new_op = operator("%s_{%s%s}" %(label, ket1, bra2), b,
                        add_to_list=False)
                else:
                    new_op = 0
                add_rule(sigmas[i]*sigmas[j], new_op)

                if bra2 == ket1:
                    new_op = operator("%s_{%s%s}" %(label, ket2, bra1), b,
                        add_to_list=False)
                else:
                    new_op = 0
                add_rule(sigmas[j]*sigmas[i], new_op)

                if bra1 == bra2:
                    new_op = operator("%s_{%s%s}" %(label, ket1, ket2), b,
                        add_to_list=False)
                else:
                    new_op = 0
                add_rule(sigmas[i]*dagger(sigmas[j]), new_op)

                if ket1 == ket2:
                    new_op = operator("%s_{%s%s}" %(label, bra1, bra2), b,
                        add_to_list=False)
                else:
                    new_op = 0
                add_rule(dagger(sigmas[i])*sigmas[j], new_op)

            else:
                bra1 = str1[-6]
                bra2 = str2[-6]
                ket1 = str1[-7]
                ket2 = str2[-7]

                if bra1 == ket2:
                    new_op = operator("%s_{%s%s,(%i)}" %(label, ket1, bra2, index),
                        b, add_to_list=False)
                else:
                    new_op = 0
                add_rule(sigmas[i]*sigmas[j], new_op)

                if bra2 == ket1:
                    new_op = operator("%s_{%s%s,(%i)}" %(label, ket2, bra1, index),
                        b, add_to_list=False)
                else:
                    new_op = 0
                add_rule(sigmas[j]*sigmas[i], new_op)

                if bra1 == bra2:
                    new_op =operator("%s_{%s%s,(%i)}" %(label, ket1, ket2, index),
                    b, add_to_list=False)
                else:
                    new_op = 0
                add_rule(sigmas[i]*dagger(sigmas[j]), new_op)

                if ket1 == ket2:
                    new_op = operator("%s_{%s%s,(%i)}" %(label, bra1, bra2, index),
                        b, add_to_list=False)
                else:
                    new_op = 0
                add_rule(dagger(sigmas[i])*sigmas[j], new_op)

    return sigmas
