from numpy import array, dot
from numpy.linalg import inv, eig, eigh

def solve(A, B):
    print "inverting"
    M = dot(inv(B), A)
    print "solving"
    w, v = eig(M)
    print "sorting the eigenvalues"

    r = []
    for i in range(len(w)):
        vec = v[:, i]
        r.append((w[i], vec))
    r.sort(key=lambda x: x[0])
    print "eigenvalues:"
    eigs = []
    for w, vec in r:
        if w > 0:
            break
        print w
        eigs.append(vec)
    return r
