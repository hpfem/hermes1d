def solve_eig_numpy(A, B):
    # in case A or B is a CooMatrix:
    if hasattr(A, "to_scipy_coo"):
        A = A.to_scipy_coo().todense()
    if hasattr(B, "to_scipy_coo"):
        B = B.to_scipy_coo().todense()

    from numpy import array, dot
    from numpy.linalg import inv, eig, eigh

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
