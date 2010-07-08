from hermes1d import Mesh

def test_mesh1():
    m = Mesh(-1, 1, 5, 2)

def test_mesh2():
    m = Mesh(-1, 1, 5, 2, 1)

def test_mesh3():
    m1 = Mesh([0, 1, 3, 4], [2, 2, 2])
    m2 = Mesh([-1, -0.5, 0, 0.5, 1], [2, 2, 2, 2])
