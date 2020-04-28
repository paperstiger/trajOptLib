# Test if autograd works
import numpy as np
from autograd import numpy
from trajoptlib import NonLinearPointConstr


class SimpleConstr(NonLinearPointConstr):
    def __init__(self):
        NonLinearPointConstr.__init__(self, 0, 4, 4, 2, 0, None, None, 'ad')

    def __ad__(self, x):
        # x is of size 1 + 4 + 2 = 7
        return numpy.array([x[3], x[4], x[5], x[6]])


def main():
    constr = SimpleConstr()
    import pdb; pdb.set_trace()
    xin = np.random.random(7)
    G = np.zeros(constr.nG)
    row, col = np.zeros((2, constr.nG), dtype=int)
    F = np.zeros(constr.nf)
    constr.__callg__(xin, F, G, row, col, True, True)
    print(F, G, row, col)


if __name__ == '__main__':
    main()
