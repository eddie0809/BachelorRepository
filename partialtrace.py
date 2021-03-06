from numpy import matmul
import prepState


###
#
# maybe use QuTiP for this? q.ptrace(sel) exists, but this also exists, i'm conflicted.
#
###

from qutip import *

# Q.ptrace(sel) Partial trace returning components selected using 'sel' parameter.


# Partial trace defined as a class because of the vector basis
class PartialTrace:
    __N = 0
    __basis = {}

    def __init__(self, n) -> None:
        self.__N = n
        for i in range(n+1):
            self.__basis['0',i] = prepState.vector0(self.__N,i)
            self.__basis['1',i] = prepState.vector1(self.__N,i)

    def __bin_set(self):
        # to improve readability
        def bin(x):
            return '{:b}'.format(x).zfill(self.__N)
        return [ bin(i) for i in range(2**self.__N) ]

    def get_first_state(self, rho):
        """Perform the partial trace over all the spins but the first one.
            ā_{x_1,x_2,...,x_N} [šāāØx_1 x_2...x_N|] Ļ [šā|x_1 x_2...x_Nā©]
        """
        result = 0
        for b in self.__bin_set():
            proj = 1
            for (index, bit) in enumerate(b):
                proj *= self.__basis[bit, index+1]
            result += matmul(proj.T, matmul(rho, proj))
        return result

    def get_last_state(self, rho):
        """Perform the partial trace over all the spins but the first one.
            ā_{x_1,x_2,...,x_N} [āØx_1 x_2...x_N|āš] Ļ [|x_1 x_2...x_Nā©āš]
        """
        result = 0
        for b in self.__bin_set():
            proj = 1
            for (index, bit) in enumerate(b):
                proj *= self.__basis[bit, index]
            result += matmul(proj.T, matmul(rho, proj))
        return result
