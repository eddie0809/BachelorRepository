from numpy import matmul
import prepState

class PartialTrace:
    __N = 0
    __basis = {}

    def __init__(self, n) -> None:
        self.__N = n
        for i in range(n+1):
            self.__basis['0',i] = prepState.vector0(self.__N,i)
            self.__basis['1',i] = prepState.vector1(self.__N,i)

    def __bin(self, x):
        return '{:b}'.format(x).zfill(self.__N)

    def __bin_set(self):
        return [ self.__bin(i) for i in range(2**self.__N) ]

    def get_first_state(self, rho):
        result = 0
        for b in self.__bin_set():
            proj = 1
            for (index, bit) in enumerate(b):
                proj *= self.__basis[bit, index+1]
            result += matmul(proj.T, matmul(rho, proj))
        return result

    def get_last_state(self, rho):
        result = 0
        for b in self.__bin_set():
            proj = 1
            for (index, bit) in enumerate(b):
                proj *= self.__basis[bit, index]
            result += matmul(proj.T, matmul(rho, proj))
        return result
