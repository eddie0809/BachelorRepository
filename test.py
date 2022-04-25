from numpy import array
import partialtrace

# this is just to test if the partial trace works

rho = array([[ 0.56, 0.1 ,  -0.12,  0.  ],
             [ 0.1 , 0.24,   0.  ,  0.12],
             [-0.12, 0.  ,   0.14, -0.1 ],
             [ 0.  , 0.12,  -0.1 ,  0.06]])

N = 1
pt = partialtrace.PartialTrace(N)
print(rho)
print()
print(pt.get_first_state(rho))
print()
print(pt.get_last_state(rho))
print()
