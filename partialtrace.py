import prepState

N=2

basis = {}
for i in range(N+1):
	basis[0,i] = prepState.state0(N,i)
	basis[1,i] = prepState.state1(N,i)

def bin(x, N):
  return [ int(d) for d in '{:b}'.format(x).zfill(N) ]

for i in range(2**N):
  proj = 1
  for k in range(1,N+1):
    print(k)
    for j in bin(i,N):
      proj *= basis[j,k]
  print(proj)