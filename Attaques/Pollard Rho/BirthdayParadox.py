from numpy.random import randint as rd
import matplotlib.pyplot as plt

N = 10000
n = 500
res = []
moy = []
delt = []
lim = [ (3.14*N/2)**0.5 for i in range(n)]

for i in range(n):
    res.append(1)
    S = [rd(0, N)]
    while not S[-1] in S[:-1]:
        S.append(rd(0, N))
        res[-1] += 1
    moy.append(sum(i for i in res)/len(res))
    delt.append(abs((3.14*N/2)**0.5 - moy[-1]))
   
# plt.plot(list(range(n)),delt)
plt.plot(list(range(n)),res,"")
plt.plot(list(range(n)),moy,"orange")
plt.plot(list(range(n)),lim,"r")
