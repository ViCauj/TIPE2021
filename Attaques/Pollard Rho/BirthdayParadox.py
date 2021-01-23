from numpy.random import randint as rd
import matplotlib.pyplot as plt


N=10000
n=500
res=[]
moy=[]
lim=[(3.14*N/2)**0.5 + 2 for i in range(n)]

for i in range(n):
    res.append(1)
    S=[rd(0,N)]
    while not S[-1] in S[:-1]:
        S.append(rd(0,N))
        res[-1]+=1
    moy.append(sum(i for i in res)/len(res))
    
plt.plot(list(range(n)),res,"grey")
plt.plot(list(range(n)),moy)
plt.plot(list(range(n)),lim,"r")
