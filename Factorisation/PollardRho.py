def prod(L):
    r = 1
    for i in L:
        r *= i
    return r

def PollardFact(N, a = 2, k = 1, lim = 10000):
    def f(x):
        return x**2 + k
    b = f(a)%N
    d = 0
    compt = 0
    while not 1<d<N :
        a = f(a)%N
        b = f(f(b))%N
        d = PGCD(b-a,N)
        compt += 1
        if compt > lim:
            return PollardFact(N, random.randint(3, N), random.randint(-10, 10), lim)
    return d

def listeD(N):
    n=N
    D=[]
    if n%2 == 0:
        n//=2
        D.append(2)
    while prod(D)!=N:
        if Prime(n):
            D.append(n)
        else:
            r=PollardFact(n)
            if Prime(r):
                 D.append(r)
            else:
                D+=listeD(r)
            n//=D[-1]
        #print("éhop", D)
    res = []
    for i in D:
        if not (i in res):
            res.append(i)
    return res
        
