def prod(L):
    r = 1
    for i in L:
        r *= i
    return r

def PollardFact(N, k = 1):
    def f(x):
        return x**2 + k
    a = 2
    b = f(a)%N
    d = 0
    while not 1 < d < N :
        a = f(a)%N
        b = f(f(b))%N
        d = PGCD(b-a,N)
    return d

def listeD(N):
    n = N
    D = []
    while prod(D) != N:
        if Prime(n):
            D.append(n)
        else:
            r = PollardFact(n)
            # if Prime(r):
            #     D.append(r)
            # else:
            #     D.append(listeD(r))
            D.append(r)
            n //= D[-1]
    return D
        
