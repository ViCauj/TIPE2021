def occurence_commune(L1, L2):
    for i in range(len(L1)):
        for j in range(len(L2)):
            if L1[i] == L2[j]:
                return (i, j)

def Shanks(A, g, p):
    s = 1 + int(np.sqrt(p))
    b = euclEtend(g**s, p)[1]%p
    L1 = [1] + [(g**i)%p for i in range(1,s+1)]
    L2 = [(A*b**i)%p for i in range(s+1)]
    x, y = occurence_commune(L1,L2)
    return x + y*s
