def occurence_commune(L1, L2):
    for i in range(len(L1)):
        for j in range(len(L2)):
            if L1[i] == L2[j]:
                return (i, j)
            
def Shanks2(p, g, h, r):
    h = h%p
    s = 1 + int(r**0.5)
    b = euclEtend(expoRapMod(g, s, p), p)[1]%p
    L1 = [1] + [expoRapMod(g, i, p) for i in range(1,s+1)]
    L2 = [h*expoRapMod(b, i, p)%p for i in range(s+1)]
    x, y = occurence_commune(L1,L2)
    return (x + y*s)%r

