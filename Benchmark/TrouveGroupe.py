from numpy.random import randint as rd

def trouveG(listeD, p):
    res = 0
    test = True
    while test:
        g = rd(1, p)
        for li in listeD:
            if expoRapMod(g, li ,p) == 1:
                res = (g, li)
                test = False
    return res
