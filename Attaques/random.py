from random import randint as rd

def random(p, g, h, r):
    expo = rd(0, r)
    res = expoRapMod(g, expo, p)
    while res != h%p:
        expo = rd(0, r)
        res = expoRapMod(g, expo, p)
    return expo
