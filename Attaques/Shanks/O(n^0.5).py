def Shanks05(p, g, h, r):
    h = h%p
    s = 1 + int(r**0.5)
    b = euclEtend(expoRapMod(g, s, p), p)[1]%p
    d1 = {1:0}
    for i in range(1, s + 1): #on crée le dictionnaire#
        d1[expoRapMod(g, i, p)] = i # les elts de L1 sont les clefs#
    for i in range(s + 1):#recherche d'une occurrence#
        k = h*expoRapMod(b, i, p)%p
        if k in d1: # commune#
            x, y = d1[k], i
            return (x + y*s)%r  #on stoppe l'algorithme #
            #quand on trouve la solution#
