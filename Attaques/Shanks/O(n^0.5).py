def Shanks_en_mieux(h, g, p, r):
    s = 1 + int(p**0.5)
    b = euclEtend(expoRapMod(g, s, p),p)[1]%p
    d1 = {1:0}
    for i in range(1, s + 1): #on crée le dictionnaire#
        d1[expoRapMod(g, i, p)] = i # les elts de L1 sont les clefs#
    for i in range(s + 1):#recherche d'une occurrence#
        k = expoRapMod(h*b, i, p)
        if k in d1: # commune#
            x, y = d1[k], i
            return (x + y*s)%r  #on stoppe l'algorithme #
            #quand on trouve la solution#
