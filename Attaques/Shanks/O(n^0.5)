def Shanks_en_mieux(A,g,p):
    s=1+int(np.sqrt(p))
    b=euclEtend(g**s,p)[1]%p
    d1={1:0}
    for i in range(1,s+1): #on crée le dictionnaire#
        d1[(g**i)%p]=i # les elts de L1 sont les clefs#
    for i in range(s+1):#recherche d'une occurrence#
        if ((A*b**i)%p) in d1: # commune#
            x,y=d1[(A*b**i)%p],i
            return x+y*s  #on stoppe l'algorithme #
            #quand on trouve la solution#
