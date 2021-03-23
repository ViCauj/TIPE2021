from numpy.random import randint as rds

def Fp(p, g):
    '''
    g dans Z/pZ,
    G = <g>
    '''
    G = [1, g]
    while G[-1] != 1:
        G.append((G[-1]*g)%p)
    G.pop()
    return (G, len(G))

def inv(a, p):
    '''inverse de a dans Z/pZ'''
    return euclEtend(a, p)[1]%p
    
##################The Pseudorandom Walk######################

def S(g, n):
    '''
    g dans Z/pZ,
    G = <g>,
    n le nombre de sous ensembles Si de G = S0 U S1 U ... U Sn-1
    '''
    return int(bin(g)[2:])%n

def DecoupeG(G, n):          #calculs + rapides avec n une puiss de 2
    GBis=[ [] for i in range(n)]
    for g in G:
        GBis[S(g, n)].append(g)
    return GBis

def Walk(trouple, cste):
    '''
    Dans Z/pZ,
    Lg la liste des gj, u la liste des uj et v la liste des vj,
    x = xi, a = ai et b = bi,
    n le nombre de subdivisions Si de G = S0 U S1 U ... U Sn-1,
    r l'ordre de G = <g>
    '''
    Lg, u, v, n, p, r = cste
    x, a, b = trouple
    j = S(x, n)
    if j == 0:
        return ((x**2)%p, (2*a)%r, (2*b)%r)
    return ((x*Lg[j])%p, (a+u[j])%r, (b+v[j])%r)
        
def PollardRho(p, g, h, r = -1, n = 32):
    '''
    g dans Z/pZ,
    G = <g>,
    n le nombre de subdivisions Si de G = S0 U S1 U ... U Sn-1,
    h = g^a avec a le truc qu'on cherche (a = logg(h)) et h dans G
    r l'ordre de g que l'on peut donner pour accelerer l'algo si on le souhaite
    '''
    if r == -1:
        G, r = Fp(p, g)
    
    ############ CHOIX DE LA MARCHE  A L E A T O I R E !!!!!!!!!!! ############
    #u, v = list(rds(0, r, n)), list(rds(0, r, n))
    u = [rds(0,r) for i in range(n)]
    v = [rds(0,r) for i in range(n)]
    Lg = [ (((g**u[i])%p)*((h**v[i])%p))%p for i in range(n)]
    cste = Lg, u, v, n, p, r
    trouple1 = (g, 1, 0)
    trouple2 = Walk(trouple1, cste)
    ########## FIN CHOIX DE LA MARCHE  A L E A T O I R E !!!!!!!!!!! ##########
    
    while trouple1[0] != trouple2[0]:
        trouple1 = Walk(trouple1, cste)
        trouple2 = Walk(Walk(trouple2, cste), cste)
    if trouple1[2]%r == trouple2[2]%r:
        return -1
    res = ((trouple2[1]-trouple1[1])*inv(trouple1[2]-trouple2[2], r))%r
    return res, (g**res)%p
