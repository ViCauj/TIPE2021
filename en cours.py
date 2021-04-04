from time import time
from random import randint as rd
import matplotlib.pyplot as plt

def expoRapMod(x, n, p):
    '''renvoit x**n % p'''
    nBin = bin(n)[2:]
    res = x
    for i in range(1,len(nBin)):
        if nBin[i] == '1':
            res = ((res**2)*x)%p
        else:
            res = (res**2)%p
    return res

##################################################################################################
############################################## POLLARD ###########################################
##################################################################################################

def euclEtend(a,b):
    if b==0:
        return (a,1,0)
    else:
        d,u,v=euclEtend(b,a%b)
        return (d,v,u-(a//b)*v)

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
        return (expoRapMod(x, 2, p), (2*a)%r, (2*b)%r)
    return ((x*Lg[j])%p, (a+u[j])%r, (b+v[j])%r)

def Pollard(p, g, h, r = -1, n = 3):
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
    u = [rd(0,r+1) for i in range(n)]
    v = [rd(0,r+1) for i in range(n)]
    Lg = [ (expoRapMod(g, u[i], p)*expoRapMod(h, v[i], p))%p for i in range(n)]
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
    return res

##################################################################################################
############################################# FIN POLLARD ########################################
##################################################################################################

##################################################################################################
############################################## SHANKS 0.5 ########################################
##################################################################################################

def Shanks05(p, g, h, r):
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

##################################################################################################
########################################### FIN SHANKS 0.5 #######################################
##################################################################################################

##################################################################################################
############################################## SHANKS 2 ##########################################
##################################################################################################

def occurence_commune(L1, L2):
    for i in range(len(L1)):
        for j in range(len(L2)):
            if L1[i] == L2[j]:
                return (i, j)

def Shanks2(p, g, h, r):
    s = 1 + int(p**0.5)
    b = euclEtend(expoRapMod(g, s, p), p)[1]%p
    L1 = [1] + [expoRapMod(g, i, p) for i in range(1,s+1)]
    L2 = [expoRapMod(h*b, i, p) for i in range(s+1)]
    x, y = occurence_commune(L1,L2)
    return (x + y*s)%r

##################################################################################################
############################################ FIN SHANKS 2 ########################################
##################################################################################################

##################################################################################################
############################################### BRUTE ############################################
##################################################################################################

def brute(p, g, h, r):
    res = 1
    expo = 0
    while res != h%p :
        res *= g
        res %= p
        expo += 1
    return expo

##################################################################################################
############################################# FIN BRUTE ##########################################
##################################################################################################

##################################################################################################
############################################### RANDOM ###########################################
##################################################################################################

def random(p, g, h, r):
    expo = rd(0, r)
    res = expoRapMod(g, expo, p)
    while res != h%p:
        expo = rd(0, r)
        res = expoRapMod(g, expo, p)
    return expo

##################################################################################################
############################################# FIN RANDOM #########################################
##################################################################################################

def benBenchmark(f, p, g, r, h):
    '''
    Dans Z/pZ 
    f la fonction à tester 
    g le générateur du sous groupe d'ordre r
    '''
    test = True
    
    while test:
        res = time()
        r = f(p, g, h, r)
        res = time() - res
        if r != -1 :
            test = False
    return res*(10**3)
    
def benchmark(L):
    '''
    L = [Liste des p, Liste des g, Liste des ordres (= r)]
    '''
    tt = []
    
    tPollard = []
    tShanks05 = []
    tShanks2 = []
    tBrute = []
    tRandom = []
    
    for i in range(len(L)):
        PGR = L[i]
        h = expoRapMod(PGR[1], rd(1, PGR[2] - 1), PGR[0])
        tt.append(PGR[0])
        
        a, b, n = "O", "-", 10
        print(tt[-1]," :")
        print(a*n*0 + b*n*5)
        tPollard.append(benBenchmark(Pollard, PGR[0], PGR[1], PGR[2], h))
        print(a*n*1 + b*n*4)
        tShanks05.append(benBenchmark(Shanks05, PGR[0], PGR[1], PGR[2], h))
        print(a*n*2 + b*n*3)
        tShanks2.append(benBenchmark(Shanks2, PGR[0], PGR[1], PGR[2], h))
        print(a*n*3 + b*n*2)
        tBrute.append(benBenchmark(brute, PGR[0], PGR[1], PGR[2], h))
        print(a*n*4 + b*n*1)
        tRandom.append(benBenchmark(random, PGR[0], PGR[1], PGR[2], h))
        print(a*n*5 + b*n*0)
        print("\n")
        
    plt.plot(tt, tPollard, label = "Pollard")
    plt.plot(tt, tShanks05, label = "Shank O(n**0.5)")
    plt.plot(tt, tShanks2, label = "Shank O(n**2)")
    plt.plot(tt, tBrute, label = "Brute")
    plt.plot(tt, tRandom, label = "Random")
    plt.legend()
    plt.show()

L = [(3, 2, 2), (5, 4, 2 ), (7, 2, 3), (11, 9, 5), (13, 9, 3), 
     (31, 2, 5), (37, 26, 3), (41, 40, 2), (47, 7, 23), (53, 49, 13),
     (59, 12, 29), (67, 40, 11), (71, 45, 7), (79, 38, 13), (83, 69, 41),
     (89, 64, 11), (499, 345, 83), (797, 256, 199), (1091, 3, 109), (1447, 684, 241),
     (2003, 484, 13), (8111, 5152, 811), (8699, 4856, 4349), (10007, 6674, 5003),
     (11587, 6615, 1931), (14143, 4660, 2357), (66109, 32788, 787), (65467463, 38859570, 1103),
     (654656818789, 171720815255, 244411)]
