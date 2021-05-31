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
    h = h%p
    s = 1 + int(r**0.5)
    b = euclEtend(expoRapMod(g, s, p), p)[1]%p
    L1 = [1] + [expoRapMod(g, i, p) for i in range(1,s+1)]
    L2 = [h*expoRapMod(b, i, p)%p for i in range(s+1)]
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
    return res
    
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
        #tShanks2.append(benBenchmark(Shanks2, PGR[0], PGR[1], PGR[2], h))
        print(a*n*3 + b*n*2)
        #tBrute.append(benBenchmark(brute, PGR[0], PGR[1], PGR[2], h))
        print(a*n*4 + b*n*1)
        #tRandom.append(benBenchmark(random, PGR[0], PGR[1], PGR[2], h))
        print(a*n*5 + b*n*0)
        print("\n")
        
    plt.plot(tt, tPollard, label = "Pollard")
    plt.plot(tt, tShanks05, label = "Shanks O(r**0.5)")
    #plt.plot(tt, tShanks2, label = "Shanks O(r**2)")
    #plt.plot(tt, tBrute, label = "Brute")
    #plt.plot(tt, tRandom, label = "Random")
    plt.legend()
    plt.xlabel('p')
    plt.ylabel("temps d'éxécution (en s)")
    
    plt.savefig('10000.png',dpi=500)
    plt.show()
    
def benchmark2(L):
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
        tt.append(PGR[2])
        
        a, b, n = "O", "-", 10
        print(tt[-1]," :")
        print(a*n*0 + b*n*5)
        tPollard.append(benBenchmark(Pollard, PGR[0], PGR[1], PGR[2], h))
        print(a*n*1 + b*n*4)
        tShanks05.append(benBenchmark(Shanks05, PGR[0], PGR[1], PGR[2], h))
        print(a*n*2 + b*n*3)
        #tShanks2.append(benBenchmark(Shanks2, PGR[0], PGR[1], PGR[2], h))
        print(a*n*3 + b*n*2)
        #tBrute.append(benBenchmark(brute, PGR[0], PGR[1], PGR[2], h))
        print(a*n*4 + b*n*1)
        #tRandom.append(benBenchmark(random, PGR[0], PGR[1], PGR[2], h))
        print(a*n*5 + b*n*0)
        print("\n")
        
    plt.plot(tt, tPollard, label = "Pollard")
    plt.plot(tt, tShanks05, label = "Shanks O(r**0.5)")
    #plt.plot(tt, tShanks2, label = "Shanks O(r**2)")
    #plt.plot(tt, tBrute, label = "Brute")
    #plt.plot(tt, tRandom, label = "Random")
    plt.legend()
    plt.xlabel('r')
    plt.ylabel("temps d'éxécution (en s)")
    
    name = "".join(chr(rd(97,122)) for i in range(10))
    plt.savefig("image\\"+ name+'.png',dpi=500)
    print(name)
    
    plt.show()
    
L = [(3, 2, 2), (5, 4, 2 ), (7, 2, 3), (11, 9, 5), (13, 9, 3), 
     (31, 2, 5), (37, 26, 3), (41, 40, 2), (47, 7, 23), (53, 49, 13),
     (59, 12, 29), (67, 40, 11), (71, 45, 7), (79, 38, 13), (83, 69, 41),
     (89, 64, 11), (499, 345, 83), (797, 256, 199), (1091, 3, 109), (1447, 684, 241),
     (2003, 484, 13), (8111, 5152, 811), (8699, 4856, 4349), (10007, 6674, 5003),
     (11587, 6615, 1931), (14143, 4660, 2357), (66109, 32788, 787), (65467463, 38859570, 1103),
     (654656818789, 171720815255, 244411)]

L2 = [(3, 2, 2),
 (5, 4, 2),
 (7, 2, 3),
 (11, 9, 5),
 (13, 9, 3),
 (31, 2, 5),
 (37, 26, 3),
 (41, 40, 2),
 (47, 7, 23),
 (53, 49, 13),
 (59, 12, 29),
 (67, 40, 11),
 (71, 45, 7),
 (79, 38, 13),
 (83, 69, 41),
 (89, 64, 11),
 (499, 345, 83),
 (743, 127, 53),
 (797, 256, 199),
 (971, 655, 97),
 (1091, 3, 109),
 (1447, 684, 241),
 (2003, 484, 13),
 (8111, 5152, 811),
 (8699, 4856, 4349),
 (10007, 6674, 5003),
 (11587, 6615, 1931),
 (14143, 4660, 2357),
 (66109, 32788, 787),
 (67157, 50100, 163),
 (67343, 20661, 3061),
 (67567, 4462, 11261),
 (65467463, 38859570, 1103),
 (66565337, 60653923, 5897),
 (66566527, 42633117, 2953),
 (66567343, 24578643, 652621),
 (586455911, 207714039, 327629),
 (665965877, 229364624, 166491469),
 (665966879, 317443264, 543203),
 (665967019, 366963373, 2707183),
 (5644648043, 1556261548, 256574911),
 (5644648781, 529414465, 21001),
 (5812646633, 1018572633, 55890833),
 (6424647101, 266865476, 64246471),
 (6424647923, 1540415202, 12457),
 (56424647837, 31784652565, 243533),
 (56424649637, 50249230461, 14106162409),
 (581264556647, 187682172571, 290632278323),
 (581264556839, 452343406884, 18979),
 (581264556883, 250754519206, 12123317),
 (581264557109, 383444230880, 7388831),
 (654656818789, 171720815255, 244411),
 (954564615737, 452626986020, 242453),
 (954564618143, 776501551948, 912585677)]

M = [(3, 2, 2),
 (7, 2, 3),
 (31, 2, 5),
 (71, 45, 7),
 (67, 40, 11),
 (53, 49, 13),
 (47, 7, 23),
 (59, 12, 29),
 (83, 69, 41),
 (743, 127, 53),
 (499, 345, 83),
 (971, 655, 97),
 (1091, 3, 109),
 (67157, 50100, 163),
 (797, 256, 199),
 (1447, 684, 241),
 (66109, 32788, 787),
 (8111, 5152, 811),
 (65467463, 38859570, 1103),
 (11587, 6615, 1931),
 (14143, 4660, 2357),
 (66566527, 42633117, 2953),
 (67343, 20661, 3061),
 (8699, 4856, 4349),
 (10007, 6674, 5003),
 (66565337, 60653923, 5897),
 (67567, 4462, 11261),
 (6424647923, 1540415202, 12457),
 (581264556839, 452343406884, 18979),
 (5644648781, 529414465, 21001),
 (954564615737, 452626986020, 242453),
 (56424647837, 31784652565, 243533),
 (654656818789, 171720815255, 244411),
 (586455911, 207714039, 327629),
 (665966879, 317443264, 543203),
 (66567343, 24578643, 652621),
 (665967019, 366963373, 2707183),
 (581264557109, 383444230880, 7388831),
 (581264556883, 250754519206, 12123317),
 (5812646633, 1018572633, 55890833),
 (6424647101, 266865476, 64246471),
 (665965877, 229364624, 166491469),
 (5644648043, 1556261548, 256574911),
 (954564618143, 776501551948, 912585677),
 (56424649637, 50249230461, 14106162409),
 (581264556647, 187682172571, 290632278323)]
