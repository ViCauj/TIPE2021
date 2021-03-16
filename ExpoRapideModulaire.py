#REM : peut être coder le passage en binaire 

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
