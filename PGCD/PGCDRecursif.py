def PGCD(a, b):
    if b == 0:
        return a
    else:
        return PGCD(b, a%b)
