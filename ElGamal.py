import numpy.random as rd

def clePub(a,p,g):
    return (p,g,g**a%p)     

def cryptCar(car,a,p,g):
    cleEph = rd.randint(1,p)
    return (g**cleEph%p,ord(car)*(clePub(a,p,g)[2]**cleEph)%p)

def deCryptCar(carChif,a,p):
    return chr((carChif[0]**(p-a-1)%p)*carChif[1]%p)

def cryptTout(message,a,p,g):
    mesCPC=list(message)
    return [cryptCar(car,a,p,g) for car in mesCPC]

def deCryptTout(messageChif,a,p):
    mes=""
    for carChif in messageChif:
        mes+=deCryptCar(carChif,a,p)
    return mes
