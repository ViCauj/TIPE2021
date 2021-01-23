def euclEtend(a, b):
    u0,v0,u1,v1=1,0,0,1
    while b:
        q,r=a//b,a%b
        a,b=b,r
        u0,v0,u1,v1=u1,v1,u0-q*u1,v0-q*v1
    return u0,v0   
