def euclEtend(a,b):
    if b==0:
        return (a,1,0)
    else:
        d,u,v=euclEtend(b,a%b)
        return (d,v,u-(a//b)*v)
