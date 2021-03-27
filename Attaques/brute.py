def brute(p, g, h, r):
    res = 1
    expo = 0
    while res != h%p :
        res *= g
        res %= p
        expo += 1
    return expo
