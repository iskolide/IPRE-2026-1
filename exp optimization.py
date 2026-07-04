
primos_init = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 
                 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 
                 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 
                 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 
                 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 
                 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 
                 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 
                 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 
                 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 
                 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 
                 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 
                 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 
                 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 
                 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 
                 947, 953, 967, 971, 977, 983, 991, 997]

#right to left binary exponentiation
def exp1(a, b, n):
    a = a % n
    r = 1
    while b > 0:
        if b % 2 == 1:
            r = r * a % n
        a = a * a % n
        b >>= 1
    return r

#left to right binary exponentiation
def exp2(a, b, n):
    a = a % n
    r = 1
    j = 1 << (b.bit_length() - 1)
    while j > 0:
        r = r * r % n
        if b & j:
            r = r * a % n
        j >>= 1
    return r

#left to right k-ary exponentiation
def expk(a, b, n, k = 5):
    a = a % n
    pre = [1, a]
    g = a
    for i in range(2**k - 2):
        g *= a
        pre.append(g)
    r = 1
    j = (2**k - 1) << (b.bit_length() // k * k)
    l = b.bit_length() // k 
    while j > 0:
        for _ in range(k):
            r = r * r % n
        r = r * pre[(b & j) >> (l * k)] % n
        l -= 1
        j >>= k
    return r

#left to right k-ary exponentiation with bit-selection optimization
def expk2(a, b, n, k = 5):
        a = a % n
        m = 2 ** k - 1
        pre = [1, a]
        g = a
        for i in range(2**k - 2):
            g *= a
            pre.append(g)
        r = 1
        j = b.bit_length() // k * k
        while j >= 0:
            for _ in range(k):
                r = r * r % n
            i = (b >> j) & m
            if i:
                r = r * pre[i] % n
            j -= k
        return r

#elementary mod exp implementation
def exp_base(a, b, n):
    r = a
    for i in range(b - 1):
        r *= a
        r %= n
    return r

#optimization over exp_base that uses prime factor exponentiation
def exp(a, b, n):
    i = 0
    l = len(primos_init)
    a %= n
    while b != 1 and i < l:
        if b % primos_init[i] == 0:
            b //= primos_init[i]
            a = a ** primos_init[i] % n
        else:
            i += 1
    return a if i < l else exp_base(a, b, n)
