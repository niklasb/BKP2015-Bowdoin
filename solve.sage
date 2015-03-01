# source: https://github.com/mimoo/RSA-and-LLL-attacks
def coppersmith_howgrave_univariate(pol, modulus, beta, mm, tt, XX):
    """
    Howgrave-Graham revisited method of Coppersmith following theorem:
    
    modulus integer of unknown factorization,
    0 < beta <= 1
    0 < epsilon <= beta / 7
    b | modulus and b >= modulus^beta,
    pol is a monic polynomial of degree dd,
    
    THEN
    
    we can find roots of pol(x) = 0 mod b with
    |root| <= |X|
    """
    #
    # init
    #
    dd = pol.degree()
    nn = dd * mm + tt

    #
    # checks
    #
    if not 0 < beta <= 1:
        raise ValueError("beta should belongs in (0, 1]")

    if not pol.is_monic():
        raise ArithmeticError("Polynomial must be monic.")

    #
    # calculate bounds and display them
    #
    """
    * we want to find g(x) such that ||g(xX)|| <= b^m / sqrt(n)

    * we know LLL will give us a short vector v such that:
    ||v|| <= 2^((n - 1)/4) * det(L)^(1/n)

    * we will use that vector as a coefficient vecotr for our g(x)
    
    * so we want to satisfy:
    2^((n - 1)/4) * det(L)^(1/n) < N^(beta*m) / sqrt(n)
    
    so we can obtain ||v|| < N^(beta*m) / sqrt(n) <= b^m / sqrt(n)
    (it's important to play with N because we might not know b)
    """
    
    # t optimized?
    print "\n# Optimized t?\n"
    print "we want X^(n-1) < N^(beta*m) so that each vector is helpful"
    cond1 = RR(XX^(nn-1))
    print "* X^(n-1) = ", cond1
    cond2 = pow(modulus, beta*mm)
    print "* N^(beta*m) = ", cond2
    print "* X^(n-1) < N^(beta*m) \n-> GOOD" if cond1 < cond2 else "* X^(n-1) >= N^(beta*m) \n-> NOT GOOD"
    
    # bound for X
    print "\n# X bound respected?\n"
    print "we want X <= N^(((2*beta*m)/(n-1)) - ((delta*m*(m+1))/(n*(n-1)))) / 2 = M"
    print "* X =", XX
    cond2 = RR(modulus^(((2*beta*mm)/(nn-1)) - ((dd*mm*(mm+1))/(nn*(nn-1)))) / 2)
    print "* M =", cond2
    print "* X <= M \n-> GOOD" if XX <= cond2 else "* X > M \n-> NOT GOOD"

    # solution possible?
    print "\n# Solutions possible?\n"
    detL = RR(modulus^(dd * mm * (mm + 1) / 2) * XX^(nn * (nn - 1) / 2))
    print "we can find a solution if 2^((n - 1)/4) * det(L)^(1/n) < N^(beta*m) / sqrt(n)"
    cond1 = RR(2^((nn - 1)/4) * detL^(1/nn))
    print "* 2^((n - 1)/4) * det(L)^(1/n) = ", cond1
    cond2 = RR(modulus^(beta*mm) / sqrt(nn))
    print "* N^(beta*m) / sqrt(n) = ", cond2
    print "* 2^((n - 1)/4) * det(L)^(1/n) < N^(beta*m) / sqrt(n) \n-> SOLUTION WILL BE FOUND" if cond1 < cond2 else "* 2^((n - 1)/4) * det(L)^(1/n) >= N^(beta*m) / sqroot(n) \n-> NO SOLUTIONS MIGHT BE FOUND (but we never know)"

    # warning about X
    print "\n# Note that no solutions will be found _for sure_ if you don't respect:\n* |root| < X \n* b >= modulus^beta\n"
    
    #
    # Coppersmith revisited algo for univariate
    #

    # change ring of pol and x
    polZ = pol.change_ring(ZZ)
    x = polZ.parent().gen()

    # compute polynomials
    gg = []
    for ii in range(mm):
        for jj in range(dd):
            gg.append((x * XX)**jj * modulus**(mm - ii) * polZ(x * XX)**ii)
    for ii in range(tt):
        gg.append((x * XX)**ii * polZ(x * XX)**mm)
    
    # construct lattice B
    BB = Matrix(ZZ, nn)

    for ii in range(nn):
        for jj in range(ii+1):
            BB[ii, jj] = gg[ii][jj]

    # LLL
    BB = BB.LLL()

    # transform shortest vector in polynomial    
    new_pol = 0
    for ii in range(nn):
        new_pol += x**ii * BB[0, ii] / XX**ii

    # factor polynomial
    potential_roots = new_pol.roots()

    # test roots on original pol
    roots = []
    for root in potential_roots:
        result = ZZ(pol(ZZ(root[0])))

        if gcd(modulus, result) >= modulus**beta:
            roots.append(root[0])

    # 
    return roots

############################################
# Test on Factoring with High Bits Known
##########################################

import os
RSATOOL = 'python2 ~/toolz/rsatool/rsatool.py'
alph = "abcdefghijklmnopqrstuvwxyz"
table = alph.upper() + alph.lower() + "0123456789+/"
rev = dict((c,i) for i,c in enumerate(table))

def decode(s):
    res = ""
    for c in s:
        res += bin(rev[c])[2:].rjust(6,'0')
    return ''.join(chr(int(res[i:i+8],2)) for i in xrange(0, len(res), 8))

def read_file(fname):
    with open(fname) as f:
        return ''.join(f.read().split())

s1 = decode(read_file('part1.txt'))
s2 = decode(read_file('part2.txt'))
n = int(s1[10:][:0x81].encode('hex'), 16)
qhi = int(s2[10:].encode('hex'),16)
XX = 1 << (8*(64-len(s2)+10))
qbar = qhi * XX

F.<x> = PolynomialRing(Zmod(n), implementation='NTL'); 
f = x - qbar;

# PLAY WITH THOSE:
beta = 0.4 # we should have q >= N^beta
dd = f.degree()
epsilon = beta / 7
mm = ceil(beta**2 / (dd * epsilon))
tt = floor(dd * mm * ((1/beta) - 1))
roots = coppersmith_howgrave_univariate(f, n, beta, mm, tt, XX)

# output
print "solutions:", roots
assert len(roots) > 0
q = qbar - roots[0]
p = n / q
assert p * q == n
print "Writing key to file privkey.pem"
os.system('%s -n %s -p %d -q %d -o privkey.pem' % (RSATOOL, n, p, q))
print "Extracting flag"
os.system('tar xf foundkeys.tar.gz.*')
print "Decrypting flag"
os.system('openssl rsautl -in distribute/flag.encrypted -inkey privkey.pem -decrypt')
