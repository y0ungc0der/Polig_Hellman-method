import time
from sage.all import *

# ------- Pohlig-Hellman method ------- #

# p - field characteristic
# F_p = 9103

F_p = 719

# x, y - point coordinates
# x_p = 4515
# y_p = 5197
# x_q = 8132
# y_q = 7165

x_p = 107
y_p = 443
x_q = 608
y_q = 427

# q - group order
# q = 9108

q = 699

# q_j - group order decomposition
#q_j = [[2, 2], [3, 2], [11, 1], [23, 1]]

q_j = [[3, 1], [233, 1]]

# a4, a6 - coefficients of the curve equation
# a_4 = 3
# a_6 = 16

a_4 = 130
a_6 = 565


def main():
    # elliptic curve defined by y^2 = x^3 + a_4*x + a_6 over finite field of size F_p
    E = EllipticCurve(GF(F_p), [a_4, a_6])
    # point P of the elliptic curve
    P = E(x_p, y_p)
    # point Q of the elliptic curve E: Q = dP 
    #Q = E(x_q, y_q)
    Q = P * 83
    # number of factors in the decomposition of the value q
    l_q = len(q_j)
    # d_j = d (mod ⁡ (q_j)^(a_j))
    d_j = []
    moduli = []
    
    for j in range(l_q):
        print (f'\n{j + 1} - iteration')
        dj = 0
        z = []
        p = q_j[j][0] # p ← q_j
        a = q_j[j][1] # a ← a_j 
        S = E(0)
        z.append(0) # z[0] = z_(-1)
        P_0 = (q // p) * P
        print ('P_0 = ', P_0)

        for k in range(a):

            if (k == 0):
                S = S + (z[k] // p) * P
                print ('S = ', S)

            else:
                S = S + z[k] * pow(p, k - 1) * P
                print ('S = ', S)

            Q_k = (q // pow(p, k + 1)) * (Q - S)
            print (f'Q_{k} = ', Q_k)
            z.append(P_0.discrete_log(Q_k)) # z[k + 1] = log_(P_0)Q_k

        for i in range(a):
            dj += z[i + 1] * pow(p, i)

        d_j.append(dj, )
        moduli.append(pow(p, a))

    # using the Chinese remainder theorem, reconstruct the value of the logarithm
    print (U'\nChinese remainder theorem:')
    for i in range(l_q):
        print(f'    d_{i + 1} ≡ {d_j[i]}(mod {moduli[i]})')
    d = CRT_list(d_j, moduli)
    
    if (P * d == Q):
        return (d)
    else:
        return (-1) 
    
if __name__ == '__main__':
    start_time = time.time()
    d = main()
    if (d != -1):
        print(f'Result: logarithm value d ≡ {d}(mod {q})')
        print(f'\nRuntime = {time.time() - start_time} seconds\n')
    else:
        print('!ERROR!')