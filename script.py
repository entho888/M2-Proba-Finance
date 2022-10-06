#Libraries
import numpy as np 
import numpy.random as rd
import matplotlib.pyplot as plt

#Variables
N = 100
sample = rd.randn(N)

# plus tard : box muller et marsaglia

#relation prix comptant / prix forward 

def zero_coupon(taux, maturite) :
    return 1/((1 + taux)**maturite)

def prix_comptant_R(prix_forward, zero_coupon_v) :
    return prix_forward*zero_coupon_v

def prix_forward_R(prix_comptant, zero_coupon_v) :
    return prix_comptant/zero_coupon_v

#relation parite call / put 

def call_parity(put, zero_coupon_v, forward, strike) :
    return put + zero_coupon_v*(forward - strike)

def put_parity(call, zero_coupon_v, forward, strike) :
    return call - zero_coupon_v*(forward - strike)

def zero_coupon_parity(call, put, forward, strike) :
    return (call - put)/(forward - strike)

def forward_parity(call, put, zero_coupon_v, strike) :
    return (call - put)/zero_coupon_v + strike

## Bornes d'arbitrage 
# Relation de Hedge


# Relation Bull-Spread


# Butterfly : convexity check
def convexity_check(call_list, strike_list) :
    return 0

# Décroissance du call en maturité


## Formule de Carr


## Cox-Ross-Rubinstein
def identity_func(x) :
    return x

def Pascal_triangle(row) :
    Row1 = []
    Row0 = [1]
    for i in range(1, row) :
        Row1.append(1)
        for j in range(1, i) :
            Row1.append(Row0[j-1] + Row0[j])
        Row1.append(1)
        Row0 = Row1
        Row1 = []

    return Row0

def CRRd1(S0,r,u,d, h=identity_func) :
    hSu = h(S0*u)
    hSd = h(S0*d)
    return [(1/((1+r)*(u-d)))*(hSd*u - hSu*d) , (hSu - hSd)/(S0*(u-d)), (1/((1+r)*(u-d)))*(hSu*(1+r-d) + hSd*(u - (1+r)))]

def CRR(S0,T,r,u,d, h=identity_func) :
    hSu = h(S0*u)
    hSd = h(S0*d)

    Binomial_Coefficients = Pascal_triangle(T)

    p = (1+r-d)/(u-d)

    V = 0

    for i in range(T) :
        V+= Binomial_Coefficients[i]*(p**i)*((1-p)**(T-i))*h(S0*(u**i)*(d**(T-i)))

    return V/((1+r)**T)

    
    
    



# Simulation de brownien 
def standard_brownian_simulation(time_list, gaussian_list) :
    Brownian_list = [gaussian_list[0]]

    for i in range(1, len(time_list)) :
        t2 = time_list[i]
        t1 = time_list[i-1]

        Brownian_list.append(Brownian_list[i-1] + np.sqrt((t2-t1))*gaussian_list[i])

    return Brownian_list
