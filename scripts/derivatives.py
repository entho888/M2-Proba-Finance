#Libraries
import numpy as np 
import numpy.random as rd
import matplotlib.pyplot as plt

#Variables
N = 10**5
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


def Pascal_triangle(row_number) :
    present_Row = []
    past_Row = [1]
    for i in range(1, row_number+1) :
        present_Row.append(1)
        for j in range(1, i) :
            present_Row.append(past_Row[j-1] + past_Row[j])
        present_Row.append(1)
        past_Row = present_Row
        present_Row = []

    return past_Row


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

