#Libraries
import numpy as np 
import scipy
import matplotlib.pyplot as plt


#Variables

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


################################################################################################################
################################################################################################################
### Modèle de BS
################################################################################################################
################################################################################################################
# Prix d'un call modèle de BS -> portefeuille dynamique
def normal_density(x) :
    return (1/np.sqrt(2*np.pi))*np.exp(-.5*(x**2))

def d_plus(x, y, T, sigma) :
    return (1/(sigma*np.sqrt(T)))*np.log(x/y) + .5*sigma*np.sqrt(T)

def d_minus(x, y, T, sigma) :
    return (1/(sigma*np.sqrt(T)))*np.log(x/y) - .5*sigma*np.sqrt(T)

def d1(St, K, t, T, sigma, r) :
    return (1/(sigma*np.sqrt(T-t)))*(np.log(St/K) + (r + .5*(sigma**2))*(T-t))

def d2(St, K, t, T, sigma, r) :
    return d1(St, K, t, T, sigma, r) - sigma*np.sqrt(T-t)

def AS_standard_gaussian_cdf_approx(x) :
    b0 = 0.2316419
    b1 = 0.319381530
    b2 = -0.356563782
    b3 = 1.781477937
    b4 = -1.821255978
    b5 = 1.330274429
    t = 1/(1 + b0*x)

    return 1 - (1/np.sqrt(2*np.pi))*np.exp(-.5*(x**2))*(b1*(t**2) + b2*(t**2) + b3*(t**2) + b4*(t**2) + b5*(t**2))

def call_BS_at0(T, S0, K, sigma, r) :
    return S0*scipy.stats.norm.cdf(d_plus(S0*np.exp(r*T), K, T, sigma)) - K*np.exp(-r*T)*scipy.stats.norm.cdf(d_minus(S0*np.exp(r*T), K, T, sigma))

def put_BS_at0(T, S0, K, sigma, r) :
    return K*np.exp(-r*T)*scipy.stats.norm.cdf(-d_minus(S0*np.exp(r*T), K, T, sigma)) - S0*scipy.stats.norm.cdf(-d_plus(S0*np.exp(r*T), K, T, sigma))
    
def call_BS_at_t(St, K, t, T, sigma, r) :
    return St*scipy.stats.norm.cdf(d1(St, K, t, T, sigma, r)) - K*np.exp(-r*(T-t))*scipy.stats.norm.cdf(d2(St, K, t, T, sigma, r))

def put_BS_at_t(St, K, t, T, sigma, r) :
    return K*np.exp(-r*(T-t))*scipy.stats.norm.cdf(-d2(St, K, t, T, sigma, r)) - St*scipy.stats.norm.cdf(-d1(St, K, t, T, sigma, r)) 

# Sensibilités (delta, vega, gamma, theta)
def delta_call_BS_at0(S0, K, T, sigma, r) : return d_plus(S0*np.exp(r*T), K, T, sigma)

def delta_put_BS_at0(S0, K, T, sigma, r) : return d_plus(S0*np.exp(r*T), K, T, sigma) - 1

def delta_call_BS_at_t(St, K, t, T, sigma, r) : return d1(St, K, t, T, sigma, r)

def delta_put_BS_at_t(St, K, t, T, sigma, r) : return d1(St, K, t, T, sigma, r) - 1

def gamma_BS_at0(S0, K, T, sigma, r) : return (1/(S0*sigma*np.sqrt(T)))*normal_density(d_plus(S0*np.exp(r*T), K, T, sigma))

def gamma_BS_at_t(St, K, t, T, sigma, r) : return (1/(S0*sigma*np.sqrt(T-t)))*normal_density(d1(St, K, t, T, sigma, r))

def delta_strike_call_BS_at0(S0, K, T, sigma, r) : return -np.exp(-r*T)*d_minus(S0*np.exp(r*T), K, T, sigma)

def delta_strike_put_BS_at0(S0, K, T, sigma, r) : return np.exp(-r*T)*(1 - d_minus(S0*np.exp(r*T), K, T, sigma) )

def vega_BS_at0(S0, K, T, sigma, r) : return S0*normal_density(d_plus(S0*np.exp(r*T), K, T, sigma))*np.sqrt(T)

def vega_BS_at_t(St, K, t, T, sigma, r) : return St*normal_density(d1(St, K, t, T, sigma, r))*np.sqrt(T-t)

def theta_Call_BS_at0(S0, K, T, sigma, r) : 
    return (S0*sigma)/(2*np.sqrt(T))*normal_density(d_plus(S0*np.exp(r*T), K, T, sigma)) + r*K*np.exp(-r*T)*scipy.stats.norm.cdf(d_minus(S0*np.exp(r*T), K, T, sigma))

def theta_Put_BS_at0(S0, K, T, sigma, r) : 
    return (S0*sigma)/(2*np.sqrt(T))*normal_density(d_plus(S0*np.exp(r*T), K, T, sigma)) - r*K*np.exp(-r*T)*(1 - scipy.stats.norm.cdf(d_minus(S0*np.exp(r*T), K, T, sigma)))

def theta_Call_BS_at_t(St, K, t, T, sigma, r)  : 
    return (St*sigma)/(2*np.sqrt(T-t))*normal_density(d1(St, K, t, T, sigma, r)) + r*K*np.exp(-r*(T-t))*scipy.stats.norm.cdf(d2(St, K, t, T, sigma, r))

def theta_Put_BS_at_t(St, K, t, T, sigma, r)  : 
    return (St*sigma)/(2*np.sqrt(T-t))*normal_density(d1(St, K, t, T, sigma, r)) - r*K*np.exp(-r*(T-t))*(1 - scipy.stats.norm.cdf(d2(St, K, t, T, sigma, r)))

def rho_Call_BS_at0(S0, K, T, sigma, r) : return T*K*np.exp(-r*T)*scipy.stats.norm.cdf(d_minus(S0*np.exp(r*T), K, T, sigma))

def rho_Put_BS_at0(S0, K, T, sigma, r) : return -T*K*np.exp(-r*T)*(1 - scipy.stats.norm.cdf(d_minus(S0*np.exp(r*T), K, T, sigma)))

def rho_Call_BS_at_t(St, K, t, T, sigma, r) :
    return (T-t)*K*np.exp(-r*(T-t))*scipy.stats.norm.cdf(d2(St, K, t, T, sigma, r))

def rho_Put_BS_at_t(St, K, t, T, sigma, r) :
    return -(T-t)*K*np.exp(-r*(T-t))*(1 - scipy.stats.norm.cdf(d2(St, K, t, T, sigma, r)))


### Vol implicite 
# Dichotomie -> vol implicite
def implied_volatility_bisection(S0, T, K, r, Market_price, max_ite, eps) :
    x, y = 0, 1
    i = 0
    error = 2*eps

    while i <= max_ite or error > eps :
        z = (x + y)/2
        call_BS = call_BS_at0(T, S0, K, z/(1-z), r)

        if Market_price > call_BS : x, y = z, y
        if Market_price < call_BS : x, y = x, z
        
        error = np.abs(Market_price - call_BS)
        i += 1

    if i > max_ite :
        print("'implied_volatility_bisection' : number of iterations too high")
        return z
    return z


# Descente de gradients -> vol implicite
def implied_volatility_NewtonRaphson(S0, T, K, r, Market_price, max_ite, eps) :
    sigma_ = np.sqrt((2/T)*np.abs(np.log(S0*np.exp(r*T)/K)))
    i = 0
    call_BS = call_BS_at0(T, S0, K, sigma_, r)
    error = np.abs(Market_price - call_BS)

    while i <= max_ite or error > eps :
        sigma_ = sigma_ - (call_BS - Market_price)/vega_BS_at0(S0, K, T, sigma_, r)
        
        call_BS = call_BS_at0(T, S0, K, sigma_, r)

        error = np.abs(Market_price - call_BS)
        i += 1

    if i > max_ite :
        print("'implied_volatility_NewtonRaphson' : number of iterations too high")
        return sigma_
    return sigma_


# Estimateur de la variance 
def volatility_approximation(asset_values_sample, times) :
    sum = 0

    for i in range(1, len(times)) :
        sum += ((np.log(asset_values_sample[i]/asset_values_sample[i-1]))**2)/(times[i]-times[i-1])

    return sum



### Modèle log décalé
# Pricer
# Sensi

### Autre modèle 
# Pricer
# Sensi

### Methode de Monte-Carlo, intervalle de confiance, etc 