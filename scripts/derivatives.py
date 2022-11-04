#Libraries
import numpy as np 
import scipy
import matplotlib.pyplot as plt


#Variables
sigma = 0.3
r = 0.004
q = 0.002

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
    """Return the 'row_number'-th row of the pascal triangle (meaning is 'row_number'=n, it returns a list of length n+1)

    Args:
        row_number (int): positiv integer

    Returns:
        list: list of binomial coefficient of 'degree' 'row_number'.
    """
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
    """Cox Rox Rubinstein with 1 period

    Args:
        S0 (float): price of the asset at t=0
        r (float): interest rate
        u (float): increasement factor
        d (float): decreasement (?) factor
        h (function, optional): payoff. Defaults to identity_func.

    Returns:
        list : [amount of money placed in cash, parts of asset bought, price of the option]
    """
    hSu = h(S0*u)
    hSd = h(S0*d)
    return [(1/((1+r)*(u-d)))*(hSd*u - hSu*d) , (hSu - hSd)/(S0*(u-d)), (1/((1+r)*(u-d)))*(hSu*(1+r-d) + hSd*(u - (1+r)))]


def CRR(S0,T,r,u,d, h=identity_func) :
    """Cox Rox Rubinstein with T periods

    Args:
        S0 (float): price of the asset at t=0
        T (int): positiv integer. Number of periods 
        r (float): interest rate
        u (float): increasing factor
        d (float): descreasing factor
        h (function, optional): payoff. Defaults to identity_func.

    Returns:
        float: price of the option.
    """
    hSu = h(S0*u)
    hSd = h(S0*d)

    Binomial_Coefficients = Pascal_triangle(T)

    p = (1+r-d)/(u-d)

    V = 0

    for i in range(T+1) :
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

def d_plus(x, y, T, sigma=sigma) :
    return (1/(sigma*np.sqrt(T)))*np.log(x/y) + .5*sigma*np.sqrt(T)

def d_minus(x, y, T, sigma=sigma) :
    return (1/(sigma*np.sqrt(T)))*np.log(x/y) - .5*sigma*np.sqrt(T)

def d1(St, K, t, T, sigma=sigma, r=r) :
    return (1/(sigma*np.sqrt(T-t)))*(np.log(St/K) + (r + .5*(sigma**2))*(T-t))

def d2(St, K, t, T, sigma=sigma, r=r) :
    return d1(St, K, t, T, sigma, r) - sigma*np.sqrt(T-t)

def AS_standard_gaussian_cdf_approx(x) :
    """Abramowitz-Stegun gaussian cdf approximation

    Args:
        x (float): number when the cdf is evaluated

    Returns:
        float: float between 0 and 1. Approximation of the gaussian cdf evaluated in x?
    """
    b0 = 0.2316419
    b1 = 0.319381530
    b2 = -0.356563782
    b3 = 1.781477937
    b4 = -1.821255978
    b5 = 1.330274429
    t = 1/(1 + b0*x)

    return 1 - (1/np.sqrt(2*np.pi))*np.exp(-.5*(x**2))*(b1*(t**2) + b2*(t**2) + b3*(t**2) + b4*(t**2) + b5*(t**2))

def call_BS_at0(S0, K, T, sigma=sigma, r=r, q=q) :
    return S0*np.exp(-q*T)*scipy.stats.norm.cdf(d_plus(S0*np.exp((r-q)*T), K, T, sigma)) - K*np.exp(-r*T)*scipy.stats.norm.cdf(d_minus(S0*np.exp((r-q)*T), K, T, sigma))

def put_BS_at0(S0, K, T, sigma=sigma, r=r, q=q) :
    return K*np.exp(-r*T)*scipy.stats.norm.cdf(-d_minus(S0*np.exp((r-q)*T), K, T, sigma)) - S0*np.exp(-q*T)*scipy.stats.norm.cdf(-d_plus(S0*np.exp(r*T), K, T, sigma))
    
def call_BS_at_t(St, K, t, T, sigma=sigma, r=r, q=q) :
    return call_BS_at0(St, K, T-t, sigma, r, q)

def put_BS_at_t(St, K, t, T, sigma=sigma, r=r, q=q) :
    return put_BS_at0(St, K, T-t, sigma, r, q)

# Sensibilités (delta, vega, gamma, theta)
def delta_call_BS_at0(S0, K, T, sigma=sigma, r=r, q=q) : return np.exp(-q*T)*scipy.stats.norm.cdf(d_plus(S0*np.exp((r-q)*T), K, T, sigma))

def delta_put_BS_at0(S0, K, T, sigma=sigma, r=r, q=q) : return np.exp(-q*T)*(scipy.stats.norm.cdf(d_plus(S0*np.exp((r-q)*T), K, T, sigma)) - 1)

def delta_call_BS_at_t(St, K, t, T, sigma=sigma, r=r, q=q) : return delta_call_BS_at0(St, K, T-t, sigma, r, q)

def delta_put_BS_at_t(St, K, t, T, sigma=sigma, r=r, q=q) : return delta_put_BS_at0(St, K, T-t, sigma, r, q)

def gamma_BS_at0(S0, K, T, sigma=sigma, r=r, q=q) : return (1/(S0*sigma*np.sqrt(T)))*np.exp(-q*T)*normal_density(d_plus(S0*np.exp((r-q)*T), K, T, sigma))

def gamma_BS_at_t(St, K, t, T, sigma=sigma, r=r, q=q) : return gamma_BS_at0(St, K, T-t, sigma, r, q)

def delta_strike_call_BS_at0(S0, K, T, sigma=sigma, r=r, q=q) : return -np.exp(-r*T)*scipy.stats.norm.cdf(d_minus(S0*np.exp((r-q)*T), K, T, sigma))

def delta_strike_put_BS_at0(S0, K, T, sigma=sigma, r=r, q=q) : return np.exp(-r*T)*(1 - scipy.stats.norm.cdf( d_minus(S0*np.exp((r-q)*T), K, T, sigma)) )

def delta_strike_call_BS_at_t(St, K, t, T, sigma=sigma, r=r, q=q) : return delta_call_BS_at0(St, K, T-t, sigma, r, q)

def delta_strike_put_BS_at_t(St, K, t, T, sigma=sigma, r=r, q=q) : return delta_put_BS_at0(St, K, T-t, sigma, r, q)

def vega_BS_at0(S0, K, T, sigma=sigma, r=r, q=q) : return S0*normal_density(d_plus(S0*np.exp((r-q)*T), K, T, sigma))*np.sqrt(T)*np.exp(-q*T)

def vega_BS_at_t(St, K, t, T, sigma=sigma, r=r, q=q) : return vega_BS_at0(St, K, T-t, sigma, r, q)

def theta_Call_BS_at0(S0, K, T, sigma=sigma, r=r, q=q) : 
    return np.exp(-q*T)*(S0*sigma)/(2*np.sqrt(T))*normal_density(d_plus(S0*np.exp((r-q)*T), K, T, sigma)) + r*K*np.exp(-r*T)*scipy.stats.norm.cdf(d_minus(S0*np.exp((r-q)*T), K, T, sigma)) -q*S0*np.exp(-q*T)*scipy.stats.norm.cdf(d_plus(S0*np.exp((r-q)*T), K, T, sigma))

def theta_Put_BS_at0(S0, K, T, sigma=sigma, r=r, q=q) : 
    return theta_Call_BS_at0(S0, K, T, sigma, r) + q*S0*np.exp(-q*T) -r*K*np.exp(-r*T)

def theta_Call_BS_at_t(St, K, t, T, sigma=sigma, r=r, q=q)  : 
    return theta_Call_BS_at0(St, K, T-t, sigma, r, q)

def theta_Put_BS_at_t(St, K, t, T, sigma=sigma, r=r, q=q)  : 
    return theta_Put_BS_at0(St, K, T-t, sigma, r, q)

def rho_Call_BS_at0(S0, K, T, sigma=sigma, r=r, q=q) : return T*K*np.exp(-r*T)*scipy.stats.norm.cdf(d_minus(S0*np.exp((r-q)*T), K, T, sigma))

def rho_Put_BS_at0(S0, K, T, sigma=sigma, r=r, q=q) : return -T*K*np.exp(-r*T)*(1 - scipy.stats.norm.cdf(d_minus(S0*np.exp((r-q)*T), K, T, sigma)))

def rho_Call_BS_at_t(St, K, t, T, sigma=sigma, r=r, q=q) :
    return rho_Call_BS_at0(St, K, T-t, sigma, r, q)

def rho_Put_BS_at_t(St, K, t, T, sigma=sigma, r=r, q=q) :
    return rho_Put_BS_at0(St, K, T-t, sigma, r, q)

### Vol implicite 
# Dichotomie -> vol implicite
def implied_volatility_bisection(S0, T, K, Market_price, max_ite, eps=10**(-5), r=r, q=q) :
    x, y = 0, 1
    i = 0
    error = 2*eps

    while i <= max_ite or error > eps :
        z = (x + y)/2
        call_BS = call_BS_at0(T, S0, K, z/(1-z), r, q)

        if Market_price > call_BS : x, y = z, y
        if Market_price < call_BS : x, y = x, z
        
        error = np.abs(Market_price - call_BS)
        i += 1

    if i > max_ite :
        print("'implied_volatility_bisection' : number of iterations too high")
        return z/(1-z), error, i
    return z/(1-z), error, i


# Descente de gradients -> vol implicite
def implied_volatility_NewtonRaphson(S0, T, K, Market_price, max_ite, eps=10**(-5), r=r, q=q) :
    sigma_ = np.sqrt((2/T)*np.abs(np.log(S0*np.exp((r-q)*T)/K))) #PAS SUR DU q ICI !!!!!!
    i = 0
    call_BS = call_BS_at0(T, S0, K, sigma_, r, q)
    error = np.abs(Market_price - call_BS)

    while i <= max_ite or error > eps :
        sigma_ = sigma_ - (call_BS - Market_price)/vega_BS_at0(S0, K, T, sigma_, r, q)
        
        call_BS = call_BS_at0(T, S0, K, sigma_, r, q)

        error = np.abs(Market_price - call_BS)
        i += 1

    if i > max_ite :
        print("'implied_volatility_NewtonRaphson' : number of iterations too high")
        return sigma_, error, i
    return sigma_, error, i


# Estimateur de la variance 
def volatility_approximation(asset_values_sample, times) :
    sum = 0

    for i in range(1, len(times)) :
        sum += ((np.log(asset_values_sample[i]/asset_values_sample[i-1]))**2)/(times[i]-times[i-1])

    return sum


##### OPTIONS BARRIERES (besoin de créer une classe black scholes parce que ça commence à être chiant)
def BinCall_at0(S0, K, T, sigma=sigma, r=r, q=q) : 
    return np.exp(-r*T)*scipy.stats.norm.cdf(d_minus(S0*np.exp((r-q)*T), K, T, sigma))

def BinPut_at0(S0, K, T, sigma=sigma, r=r, q=q) : 
    return np.exp(-r*T)*scipy.stats.norm.cdf(-d_minus(S0*np.exp((r-q)*T), K, T, sigma))

def BinCall_at_t(St, K, t, T, sigma=sigma, r=r, q=q) : return BinCall_at0(St, K, T-t, sigma, r, q)

def BinPut_at_t(St, K, t, T, sigma=sigma, r=r, q=q) : return BinPut_at0(St, K, T-t, sigma, r, q)

def DIC_regular(St, D, K, t, T, sigma=sigma, r=r, q=q) :
    if D > K :
        print('Erreur dans DIC_regular : D > K pas regular')
        return
    nu = r - q
    gamma = 1 - (2*nu)/(sigma**2)

    if St > D :
        return ((St/D)**(gamma - 1))*call_BS_at_t(D, K*St/D, t, T, sigma=sigma, r=r, q=q)
    else : 
        return call_BS_at_t(St, K, t, T, sigma=sigma, r=r, q=q)

def BinDIC_regular(St, D, K, t, T, sigma=sigma, r=r, q=q) :
    if D > K :
        print('Erreur dans BinDIC_regular : D > K pas regular')
        return
    nu = r - q
    gamma = 1 - (2*nu)/(sigma**2)

    if St > D :
        return ((St/D)**(gamma))*BinCall_at_t(D, K*St/D, t, T, sigma=sigma, r=r, q=q)
    else : 
        return BinCall_at_t(St, K, t, T, sigma=sigma, r=r, q=q)

def BinDIC(St, D, K, t, T, sigma=sigma, r=r, q=q) :
    if St <= D : return BinCall_at_t(St, K, t, T, sigma=sigma, r=r, q=q)
    if D <= K : 
        return BinDIC_regular(St, D, K, t, T, sigma=sigma, r=r, q=q)

    else :
        return BinCall_at_t(St, K, t, T, sigma=sigma, r=r, q=q) - np.exp(-r*(T-t)) + BinPut_at_t(St, D, t, T, sigma=sigma, r=r, q=q) + BinDIC_regular(St, D, D, t, T, sigma=sigma, r=r, q=q)

def DIC(St, D, K, t, T, sigma=sigma, r=r, q=q) :
    if St <= D : return call_BS_at_t(St, K, t, T, sigma=sigma, r=r, q=q)
    if D <= K : return DIC_regular(St, D, K, t, T, sigma=sigma, r=r, q=q)
    return call_BS_at_t(St, K, t, T, sigma=sigma, r=r, q=q) - call_BS_at_t(St, D, t, T, sigma=sigma, r=r, q=q) - (D-K)*BinCall_at_t(St, D, t, T, sigma=sigma, r=r, q=q) + DIC_regular(St, D, D, t, T, sigma=sigma, r=r, q=q) + (D-K)*BinDIC_regular(St, D, D, t, T, sigma=sigma, r=r, q=q)

def DOC(St, D, K, t, T, sigma=sigma, r=r, q=q) :
    return call_BS_at_t(St, K, t, T, sigma=sigma, r=r, q=q) - DIC(St, D, K, t, T, sigma=sigma, r=r, q=q)

def BinDOC(St, D, K, t, T, sigma=sigma, r=r, q=q) :
    return BinCall_at_t(St, K, t, T, sigma=sigma, r=r, q=q) - BinDIC(St, D, K, t, T, sigma=sigma, r=r, q=q)

def call_power_asset(St, K, t, T, n, sigma=sigma, r=r, q=q) :
    d1 = (np.log(K/(St**n)) - (T-t)*n*(r-q-(sigma**2)/2))/(sigma*n*np.sqrt(T-t))
    d2 = d1 - sigma*n*np.sqrt(T-t)
    return (St**n)*np.exp((T-t)*(-r+r*n-q*n-((sigma**2)*(n/2)) + (sigma**2)*(n**2)/2 ))*(1 - scipy.stats.norm.cdf(d2)) - K*np.exp(-r*(T-t))*(1 - scipy.stats.norm.cdf(d1))


#Test 
"""
K = 90
T = 4
S0 = 100
sigma = 0.3
r = 0.004
q = 0.002
D = 85

L = []
L.append(BinPut_at0(S0, K, T))
L.append(BinCall_at0(S0, K, T))
L.append(DIC(S0, D, K, 0, T))
L.append(BinDIC(S0, D, K, 0, T))
L.append(DOC(S0, D, K, 0, T))
print(L)

D = 95

L = []
L.append(BinPut_at0(S0, K, T))
L.append(BinCall_at0(S0, K, T))
L.append(DIC(S0, D, K, 0, T))
L.append(BinDIC(S0, D, K, 0, T))
L.append(DOC(S0, D, K, 0, T))
print(L)
"""

### Modèle log décalé
# Pricer
# Sensi

### Autre modèle 
# Pricer
# Sensi

### Methode de Monte-Carlo, intervalle de confiance, etc 