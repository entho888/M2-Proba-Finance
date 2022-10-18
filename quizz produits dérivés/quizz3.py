### QUIZZ 2 : 18/10/2022 ###
# Note : 8.5 
# Note si j'étais pas un bouffon : 17 ... 
# Note si j'étais réveillé et en forme : 20 .......................;

# On importe les bonnes librairies
import sys
sys.path.append('/workspace/M2-Proba-Finance/scripts')
import derivatives as deriv
import numpy as np
# import simulation as sim  # pas nécessaire 

##########################################################################################################################
# TIPS : double clique sur les trucs écrits en latex du quizz ouvre un page où on peut selectionner ces trucs et les c/c.
# Utile pour les données qui sont souvent en latex 
##########################################################################################################################

# Question 1
print('Q1 answer : ', 1000 - 2*(100 +80 + 8*8)) # J'aurais du le faire sur python car je me suis trompé sur le papier ... erreur de calcul

# Question 4
S0 = 94.62
T = 0.612
K = 98.62
r = 0
sigma = 0.2286
a = 40.51
Market_price = deriv.call_BS_at0(T, S0+a, K+a, sigma, r)

# J'ai pas vérifier l'erreur et j'avais mal implémenté l'algo bisection + pas croisé les sources avec Newton-Raphson ...
vol_imp1, error1, n_iteration1 = deriv.implied_volatility_bisection(S0, T, K, r, Market_price, max_ite=10**3, eps=10**(-7)) 
print('Q4 answers with bisection : ')
print('implied volatility = ', np.round(vol_imp1, decimals = 4))
print('error = ', error1)
vol_imp2, error2, n_iteration2 = deriv.implied_volatility_NewtonRaphson(S0, T, K, r, Market_price, max_ite=10**3, eps=10**(-6))
print('Q4 answers with Newton-Raphson : ')
print('implied volatility = ', np.round(vol_imp2, decimals = 4))
print('error = ', error2)


# Question 5
S0 = 91.6
T = 2.495
K = 92.68
r = 0.082

sigma_0 = np.sqrt((2/T)*np.abs(np.log(S0*np.exp(r*T)/K)))
print('Q5 answer : ', np.round(sigma_0, decimals=4))

# Question 6 
S0 = 99.87
T = 0.389
K = 90.43
r = 0.079
sigma = 0.2157  #dans le test j'ai mis 21.57 .............

call_BS = deriv.call_BS_at0(T, S0, K, sigma, r)
print('Q6 answer : ', np.round(call_BS, decimals=3))

# Question 7 (pas faite pendant le quizz)
print('Q7 answer : ', np.round(((6.45)**2)+ ((2.72)**3)/3, decimals=2))


# Question 8
S0 = 99.11
T = 0.963
K = 96.71
r = 0.007
sigma = 0.2771

put_BS = deriv.put_BS_at0(T, S0, K, sigma, r) 
print('Q8 answer : ', np.round(put_BS, decimals=3))

































































































































