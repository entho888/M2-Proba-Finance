### QUIZZ  4 : 27/10/2022 ###
# Note : 

# On importe les bonnes librairies
import sys
sys.path.append('/workspace/M2-Proba-Finance/scripts')
import derivatives as deriv
import simulation as sim  
import numpy as np

##########################################################################################################################
# TIPS : double clique sur les trucs écrits en latex du quizz ouvre un page où on peut selectionner ces trucs et les c/c.
# Utile pour les données qui sont souvent en latex 
##########################################################################################################################

answers = open("/workspace/M2-Proba-Finance/quizz produits dérivés/quizz_4_answers.txt", "w")

# Question 1

t=0
n=5
S_0=1.080000
T=0.718000
K=1.060000
r=0.040000
q=0.091000
sigma=0.132600

A = deriv.call_power_asset(S_0, K, 0, T, n, sigma=sigma, r=r, q=q)

print('Q1 answer : ', np.round(  A , decimals = 4  ))
answers.write('Q1 answer : ' + str( np.round( A  , decimals = 4  )))


# Question 2

S_0=90.610000
T=0.600000
K=102.440000
r=0.071000
q=0.086000
sigma=0.230800

A = deriv.BinCall_at0(S_0, K, T, sigma = sigma, r=r, q=q)

print('Q2 answer : ', np.round( A  , decimals = 4 ))
answers.write('\nQ2 answer : ' + str( np.round( A  , decimals =  4 )))


# Question 3

S_0=99.300000
K=100.130000
r=0.003000
p=0.600000
x=2.590000
y=1.620000
u=1 + (x/100)
d=1 - (y/100)

def payoff1(x) :
    if x - K >= 0 : return x - K
    else : return 0

A = deriv.CRRd1(S_0, r, u, d, h=payoff1)

print('Q3 answer : ', np.round( A  , decimals = 3 ))
answers.write('\nQ3 answer : ' + str( np.round(  A , decimals =  3 )))


# Question 4

S_0=92.910000
T=0.452000
K=90.970000
r=0.009000
q=0.025000
sigma=0.116200

A = deriv.vega_BS_at0(S_0, K, T, sigma=sigma, r=r, q=q)

print('Q4 answer : ', np.round(  A , decimals =  2))
answers.write('\nQ4 answer : ' + str( np.round(  A , decimals =  2 )))


# Question 5

S_0=101.230000
T=0.726000
K=103.670000
D=98.860000
r=0.016000
q=0.045000
sigma=0.172800

A = deriv.DIC(S_0, D, K, 0, T, sigma=sigma, r=r, q=q)

print('Q5 answer : ', np.round(  A , decimals = 2 ))
answers.write('\nQ5 answer : ' + str( np.round( A  , decimals =  2 )))


# Question 6

"""
print('Q6 answer : ', np.round(   , decimals =  ))
answers.write('\nQ6 answer : ' + str( np.round(   , decimals =   )))
"""

# Question 7

"""
print('Q7 answer : ', np.round(   , decimals =  ))
answers.write('\nQ7 answer : ' + str( np.round(   , decimals =   )))
"""

# Question 8

"""
print('Q8 answer : ', np.round(   , decimals =  ))
answers.write('\nQ8 answer : ' + str( np.round(   , decimals =   )))
"""

# Question 9

"""
print('Q9 answer : ', np.round(   , decimals =  ))
answers.write('\nQ9 answer : ' + str( np.round(   , decimals =   )))
"""

# Question 10

"""
print('Q10 answer : ', np.round(   , decimals =  ))
answers.write('\nQ10 answer : ' + str( np.round(   , decimals =   )))
"""

# Question 11

"""
print('Q11 answer : ', np.round(   , decimals =  ))
answers.write('\nQ11 answer : ' + str( np.round(   , decimals =   )))
"""

# Question 12

"""
print('Q12 answer : ', np.round(   , decimals =  ))
answers.write('\nQ12 answer : ' + str( np.round(   , decimals =   )))
"""

answers.close()