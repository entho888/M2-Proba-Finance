### QUIZZ 2 : 6/10/2022 ###
# Note : 6.21
# N'est pas le document utilisé ce jour là, contient des modifications apportées le 7/10/2022 après avoir vu la note. 

# On importe les bonnes librairies
import sys
sys.path.append('/workspace/Derivatives/scripts')
from derivatives import *

# Variables
S0 = 36.43
Su = 62.88
Sd = 20.34
u = Su/S0
d = Sd/S0
a =0.8
b = 7.15

def ATM_call(x):
    if x - S0 >= 0 : return x- S0
    else : return 0

# Résultats
print(CRRd1(S0, 0, u, d, ATM_call))
print((1 - Pascal_triangle(7)[4]/(7**4))*100)
print(a*(b-a) + 3*(a**2))


# Test personnels, rien à voir avec le Quizz
"""
print(np.cos(2)) # En important derivatives normalement avec un alias : ERREUR
print(deriv.cosinus(2)) # Fonctionne. Cosinus est une fonction définit dans derivatives qui fait x : np.cos(x).
print(deriv.np.cos(2)) # Fonctionne. Techniquement pas besoin d'importer numpy une seconde fois.
"""
print(np.cos(2)) # En important derivatives avec * : fonctionne