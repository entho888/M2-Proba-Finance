#%%
#Libraries
import numpy as np 
import numpy.random as rd
import matplotlib.pyplot as plt

def uniform(a,b,N) :
    return "pas fini, mais je sais pas si je fais"

def standard_brownian_simulation(time_list, gaussian_list) :
    Brownian_list = [gaussian_list[0]]

    for i in range(1, len(time_list)) :
        t2 = time_list[i]
        t1 = time_list[i-1]

        Brownian_list.append(Brownian_list[i-1] + np.sqrt((t2-t1))*gaussian_list[i])

    return Brownian_list