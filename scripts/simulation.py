#%%
#Libraries
import numpy as np 
import matplotlib.pyplot as plt

rng = np.random.default_rng()

def uniform(a,b,N) :
    """uniform _summary_

    _extended_summary_

    Args:
        a (_type_): _description_
        b (_type_): _description_
        N (_type_): _description_
    """
    return "on verra"


def standard_brownian_simulation(time_list, gaussian_list) :
    Brownian_list = [gaussian_list[0]]

    for i in range(1, len(time_list)) :
        t2 = time_list[i]
        t1 = time_list[i-1]

        Brownian_list.append(Brownian_list[i-1] + np.sqrt((t2-t1))*gaussian_list[i])

    return Brownian_list