#%%
#Libraries
import numpy as np 
import matplotlib.pyplot as plt

rng = np.random.default_rng()

def uniform(low=float, high=float, size=None) :
    """Return a sample of size "size" of uniform distribution over the interval [low, high)

    Args:
        low (float): lower bound of the interval the generated numbers are in
        high (float): upper bound of the interval the generated numbers are in
        size (tupple): output shape. if the given shape is (m*n*k) then m*n*k samples are generated. 
            Default is None, in the case the function returns a single value.

    Returns:
        Array of random floats of shape size.
    """
    if low > high : 
        print("Error from 'uniform' function of 'simulation' : lower bound of interval greater than upper bound")
        return 
    return (high - low)*rng.random(size=size) + low

def standard_brownian_simulation(time_list, gaussian_list) :
    Brownian_list = [gaussian_list[0]]

    for i in range(1, len(time_list)) :
        t2 = time_list[i]
        t1 = time_list[i-1]

        Brownian_list.append(Brownian_list[i-1] + np.sqrt((t2-t1))*gaussian_list[i])

    return Brownian_list

"""
times = np.linspace(0,100,N)
Brownian = standard_brownian_simulation(times, sample)
plt.plot(times, Brownian)
plt.show()
"""
# %%
