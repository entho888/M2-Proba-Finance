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


def box_muller_method(size: int=1) :
    """Implementation of the Box-Muller method to draw independent samples of standard normal distribution 
    (centered and variance is identity). 

    Args:
        size (non negativ integer): number of samples draw
    """
    
    try :
        iteration_number = int(np.ceil(size/2)) #number of iteration of the method needed
        uniform_sample = rng.random((iteration_number,2))
        normal_sample = []

        for i in range(iteration_number) :
            u1, u2 = uniform_sample[i,]
            exponential = (-2*np.log(u1))**(.5)
            uniform_angle = 2*np.pi*u2

            normal_sample.extend((exponential*np.cos(uniform_angle), exponential*np.sin(uniform_angle)))

        return normal_sample[:size]

    except TypeError :
        print("Enter an int")
        raise
    except ValueError :
        print("Enter a non negativ int")
        raise


def marsaglia_method(size: int=1) :
    """Implementation of the Marsaglia's polar method to draw independent samples of standard normal distribution 
    (centered and variance is identity). 

    Args:
        size (non negativ integer): number of samples draw
    """

    try :
        iteration_number = int(np.ceil(size/2)) #number of iteration of the method needed
        uniform_sample = uniform(-1,1,(iteration_number,2))
        normal_sample = []

        for i in range(iteration_number) :
            x, y = uniform_sample[i,]
            r = x**2 + y**2
            temp = (-2*np.log(r)/r)**(.5)

            normal_sample.extend((temp*x, temp*y))

        return normal_sample[:size]

    except TypeError :
        print("Enter an int")
        raise
    except ValueError :
        print("Enter a non negativ int")
        raise


import time

start_time = time.time()
box_muller_method(10**7)
print("Box-Muller a mis %s secondes à tourner ---" % (time.time() - start_time))

start_time = time.time()
marsaglia_method(10**7)
print("Marsaglia a mis %s secondes à tourner ---" % (time.time() - start_time))







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
