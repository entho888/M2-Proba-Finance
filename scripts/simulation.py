#%%
#Libraries
import numpy as np 
import matplotlib.pyplot as plt

rng = np.random.default_rng()

def uniform_distribution(low=float, high=float, size=None) :
    """Return a sample of size "size" of uniform distribution over the interval [low, high)

    Args:
        low (float): lower bound of the interval the generated numbers are in
        high (float): upper bound of the interval the generated numbers are in
        size (tupple): output shape. if the given shape is (m,n,k) then m*n*k samples are generated. 
            Default is None, in the case the function returns a single value.

    Returns:
        Array of random floats of shape "size".
    """
    if low > high : 
        print("Error from 'uniform_distribution' function of 'simulation' : lower bound of interval greater than upper bound")
        return 
    return (high - low)*rng.random(size=size) + low


def exponential_distribution(lamda=float, size=None) :
    """exponential_distribution : Return of shape "size" of exponential distribution 
    of paramter lambda = lamda

    Args:
        lamda (float, optional): non negativ float, parameter of the exponential distribution.
         Defaults to float.
        size (tupple, optional): output shape. if given shape is (m,,n,k) then m*n*k samples are generated.
         Defaults to None.

    Returns :   
        Array of random floats of shape "size".
    """
    if lamda <= 0 : 
        print("Error from 'exponential_distribution' function of 'simulation' : parameter lambda is negativ")
        return
    return (-1/lamda)*np.log(rng.random(size=size))


def bernoulli(p=float, size=None) :
    def jump_function(x) :
        if x < p : return 1
        else: return 0
    jump_function_vectorized = np.vectorize(jump_function)
    return jump_function_vectorized(rng.random(size=size))


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


def binomial(n=int, p=float, size=None) :
    binomial_coefficients = Pascal_triangle(n)
    cummulated_probability_list = [0]

    Sum = 0
    for k in range(n) :
        Sum += binomial_coefficients[k]*(p**(k))*((1-p)**(n-k))
        cummulated_probability_list.append(Sum)
    cummulated_probability_list.append(1)

    def multiple_jumps_function(x) :
        for k in range(len(cummulated_probability_list)-1) :
            if x >= cummulated_probability_list[k] and x < cummulated_probability_list[k+1] :
                return k

    multiple_jumps_function_vect = np.vectorize(multiple_jumps_function)

    return multiple_jumps_function_vect(rng.random(size=size))
    

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

"""
import time

start_time = time.time()
box_muller_method(10**7)
print("Box-Muller a mis %s secondes à tourner ---" % (time.time() - start_time))

start_time = time.time()
marsaglia_method(10**7)
print("Marsaglia a mis %s secondes à tourner ---" % (time.time() - start_time))
"""






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
