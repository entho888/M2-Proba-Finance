#%%
#Libraries
import numpy as np 
import matplotlib.pyplot as plt
import scipy


sq = np.random.SeedSequence()
seed = sq.entropy
print('seed = ', seed)
rng = np.random.default_rng(sq)

def uniform_distribution(low=float, high=float, size=None, random_state: np.random.Generator=rng) :
    """Return an independent sample of size "size" of uniform distribution over 
    the interval [low, high).

    Args:
        low (float): lower bound of the interval the generated numbers are in
        high (float): upper bound of the interval the generated numbers are in
        size (tuple): output shape. if the given shape is (m,n,k) then m*n*k samples are generated. 
            Default is None, in the case the function returns a single value.
        random_state : 'np.random.Generator' used for simulation. Defaults to rng.

    Returns:
        Array of random floats of shape "size".
    """
    if low > high : 
        print("Error from 'uniform_distribution' function of 'simulation' : lower bound of interval greater than upper bound")
        return 

    return (high - low)*random_state.random(size=size) + low


def exponential_distribution(lamda=float, size=None, random_state: np.random.Generator=rng) :
    """exponential_distribution : Return an independent sample of shape "size" of 
    exponential distribution of parameter lambda = lamda. 
    Density is lamda*np.exp(-lamda*x).
    Generated by the inverse distribution function method.

    Args:
        lamda (float, optional): non negativ float, parameter of the exponential distribution.
            Defaults to float.
        size (tuple, optional): output shape. if given shape is (m,,n,k) then m*n*k samples are generated.
            Defaults to None.
        random_state : 'np.random.Generator' used for simulation. Defaults to rng.

    Returns :   
        Array of random floats of shape "size".
    """

    if lamda <= 0 : 
        print("Error from 'exponential_distribution' function of 'simulation' : parameter lambda is negativ")
        return

    return (-1/lamda)*np.log(random_state.random(size=size))


def bernoulli(p=float, size=None, random_state: np.random.Generator=rng) :
    """Bernoulli distribution

    Return an array of shape 'size' of independent random variables of law bernoulli('p').

    Args:
        p (float, optional): Float between 0 and 1. Parameter of the law. Defaults to float.
        size (tuple, optional): shape of the output. Defaults to None.
        random_state : 'np.random.Generator' used for simulation. Defaults to rng.
    """
    if (p < 0) or (p >1) :
        print("Error in 'binomial' function of 'simulation' : p must be include in [0,1]")
        return

    def jump_function(x) :
        if x < p : return 1
        else: return 0

    jump_function_vectorized = np.vectorize(jump_function)

    return jump_function_vectorized(random_state.random(size=size))


def Pascal_triangle(row_number) :
    """Pascal_triangle : Return a list of all binomial coefficients of 'degree' equal to
    'row_number'. Let's say row_number is equal to 2 then it returns [1,2,1].
    It uses the pascal triangle relation.

    Args:
        row_number (int): positiv integer, is equal to the 'degree' in the 
        binomial relation of Newton

    Returns:
        List : list of the binomial coefficents of 'degree' equal to 'row_number'
    """

    if row_number < 0 : 
        print("Error from 'Pascal_triangle' function of 'simulation' : row_number is negativ")
        return

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


def binomial(n=int, p=float, size = None, random_state: np.random.Generator=rng) :
    """binomial : Return an independent sample of shape 'size' of binomial 
    distribution of parameter 'n' and 'p'. 
    Generated by the inverse distribution function method.

    Args:
        n (int, optional): non negativ integer, maximum number the distribution may take.
         Defaults to int.
        p (float, optional): float between 0 and 1, parameter of the law.
         Defaults to float.
        size (tuple, optional): output shape. if given shape is (m,,n,k) then m*n*k samples are generated.
         Defaults to None.
        random_state : 'np.random.Generator' used for simulation. Defaults to rng.

    Returns:
        Array: array of random positiv integer lower than 'n' following the binomial(n,p) distribution.
    """

    if (n <= 0) :
        print("Error in 'binomial' function of 'simulation' : n is negativ")
        return
    if (p < 0) or (p >1) :
        print("Error in 'binomial' function of 'simulation' : p must be include in [0,1]")
        return

    binomial_coefficients = Pascal_triangle(n)
    cumulated_probability_list = [0]

    Sum = 0
    for k in range(n) :
        Sum += binomial_coefficients[k]*(p**(k))*((1-p)**(n-k))
        cumulated_probability_list.append(Sum)
    cumulated_probability_list.append(1)

    def multiple_jumps_function(x) :
        for k in range(len(cumulated_probability_list)-1) :
            if x >= cumulated_probability_list[k] and x < cumulated_probability_list[k+1] :
                return k

    multiple_jumps_function_vect = np.vectorize(multiple_jumps_function)

    return multiple_jumps_function_vect(random_state.random(size=size))
    

def box_muller_method(size = 1, random_state: np.random.Generator=rng) :
    """Implementation of the Box-Muller method to draw independent samples of standard normal distribution 
    (centered and variance is identity). 

    Args:
        size (tuple): shape of the output
        random_state : 'np.random.Generator' used for simulation. Defaults to rng.

    Returns:
        Array of shape 'size' of independent draws of standard normal distribution.
    """

    try :

        iteration_number = int(np.ceil(np.prod(size)/2)) #number of iteration of the method needed
        uniform_sample = random_state.random((iteration_number,2))
        normal_sample = []

        for i in range(iteration_number) :
            u1, u2 = uniform_sample[i,]
            exponential = (-2*np.log(u1))**(.5)
            uniform_angle = 2*np.pi*u2

            normal_sample.extend((exponential*np.cos(uniform_angle), exponential*np.sin(uniform_angle)))

        return (np.array(normal_sample[:np.prod(size)])).reshape(size)

    except TypeError :
        print("Enter an int")
        raise
    except ValueError :
        print("Enter a non negativ int")
        raise


def uniform_unit_ball(dimension: int=1, size = 1, random_state: np.random.Generator=rng) :
    """uniform_unit_ball : Return a list of length 'size' of independent realisations of
    uniform distribution over the unit ball of the real space to the power of 'dimension'.
    Generated using the acceptance-reject method.
    It's the unit ball for the L2 norm.

    Args:
        dimension (int, optional): dimension of the space of the ball. Defaults to 1.
        size (tuple, optional): 'shape of the output'. The function draws a sample of shape 'size'
            following the uniform unit ball distribution of dimension 'dimension'.
            Defaults to 1.
        random_state : 'np.random.Generator' used for simulation. Defaults to rng.
        
    Returns :
        Array of shape 'dimension, size' containing independent draws of uniform distribution over the unit ball of R**(dimension).
    """

    iteration_number = np.prod(size)
    U = np.empty((dimension, iteration_number))

    for i in range(iteration_number)  :
        u = 2*np.ones(dimension)

        while np.linalg.norm(u) > 1 :
            u = np.array(uniform_distribution(-1,1,dimension, random_state = random_state))
        
        U[:,i] = u

    if type(size) == int :
        return np.array(U).reshape((dimension, size))

    else : 
        return np.array(U).reshape((dimension,) + size)


def marsaglia_method(size = 1, random_state: np.random.Generator=rng) :
    """Implementation of the Marsaglia's polar method to draw independent draws of standard normal distribution 
    (centered and variance is identity). 

    Improvements to make : vectorize more, don't use a loop and generate a fixed number of uniform[0,1] chose so
    we accept at least 'size' points with high probability (cf : https://quantgirl.blog/comparing-box-muller-and-marsaglia-bray/)

    Args:
        size (tuple): shape of the output
        random_state : 'np.random.Generator' used for simulation. Defaults to rng.

    Returns:
        Array of shape 'size' of independent draws of standard normal distribution.
    """

    try :
        iteration_number = int(np.ceil(np.prod(size)/2)) #number of iteration of the method needed
        normal_sample = []

        for i in range(iteration_number) :
            r=2
            while r > 1 :
                x, y = 2*random_state.random(2) - 1
                r = x**2 + y**2
            
            temp = (-2*np.log(r)/r)**(.5)

            normal_sample.extend((temp*x, temp*y))

        return (np.array(normal_sample[:np.prod(size)])).reshape(size)

    except TypeError :
        print("Enter an int")
        raise
    except ValueError :
        print("Enter a non negativ int")
        raise


def gaussian_vector_distribution(size =1, mu = 0, sigma = 1, cholesky: bool=True, 
                                method = 'default',
                                random_state: np.random.Generator=rng) :
    """gaussian_distribution : Return a gaussian vector of mean 'mu' and covariance matrix
    'sigma'. 
    Uses Box-Muller ...

    Args:
        size (tuple) : 'shape of the output'. The function draws a sample of shape 'size'
            following the uniform unit ball distribution of dimension 'dimension'.
            Defaults to 1.
        mu (array_like): 1d array, mean vector of the gaussian vector
        sigma (array_like): 2d array of shape (n,n), where n is the shape of mu. Sigma
            is the matrix of the covariances between the variables within our gaussian vector.
            The function assumes that sigma really is a covariance matrix.
        cholesky (Boolean): if equals to 'True' then the algorithm will use the Cholesky Decomposition.
            Otherwise, it uses scipy.linalg.sqrtm. 
        method (string): method is either 'default', 'bm', 'marsaglia'. If method is 'default', the standard gaussians
            will be generated by the function of numpy. If methode is 'bm', it will use the function 'box_muller_method'
            and if its 'marsaglia' it will use 'marsaglia_method'. Defaults to 'default'
        random_state : 'np.random.Generator' used for simulation. Defaults to rng.

    Returns:
        1 array_like with same shape as mu. 
    """

    try : 
        mu, sigma = np.array(mu), np.array(sigma)
        vector_dimension = mu.size
        
        # Dimensions check
        if sigma.size == 1 and vector_dimension == 1 :
            if method == 'default' :
                standard_gaussian = random_state.standard_normal(size)
            elif method == 'bm' :
                standard_gaussian = box_muller_method(size, random_state=random_state)
            elif method == 'marsaglia' :
                standard_gaussian = marsaglia_method(size, random_state=random_state)
            else :
                print("Error in function 'gaussian_vector_distribution' of 'simulation' : method argument is not valid")
                return
            return np.sqrt(sigma)*standard_gaussian + mu*np.ones(size)
        else :
            if sigma.size == 1 and vector_dimension != 1 :
                sigma = sigma*np.eye( vector_dimension )
            else : 
                if sigma.size !=1 and vector_dimension == 1 :
                    mu = mu*np.ones( sigma.shape[0] )
                    vector_dimension = mu.size
                else :
                    if sigma.shape != (vector_dimension,vector_dimension) :
                        print("Error in 'gaussian_vector' of 'simulation' : shapes of parameters aren't compatible.")
                        return

        # Reshaping mu and drawing the gaussians
        if type(size) == int :
            mu = (np.repeat(mu, size, axis =0)).reshape((vector_dimension, size))
            if method == 'default' :
                standard_gaussian = random_state.standard_normal((vector_dimension, size))
            elif method == 'bm' :
                standard_gaussian = box_muller_method((vector_dimension, size), random_state=random_state)
            elif method == 'marsaglia' :
                standard_gaussian = marsaglia_method((vector_dimension, size), random_state=random_state)
            else :
                print("Error in function 'gaussian_vector_distribution' of 'simulation' : method argument is not valid")
                return
        else :
            mu = (np.repeat(mu, np.prod(size), axis=0)).reshape((vector_dimension,) + size)
            if method == 'default' :
                standard_gaussian = random_state.standard_normal((vector_dimension,) + size)
            elif method == 'bm' :
                standard_gaussian = box_muller_method((vector_dimension,) + size, random_state=random_state)
            elif method == 'marsaglia' :
                standard_gaussian = marsaglia_method((vector_dimension,) + size, random_state=random_state)
            else :
                print("Error in function 'gaussian_vector_distribution' of 'simulation' : method argument is not valid")
                return

        # Easy case : sigma is identity
        if np.allclose( sigma, np.eye(vector_dimension) ) :
            return standard_gaussian + mu

        # Computing square root of the covariance matrix
        if cholesky == True :
            covariance_matrix_sqrt = np.linalg.cholesky(sigma)
        else : 
            covariance_matrix_sqrt = (scipy.linalg.sqrtm(sigma))[0]
            

        return np.einsum('ij,j...->i...', covariance_matrix_sqrt, standard_gaussian) + mu

    except TypeError :
        print("Error in 'gaussian_vector' of 'simulation' : enter array-like parameters")
        raise


def standard_brownian_motion_1d_timeList(time_list, n_paths: int = 1, 
                                    increments : bool = False,
                                    random_state: np.random.Generator=rng) :
    """Standard brownian motion (1 dimension) generator

    Draw 'n_paths' standard brownian motion using the independant gaussians generated with numpy.
    Generated using independent increasing increments.

    Args:
        time_list (array_like): non decreasing sequence of non negativ float representing time. 
            Tis function doesn't check if 'non-decreasing' caracteristic so it won't return a brownian motion
            if time_list is not non decreasing.
        n_paths (int): number of paths to generate
        increments (bool) : Determine if we return the increments array or not.
        random_state : 'np.random.Generator' used for simulation. Defaults to rng.

    Returns:
        Array of float representing the values of the brownian motion at times in the time_list.
    """

    n_times = len(time_list) #number of gaussian needed
    gaussian_vector = random_state.normal(loc=0, scale=1, size=(n_times, n_paths))

    Diagonal = np.sqrt( time_list - (np.concatenate([0], time_list[:n-1]))  )
    increments_array = np.diag(D).dot(gaussian_vector)

    if increments :
        return increments_array

    return np.cumsum(increments_array, axis = 0)

def standard_brownian_motion_1d_timeParameters(n_times: int, n_paths : int,
                                            final_time: float = 1.0,
                                            increments : bool = False,
                                            random_state: np.random.Generator=rng) :
    """Standard brownian motion (1 dimension) generator

    Draw 'n_paths' standard brownian motion using the independant gaussians generated with numpy.
    Generated using independent increasing increments.

    Args:
        n_times (int): number of timesteps
        n_paths (int): number of paths simulated
        final_time (float, optional): Final time of simulation. Defaults to 1.0.
        random_state : 'np.random.Generator' used for simulation. Defaults to rng.

    Returns:
        `np.array` of shape `(n_times+1, n_paths)` containing the paths.
    """

    increments_list = np.sqrt(final_time / n_times) * random_state.standard_normal((n_times, n_paths))
    if increments:
        return increments_list
    else:
        brownian = np.zeros((n_times+1, n_paths))
        brownian[1:] = np.cumsum(increments_list, axis=0)
        return brownian






#####################################################
#####################################################
#####################################################
#####################################################
# time_list == None pose problème, suffirait pas de check le type (j'ai pas le temps là) ?

def refine_brownian_motion_1d(paths, time_list = None, 
                            constant_step : bool = True, 
                            random_state: np.random.Generator=rng) :
    """Refine a given brownian motion (1 dimension):

    For every 2 succesives times, add a point at the middle time on every path. 
    New points are determined using a formula (Lévy) that gives the conditionnal law of a Brownian knowing a point in the past 
    and a point in the future.

    Args:
        paths (2d array): First dimmension is the trajectories, second is the number of trajactories. The brownian motion to refine.
        time_list (1d array like, optional): list of the times associated with the paths. Defaults to None.
        constant_step (bool, optional): determines if the step of the time list is constant or not. Defaults to True.
        random_state : 'np.random.Generator' used for simulation. Defaults to rng.

    Returns:
        If time_list == None : 2d array representing the refined brownian motion.
        Else : 2d array representing the refined brownian motion and its associated new time_list.
    """
    n_times, n_paths = paths.shape
    if constant_step == True :
            # We use the expression of the brownian conditionally to his future and past
            temp = (paths.repeat(2, axis=0))
            non_random_part_paths = ( temp[1:,:] + temp[:-1,:] )/2 #creates an offset then calculates the mean
            
            step = 1/n_times
            random_part_paths = np.zeros(non_random_part_paths.shape)
            random_part_paths[1::2] = np.sqrt(step/4)*random_state.standard_normal((n_times - 1,n_paths)) #We add the random part to new points

            refined_paths = random_part_paths + non_random_part_paths

            if time_list == None : 
                return refined_paths

            else :
                return refined_paths, np.linspace(time_list[0], time_list[-1], step/2, endpoint=True)

    else :
        if time_list == None : 
            print("Error in 'refine_brownian_motion' in 'simulation' : arguments are missing")
            return
        else :
            temp = (paths.repeat(2, axis=0))
            non_random_part_paths = ( temp[1:,:] + temp[:-1,:] )/2

            random_part_paths = np.zeros(non_random_part_paths.shape)
            temp = (time_list.repeat(2, axis=0))
            sqrt_covariance_matrix = np.sqrt(( temp[1:] - temp[:-1] ))/2
            random_part_paths[1::2] = (np.diag(sqrt_covariance_matrix)).dot(random_state.standard_normal((n_times - 1,n_paths)))

            refined_time_list = ( (time_list.repeat(2, axis=0))[1:,:] + (time_list.repeat(2, axis=0))[:-1,:] )/2

            return non_random_part_paths + random_part_paths, refined_time_list

"""
B = standard_brownian_motion_timeParameters(50, 3, final_time = 2.0)
time_list = np.linspace(0, 2, 50, endpoint = True)
refined_B, refined_time_list = refine_brownian_motion(B, time_list)

plt.plot(B, time_list)
plt.show()
plt.plot(refined_B, refined_time_list)
plt.show()
"""

def fractionnal_brownian_motion() :
    return


def brownian_bridge() :
    return


def geometric_brownian_motion() :
    return


###
"""
N = 10**3
times = np.linspace(0,10,N)
Brownian = standard_brownian_simulation(times)
plt.plot(times, Brownian)
plt.show()
"""

# %%
