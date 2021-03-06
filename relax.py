"""
Adaptation of relax.py (github.com/fraser-lab/relax) licensed via BSD 3-clause license.

Simplified to behave well in python 2 without signature or warnings regarding poor fits.
(probably put this back in if you are working with new data!)

Adds the ability to use bootstrap-based parameter error estimation (Efron B 1979 doi:10.1214/aos/1176344552)
"""
from scipy.optimize import curve_fit, leastsq, least_squares
import numpy as np



def single_step_relaxation(x,a,b,c):
    """
    Single step relaxation function based on 3 parameters a, b, c

    :param float x: Time variable (unit agnostic)
    :param float a: Total change in signal after relaxation
    :param float b: Observed rate of change relative to time variable
    :param float c: Initial signal offset prior to relaxation
    :return: Calculated signal at time x given relaxation parameters a, b, c
    """
    # A is the asymptote
    # B is the kobs
    # C is the offset
    return a*(1-np.exp(-b*x))+c

def expand_relax_single(x, p):
    """
    Utility function to pass a single argument tuple to single step relaxation, 
    for the purposes of bootstrapping with explicit residuals

    :param float x: Time variable
    :param array-like p: tu
    :return: Calculated signal at time x given parameters p
    """
    return single_step_relaxation(np.asarray(x), *p)


def two_step_relaxation(x,a,b,c,d,e):
    """
    Two-step relaxation function based on 5 parameters a, b, c, d, e
    Recommend using only if single-step relaxation is unable to converge effectively.

    :param float x: Time variable
    :param float a: Total change in signal from first relaxation step
    :param float b: Observed rate of first relaxation step
    :param float c: Total change in signal from second relaxation step
    :param float d: Observed rate of second relaxation step
    :param float e: Initial signal offset prior to relaxation
    :return: Calculated signal at time x given relaxation parameters a, b, c, d, e
    """
    return a*(1-np.exp(-b*x))+c*(1-np.exp(-d*x))+e


def expand_relax_two(x, p):
    """
    Utility function to pass a single argument tuple to two-step relaxation, 
    for the purposes of bootstrapping with explicit residuals

    :param float x: Time variable
    :param array-like p: 
    :return: Calculated signal at time x given parameters p
    """
    return two_step_relaxation(np.asarray(x),*p)


def three_step_relaxation(x,a,b,c,d,e,f,g):
    """
     Three-step relaxation function based on 7 parameters a, b, c, d, e, f, g
     Not advisable for any but very highly observed datasets due to high risk of overfitting

     :param float x: Time variable
     :param float a: Total change in signal from first relaxation step
     :param float b: Observed rate of first relaxation step
     :param float c: Total change in signal from second relaxation step
     :param float d: Observed rate of second relaxation step
     :param float e: Total change in signal from third relaxation step
     :param float f: Observed rate of change in second relaxation step
     :param float g: Initial signal offset prior to relaxation
     :return: Calculated signal at time x given relaxation parameters a, b, c, d, e, f, g
     """
    return a*(1-np.exp(-b*x))+c*(1-np.exp(-d*x))+e*(1-np.exp(-f*x))+g


def relaxation_fit(x, y, relaxation_function=single_step_relaxation, initial_guess=(1, 1, 1), maxfev=5000, sigma=None, 
    absolute_sigma=True):
    """
    Function to fit relaxation to observed signals y over times x using arbitrary functions.

    :param array-like x: times of observations
    :param array-like y: observed signal as a function of time
    :param function relaxation_function: relaxation function used to fit observed signal
    :param array-like initial_guess: initial guess for relaxation function parameters (default 1,1,1)
        Must match number of parameters in function
    :param int maxfev: Maximum cycles used for curve-fitting (default 5000)
    :param array-like sigma: array of standard errors
    :param bool absolute_sigma: absolute vs relative errors
    
    :return: tuple of calculated parameters and parameter covariance matrix and calculated y values for the parameters given data_x
    :rtype: tuple: array-like, matrix-like, array-like
    """

    # Signature.parameters gets the number of arguments in a function - that is 1 + the number of parameters.
    #assert len(initial_guess) == len(signature(relaxation_function).parameters)-1

    # Curve fit does the majority of the work.
    parameters, covariance = curve_fit(relaxation_function, x, y, p0=initial_guess, maxfev=maxfev, method='lm', sigma=sigma, 
        absolute_sigma = absolute_sigma)
    
    # Check that the parameters converged reasonably well and warn if they don't.
    for index, value in enumerate(parameters):
        standard_dev = np.sqrt(covariance[index,index])
        if np.abs(value) < np.abs(standard_dev):
            parameter_letter = "abcdefghijklmnopqrstuvwxyz"[index]
    #       warn(f"Parameter {parameter_letter} has standard deviation ({standard_dev}) larger than its value({value})")
    
    y_calc = [relaxation_function(i,*parameters) for i in x]
    return parameters, covariance, y_calc


def fit_bootstrap_two(data_y,data_x,relaxation_function=two_step_relaxation, initial_guess=(1,1,1,1), maxfev=5000, 
    sigma=None, absolute_sigma=True, bootstrap_count = 10000):
    """
    Function to fit relaxation to observed signals y over times x, using bootstrap parameter error estimation and two-step relaxation.

    :param array-like x: times of observations
    :param array-like y: observed signal as a function of time
    :param function relaxation_function: relaxation function used to fit observed signal
    :param array-like initial_guess: initial guess for relaxation function parameters (default 1,1,1,1)
        Must match number of parameters in function
    :param int maxfev: Maximum cycles used for curve-fitting (default 5000)
    :param array-like sigma: array of standard errors
    :param bool absolute_sigma: absolute vs relative errors
    :param int bootstrap_count: number of cycles to run bootstrap (default 10000)
    
    :return: tuple of calculated parameters and parameter covariance matrix and calculated y values for the parameters given data_x
    :rtype: tuple: array-like, matrix-like, array-like
    """
    
    # Errfunc returns residuals instead of y, which is how leastsq works.
    errfunc = lambda p,x,y: y-expand_relax_two(x,p)

    # First start with an initial error model.
    p_fit, p_err = leastsq(errfunc, x0=initial_guess, args=(data_x,data_y), full_output=0, maxfev=30000)
    residuals = errfunc(p_fit, data_x, data_y)
    sigma_res = np.std(residuals)
    sigma_err_total = sigma_res
    print
    ps = []
    
    # Bootstrap Loop!
    for i in range(bootstrap_count):
        random_delta = np.random.normal(0., sigma_err_total, len(data_y))
        rand_data_y = data_y + random_delta
        # if i < 5:
        #     print (rand_data_y)
        result = least_squares(errfunc, x0=p_fit, bounds=((-10, 10**-8, -10, 10**-8, -np.inf),(10,1,10,1,np.inf)), 
            args=(data_x,rand_data_y), method='dogbox',   max_nfev=200)
        rand_fit=result['x']
        ps.append(rand_fit)

    Determine errors based on results of bootstrap.
    ps = np.array(ps)
    mean_pfit = np.median(ps, 0)
    # print (mean_pfit)

    NSigma = 1
    err_pfit = NSigma * np.std(ps,0)

    y_calc =  expand_relax_two(data_x, mean_pfit)


    return mean_pfit, err_pfit, y_calc


def fit_bootstrap_single(data_y,data_x,relaxation_function=single_step_relaxation, initial_guess=(1,1,1), maxfev=5000, 
    sigma=None, absolute_sigma=True, bootstrap_count=10000):
    """
    Function to fit relaxation to observed signals y over times x, using bootstrap parameter error estimation and two-step relaxation.

    :param array-like x: times of observations
    :param array-like y: observed signal as a function of time
    :param function relaxation_function: relaxation function used to fit observed signal
    :param array-like initial_guess: initial guess for relaxation function parameters (default 1,1,1,1)
        Must match number of parameters in function
    :param int maxfev: Maximum cycles used for curve-fitting (default 5000)
    :param array-like sigma: array of standard errors
    :param bool absolute_sigma: absolute vs relative errors
    :param int bootstrap_count: number of cycles to run bootstrap (default 10000)
    
    :return: tuple of calculated parameters and parameter covariance matrix and calculated y values for the parameters given data_x
    :rtype: tuple: array-like, matrix-like, array-like
    """

    # Errfunc returns residuals instead of y, which is how leastsq works.
    errfunc = lambda p,x,y: y-expand_relax_single(x,p)

    # First start with an initial error model.
    p_fit, p_err = leastsq(errfunc, x0=initial_guess, args=(data_x,data_y), full_output=0, maxfev=30000)
    residuals = errfunc(p_fit, data_x, data_y)
    sigma_res = np.std(residuals)
    sigma_err_total = sigma_res
    print
    ps = []

    # Bootstrap Loop! 
    for i in range(bootstrap_count):
        random_delta = np.random.normal(0., sigma_err_total, len(data_y))
        rand_data_y = data_y + random_delta
        # if i < 5:
        #     print (rand_data_y)
        result = least_squares(errfunc, x0=p_fit, args=(data_x,rand_data_y), method='dogbox',   max_nfev=200)
        rand_fit=result['x']
        ps.append(rand_fit)

    # Determine errors based on results of bootstrap.
    ps = np.array(ps)
    mean_pfit = np.median(ps, 0)
    # print (mean_pfit)

    NSigma = 1
    err_pfit = NSigma * np.std(ps,0)

    y_calc =  expand_relax_single(data_x, mean_pfit)


    return mean_pfit, err_pfit, y_calc






    
