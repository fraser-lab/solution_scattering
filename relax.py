#from inspect import signature
#from warnings import warn

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


def relaxation_fit(x, y, relaxation_function=single_step_relaxation, initial_guess=(1, 1, 1), maxfev=5000, sigma=None, absolute_sigma=True):
    """
    Function to fit relaxation to observed signals y over times x

    :param x: times of observations
    :param y: observed signal as a function of time
    :param relaxation_function: relaxation function used to fit observed signal
    :param initial_guess: initial guess for relaxation function parameters
        Must match number of parameters in function
    :param maxfev: Maximum cycles used for curve-fitting (default 5000)
    :type x: array-like
    :type y: array-like
    :type relaxation_function: function
    :type initial_guess: array-like
    :type maxfev: int
    :return: tuple of calculated parameters and parameter covariance matrix
    :rtype: tuple: array-like, matrix-like
    """

    # Signature.parameters gets the number of arguments in a function - that is 1 + the number of parameters.
    #assert len(initial_guess) == len(signature(relaxation_function).parameters)-1
    parameters, covariance = curve_fit(relaxation_function, x, y, p0=initial_guess, maxfev=maxfev, method='lm', sigma=sigma, absolute_sigma = absolute_sigma)
    for index, value in enumerate(parameters):
        standard_dev = np.sqrt(covariance[index,index])
        if np.abs(value) < np.abs(standard_dev):
            parameter_letter = "abcdefghijklmnopqrstuvwxyz"[index]
    #       warn(f"Parameter {parameter_letter} has standard deviation ({standard_dev}) larger than its value({value})")
    y_calc = [relaxation_function(i,*parameters) for i in x]
    return parameters, covariance, y_calc


def fit_bootstrap_two(data_y,data_x,relaxation_function=two_step_relaxation, initial_guess=(1,1,1), maxfev=5000, sigma=None, absolute_sigma=True):
    errfunc = lambda p,x,y: y-expand_relax_two(x,p)
    p_fit, p_err = leastsq(errfunc, x0=initial_guess, args=(data_x,data_y), full_output=0, maxfev=30000)
    residuals = errfunc(p_fit, data_x, data_y)
    sigma_res = np.std(residuals)
    sigma_err_total = sigma_res
    print
    ps = []
    for i in range(10000):
        random_delta = np.random.normal(0., sigma_err_total, len(data_y))
        rand_data_y = data_y + random_delta
        # if i < 5:
        #     print (rand_data_y)
        result = least_squares(errfunc, x0=p_fit, bounds=((-10, 10**-8, -10, 10**-8, -np.inf),(10,1,10,1,np.inf)), args=(data_x,rand_data_y), method='dogbox',   max_nfev=200)
        rand_fit=result['x']
        ps.append(rand_fit)


    ps = np.array(ps)
    mean_pfit = np.median(ps, 0)
    # print (mean_pfit)

    NSigma = 1
    err_pfit = NSigma * np.std(ps,0)

    y_calc =  expand_relax_two(data_x, mean_pfit)


    return mean_pfit, err_pfit, y_calc


def fit_bootstrap_single(data_y,data_x,relaxation_function=single_step_relaxation, initial_guess=(1,1,1), maxfev=5000, sigma=None, absolute_sigma=True):
    errfunc = lambda p,x,y: y-expand_relax_single(x,p)
    p_fit, p_err = leastsq(errfunc, x0=initial_guess, args=(data_x,data_y), full_output=0, maxfev=30000)
    residuals = errfunc(p_fit, data_x, data_y)
    sigma_res = np.std(residuals)
    sigma_err_total = sigma_res
    print
    ps = []
    for i in range(10000):
        random_delta = np.random.normal(0., sigma_err_total, len(data_y))
        rand_data_y = data_y + random_delta
        # if i < 5:
        #     print (rand_data_y)
        result = least_squares(errfunc, x0=p_fit, args=(data_x,rand_data_y), method='dogbox',   max_nfev=200)
        rand_fit=result['x']
        ps.append(rand_fit)


    ps = np.array(ps)
    mean_pfit = np.median(ps, 0)
    # print (mean_pfit)

    NSigma = 1
    err_pfit = NSigma * np.std(ps,0)

    y_calc =  expand_relax_single(data_x, mean_pfit)


    return mean_pfit, err_pfit, y_calc






    
