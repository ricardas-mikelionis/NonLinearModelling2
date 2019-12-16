import math

import Constants

from numpy.random import random

import Functions


def first_test(x, t):
    u_j = Functions.u_function(x, t)
    u_j_abs_sqr = Functions.u_abs_sqr(x, t)
    u_j_p_1 = Functions.u_function(x + Constants.H, t)
    u_j_m_1 = Functions.u_function(x - Constants.H, t)
    u_r_j = Functions.u_function(x, t + Constants.T_INTERVAL)
    u_r_j_abs_sqr = Functions.u_abs_sqr(x, t + Constants.T_INTERVAL)
    u_r_j_p_1 = Functions.u_function(x + Constants.H, t + Constants.T_INTERVAL)
    u_r_j_m_1 = Functions.u_function(x - Constants.H, t + Constants.T_INTERVAL)
    f_j = Functions.f_function(x, t, Constants.A, Constants.C, Constants.D)
    f_r_j = Functions.f_function(x, t + Constants.T_INTERVAL, Constants.A, Constants.C, Constants.D)

    u_t_aproximation = (u_r_j - u_j)/Constants.T_INTERVAL
    u_nlm_approximation = 0.5 * ((u_r_j_p_1 - (2 * u_r_j) + u_r_j_m_1 + u_j_p_1 - (2 * u_j)
                                  + u_j_m_1) / math.pow(Constants.H, 2))
    u_x_aproximation = 0.5 * (u_r_j + u_j)
    u_x_abs_approximation = 0.5 * (u_r_j_abs_sqr * u_r_j + u_j_abs_sqr * u_j)
    f_approximation = 0.5 * (f_r_j + f_j)
    right_side = complex(math.pow(Constants.A, 2), 1) * u_nlm_approximation + complex(0, Constants.C)\
                 * u_x_aproximation + complex(0, Constants.D) * u_x_abs_approximation + f_approximation

    residual = abs(u_t_aproximation - right_side)

    return residual


def second_test(x, t):
    u_j = Functions.u_function(x, t)
    u_j_p_1 = Functions.u_function(x + Constants.H, t)
    u_j_m_1 = Functions.u_function(x - Constants.H, t)

    u_r_j = Functions.u_function(x, t + Constants.T_INTERVAL)
    u_r_j_p_1 = Functions.u_function(x + Constants.H, t + Constants.T_INTERVAL)
    u_r_j_m_1 = Functions.u_function(x - Constants.H, t + Constants.T_INTERVAL)

    u_j_abs = Functions.u_abs_sqr(x, t)
    u_r_j_abs = Functions.u_abs_sqr(x, t + Constants.T_INTERVAL)

    c_big = Functions.c_function(Constants.H, Constants.T_INTERVAL, Constants.A)
    f_small = Functions.f_function(x, t, Constants.A, Constants.C, Constants.D)
    f_r_small = Functions.f_function(x, t + Constants.T_INTERVAL, Constants.A, Constants.C, Constants.D)
    f_big = Functions.f_big_function(u_j, u_j_p_1, u_j_m_1, u_r_j, u_j_abs, u_r_j_abs, f_small, f_r_small,
                                     Constants.A, Constants.C, Constants.D, Constants.H, Constants.T_INTERVAL)
    residual = abs(u_r_j_m_1 - (c_big * u_r_j) + u_r_j_p_1 + f_big)

    return residual


def third_test(n, t_interval, a):
    print '\nThird test'

    h = 1.0 / n
    y = [complex(0, 0) for i in range(n + 1)]
    f_big = [0 for i in range(n+1)]
    c_big = Functions.c_function(h, t_interval, a)

    for i in range(1, n):
        y[i] = complex(random(), random())

    for j in range(1, n):
        f_big[j] = c_big * y[j] - y[j - 1] - y[j + 1]

    y_new = Functions.thomas_algorithm(n, t_interval, a, f_big)

    max_delta = Functions.max_difference(y_new, y)

    print max_delta


def global_test(t, n, t_interval, a, c, d, delta):
    print '\nGlobal Test'

    print "Pradines reiksmes"
    print "tau = ", t_interval, ", h = ", (1.0 / n)

    print "Paklaida:"
    max_delta = Functions.full_algorithm(t_interval, t, n, a, c, d, delta)
    print max_delta

    print"h, bei t zingsniai sumazinami 10 kartu"
    n = n * 10
    t_interval = t_interval/10
    print "tau = ", t_interval, ", h = ", (1.0 / n)

    print "Paklaida:"
    max_delta_2 = Functions.full_algorithm(t_interval, t, n, a, c, d, delta)
    print max_delta_2

    print "Paklaida ", max_delta/max_delta_2, " kartu"




