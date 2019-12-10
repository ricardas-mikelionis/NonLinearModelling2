import math

import Functions


def first_test(x, t, h, t_interval, a, c, d):
    u_j = Functions.u_function(x, t)
    u_j_abs_sqr = Functions.u_abs_sqr(x, t)
    u_j_p_1 = Functions.u_function(x + h, t)
    u_j_m_1 = Functions.u_function(x - h, t)
    u_r_j = Functions.u_function(x, t + t_interval)
    u_r_j_abs_sqr = Functions.u_abs_sqr(x, t + t_interval)
    u_r_j_p_1 = Functions.u_function(x + h, t + t_interval)
    u_r_j_m_1 = Functions.u_function(x - h, t + t_interval)
    f_j = Functions.f_function(x, t, a, c, d)
    f_r_j = Functions.f_function(x, t + t_interval, a, c, d)

    u_t_aproximation = (u_r_j - u_j)/t_interval
    u_nlm_approximation = 0.5 * ((u_r_j_p_1 - (2 * u_r_j) + u_r_j_m_1 + u_j_p_1 - (2 * u_j) + u_j_m_1)/math.pow(h, 2))
    u_x_aproximation = 0.5 * (u_r_j + u_j)
    u_x_abs_approximation = 0.5 * (u_r_j_abs_sqr * u_r_j + u_j_abs_sqr * u_j)
    f_approximation = 0.5 * (f_r_j + f_j)
    right_side = complex(pow(a, 2), 1) * u_nlm_approximation + complex(0, c) * u_x_aproximation \
                + complex(0, d) * u_x_abs_approximation + f_approximation

    residual = abs(u_t_aproximation - right_side)

    return residual


def second_test(x, t, h, t_interval, a, c, d):
    u_j = Functions.u_function(x, t)
    u_j_p_1 = Functions.u_function(x + h, t)
    u_j_m_1 = Functions.u_function(x - h, t)
    u_r_j = Functions.u_function(x, t + t_interval)
    u_r_j_p_1 = Functions.u_function(x + h, t + t_interval)
    u_r_j_m_1 = Functions.u_function(x - h, t + t_interval)
    u_j_abs = Functions.u_abs_sqr(x, t)
    u_r_j_abs = Functions.u_abs_sqr(x, t + t_interval)
    c_big = Functions.c_function(h, t_interval, a)
    f_small = Functions.f_function(x, t, a, c, d)
    f_r_small = Functions.f_function(x, t + t_interval, a, c, d)
    f_big = Functions.f_big_function(u_j, u_j_p_1, u_j_m_1, u_r_j,
                                     u_j_abs, u_r_j_abs, f_small, f_r_small,
                                     a, c, d, h, t_interval)
    residual = abs(u_r_j_p_1 - (c_big * u_r_j) + u_r_j_m_1 + f_big)

    return residual


def third_test(t, t_interval, n, a, c, d, delta):
    print 'Third Test'

    print "Pradines reiksmes"
    u_current = Functions.u_precise(t, n)
    print u_current

    print "Tikimasi reiksmiu:"
    u_new_percise = Functions.u_precise(t + t_interval, n)
    print u_new_percise

    print "Gauta:"
    u_new_calc = Functions.u_new_function(n, t, t_interval, a, c, d, Functions.u_precise(0, n), u_current, delta)
    print u_new_calc

    max_delta = 0.0

    for i in range (0, n+1):
        delta_u = abs(u_new_calc[i]) - abs(u_new_percise[i])
        if delta_u.real > max_delta:
            max_delta = delta_u.real

    print "Netiktis:"
    print max_delta


def fourth_test(t, n, t_interval, a, c, d, delta):
    print 'Fourth Test'

    print "Pradines reiksmes"
    u_current = Functions.u_precise(t, n)
    print u_current

    print "Tikimasi reiksmiu:"
    u_new_percise = Functions.u_precise(t + t_interval, n)
    print u_new_percise

    print "Gauta:"
    u_new_calc = Functions.u_new_function(n, t, t_interval, a, c, d, Functions.u_precise(0, n), u_current, delta)
    print u_new_calc

    max_delta = 0.0

    for i in range (0, n+1):
        delta_u = abs(u_new_calc[i]) - abs(u_new_percise[i])
        if delta_u.real > max_delta:
            max_delta = delta_u.real

    print "Netiktis:"
    print max_delta

    print"h, bei t zingsniai sumazinami 10 kartu"
    n = n*10
    t_interval = t_interval/10

    u_current = Functions.u_precise(t, n)
    u_new_percise = Functions.u_precise(t + t_interval, n)
    u_new_calc = Functions.u_new_function(n, t, t_interval, a, c, d, Functions.u_precise(0, n), u_current, delta)

    max_delta_2 = 0.0

    for i in range (0, n+1):
        delta_u = abs(u_new_calc[i]) - abs(u_new_percise[i])
        if delta_u.real > max_delta_2:
            max_delta_2 = delta_u.real

    print "Netiktis:"
    print max_delta_2

    print "Netiktis sumazejo ", max_delta/max_delta_2, " kartu"

def global_test():


