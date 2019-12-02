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
    right_side = complex(a^2, 1) * u_nlm_approximation + complex(0, c) * u_x_aproximation \
                + complex(0, d) * u_x_abs_approximation + f_approximation

    residual = abs(u_t_aproximation - right_side)

    return residual


def second_test(t_total, t_step, x, n, h, a, c, d, delta):
    u_j = Functions.u_function(x, 0)
    u_j_p_1 = Functions.u_function(x + h, 0)
    u_j_m_1 = Functions.u_function(x - h, 0)
    u_r_j = Functions.u_function(x, 0)
    u_r_j_p_1 = Functions.u_function(x + h, 0)
    u_r_j_m_1 = Functions.u_function(x - h, 0)

    t_in_loop = 0 + t_step

    while t_in_loop < t_total:
        #do stuff
        t_in_loop = t_in_loop + t_step

