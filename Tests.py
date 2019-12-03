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


def third_test(n, h, t_interval, a, c, d, delta):
    alpha_array = [0 for i in range(n + 1)]
    beta_array = [0 for i in range(n + 1)]
    t = 0
    cont = True

    u = [0 for i in range(n + 1)]
    u_new = [0 for i in range(n + 1)]
    u_alph_sqr = [0 for i in range(n + 1)]
    f_small = [0 for i in range(n + 1)]
    new_old_deltas = [0 for i in range(n + 1)]
    c_big = Functions.c_function(h, t_interval, a)

    for j in range(0, n + 1):
        u[j] = Functions.u_function(h * j, 0)
        u_alph_sqr[j] = Functions.u_abs_sqr(h * j, 0)
        f_small[j] = Functions.f_function(h * j, 0, a, c, d)

    u_r = u.copy()

    for i in range(1, n):
        alpha_array[i] = 1 / (c_big * alpha_array[i - 1])

    while cont:
        for i in range(1, n):
            u_r_abs_sqr = Functions.u_abs_sqr(i * h, t)
            f_r_small = Functions.f_function(h*i, t, a, c, d)
            f_big = Functions.f_big_function(u[i], u[i+1], u[i-1], u_r[i], u_alph_sqr[i], u_r_abs_sqr,
                                             f_small[i], f_r_small, a, c, d, h, t_interval)
            beta_array[i] = (f_big + beta_array[i - 1]) / (c_big * alpha_array[i - 1])

        for k in range(1, n):
            u_new[n+1-k] = (alpha_array[n-k] * u_new[n-k]) + beta_array[n-k]
            new_old_deltas[n+1-k] = u_new[n + 1 - k] - u_r[n + 1 - k]

        if max(new_old_deltas) < delta:
            cont = False
        else:
            t = t + t_interval

    return u_new[5]
