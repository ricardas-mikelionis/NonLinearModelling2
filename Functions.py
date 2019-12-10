import math


def f_function(x_f, t_f, a, c, d):
    com_a = complex(pow(a, 2), 1)
    com_c = complex(0, c)
    com_d = complex(0, d)
    f = u_deriv_t(x_f) - (com_a * u_sqr_deriv_x(x_f, t_f)) \
        - (com_c * u_function(x_f, t_f)) \
        - (com_d * u_abs_sqr(x_f, t_f) * u_function(x_f, t_f))
    return f


def u_function(x_f, t_f):
    com = complex(1, t_f)
    pi_x = math.pi * x_f
    cos = math.cos(math.pi/2 + pi_x)
    u = com * cos
    return u


def u_deriv_t(x_f):
    com = complex(0, 1)
    pi_x = math.pi * x_f
    sin = math.sin(pi_x)
    u_derived_t = -1 * com * sin
    return u_derived_t


def u_deriv_x(x_f, t_f):
    com = complex(1, t_f)
    pi_x = math.pi * x_f
    cos = math.cos(pi_x)
    u_derived_x = -1 * math.pi * com * cos
    return u_derived_x


def u_sqr_deriv_x(x_f, t_f):
    com = complex(1, t_f)
    pi_x = math.pi * x_f
    sin = math.sin(pi_x)
    u_sqr_derived_x = math.pow(math.pi, 2) * com * sin
    return u_sqr_derived_x


def u_abs_sqr(x_f, t_f):
    return (math.pow(t_f, 2) + 1) * (math.pow(math.sin(math.pi * x_f), 2))


def c_function(h, t_interval, a):
    com = complex(math.pow(a, 2), 1)
    c = 2 + ((2 * math.pow(h, 2))/(t_interval * com))
    return c


def f_big_function(u_j, u_j_p_1, u_j_m_1, u_r_j, u_j_abs_sqr, u_r_j_abs_sqr, f, f_r, a, c, d, h, t_interval):
    com_d = complex(0, d)
    com_a = complex(math.pow(a, 2), 1)
    com_c = complex(0, c)

    f_part = f + f_r
    id_part = u_r_j_abs_sqr * u_r_j + u_j_abs_sqr * u_j
    ic_part = u_r_j + u_j
    h_a = math.pow(h, 2) / com_a
    two_tau = 2 / t_interval

    f_big = u_j_p_1 - (2 * u_j) + u_j_m_1 + (h_a * ((two_tau * u_j) + (com_c * ic_part) + (com_d * id_part) + f_part))

    return f_big


def u_precise(t, n):
    u = [0 for i in range(n + 1)]
    h = 1.0 / n

    for j in range(0, n + 1):
        u[j] = u_function(h * j, t)

    return u


def alpha_array(n, t_interval, a):
    h = 1 / n
    c_big = c_function(h, t_interval, a)
    alpha_array = [0 for i in range(n + 1)]

    for i in range(1, n):
        alpha_array[i] = 1 / (c_big - alpha_array[i - 1])
    return alpha_array


def u_new_function(n, t_current, t_interval, a, c, d, u, u_r, delta):
    h = 1 / n
    beta_array = [0 for i in range(n + 1)]
    u_new = [complex(0, 0) for i in range(n+1)]
    c_big = c_function(h, t_interval, a)
    aph_array = alpha_array(n, t_interval, a)
    new_old_deltas = [0 for i in range(n + 1)]
    cont = True

    while cont:
        for i in range(1, n+1):
            u_r_abs_sqr = math.pow(u_r[i-1].real, 2) * math.pow(u_r[i-1].imag, 2)
            u_abs = u_abs_sqr(h * i-1, 0)
            f_r_small = f_function(h*i-1, t_current, a, c, d)
            f_big = f_big_function(u[i-1], u[i], u[i-2], u_r[i-1], u_abs, u_r_abs_sqr,
                                   f_function(h * i-1, 0, a, c, d), f_r_small, a, c, d, h, t_interval)
            beta_array[i] = (f_big + beta_array[i - 1]) / (c_big - aph_array[i - 1])

        for k in reversed(range(1, n)):
            u_new[k] = (aph_array[k+1] * u_new[k+1]) + beta_array[k+1]
            delta_new_old = u_new[k] - u_r[k]
            new_old_deltas[k] = delta_new_old.real

        if max(new_old_deltas) < delta:
            cont = False
        else:
            u_r = list(u_new)

    return u_new


def full_algorithm(t_interval, t_total, n, a, c, d, delta):
    t = 0
    u = u_precise(t, n)
    while t < t_total:
        u_r = list(u)
        u_new = u_new_function(n, t, t_interval, a, c, d, u, u_r, delta)
        u = list(u_new)
        t = t + t_interval

    return u
