import math
import Constants


def f_function(x_f, t_f):
    com_a = complex(pow(Constants.A, 2), 1)
    com_c = complex(0, Constants.C)
    com_d = complex(0, Constants.D)
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


def c_function():
    com = complex(math.pow(Constants.A, 2), 1)
    bott = Constants.T_INTERVAL * com
    top = 2 * math.pow(Constants.H, 2)
    c = 2 + (top/bott)

    return c


def f_big_function(u_j, u_j_p_1, u_j_m_1, u_r_j, u_j_abs_sqr, u_r_j_abs_sqr, f, f_r):
    com_d = complex(0, Constants.D)
    com_a = complex(math.pow(Constants.A, 2), 1)
    com_c = complex(0, Constants.C)

    f_part = f + f_r
    id_part = u_r_j_abs_sqr * u_r_j + u_j_abs_sqr * u_j
    ic_part = u_r_j + u_j
    h_a = math.pow(Constants.H, 2) / com_a
    two_tau = 2 / Constants.T_INTERVAL

    f_big = u_j_p_1 - (2 * u_j) + u_j_m_1 + (h_a * ((two_tau * u_j) + (com_c * ic_part) + (com_d * id_part) + f_part))

    return f_big


def u_precise(t):
    u = []

    for j in range(0, Constants.N + 1):
        u.append(u_function(Constants.H * j, t))

    return u


def alpha_array():
    c_big = c_function()
    aph_array = [0, 0] #Alpha_0 neegzistuoja, o Alpha_1 = 0, del I K.S. palikta del kodo skaitomumo

    for i in range(2, Constants.N):
        aph_array.append(1 / (c_big - aph_array[i - 1]))

    aph_array.append(0)  # Alpha_N = 0 del I K.S
    return aph_array


def f_big_list(t_current, u, u_r):
    f_big_array = [0j]  # Pirmas narys bus nenaudojamas(palikta del kodo skaitomumo

    for i in range(1, Constants.N):
        u_r_abs_sqr = math.pow(u_r[i].real, 2) + math.pow(u_r[i].imag, 2)
        u_abs = math.pow(u[i].real, 2) + math.pow(u[i].imag, 2)
        f_small = f_function(Constants.H * i, t_current)
        f_r_small = f_function(Constants.H * i, t_current + Constants.T_INTERVAL)
        f_big_array.append(f_big_function(u[i], u[i + 1], u[i - 1], u_r[i], u_abs, u_r_abs_sqr, f_small, f_r_small))

    return f_big_array


def u_new_function(u_r, f_big):
    u_new = []
    cont = True

    while cont:
        u_new = thomas_algorithm(f_big)
        max_diff = max_difference(u_new, u_r)

        if max_diff < Constants.DELTA:
            cont = False
        else:
            u_r = list(u_new)

    return u_new


def thomas_algorithm(f_big):
    aph_array = alpha_array()
    c_big = c_function()
    beta_array = [0, 0]  # beta_0 - neegzistuoja, o beta_1 = 0 del I K.S.
    u_new = [complex(0, 0) for i in range(Constants.N + 1)]

    for i in range(2, Constants.N + 1):
        top = f_big[i - 1] + beta_array[i - 1]
        bottom = c_big - aph_array[i - 1]
        beta_array.append(top / bottom)

    for k in reversed(range(1, Constants.N)):
        u_new[k] = (aph_array[k + 1] * u_new[k + 1]) + beta_array[k + 1]

    return u_new


def full_algorithm(t_total):
    t = 0
    u = u_precise(t)
    max_difference_array = []
    u_r = list(u)

    while t < t_total:
        f_big = f_big_list(t, u, u_r)
        u_new = u_new_function(u_r, f_big)
        t = t + Constants.T_INTERVAL
        u = u_precise(t)
        max_difference_array.append(max_difference(u_new, u))
        u_r = list(u_new)

    return max(max_difference_array)


def max_difference(u_new, u_old):
    max_diff = 0

    for i in range(0, len(u_new)):
        difference = abs(u_old[i] - u_new[i])
        if max_diff < difference:
            max_diff = difference

    return max_diff


def first_test(x, t):
    u_j = u_function(x, t)
    u_j_abs_sqr = u_abs_sqr(x, t)
    u_j_p_1 = u_function(x + Constants.H, t)
    u_j_m_1 = u_function(x - Constants.H, t)
    u_r_j = u_function(x, t + Constants.T_INTERVAL)
    u_r_j_abs_sqr = u_abs_sqr(x, t + Constants.T_INTERVAL)
    u_r_j_p_1 = u_function(x + Constants.H, t + Constants.T_INTERVAL)
    u_r_j_m_1 = u_function(x - Constants.H, t + Constants.T_INTERVAL)
    f_j = f_function(x, t)
    f_r_j = f_function(x, t + Constants.T_INTERVAL)

    u_t_approximation = (u_r_j - u_j)/Constants.T_INTERVAL
    u_nlm_approximation = 0.5 * ((u_r_j_p_1 - (2 * u_r_j) + u_r_j_m_1 + u_j_p_1 - (2 * u_j)
                                  + u_j_m_1) / math.pow(Constants.H, 2))
    u_x_approximation = 0.5 * (u_r_j + u_j)
    u_x_abs_approximation = 0.5 * (u_r_j_abs_sqr * u_r_j + u_j_abs_sqr * u_j)
    f_approximation = 0.5 * (f_r_j + f_j)
    right_side = complex(math.pow(Constants.A, 2), 1) * u_nlm_approximation\
                 + complex(0, Constants.C) * u_x_approximation + complex(0, Constants.D)\
                 * u_x_abs_approximation + f_approximation

    residual = abs(u_t_approximation - right_side)

    return residual


def second_test(x, t):
    u_j = u_function(x, t)
    u_j_p_1 = u_function(x + Constants.H, t)
    u_j_m_1 = u_function(x - Constants.H, t)

    u_r_j = u_function(x, t + Constants.T_INTERVAL)
    u_r_j_p_1 = u_function(x + Constants.H, t + Constants.T_INTERVAL)
    u_r_j_m_1 = u_function(x - Constants.H, t + Constants.T_INTERVAL)

    u_j_abs = u_abs_sqr(x, t)
    u_r_j_abs = u_abs_sqr(x, t + Constants.T_INTERVAL)

    c_big = c_function()
    f_small = f_function(x, t)
    f_r_small = f_function(x, t + Constants.T_INTERVAL)
    f_big = f_big_function(u_j, u_j_p_1, u_j_m_1, u_r_j, u_j_abs, u_r_j_abs, f_small, f_r_small)
    residual = abs(u_r_j_m_1 - (c_big * u_r_j) + u_r_j_p_1 + f_big)

    return residual
