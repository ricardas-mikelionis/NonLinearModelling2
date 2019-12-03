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
