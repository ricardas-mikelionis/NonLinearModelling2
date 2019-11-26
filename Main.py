# Ricardas Mikelionis Kompiuterinis Modeliavimas 2019
# Uzduotis 23
# Lygtis F, algoritmas F2, kraštinės sąlygos I

import cmath as cm
import math


def f_function(x_f, t_f):
    f = 0
    return f


def u_function(x_f, t_f):
    com = complex(1, t_f)
    pi_x = math.pi * x_f
    sin = math.pow(math.sin(pi_x), 2)
    u = com * sin
    return u


# Apsirašomos konstantos:
a = 0
c = 0
d = 0
k = 0
alpha = 1.5
beta = 0
gamma = 0
theta = 0
delta = math.pow(10, -8)  # vartotojo pasirenkamas tikslumas

n = 100  # somebignumber
h = 1 / n  # interval, x-step (1/N)
t_interval = 0.01  # t step
time = 2  # time in seconds
x_min = 0  # 0 < x < 1

t = 0.0  # initial
count = 0  # initial
u_matrix = [[0 for x in range(n)] for y in range(int(time // t_interval)+1)]

while t < time:
    x = x_min
    for j in range(0, n):
        u_matrix[count][j] = u_function(x, t)
        x = x + h
    count = count + 1
    t = t + t_interval

"""
x=x_min
for j in range(0, n):
    u_matrix[0][j] = u_function(x, t)
    x = x + h
    print(x, u_matrix[0][j])
"""