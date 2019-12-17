from random import random
import Constants
# from numpy.random import random
import Functions


def first_test(x, t):
    print 'First Test'
    print "Pradines reiksmes: t_interval = ", Constants.T_INTERVAL, ", h = ", Constants.H
    residual = Functions.first_test(x, t)
    Constants.N = Constants.N * 10
    Constants.H = Constants.H / 10
    Constants.T_INTERVAL = Constants.T_INTERVAL / 10
    residual_2 = Functions.first_test(x, t)
    print "Zingsni sumazinus 10 kart: t_interval = ", Constants.T_INTERVAL, ", h = ", Constants.H
    print "Netiktis sumazeja: ", residual / residual_2

    # Reiksmes atstatomos:
    Constants.N = Constants.N / 10
    Constants.H = Constants.H * 10
    Constants.T_INTERVAL = Constants.T_INTERVAL * 10


def second_test(x, t):
    print '\nSecond Test'
    print "Pradines reiksmes: t_interval = ", Constants.T_INTERVAL, ", h = ", Constants.H
    residual = Functions.second_test(x, t)
    Constants.N = Constants.N * 10
    Constants.H = Constants.H / 10
    Constants.T_INTERVAL = Constants.T_INTERVAL / 10
    residual_2 = Functions.second_test(x, t)
    print "Zingsni sumazinus 10 kart: t_interval = ", Constants.T_INTERVAL, ", h = ", Constants.H
    print "Netiktis sumazeja: ", residual / residual_2

    # Reiksmes atstatomos:
    Constants.N = Constants.N /10
    Constants.H = Constants.H * 10
    Constants.T_INTERVAL = Constants.T_INTERVAL * 10


def third_test():
    print '\nThird test'

    y = [0j]  # Del I K.S. u_0 = 0
    f_big = [0 for i in range(Constants.N + 1)]
    c_big = Functions.c_function()

    for i in range(1, Constants.N):
        y.append(complex(random(), random()))

    y.append(0j)  # Del I K.S. u_n = 0

    for j in range(1, Constants.N):
        f_big[j] = c_big * y[j] - y[j - 1] - y[j + 1]

    y_new = Functions.thomas_algorithm(f_big)

    max_delta = Functions.max_difference(y_new, y)

    print "Sugeneruoto atsitiktiniu reiksmiu masyvo Y paskaiciuotos naujos " \
          "reiksmes naudojant Thomas Algoritma su paklaida:", max_delta


def global_test(t):
    print '\nGlobal Test'

    print "Pradines reiksmes:"
    print "tau = ", Constants.T_INTERVAL, ", h = ", Constants.H

    print "Paklaida:"
    max_delta = Functions.full_algorithm(t)
    print max_delta

    print"Zingsniai h, bei t sumazinami 10 kartu:"
    Constants.N = Constants.N * 10
    Constants.T_INTERVAL = Constants.T_INTERVAL / 10
    Constants.H = Constants.H / 10
    print "tau = ", Constants.T_INTERVAL, ", h = ", Constants.H

    print "Paklaida:"
    max_delta_2 = Functions.full_algorithm(t)
    print max_delta_2

    print "Paklaida sumazeja ", max_delta/max_delta_2, " kartu"




