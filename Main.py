# Ricardas Mikelionis Kompiuterinis Modeliavimas 2019
# Uzduotis 23
# Lygtis F, algoritmas F2, krastines salygos I
import math
import Tests
import Functions


# Apsirasomos konstantos:
a = -0.7
c = 11.5
d = 8.3
delta = math.pow(10, -8)  # vartotojo pasirenkamas tikslumas

n = 10  # somebignumber
h = 1.0 / n  # interval, x-step (1/N)
t_interval = 0.1  # t step
time = 2  # total time in seconds

x = 0.39  # 0 < x < 1
t = 1.17 # initial place in time
'''
print 'First Test'
print Tests.first_test(x, t, h, t_interval, a, c, d) / Tests.first_test(x, t, h/10, t_interval/10, a, c, d)
print Tests.first_test(x, t, h/10, t_interval/10, a, c, d) / Tests.first_test(x, t, h/100, t_interval/100, a, c, d)

print 'Second Test'
print Tests.second_test(x, t, h, t_interval, a, c, d) / Tests.second_test(x, t, h/10, t_interval/10, a, c, d)
print Tests.second_test(x, t, h/10, t_interval/10, a, c, d) / Tests.second_test(x, t, h/100, t_interval/100, a, c, d)
'''
Tests.third_test(n, t_interval, a)

#Tests.global_test(t, n, t_interval, a, c, d, delta)

