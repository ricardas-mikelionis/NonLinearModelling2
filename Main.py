# Ricardas Mikelionis Kompiuterinis Modeliavimas 2019
# Uzduotis 23
# Lygtis F, algoritmas F2, krastines salygos I
import math
import Tests
import Functions


x = 0.39  # 0 < x < 1
t = 1.17  # place in time

print 'First Test'
print Tests.first_test(x, t) / Tests.first_test(x, t) #reikia pamazinti global konstanta h ir t interval
print Tests.first_test(x, t) / Tests.first_test(x, t)

print '\nSecond Test'
print Tests.second_test(x, t) / Tests.second_test(x, t)#reikia pamazinti global konstanta h ir t interval
print Tests.second_test(x, t) / Tests.second_test(x, t)

'''
Tests.third_test(n, t_interval, a)

Tests.global_test(t, n, t_interval, a, c, d, delta)
'''
