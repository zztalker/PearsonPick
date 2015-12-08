# -*- coding: utf-8 -*-
import sys
import scipy.integrate as spint
import scipy.stats as spstat
import numpy as np
from sympy import integrate

import pearson
import matplotlib.pyplot as plt
from matplotlib import rc
import csv

rc('font', **{'family': 'verdana'})
rc('text.latex', unicode=True)
rc('text.latex', preamble='sepackage[utf8]{inputenc}')
rc('text.latex', preamble='sepackage[russian]{babel}')

print("введите порядковый номер файла (1-7)")
data_array = int(input())

# загрузка данных из файла
x, y = np.loadtxt("data{0}.csv".format(data_array), delimiter=",")
print("Test data")
p = y / y.sum()
p_cal = p.copy()

h = 1.0
# проверка на ошибки размерности
if round(p.sum(), 15) != 1.:
    print("Ошибка массива Y summ <> 1 ", p.sum())
    sys.exit(0)

if x.shape != y.shape:
    print("Ошибка размерности!")
    print("X", x.shape)
    print("Y", y.shape)
    sys.exit(0)
# конец проверок ошибок на размерность

# x - массив значений случайной величины, диапазон смещен
# y - массив частот попадания
# p - массив вероятностей
pears = pearson.pearson(x, p)
# x = _x[i] + self.c - self.Xa
# x_= (x[i] + self.c - self.Xa) / self.c
print("X", "N(p)", "X-Xa", "(X-Xa)/c", sep="\t")
for i in range(0, x.shape[0]):
    print(x[i], y[i], pears.x_[i], pears.x__[i], sep="\t")

print(pears)

# Calc function points
lB = np.min(pears.x__)
rB = np.max(pears.x__)
x1 = np.linspace(lB, rB, 1000)
y1 = np.ndarray(x1.shape)
# расчет y_i    ####косяяяяяяк
print('y1 = f.fun(x[i])\n')
i = 1
for i in x1:
    y1[i] = pears.f.fun(i)
    print(str(y1[i]) + ' ')
    i += 1

# Claculation for BAR - there is an ERROR
leftBorder = lB
ii = 0
for i in x:
    rightBorder = i - pears.Mx
    p_cal[ii] = spint.quad(pears.f.fun, leftBorder, rightBorder)[0]
    ii += 1
leftBorder = rightBorder
print(p_cal)
#
# # Original Datas
# plt.bar(pears.x__, y, color="white")  # width=pears.c,
# plt.ylabel("N")
# plt.xlabel(r"$x = (X - X_a)/c$")
# plt.title(r'Histogram of $x$')
# plt.savefig('data_bar_%d.png' % data_array)
# plt.clf()
#
# # Curve
# plt.plot(x1, y1)
# plt.ylabel("P")
# plt.xlabel("X")
# plt.title('Pearson-Curve TYPE: {0}'.format(pears.type))
# plt.savefig('data_curve_%d.png' % data_array)
# plt.clf()
#
# # Curve and Bar
# plt.bar(pears.x__, p, color="white")  # width=pears.c,
# plt.plot(x1, y1)
# plt.ylabel("P")
# plt.xlabel("X")
# plt.title('Compare Bar and Pearson-Curve TYPE: {0}'.format(pears.type))
# plt.savefig('data_curve_bar_%d.png' % data_array)
#
# plt.subplot(122)

# посчитаем квантиль распределения хи квадрат
chi2 = 0.0
ii = 0
for p0 in p_cal:
    chi2 += ((p0 - p[ii]) ** 2) / (p0)
    ii += 1
i_chi2 = 0
b_chi2 = False
for chi2critical in spstat.chi2.ppf([0.01, 0.05, 0.1], ii - 1):
    i_chi2 += 1
    if chi2 <= chi2critical:
        b_chi2 = True
        print(
            "Pearson criteria SUCCESSFUL for alpha=", 0.01 * i_chi2, " chi2_emp = ", chi2, "<chi2(", 0.01 * i_chi2,
            ",",
                                                      ii - 1, ")=", chi2critical)
        break

if not b_chi2:
    print("Pearson criteria UNSUCCESSFUL for alpha<0.1, gamma>0.9")

for i in range(0, x.shape[0]):
    print(x[i])
test = []
for i in pears.x_:
    test.append(str(i))

with open('x_.csv', 'w') as csvfile:
    csv.writer(csvfile).writerow(test)

test.clear()

for i in pears.x__:
    test.append(i)
with open('data/x__.csv', 'w') as csvfile:
    csv.writer(csvfile).writerow(test)

test.clear()

for i in p:
    test.append(str(i))
with open('data/p.csv', 'w') as csvfile:
    csv.writer(csvfile).writerow(test)

test.clear()

for i in y1:
    test.append(str(i))
with open('data/y.csv', 'w') as csvfile:
    csv.writer(csvfile).writerow(test)

test.clear()

print('\nTEST\n')
print("C=",pears.c,'\n')
for i in x1:
    print (i)
    test.append(str(spint.quad(pears.f.fun, i - 0.5, i + 0.5)))
with open('data/y_.csv', 'w') as csvfile:
    csv.writer(csvfile).writerow(test)
