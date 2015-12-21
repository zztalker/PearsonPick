# -*- coding: utf-8 -*-
import sys
import scipy.integrate as spint
import scipy.stats as spstat
import numpy as np

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
# x_ = _x[i] + self.c - self.Xa                  - смещенная относительно центра сл.вел
# x__= (x[i] + self.c - self.Xa) / self.c        - смещенная и нормированная сл.вел

with open('data/data_calc.csv', 'w') as csvfile:
    for i in range(0, x.shape[0]):
        test = [x[i],y[i],p[i],pears.x_[i],pears.x__[i]]
        csv.writer(csvfile).writerow(test)

with open('data/report.txt','w') as reportfile:
    reportfile.write(pears.__str__())

# Calc function points
lB = np.min(pears.x__)-2 # левая граница, как миниму смещенного и центрированного массива
rB = np.max(pears.x__) # правая как максимум
x1 = np.linspace(lB, rB, 1000) # 1000 точек
y1 = np.ndarray(x1.shape) # значения функции той же размерности

print('y1 = f.fun(x[i])\n')
i = 0
for xi in x1:
    y1[i] = pears.f.fun(xi)
    i += 1

with open('data/test_function.csv','w') as csvfile:
    i = 0
    for xi in x1:
        csvfile.write("{0},{1}\n".format(xi,y1[i]))
        i += 1


# Claculation for BAR - there is an ERROR
with open('data/test_int.csv','w') as csvfile:
    ii = 1
    leftBorder = np.min(pears.x__)
    p_cal[0] = spint.quad(pears.f.fun, -4.111, leftBorder)[0] # ЗДЕСЬ 4.111 - это левая граница для функции
    csvfile.write("{0},{1},{2},{3}\n".format(leftBorder,-np.inf,leftBorder,p_cal[0]))
    rightBorder = np.max(pears.x__)
    count = len(pears.x__)
    delta = (rightBorder-leftBorder)/count
    print("Tets integrall:",spint.quad(pears.f.fun,-4.111,9)) # ЗДЕСЬ 4.111 - это левая граница для функции, 9 - правая
    for i in pears.x__[1:]:
        rightBorder = i
        p_cal[ii] = spint.quad(pears.f.fun, leftBorder, rightBorder)[0]
        csvfile.write("{0},{1},{2},{3}\n".format(i,leftBorder,rightBorder,p_cal[ii]))
        leftBorder = rightBorder
        ii += 1
    p_cal_ = spint.quad(pears.f.fun, leftBorder,9)[0] # ЗДЕСЬ 9 - это правая граница для функции
    csvfile.write("{0},{1},{2},{3}\n".format(leftBorder,leftBorder,np.infty,p_cal_))

# # Original Datas
plt.bar(pears.x__, y, color="white")  # width=pears.c,
plt.ylabel("N")
plt.xlabel(r"$x = (X - X_a)/c$")
plt.title(r'Histogram of $x$')
plt.savefig('data/data_bar_%d.png' % data_array)
plt.clf()
#
# # Curve
x_bar = pears.x__.copy()-delta
plt.plot(x1, y1)
plt.ylabel("P")
plt.xlabel("X")
plt.title('Pearson-Curve TYPE: {0}'.format(pears.type))
plt.bar(x_bar,p,  color="white")
plt.bar(x_bar,p_cal, width=0.4, color="0.8")
plt.savefig('data/data_curve_%d.png' % data_array)
plt.clf()

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
