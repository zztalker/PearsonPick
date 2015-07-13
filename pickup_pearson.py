# -*- coding: utf-8 -*-
import sys
import scipy.integrate as spint
import scipy.stats as spstat
import numpy as np
import pearson
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'verdana'})
rc('text.latex',unicode=True)
rc('text.latex',preamble='sepackage[utf8]{inputenc}')
rc('text.latex',preamble='sepackage[russian]{babel}')


#x = np.arange(1.5,21,1)
#y = np.array([0.11,0.69,3.67,4.93,2.23,1.1,1.8,2.62,4.69,2.41,1.14,7.26,9.34,5.70,8.15,12.25,19.64,8.42,3.79,0.51]);

#x = np.array([20,25,30,35,40,45,50,55,60,65,70,75,80],'float')
#x = x + 2.5
#y = np.array([11.,93.,162.,178.,176.,132.,101.,67.,40.,24.,12.,3.,1.])


#x = np.array([-0.5,	63.7,	76.09,	86.58,	101.47,	122.12,	219.2,	248.24],'float')
#y = np.array([0.5,	1	,2.5	,3.5,	10,	25,	30,	27.5])

#x = np.array([0.5,	13.7,	36.09,	56.58,	61.47,	80.12,	90.2,	108.24],'float')
#y = np.array([30,	20,	10,	16,	5,	4,	3,	2])

x = np.array([0,	1,	2,	3,	4,	5,	6,	7, 	8,	 9],'float')
y = np.array([29.5,	15,	12,	11,	10,	6.5,	6,	5,	4,	3],'float')


#y = y/100
#x = np.array([3,4,5,6,7,8],'float')
#y = np.array([2,8,20,40,20,10],'float')
p = y / y.sum()
p_cal = p.copy()

h = 1.0

if round(p.sum(),15)!=1.:
    print("Ошибка массива Y summ <> 1 ",p.sum())
    sys.exit(0)

if x.shape != y.shape:
    print( "Ошибка размерности!")
    print( "X", x.shape)
    print( "Y", y.shape)
    sys.exit(0)

# x - массив значений случайной величины, диапазон смещен
# y - массив частот попадания
# p - массив вероятностей

pears = pearson.pearson(x,p)
print(pears)

x1 = np.linspace(-pears.f.a[1],pears.f.a[2],100)
y1 = np.ndarray(x1.shape)
for i in range(0,100):
    y1[i] = pears.f.fun(x1[i])

leftBorder = -pears.f.a[1]
ii = 0
for i in x:
    rightBorder = i-pears.Mx
    p_cal[ii] = spint.quad(pears.f.fun,leftBorder,rightBorder)[0]
    ii+=1
    leftBorder = rightBorder

#plt.subplot(121)
plt.bar(x,p, width=20.0, color="red")
plt.bar(x,p_cal,width=10., color="green")
plt.ylabel("% состава смеси")
plt.xlabel("T (C) температура кипения")
plt.show()

plt.bar(x,p, width=20.0, color="red")
#plt.bar(x,p_cal,width=10., color="green")
plt.ylabel("% состава смеси")
plt.xlabel("T (C) температура кипения")
plt.show()

#plt.subplot(122)
plt.plot(x1,y1)
plt.ylabel("P")
plt.xlabel("T (C) температура кипения")
plt.show()

# посчитаем квантиль распределения хи квадрат
chi2 = 0.0
ii = 0
for p0 in p_cal:
    chi2 += ((p0-p[ii])**2)/(p0)
    ii += 1
i_chi2 = 0
b_chi2 = False
for chi2critical in spstat.chi2.ppf([0.01,0.05,0.1],ii-1):
    i_chi2 += 1 
    if chi2<=chi2critical:
        b_chi2 = True
        print( "Pearson criteria SUCCESSFUL for alpha=",0.01*i_chi2," chi2_emp = ",chi2,"<chi2(",0.01*i_chi2,",",ii-1,")=",chi2critical)
        break

if not b_chi2:
    print( "Pearson criteria UNSUCCESSFUL for alpha<0.1, gamma>0.9")
    
print( "Проверка функции плотности распредлеения, интеграл должен быть равен 1: ", spint.quad(pears.f.fun,-pears.f.a[1],pears.f.a[2]))
