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

# v1 mu3 = 0

raw_data = np.loadtxt('data.csv',dtype=float)
x = raw_data[0]
y = raw_data[1]

# Проверим данные.
# Предполагается, что два массива одинаковы по размерности
if (x.shape!=y.shape):
    print("Некорректные размерности в массиве данных. x -",x.shape,"y",y.shape)
    exit()
# Проверим, что Y - Это частота.
if not((y.sum()==1)):
    print("Y - это не частота, преобразуем")
    p = y / y.sum()
else:
    p = y

# p - частота, x

res = pearson.pearson(x,p)
print(res)

