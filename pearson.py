# -*- coding: utf-8 -*-
import numpy as np
import scipy.special as sp
import math

# прото  класс для всех будущих функций
# TODO: свойства l - границы, p - параметры, т.е. все параметры упаковываем в массив или словарь
class func_pearson:
    l = ["border",0,0]
    p = {'a':["a params",0,0]} #это для примера
    type = "Unknown"
    def __str__(self):
        param = ""
        for (key,value) in self.p:
            param = "{0} {1}, ".format(key,value)
        return "="*10+"""
        Pearson curve type {type}
        Main params: {param}
        """.format(type=self.type,param=param)+"\n="*10
    def fun(self, x):
        return 0

#TODO: убрать _str_ - сделать, что бы работал от наследуемого класса
class fI(func_pearson):
    def __init__(self, beta, r, mu, M, rk, sigma, k):
        # посчитаем необходимые коэффициенты
        sigma = math.sqrt(mu[2])
        p = np.poly1d([1, -(r - 2), -(r - 1) + (4 * r * r * (r + 1)) / (beta[1] * (r + 2) * (r + 2) + 16 * (r + 1))])
        m = p.r
        self.m = ["m coefficient", 0., 0.]
        if mu[3] > 0:
            self.m[2] = max(m)
            self.m[1] = min(m)
        else:
            self.m[1] = max(m)
            self.m[2] = min(m)
        v = sigma / (2 * (r - 2)) * math.sqrt(beta[1] * (r + 2) * (r + 2) + 16 * (r + 1))
        self.a = ["a coefficient", 0., 0.]
        self.a[1] = v * self.m[1] - M
        self.a[2] = v * self.m[2] + M
        y0_1 = pow(self.a[1], self.m[1]) * pow(self.a[2], self.m[2])
        y0_1 /= pow(self.a[1] + self.a[2], self.m[1] + self.m[2] + 1)
        y0_2 = sp.gamma(self.m[1] + self.m[2] + 2) / (sp.gamma(self.m[1] + 1) * sp.gamma(self.m[2] + 1))
        self.y0 = y0_1 * y0_2
        # считаем шраницы
        # нужные коэф
        S = (6 * (rk[4] - rk[3] ** 2 - 1)) / (3 * rk[3] ** 2 - 2 * rk[4] + 6)
        print('S = ', S)
        print(rk)
        print((S + 1) * (1 - k))
        print('k = ', k)
        t = math.sqrt(rk[3] ** 2 * (S + 2) ** 2 + 16 * (S + 1))
        # t = 4 * math.sqrt((S + 1) * (1 - k))
        q = ['q =', 0., 0.]
        q[1] = 0.5 * ((S - 2) + S * (S + 2) * rk[3] / t)
        q[2] = 0.5 * ((S - 2) - S * (S + 2) * rk[3] / t)
        print(q)
        L = sigma * t / 2
        print('L =', L)
        # сами границы
        self.l = ['border = ', -self.a[1], self.a[2]]
        #self.l[1] = q[1] * L / (S - 2)
        #self.l[2] =- q[2] * L / (S - 2)
        print(self.l)
        print(123456)

        return

    def __str__(self):
        s = "Type I function params:\n"
        s += "\n    a1,a2: " + str(self.a)
        s += "\n    m1,m2: " + str(self.m)
        s += "\n    y0: " + str(self.y0)
        s += "\n    l:  " + str(self.l)
        s += "\n"
        return s

    def fun(self, x):
        return self.y0 * pow(1 + x / self.a[1], self.m[1]) * pow(1 - x / self.a[2], self.m[2])


class fII(func_pearson):
    def __init__(self, beta, mu):
        self.m = (5 * beta[2] - 9) / (2 * (3 - beta[2]))
        self.a = math.sqrt(2 * mu[2] * beta[2] / (3 - beta[2]))
        self.y0 = sp.gamma((2 * self.m + 1) / 2) / (self.a * math.sqrt(math.pi) * sp.gamma(self.m + 1))
        self.l[1] = -self.a
        self.l[2] = self.a
        return

    def __str__(self):
        s = """ Type II function params:
            m: {0}
            a: {1}
            y0 {2}:""".format(self.m, self.a, self.y0)
        return s

    def fun(self, x):
        return self.y0 * pow(1 - x * x / (self.a * self.a), self.m)

#TODO: разобраться с реализацией функции и границами
class fIII(func_pearson):
    def __init__(self, r, sigma, summ):
        self.p = 4 / pow(r[3], 2) - 1
        l = ['granica l', 0., 0.]
        self.l[2] = sigma * (2 / r[3] - r[3] / 2)
        self.l[1] = -self.l[2]
        self.n0 = (summ / self.l[2]) * ((self.p + 1) / (pow(math.e, self.p) * sp.gamma(self.p + 1)))

        return

    def __str__(self):
        s = """ Type III function params:
         l : {0}
         p : {1}
         n0 : {2}""".format(self.l, self.p, self.n0)
        return s

    def fun(self, x):
        return self.n0 * pow(1 + x / self.l[2], self.p) * pow(math.e, -self.p * x / self.l[2])


#TODO: разобраться с реализацией функции и границами
class fIV(func_pearson):
    def __init__(self, r, sigma, summ):
        R = (6 * (r[4] - r[3] ** 2 - 1)) / (2 * r[4] - 3 * r[3] ** 2 - 6)
        print(R)
        self.q = (R + 2) / 2
        self.v = -(R * (R - 2) * r[3]) / math.sqrt(16 * (R - 1) - pow(r[3], 2) * pow(R - 2, 2))
        print('q=', self.q)
        print('v=', self.v)
        # print('sqrt=', 16 * (R - 1) - pow(r[3], 2) * pow(R - 2, 2))
        #self.l = ['granica l', 0., 0.]
        self.l[2] = (sigma / 4) * math.sqrt(16 * (R - 1) - pow(r[3], 2) * pow(R - 2, 2))
        self.l[1] = -self.l[2]
        # print(R, ' ', self.v)
        self.n0 = (summ / self.l[2]) * (1 / 1.80753)
        return

    def fun(self, x):
        return self.n0 * pow(1 + x ** 2 / self.l[2] ** 2, -self.q) * pow(math.e, -self.v * math.atan(x / self.l[2]))

    def __str__(self):
        s = "Type IV function params:\n"
        s += "\n    l: " + str(self.l)
        s += "\n    q: " + str(self.q)
        s += "\n    v: " + str(self.v)
        s += "\n"
        return s


#TODO: разобраться с реализацией функции и границами
class fV(func_pearson):
    def __init__(self, r, sigma, summ):
        self.l[1] = 0.
        self.l[2] = np.infty
        self.p = 4 + math.sqrt((8 + 4 * math.sqrt(4 + r[3] ** 2)) / r[3] ** 2)
        self.v = sigma * (self.p - 2) * math.sqrt(self.p - 3)
        self.n0 = (summ * pow(self.v, self.p - 1)) / (sp.gamma(self.p - 1))
        return

    def __str__(self):
        s = "Type V function params:\n"
        s += "\n    v: " + str(self.v)
        s += "\n    p: " + str(self.p)
        s += "\n    n0: " + str(self.n0)
        s += "\n"
        return s

    def fun(self, x):
        return self.n0 * pow(x, - self.p) * pow(math.e, - self.v / x)


class pearson:
    def __init__(self, x, p):
        mu = self.moments(x, p)
        self.l = ['l border: ',0,0]
        self.sigma = math.sqrt(mu[2])
        self.beta = ["beta array", 0., 0.]
        self.beta[1] = pow(mu[3], 2) / pow(mu[2], 3)
        self.beta[2] = mu[4] / pow(mu[2], 2)
        beta = self.beta
        self.r = 6. * (beta[2] - beta[1] - 1.) / (3. * beta[1] - 2. * beta[2] + 6.)  # похож на S коэф
        # r = self.r
        self.rk = ['r-koef', 0., 0., 0., 0.]
        self.rk[3] = mu[3] / pow(self.sigma, 3)
        self.rk[4] = mu[4] / pow(self.sigma, 4)

        self.b = [0., 0., 0.]
        b = self.b
        b[0] = mu[2] * (self.r + 1.) / (self.r - 2.)
        b[1] = mu[3] * (self.r + 2.) / (2. * mu[2] * (self.r - 2.))
        b[2] = -1. / (self.r - 2.)
        self.M = -b[1]
        self.k = pow(b[1], 2) / (4 * b[2] * b[0])
        if round(self.k, 1) == 0:
            self.type = "II"
            self.f = fII(beta, mu)
        elif round(self.k, 0) == 1:
            self.type = "V"
            self.f = fV(self.rk, self.sigma, self.Mx)
        elif (self.k < -1):
            self.f = fIII(self.rk, self.sigma, self.Mx)
            self.type = "III"
        elif self.k < 0:
            self.type = "I"
            self.f = fI(beta, self.r, mu, self.M, self.rk, self.sigma, self.k)
        elif (self.k < 1):
            self.type = "IV"
            self.f = fIV(self.rk, self.sigma, self.Mx)
        elif (self.k > 1):
            self.f = fIII(self.rk, self.sigma, self.Mx)
            self.type = "III"
        else:
            self.f = func_pearson()
            self.type = "unknown! " + str(self.k)

    def __str__(self):
        s = "\nPearson params: \n"
        s += "\nXmin: " + str(self.Xmin)
        s += "\nXmax: " + str(self.Xmax)
        s += "\nc: " + str(self.c)
        s += "\nMx: " + str(self.Mx)
        s += "\nXa: " + str(self.Xa)
        s += "\nmoments: "
        s += "\nN\tnormal\tcentral"
        for i in range(1, 5):
            s += "\n{0}\t{1}\t{2}".format(i, self.m[i], self.mu[i])
        s += "\n\nsigma: " + str(self.sigma)
        s += "\nbeta[1,2]: " + str(self.beta[1:])
        s += "\nr: " + str(self.r)
        s += "\nb[0,1,2]: " + str(self.b)
        s += "\nM: " + str(self.M)
        s += "\nk: " + str(self.k)
        s += "\ntype of curve:" + str(self.type)
        s += "\n\nCurve params:"
        s += str(self.f)
        s += "\n"
        return s

    def moments(self, x, p):
        self.m = ["momenth array", 0, 0, 0, 0]
        self.mu = ["central momenth array", 0, 0, 0, 0]
        self.Mx = 0
        self.Xmax = x[0]
        self.Xmin = x[0]
        for i in range(len(x)):
            if self.Xmax < x[i]: self.Xmax = x[i]
            if self.Xmin > x[i]: self.Xmin = x[i]
            self.Mx += p[i] * x[i]
        self.c = (self.Xmax - self.Xmin) / (x.shape[0] - 1)
        # Get an Xa - there is X from middle of array
        tSum = 0
        for i in range(len(x)):
            tSum += p[i]
            if tSum >= 0.5:
                # Get a middle! 
                self.Xa = x[i] + self.c / 2
                break

        self.x_ = x.copy()
        self.x__ = x.copy()

        for i in range(len(x)):
            self.x_[i] = x[i] + self.c/2 - self.Xa
            self.x__[i] = (x[i] + self.c/2 - self.Xa) / self.c
            self.m[1] += p[i] * (((x[i] + self.c / 2) - self.Xa) / self.c)
            self.m[2] += p[i] * pow(((x[i] + self.c / 2) - self.Xa) / self.c, 2)
            self.m[3] += p[i] * pow(((x[i] + self.c / 2) - self.Xa) / self.c, 3)
            self.m[4] += p[i] * pow(((x[i] + self.c / 2) - self.Xa) / self.c, 4)

        self.mu[1] = 0
        self.mu[2] = self.m[2] - self.m[1] ** 2
        self.mu[3] = self.m[3] - 3 * self.m[2] * self.m[1] + 2 * pow(self.m[1], 3)
        self.mu[4] = self.m[4] - 4 * self.m[3] * self.m[1] + 6 * self.m[2] * pow(self.m[1], 2) - 3 * pow(self.m[1], 4)

        return self.mu
