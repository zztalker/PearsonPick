# -*- coding: utf-8 -*-
import numpy as np
import scipy.special as sp
import math 

class func_pearson:
    def __str__(self):
        return "Unknown function"
    def fun(self,x):
        return 0

class fI:
    def __init__(self,beta,r,mu,M):
        #посчитаем необходимые коэффициенты
        sigma = math.sqrt(mu[2])
        p = np.poly1d([1,-(r-2),-(r-1)+(4*r*r*(r+1))/(beta[1]*(r+2)*(r+2)+16*(r+1))])
        m = p.r
        self.m = ["m coefficient",0.,0.]
        if mu[3]>0:
            self.m[2] = max(m)
            self.m[1] = min(m)
        else:
            self.m[1] = max(m)
            self.m[2] = min(m)
        v = sigma/(2*(r-2))*math.sqrt(beta[1]*(r+2)*(r+2)+16*(r+1))
        self.a = ["a coefficient",0.,0.]
        self.a[1] = v*self.m[1]-M
        self.a[2] = v*self.m[2]+M
        
        y0_1 = pow(self.a[1],self.m[1])*pow(self.a[2],self.m[2])
        y0_1 = y0_1/pow(self.a[1]+self.a[2],self.m[1]+self.m[2]+1)
        
        y0_2 = sp.gamma(self.m[1]+self.m[2]+2)/(sp.gamma(self.m[1]+1)*sp.gamma(self.m[2]+1))
        
        self.y0 = y0_1*y0_2        
        return
    def __str__(self):
        s = "Type I function params:\n"
        s += "\n    a1,a2: "+str(self.a)
        s += "\n    m1,m2: "+str(self.m)
        s += "\n    y0: "+str(self.y0)
        s += "\n"
        return s
    def fun(self,x):		
        return self.y0*pow(1+x/self.a[1],self.m[1])*pow(1-x/self.a[2],self.m[2])

class fIplot(fI):
    def __init__(self,y0,a1,a2,m1,m2):
        self.m = ["m coefficient",0.,0.]
        self.a = ["a coefficient",0.,0.]
        self.a[1]=a1
        self.a[2]=a2
        self.m[1]=m1
        self.m[2]=m2
        self.y0=y0


class fII:
    def __init__(self,beta,mu):
        self.m = (5*beta[2]-9)/(2*(3-beta[2]))
        self.a = math.sqrt(2*mu[2]*beta[2]/(3-beta[2]))
        self.y0 = sp.gamma((2*self.m+1)/2)/(self.a*math.sqrt(math.pi)*sp.gamma(self.m+1))
        return
    def __str__(self):
        s = """ Type II function params:
            m: {0}
            a: {1}
            y0 {2}:""".format(self.m,self.a,self.y0)
        return s
    def fun(self,x):
        return self.y0*pow(1-x*x/(self.a*self.a),self.m)

class fIIplot(fII):
    def __init__(self,y0,a,m):
        self.a=a
        self.m=m
        self.y0=y0


class fIV:
    def __init__(self,beta,r,mu,M):
        #посчитаем необходимые коэффициенты
        self.a = []
        self.m = []
        self.y0 = 0
        return
        sigma = math.sqrt(mu[2])
        p = np.poly1d([1,-(r-2),-(r-1)+(4*r*r*(r+1))/(beta[1]*(r+2)*(r+2)+16*(r+1))])
        m = p.r
        self.m = [1,0.,0.]
        if mu[3]>0:
            self.m[2] = max(m)
            self.m[1] = min(m)
        else:
            self.m[1] = max(m)
            self.m[2] = min(m)
        print( beta[1]*(r+2)*(r+2)+16*(r+1))
        v = sigma/(2*(r-2))*math.sqrt(beta[1]*(r+2)*(r+2)+16*(r+1))
        self.a = [1,0.,0.]
        self.a[1] = v*self.m[1]-M
        self.a[2] = v*self.m[2]+M
        
        y0_1 = pow(self.a[1],self.m[1])*pow(self.a[2],self.m[2])
        y0_1 = y0_1/pow(self.a[1]+self.a[2],self.m[1]+self.m[2]+1)
        
        y0_2 = sp.gamma(self.m[1]+self.m[2]+2)/(sp.gamma(self.m[1]+1)*sp.gamma(self.m[2]+1))
        
        self.y0 = y0_1*y0_2        
        return
    def __str__(self):
        s = "Type I function params:\n"
        s += "\n    a1,a2: "+str(self.a)
        s += "\n    m1,m2: "+str(self.m)
        s += "\n    y0: "+str(self.y0)
        s += "\n"
        return s
    def fun(self,x):		
        return self.y0*pow(1+x/self.a[1],self.m[1])*pow(1-x/self.a[2],self.m[2])
    

class pearson:
    def __init__(self,x,p):
        mu = self.moments(x,p)
        self.sigma = math.sqrt(mu[2])		
        self.beta = ["beta array",0.,0.]
        self.beta[1] = pow(mu[3],2)/pow(mu[2],3)
        self.beta[2] = mu[4]/pow(mu[2],2)
        beta = self.beta		
        self.r = 6.*(beta[2]-beta[1]-1.)/(3.*beta[1]-2.*beta[2]+6.)
        r = self.r
        self.b = [0.,0.,0.]
        b = self.b
        b[0] = mu[2]*(r+1.)/(r-2.)
        b[1] = mu[3]*(r+2.)/(2.*mu[2]*(r-2.))
        b[2] = -1./(r-2.)
        self.M = -b[1]
        self.k = pow(b[1],2)/(4*b[2]*b[0])
        if self.k<0:
            self.type = "I"
            self.f = fI(beta,r,mu,self.M)
        elif (self.k==0):
            self.type = "?"
            self.f = func_pearson()
        elif (self.k<1):
            self.type = "IV"
            self.f = fIV(beta,r,mu,self.M)
        elif (self.k==1):
            self.f = func_pearson()
            self.type = "?"
        elif (self.k>1):
            self.f = func_pearson()
            self.type = "?"
        else:
            self.f = func_pearson()
            self.type = "unknown! "+str(self.k)
    def __str__(self):
        s = "Pearson params: \n"
        s += "\n	Mx: "+str(self.Mx)
        s += "\n	central moments: "+str(self.mu)
        s += "\n	sigma: "+str(self.sigma)
        s += "\n	beta[1,2]: "+str(self.beta[1:])
        s += "\n	r: "+str(self.r)
        s += "\n	b[0,1,2]: "+str(self.b)
        s += "\n	M: "+str(self.M)
        s += "\n	k: "+str(self.k)
        s += "\n    type of curve:"+str(self.type)
        s += "\n\n  Curve params:"
        s += str(self.f)
        s += "\n"
        return s
    def moments(self,x,p):
        self.mu = ["central momenth array",0,0,0,0]
        self.Mx = 0
        for i in range(len(x)):
            self.Mx += p[i]*x[i]
        for i in range(len(x)):
            self.mu[1] += p[i]*((x[i]-self.Mx))
            self.mu[2] += p[i]*pow(x[i]-self.Mx,2)
            self.mu[3] += p[i]*pow(x[i]-self.Mx,3)
            self.mu[4] += p[i]*pow(x[i]-self.Mx,4)
        return self.mu
