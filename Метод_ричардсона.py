import math
import numpy as np
import sys
import copy
from numpy import array,identity,diagonal
from math import sqrt

import numpy as np

def jacobi(ain,tol = 1.0e-9): 
 
    def maxElem(a): # поиск самого большого не диагонального элемента a[k,l]
        n = len(a)
        aMax = 0.0
        for i in range(n-1):
            for j in range(i+1,n):
                if abs(a[i,j]) >= aMax:
                    aMax = abs(a[i,j])
                    k = i; l = j
        return aMax,k,l
 
    def rotate(a, k,l): # ротрация чтобы a[k,l] = 0
        n = len(a)
        aDiff = a[l,l] - a[k,k]
        if abs(a[k,l]) < abs(aDiff)*1.0e-36: t = a[k,l]/aDiff
        else:
            phi = aDiff/(2.0*a[k,l])
            t = 1.0/(abs(phi) + sqrt(phi**2 + 1.0))
            if phi < 0.0: t = -t
        c = 1.0/sqrt(t**2 + 1.0); s = t*c
        tau = s/(1.0 + c)
        temp = a[k,l]
        a[k,l] = 0.0
        a[k,k] = a[k,k] - t*temp
        a[l,l] = a[l,l] + t*temp
        for i in range(k):      # Case of i < k
            temp = a[i,k]
            a[i,k] = temp - s*(a[i,l] + tau*temp)
            a[i,l] = a[i,l] + s*(temp - tau*a[i,l])
        for i in range(k+1,l):  # Case of k < i < l
            temp = a[k,i]
            a[k,i] = temp - s*(a[i,l] + tau*a[k,i])
            a[i,l] = a[i,l] + s*(temp - tau*a[i,l])
        for i in range(l+1,n):  # Case of i > l
            temp = a[k,i]
            a[k,i] = temp - s*(a[l,i] + tau*temp)
            a[l,i] = a[l,i] + s*(temp - tau*a[l,i])
      
 
    a=np.copy(ain)
    n = len(a)
    maxRot = 5*(n**2)       
                            
    for i in range(maxRot): 
        aMax,k,l = maxElem(a)
        if aMax < tol: return diagonal(a)
        rotate(a,k,l)

def f(n, a, b, e, xt, w):
    x = np.zeros(n)
    t=False
    k=0
    l=jacobi(a)
    lmax=max(l)
    lmin=min(l)
    r0=2/(lmin+lmax)
    ro=(lmax-lmin)/(lmax+lmin)
    while(t==False):
        v=math.cos(math.pi*(2*k-1)/(2*n))
        w=r0/(1+ro*v)
        xa=[0]*n
        for j in range(n):
            xa[j]=x[j]
        for j in range(n):
            x[j]=(1-w)*x[j]*a[j][j]
            x[j]+=w*b[j]
            for l in range(j):
                x[j]-=w*x[l]*a[j][l]
            for l in range(j+1, n):
                x[j]-=w*x[l]*a[j][l]
            x[j]=x[j]/a[j][j]
        k+=1
        print(x)
        t=True
        for i in range(n):
            ee=abs(xa[i]-x[i])
            if (ee>=e):
                t=False
    print(k, " итераций, для достижения точности e =", e)    
    return x


print("z")

w=1
e=0.001

a= [[10.9, 1.2, 2.1, 0.9],
[1.2, 11.2, 1.5, 2.5],
[2.1, 1.5, 9.8, 1.3],
[0.9, 2.5, 1.3, 12.1]]




b=[-7.0,
5.3,
10.3,
24.6]

n=len(b)

xt=np.linalg.solve(a, b)
f(n, np.array(a), np.array(b), e, xt, w)

print(np.linalg.solve(a, b))
