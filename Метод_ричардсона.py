import math
import numpy as np
import sys
import copy
from numpy import array,identity,diagonal
from math import sqrt

import numpy as np

def jacobi(myA, k_max):
    A = myA
    n = A.shape[0]
    C = np.copy(A)
    k = 1
    lambda_ = np.sqrt(np.max(C.diagonal()))*10.0**(-k)
    triu = np.abs(np.triu(C))
    for i in range(triu.shape[0]): triu[i,i] = 0.0
    index_matrix_max = np.unravel_index(np.argmax(triu), C.shape)
    matrix_max = triu[index_matrix_max]
    count = 0
    while (k != k_max and count < 10):
        i = index_matrix_max[0]
        j = index_matrix_max[1]
        d = np.sqrt((A[i,i] - A[j,j])**2+4*A[i,j]**2)
        c = np.sqrt(0.5*(1+(np.abs(A[i,i]-A[j,j]))/d))
        s = np.sign(A[i,j]*(A[i,i] - A[j,j]))*np.sqrt(0.5*(1-np.abs(A[i,i]-A[j,j])/d))
        for k in range(n):
            for l in range(n):
                if (k != i and k != j and l != i and l != j):
                    C[k,l] = A[k,l]
                elif(k != i and k != j):
                    C[k,i] = c*A[k,i] + s*A[k,j]
                    C[i,k] = C[k,i]
                    C[k,j] = -s*A[k,i]+c*A[k,j]
                    C[j,k] = C[k,j]
        C[i,i] = c**2*A[i,i]+2*c*s*A[i,j]+s**2*A[j,j]
        C[j,j] = s**2*A[i,i]-2*c*s*A[i,j]+c**2*A[j,j]
        C[i,j] = 0
        C[j,i] = 0
        A = np.copy(C)
        lambda_ = np.sqrt(np.max(C.diagonal()))*10**(-k)
        triu = np.abs(np.triu(C))
        print(triu)
        for i in range(triu.shape[0]): triu[i,i] = 0.0
        index_matrix_max = np.unravel_index(np.argmax(triu), C.shape)
        matrix_max = triu[index_matrix_max]
        while(k != k_max and matrix_max < lambda_):
            k += 1
            lambda_ = np.sqrt(np.max(C.diagonal()))*10**(-k)
        count += 1

    return(A.diagonal())
 
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
    l=jacobi(a, 10000)
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
