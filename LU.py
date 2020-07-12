# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 10:59:53 2019

@author: Widmark Kauêr
VPL - Invesão de Matrizes
"""
import numpy as np
AM = np.array([[ 20, 80, 50] ,
               [ 16, 50, 60] ,
               [  8, 20, 16] ], dtype = np.float32 )
b = np.array([250, 100, 100], dtype = np.float32)
def LUdecomposition(A):
    Nrow , Ncol = A.shape
    assert Nrow == Ncol, "LUdecomposition: O método não pode ser utilizado para matrizes não quadradas"
    L = np.zeros([Nrow,Ncol])
    U = np.zeros([Nrow,Ncol])
    for i in range (0,Nrow):
        U[i,i] = 1
        for j in range (0,Ncol):
            if (i>=j):
                soma1 = 0
                for k in range (0,j):
                    soma1= soma1 + L[i,k]*U[k,j]
                L[i,j] = A[i,j] -  soma1
            else:
                assert L[i,i] != 0, "LUdecomposition: Esta matriz não pode ser fatorada pelo método de CROUT sem pivotamento"
                soma2 = 0
                for k in range (0,i):
                    soma2 = soma2 + L[i,k]*U[k,j]
                U[i,j] =( A[i,j] - soma2)/L[i,i]
    return L, U

def subreg(A,b):
    Nrow, Ncol = A.shape
    assert Nrow == Ncol, "subreg: A matriz de coeficientes não é quadrada."
    n = len(b)
    assert Nrow == n, "subreg: Os tamanhos de A e b são incompatíveis."
    x = np.zeros(Nrow)
    for i in range(Nrow-1,-1,-1):
        soma = 0
        for j in range(i+1,Nrow):
            soma =soma + A[i,j]*x[j]            
        x[i] = (b[i] - soma)/A[i,i]
    return x

def subpro(A,b):
    Nrow, Ncol = A.shape
    assert Nrow == Ncol, "subpro: A matriz de coeficientes não é quadrada."
    n = len(b)
    assert Nrow == n, "subpro: Os tamanhos de A e b são incompatíveis."
    x = np.zeros(Nrow)
    for i in range(0,Nrow):
        soma = 0
        for j in range(0,i):
            soma =soma + A[i,j]*x[j]            
        x[i] = (b[i] - soma)/A[i,i]
    return x

    
    
    
    
