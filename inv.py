# -*- coding: utf-8 -*-
"""
Created on Sat Oct  5 16:47:46 2019

@author: Widmark Kauê
"""
import numpy as np
import LU

AM = np.array([[ 1, 1, 0, 3] ,
               [ 2, 1,-1, 1] ,
               [ 3,-1,-1, 2] ,
               [-1, 2, 3,-1] ], dtype = np.float32 )

def inv(A):
    Nrow, Ncol = A.shape
    det = np.linalg.det(A) 
    assert Nrow == Ncol,  "inv: A matriz não é quadrada"
    assert det != 0, "inv: A matriz tem determinante nulo e não pode ser invertida" 
    B = np.zeros([Nrow, Ncol])
    L, U = LU.LUdecomposition(A)
    for j in range (0, Ncol ):
        b = np.zeros(Nrow)
        b[j] = 1
        x = LU.subpro(L, b)
        B[:,j] = LU.subreg(U, x)
    return B
     