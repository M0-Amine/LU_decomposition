import numpy as np

# Multiplication of two matrices 'x'
def matrixProd(A, B):
    C = [[0 for i in range(len(B[0]))] for i in range(len(A))]
    if len(A[0]) == len(B):
        for i in range(len(A)):
            for j in range(len(B[0])):
                for k in range(len(B)):
                    C[i][j] += A[i][k] * B[k][j]
    return C
         
# Addition of two matrices '+'
def matrixAdd(A, B):
    C = [[0 for i in range(len(B[0]))] for i in range(len(A))]
    if (len(A[0]) == len(B[0]) and len(A) == len(B)):
        for i in range(len(A)):
            for j in range(len(A[0])):
                C[i][j] = A[i][j] + B[i][j]
    return C

# Scalar multiplication of a matrix '.'
def matrixScalar(l, A):
    C = [[0 for i in range(len(A[0]))] for i in range(len(A))]
    for i in range(len(A)):
        for j in range(len(A[0])):
            C[i][j] = l*A[i][j]
    return C

# Identity matrix
def I(n):
    return [[1 if i == j else 0 for i in range(n)] for j in range(n)]
    
# matrice dont je ne sais pas le nom
def E(n,i,j):
    E = [[1 if (y == i and x == j) else 0 for x in range(n)] for y in range(n)]
    return E

# Transvection matrix 
def T(n,i,j,y):
    return matrixAdd( I(n), matrixScalar(y, E(n,i,j)) )

# Dilatation matrix
def D(n,i,y):
    return matrixAdd( I(n), matrixScalar(y-1, E(n,i,i)) )
   
# Permutation matrix    
def P(n,i,j):
    l = [E(n,i,j), E(n,j,i), I(n), matrixScalar(-1, E(n,i,i)), matrixScalar(-1, E(n,j,j))]
    C = [[0 for i in range(n)] for j in range(n)]
    for M in l:
        C = matrixAdd(C, M)
    return C

# Displays a matrix
def printMatrix(M):
    print(np.array(M))

# Gauss elimination
def echelonForm(M):
    n = len(M)
    if M[0][0] != 0:
        for k in range(n):
            a = M[k][k]
            M = matrixProd(D(n, k, 1/a), M)
            for i in range(k+1, n):
                b = M[i][k]
                M = matrixProd(T(n, i, k, -b), M)
    return M 

# The inverse of a matrix M using the Gauss-Jordan elimination
def inverseOf(M):
    n = len(M)
    Id = I(n)
    if M[0][0] != 0:
        for k in range(n):
                a = M[k][k]
                M = matrixProd(D(n, k, 1/a), M)
                Id = matrixProd(D(n, k, 1/a), Id)
                for i in range(n):
                    if i != k:
                        b = M[i][k]
                        M = matrixProd(T(n, i, k, -b), M)
                        Id = matrixProd(T(n, i, k, -b), Id)
                        
    else:
        for i in range(1, len(M)):
            if M[i][0] != 0:
                M = matrixProd(P(n, 0, i), M)
                for k in range(n):
                    a = M[k][k]
                    M = matrixProd(D(n, k, 1/a), M)
                    Id = matrixProd(D(n, k, 1/a), Id)
                    for i in range(n):
                        if i != k:
                            b = M[i][k]
                            M = matrixProd(T(n, i, k, -b), M)
                            Id = matrixProd(T(n, i, k, -b), Id)
            break  
    return Id, M

# LU decomposiion
def LU(M):
    n = len(M)
    L = I(n)
    U = M
    if U[0][0] != 0:
        for k in range(n):
            a = U[k][k]
            for i in range(k+1, n):
                b = U[i][k]
                U = matrixProd(T(n, i, k, -b/a), U)
                L = matrixProd(L, T(n, i, k, b/a))
    return L, U 

# Example
A = [[1, 2],
     [2, 1]]

B = [[1, 2, 2],
     [2, 1, 2],
     [2, 2, 1]]

La, Ua = LU(A)
Lb, Ub = LU(B)
print("LU decomposition of A:")
print("A = La.Ua ="); printMatrix(La); printMatrix(Ua)
print("\nLU decomposition of B:")
print("B = Lb.Ub ="); printMatrix(Lb); printMatrix(Ub)


# verification
print("Verification:")
print("\nLa.Ua = ")
printMatrix(matrixProd(La, Ua))

print("\nLb.Ub = ")
printMatrix(matrixProd(Lb, Ub))
