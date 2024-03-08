# The square-root form cookbook

Algorithms in square-root form have several advantages in terms of numerical stability. This is especially pertinent for algorithms that require matrices to stay positive-semidefinite even though eigenvalues might approach zero. 

## Notation

 - $A, B, C, ...$ Positive semidefinite matrices
 - $G$ General matrix
 - $I$ Identity matrix
 - $Q$ Unitary matrix 

## Basic definitions

Matrix square root.

$$
A = L_A L_A^T = U_A^T U_A
$$

### Cholesky

Cholesky (lower) decomposition is defined for positive-semidefinite matrices.

$$
L_A = \text{chol}(A)
$$

### QR decomposition 

QR decomposition is defined for real rectangular matrices.

$$
Q_A U_A = \text{qr}(AA)
$$

which satisfies our square-root form as,

$$
Q^T Q = I
$$
$$
(Q_A U_A)^T (Q_A U_A) = U_A^T Q^T Q U_A = U_A^T I U_A =  U_A^T  U_A  = AA
$$

Relation to Cholesky decomposition for PSD matrices.

$$
L_A = U_A^T
$$



## Square-formulae

### Outer products

$$
A = G B G^T 
$$

$$
L_A L_A^T= G L_B (G L_B)^T = G L_B L_B^T G^T = G B G^T
$$

$$
\implies L_A = G L_B
$$



### Inner products

$$
A = G^T B G 
$$

$$
L_A L_A^T= G^T L_B (G^T L_B)^T = G^T L_B L_B^T G = G^T B G
$$

$$
\implies L_A = G^T L_B
$$



### Multiple terms with QR trick

$$
A = B + G C G^T 
$$

Notice that,

$$
\begin{bmatrix} L_B & G L_C \end{bmatrix} \begin{bmatrix} L_B^T \\ (G L_C)^T \end{bmatrix} = L_B L_B^T + G L_C (G L_C)^T = B + G C G^T = A
$$

then, 

$$
Q_A, U_A = \text{qr}\left(\begin{bmatrix} L_B^T \\ (G L_C)^T \end{bmatrix}\right) 
$$

$$
(Q_A, U_A)^T (Q_A, U_A) = U_A^T Q^T Q U_A = U_A^T U_A =  A
$$

So $U_A$ is the required upper-triangular form, and $L_A = U_A^T$.


### Inverses

$$
A^{-1} = (L_A L_A^T)^{-1} = L_A^{-1} L_A^{-T}
$$

$$
\implies L_A^{-1} = \text{chol}(A^{-1})
$$

### Linear systems

The linear system,

$$
AX = B
$$

$$
(L_A L_A^T) X = L_A (L_A^T X) = B
$$

Now make use of the triangular structure of $L_A$ to solve in two parts. First using a lower triangular solve,

$$
Y = \text{Lsolve}(L_A, B)
$$

Then an upper triangular solve.

$$
X = \text{Usolve}(L_A^T, Y)
$$

Note that the system above will often be seen in the form,

$$
X = A^{-1} B
$$

But the identities above remain. 

### Inverse Cholesky

Another special case of the above formula is when $B=I$ i.e.

$$
X = A^{-1}
$$

Plugging this into the above,

$$
Y = \text{Lsolve}(L_A, I)
$$

$$
X = A^{-1} = L_A^{-1} L_A^{-T} = \text{Usolve}(L_A^T, Y)
$$

$$
L_A^T L_A^{-1} L_A^{-T} = Y = L_A^{-1}
$$

$$
\implies \text{Lsolve}(L_A, I) =  L_A^{-1}
$$

