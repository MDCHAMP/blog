## Polynomials: fractional and/or rational
---

Hello and welcome to another instalment of 'the derivation reviewer 2 rejected'. On the menu today is **your dads** favourite frequency domain modal analysis curve fitter: *Rational polynomial decomposition (RFP)*. This particular curve-fitting tool was featured in the very first IMAC all the way back in 1983 in Orlando (I wonder if academic spam and massive paywalls were already a thing?) Today we will follow the derivation of [that very paper](https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=b263f0a47cf7636385643893c2c8e5debab409b2). Along the way we will find some nice tricks with the Hermitian transpose and more than a couple of flagrant errors. 

    Richardson, Mark H., and David L. Formenti. "Parameter estimation from frequency response measurements using rational fraction polynomials." Proceedings of the 1st international modal analysis conference. Vol. 1. Union College Schenectady, NY, 1982.


Very much the great-grandaddy of modern-day paywalled frequency domain tools such as the snappily named poly-reference least-squares complex frequency domain estimator that is often known by a different name (rhymes with jolly-flax), RFP attacks the modal transfer function by expanding the sum over the modes into a ratio of two polynomials. Recalling the general form of the modal transfer function from excitation at location $i$ to response at location $j$.

$$
H_{ij}(\bm{\omega}) = \sum_{r=1}^n \frac{\phi_i^{(r)}\phi_j^{(r)}}{\omega_r^2 - 2j\zeta_r\omega_r\bm{\omega} + \bm{\omega}^2}
$$

The core idea in RFP (and many other frequency domain curve fitters) is to reparametrize the above as a ratio of two polynomials,

$$
H(\bm{\omega}) = \frac{\sum_{k=0}^m a_k (j\bm{\omega})^k}{\sum_{k=0}^{n} b_k (j\bm{\omega})^k}
$$

Equivalently in matrix from this can be written,

$$	
H = (\Phi_{\bm{a}} \bm{a}) (\Phi_{\bm{b}} \bm{b})^{-1}
$$

where the matrices $\Phi_{\bm{a}}$ and $\Phi_{\bm{b}}$ are polynomial bases of order $n$ and $m$ respectively for the vector of frequency lines $\bm{\omega}$.

In order for this formulation to be well behaved, these polynomials should have the following properties:

- The coefficients ($a_k$ and $b_k$) should be real.
- The order of the denominator polynomial $n$ should be twice the number of modes we wish to locate.
- The order of the numerator polynomial $m$ should be the same size or greater than the denominator.

The trick here is that we want to find a least-squares solution for the coefficient vectors $\bm{a}$ and $\bm{b}$ in one hit. The trouble is (much like trying to trying to unlock my bike after a few in the deer) there is an identifiability problem; we could happily extract any constant from the ratio and leave the value unchanged. To overcome this, $b_n$ is given an arbitrary value (in the paper its unity but we can use any value). Extracting this term from the denominator basis (and sloppily redefining $\bm{b}$ and $\Phi_{\bm{b}}$ as the coefficients and  basis of a polynomial with order $n-1$),

$$
b_n (j\bm{\omega})^n + \Phi_{\bm{b}} \bm{b} := \Phi_{\bm{b}} \bm{b}
$$

(notation purists look away now), the model is now,

$$
H(\bm{\omega}) = (\Phi_{\bm{a}} \bm{a}) (b_n (j\bm{\omega})^n + \Phi_{\bm{b}} \bm{b})^{-1}
$$


Least squares solutions to problems require one thing... and it isn't the the Cholesky decomposition... yet! We need to form a quadratic problem. Lets consider the error between the model above and some observed data at some observed frequency lines $\hat{H}(\bm{\omega})$. The error is,

$$
\bm{e} = H(\bm{\omega}) - \hat{H}(\bm{\omega}) = (\Phi_{\bm{a}} \bm{a}) (b_n (j\bm{\omega})^n + \Phi_{\bm{b}} \bm{b})^{-1} - \hat{H}(\bm{\omega})
$$

$$
\bm{e} = [\Phi_{\bm{a}}\bm{a} - ((j\bm{\omega})^n + \bm{b} \Phi_{\bm{b}})\hat{H}(\bm{\omega})](b_n (j\bm{\omega})^n + \Phi_{\bm{b}} \bm{b})^{-1}
$$

$$
\bm{e} \propto \Phi_{\bm{a}}\bm{a} - ((j\bm{\omega})^n + \bm{b} \Phi_{\bm{b}})\hat{H}(\bm{\omega})
$$

Eagle-eyed viewers might have spotted some... *ahem*... issues with the reasoning above, but we will return to these later.
 
Lets now make some notational changes to keep in sync with the original paper. 

$$
P = \Phi_{\bm{a}}
$$

$$
T = \Phi_{\bm{b}} \hat{H}(\bm{\omega})
$$

$$
\bm{w} = b_n (j\bm{\omega})^n	\hat{H}(\bm{\omega})
$$

So now the error is defined as 

$$
J(\bm{a},\bm{b} | \hat{H}) = \bm{e}^\dag \bm{e} = (P\bm{a} - (\bm{w} + T\bm{b}))^\dag (P\bm{a} - (\bm{w} + T\bm{b}))
$$

where the $\dag$ refers to the Hermitian transpose (i.e the transpose and conjugation operations in one). Hopefully it is clear that this expression is a positive quadratic in both $\bm{a}$ and $\bm{b}$. Quadratic expressions such as these have a single optima. In order to derive it, we must first expand out the squared error term.  

Expanding out we have several terms to consider, lets look at them one-by-one:

First up the quadratic term in $\bm{a}$,

$$
(P\bm{a})^\dag(P\bm{a}) = \bm{a}^TP^\dag P\bm{a}
$$

Where we have used the property of the transpose $(AB)^T = B^TA^T$ that also applies to Hermitian transposes. Easy enough so far...

Next up the quadratic term in $\bm{b}$,

$$
(\bm{w} + T\bm{b})^\dag (\bm{w} + T\bm{b})
$$

$$
\bm{w}\dag\bm{w} + \bm{b}^T T^\dag T\bm{b} + \bm{w}^\dag T \bm{b} + (T\bm{b})^\dag \bm{w}
$$

Notice that we can rewrite the last two terms in the above as the Hermitian transpose of each other.

$$
\bm{w}^\dag T \bm{b} + (T\bm{b})^\dag \bm{w} = ((T\bm{b})^\dag \bm{w}) + ((T\bm{b})^\dag \bm{w})^\dag
$$

Because we know that these expressions are scalar (they must be because our objective function $J$ is scalar-valued) we know that the Hermitian transpose is equivalent to the conjugate. A complex number added to the conjugate of itself is simply: $a + a^* = 2\mathbb{Re}(a)$. We can therefore rewrite our second term as,

$$
\bm{w}^\dag\bm{w} + \bm{b}^T T^\dag T\bm{b} + 2\mathbb{Re}(\bm{b} T^\dag \bm{w})
$$

Finally now considering the cross-quadratic terms

$$
-(P\bm{a})^\dag (\bm{w} + T\bm{b}) - (\bm{w} + T\bm{b})^\dag (P\bm{a})
$$

Using the same trick as above we can simplify this to,

$$
 - \bm{a}^T P^\dag \bm{w} - \bm{w}^\dag P a
 - \bm{a}^T P^\dag T \bm{b} - \bm{b}^\dag  T^\dag  P\bm{a} 
$$

$$
 - (\bm{a}^T P^\dag \bm{w}) - (\bm{a}^T P^\dag \bm{w})^\dag 
 - (\bm{a}^T P^\dag T \bm{b}) - (\bm{a}^T P^\dag T  \bm{b})^\dag 
$$

$$
 - 2\mathbb{Re}(\bm{a}^T P^\dag \bm{w})
 - 2\mathbb{Re}(\bm{a}^T P^\dag T \bm{b})
$$

Overall we can now write:

$$
J(\bm{a},\bm{b} | \hat{H}) =  
\bm{a}^TP^\dag P\bm{a} 
+ \bm{w}\dag\bm{w} + \bm{b}^T T^\dag T\bm{b} \\ 
+ 2\mathbb{Re}(\bm{b} T^\dag \bm{w}) 
- 2\mathbb{Re}(\bm{a}^T P^\dag \bm{w})
- 2\mathbb{Re}(\bm{a}^T P^\dag T \bm{b})
$$

Comparing this to the result in the IMAC paper, we see that the $2\mathbb{Re}(\bm{b} T^\dag \bm{w})$ term has opposite sign in our derivation... [This is in fact a mistake in the original paper](https://www.tandfonline.com/doi/epdf/10.1080/10652460701511301?needAccess=true&role=button) :/

The conditions for the minimum error (literally the least-squared-error) in our objective function are,

$$
\frac{\partial J(\bm{a}, \bm{b} | \hat{H})}{\partial \bm{a}} \bigg|_{a=a^*,\  b=b^*}=0
$$
and,
$$
\frac{\partial J(\bm{a}, \bm{b} | \hat{H})}{\partial \bm{b}} \bigg|_{a=a^*,\  b=b^*}=0
$$

Evaluating these derivatives is highly tedious but largely straightforward. Starting with the former condition,

$$
\frac{\partial J}{\partial \bm{a}} = 
2\mathbb{Re}(P^\dag P) \bm{a}
- 2\mathbb{Re}(P^\dag T) \bm{b}
- 2\mathbb{Re}(P^\dag \bm{w}) 
$$

Note that the coefficient in $\bm{a}$ also differs from the original paper... another mistake in their derivation I'm afraid :/

Now the second condition,

$$
\frac{\partial J}{\partial \bm{b}} = 
- 2\mathbb{Re}(T^\dag P) \bm{a}
+ 2\mathbb{Re}(T^\dag T) \bm{b}
+ 2\mathbb{Re}(T^\dag \bm{w}) 
$$


Now applying the condition of zero gradient at the optimum, we may write the simultaneous system,

$$
2\mathbb{Re}(P^\dag P) \bm{a}
- 2\mathbb{Re}(P^\dag T) \bm{b}
= 2\mathbb{Re}(P^\dag \bm{w}) 
$$

$$
- 2\mathbb{Re}(T^\dag P) \bm{a}
+ 2\mathbb{Re}(T^\dag T) \bm{b}
= - 2\mathbb{Re}(T^\dag \bm{w}) 
$$

The above block-linear system may now be succinctly be represented as,

$$
\mathbb{Re}\left(\begin{bmatrix} P^\dag P & -P^\dag T \\ - T^\dag P & T^\dag T \end{bmatrix}\right) \begin{bmatrix}\bm{a}\\ \bm{b}\end{bmatrix} = \mathbb{Re}\left(\begin{bmatrix} P^\dag \bm{w} \\ -T^\dag \bm{w} \end{bmatrix}\right)
$$

Solving this system is as easy (or as hard) as you and your data would like it to be! All of our favourite numerical tricks and recipes can now be applied and we can lose countless days in a happy trance of matrix conditioning and square-root forms and Cholesky factors and rank one updates and where was I... 

So anyway, modal analysis. With the coefficients of the numerator and denominator polynomials computed, we can recover the modal properties directly. The natural frequencies are simply the roots of the denominator polynomial (remembering to replace the coefficient of $b_n$ that we arbitrarily chose earlier!). The result is a polynomial with $n$ roots appearing in $n/2$ complex-conjugate pairs. The damping ratios can then be recovered as the negative real parts of the roots divided by the natural frequencies. Obtaining the modeshapes is more complex and is left as an exercise for the reader ;p

A python (requiring `numpy`) implementation of the RFP model is included below.

```python
# Hermitian transpose
def HT(a):
    return a.conj().T

# Rational fraction polynomial model
def RFP(H, w, n_modes, oob_terms=0):
    # Specify the orders of our approximation
    m = 1 + (n_modes * 2) + oob_terms # number of coefficients in the numerator polynomial
    n = 1 + (n_modes * 2)  # number of coefficients in the denominator polynomial

    # Build monomial basis matricies
    Phi_a = (1j*w[:, None]) ** np.arange(m) 
    Phi_b = (1j*w[:, None]) ** np.arange(n) 

    P = Phi_a 
    T = Phi_b[:,:-1] * H[:, None]
    W = Phi_b[:, -1] * H
    PT = -HT(P)@T
    
    # form the block matricies
    M = np.block([[HT(P)@P, PT], [PT.T, HT(T)@T]])
    x = np.block([HT(P)@W, -HT(T)@W])

    # Solve and extract the coefficients of the polynomials
    AB = np.linalg.solve(np.real(M), np.real(x))
    a = AB[:m, None]
    b = np.append(AB[m:], 1)[:, None]

    # Generate the predicted FRF
    H_pred = (Phi_a @ a) / (Phi_b @ b)

    # Pull out the modal porperties
    roots_b = sorted(np.roots(np.flip(b[:, 0])))[::-2]  # remove every other because they are conj pairs   
    wns = np.abs(roots_b)
    zetas = -np.real(roots_b) / wns
    return H_pred, wns, zetas
```

So far we have ignored a rather important limitation of the approach. Earlier, when defining our error criteria, we chose,

$$
\bm{e} = \Phi_{\bm{a}}\bm{a} - ((j\bm{\omega})^n + \bm{b} \Phi_{\bm{b}})\hat{H}(\bm{\omega})
$$

In fact, this expression is only proportional to the true difference between our model and the data,

$$
\bm{e} = [\Phi_{\bm{a}}\bm{a} - ((j\bm{\omega})^n + \bm{b} \Phi_{\bm{b}})\hat{H}(\bm{\omega})](b_n (j\bm{\omega})^n + \Phi_{\bm{b}} \bm{b})^{-1}
$$

This is what those in the business call '[a bit of a yikes](https://media.tenor.com/HTz-RXN1qogAAAAC/tacos-de.gif)'. Not only are we transferring output noise onto the input (bad), we are also heavily biasing the higher frequency lines in the least squares formulation (v. bad). In practice this can lead to good fits to experimental data in the highest frequencies and poor fits in the low frequency where all the interesting modal dynamics are... 

There are several ways around these issues, but that is a story for another time...