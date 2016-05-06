
# autodiff
Automatic differentiation library for matrix operations (so far just an exercise)

##Status
- non matrix is ready
- matrix is not compiling

# Basic Principle
Automatic differentiation (AD) allows for efficiently and precisely compute the jacobian of a function Y = f(X) without incurring in the limitations of both symbolic and numerical differentiation. The basic principle of AD is the chain rule: given f(y(x)) we have: df/dx = df/dy dy/dx that can be generalized to Matrix Calculus. The starting point is the creation of an expression graph obtained by the evaluation trace of the f(y(x)) function. The interesting point of AD is that, due to its efficiency, can deal with functions that are defined at run-time.

Two methods do exist for computing the jacobian, that is the partial derivative for every relevant input variable: forward and reverse modes. In the forward mode the process runs from the inputs to the outputs, while in the reverse one the opposite. The advantage of the reverse mode is that with a single expression sweep it is possible to compute all the partial derivatives, while the forward mode has the advantage that the evaluation of the derivative is carried out together with the function evaluation. More formally: given f as a function with n inputs and m outputs, we are interested in the Jacobian as a matrix m by n. The forward mode iterates for every input n moving from the leaves to the output, while the reverse mode iterates for every output m from the output to the leaves. Note that the forward mode needs to propagate larger entities.

Things become interesting when the n inputs are matrices and the operations in the function are matrix operations: the resulting Jacobian can be expresse di partitioned (stacked) for with a total size J_n = sum(numel(input_i)). A different problem arises when the function has a vectorial or matrix output. This case is typical of mechanical problem during the computation of the velocity, while in machine learning we want to obtain the likelihood that is a single scalar value. 

#Usage

Replace your double arguments in the function with the sym class that will keep track of the expressions.

Then use jaconnum passing sym objects with the name of the variables interested 

#Implementation

The implementation follows the Operator Overloading approach together with Reverse-mode AD. A hierarchy of classes implements the operators, and then two jacobian (symbolic and numeric) are provided

#Example Problem

Tikhonov regularized maximum log likelihood for learning covariance Sigma given empirical S:

	f(Sigma) = -log det(Sigma) - trace(S inv Sigma) - ||inv Sigma||_frob^2

#Related Libraries

There are several open source libraries around and for this reason this is only a simple exercise. The most promising is Adept2:

- http://www.met.reading.ac.uk/clouds/adept2/
- https://github.com/stan-dev/math


#Possible Improvements

- matrix or tensor values (reverse-mode is optimal for any to scalar functions as in machine learning) keeping type-safe expression trees with/without fixed dimensions
- contextual evaluation
- make an exercise using expression templates
- higher order (hessian)
- lie group specialties

# Notes for the Future

## Extended Jacobian
When the expression involves matrix operations there are better solutions to the reverse mode that remove the need to iterate for every output. The most promising approach is called Extended Jacobian (EJ) as discussed in (https://dspace.lib.cranfield.ac.uk/bitstream/1826/4356/1/ForthICCS2010.pdf), applicable also to sparse systems.

The EJ approach relies on introducing the intermediate variables into the Jacobian computation, filling in a specific matrix C during recursion and finally solving for the target Jacobian with an algebrical operation:  (C - I_N) J = -P 
where:
- N=n+m+v
- P=[I_n ; 0_{mxn}] that is (m+n) by n, 
- J is (m+n+v) by n
- C (m+n) by N

J = [Jx Jv Jy]

Where:
	-In Jx = -In
	A Jx + B Jv = 0
	C Jx + D Jv + E Jy = 0
Note:
	A contains the use of x variables by some intermediate v
	E should be -Im

That gives: Jy = inv E (D inv B A - C)

Note that with matrix operations this method still holds with the side effect of increasing the number of intermediate values v.

##AD of Matrix Functions
See: Olsen, Rennie, 2012 and Giles 2008 An extended collection of matrix derivative results for forward and reverse mode algorithmic differentiation. The second paper is limited basic operations while the first deals with general functions.

The first paper supports matrix derivatives via box product, obtaining the identities:

- F(X) = X  => dF(X) = Imn
- F(X) = X' => dF(X) = Im box In
- F(X) = G(H(X)) => dF(X) = dG(H(X)) H'(X)
- F(X) = G(X) H(X) => dF(X) = ...
- F(X) = X^2 => dF(X) = Im kr X' + X kr Imn
- F(X) = inv(X) => - inv(X) kr inv(X')
- F(X) = inv(X') => ...
- F(X) = X^1/2
- F(X) = X^k
- F(X) = trace(AX) = A'
- F(X) = trace(AX') = A
- F(X) = log det (X) = X^{-1'}



## Details How to deal with matrices? 
(1) expand symbolic expressions to matrices and deal with them as single scalars 
(2) recognize special nature of matrix differentiation AND identify matrix structures (storage and semantics)

Example: (x-mu)' inv(S) (x-mu) 
Evaluation is easy, differentiation uses the special rules of matrix calculus, following the adjoint form of the recursive scheme
	
Properties:
- D(det(X)) = det(X) Tr(inv(X) DX)
- D(ln(det(X)) = Tr(inv(X) DX )
- chain rule
	U = f(X) 
	Dg(U)/DX = D g(f(X))/DX
	Dg(U)/DX = sum sum D g(U)/u  du/dx
- trace
- frobenius norm

References on Matrix Calculus
- GOOD: http://www.ee.ic.ac.uk/hp/staff/dmb/matrix/calculus.html

