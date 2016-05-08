
# Notes for the Future

##AD of Matrix Functions
- Paper: "Second Order Methods for Optimizing Convex
Matrix Functions and Sparse Covariance Clustering", Chin, 2013 http://www.stevenrennie.com/papers/jfd.pdf
- Talk: Efficient Automatic Differentiation of Matrix Functions,Olsen 2012 http://www.stevenrennie.com/papers/steve_rennie_matdiff_7_19_2012.pdf 
- Giles: An extended collection of matrix derivative results for forward and reverse mode algorithmic differentiation https://people.maths.ox.ac.uk/gilesm/files/NA-08-01.pdf

See: Olsen, Rennie, 2012 and Giles 2008 An extended collection of matrix derivative results for forward and reverse mode algorithmic differentiation. The second paper is limited basic operations while the first deals with general functions.

The first paper allows to cover:
- AXB AX'B
- inv(X) inv(X')
- logdet(AX)
- trace(AX)

Then moving to functional rules:
- chain rule: F(G(X)
- product rule (interesting):
	(Ik kron G'(X)) d(H,X) + (H(X) kr I_l) d(G,X)
- trace(G(X)H(X))
- trace(A inv(F) H(X))
- log(det(F(X)))

Using also scalar-matrix derivative: vec' (d(f,X)')

##Alternative
Why not use Tensor explicitly instead of using the Kronecker and the Box product? 

d(F,X) = [k,l,m,n] vs [kl,mn] whose definition is d(vec(F'),vec'(X'))
dd(F,X) = [k,l,m,n,m,n]

F(X) = trace(AX)B  becomes d(vec(B')vec'(A))
F(X) = AXB         becomes d(A kron B')
F(X) = AX'B        becomes d(A box  B')

example o RAD: trace(((I+X)inv X')X)

outer product as tensor: classic
kron  product as tensor: a[il]b[jk]
box   product as tensor: a[ik]b[jl]

Issues: meaning of trace and det in square forms [m1,m1,m2,m2] ~= [m1 m2,m1 m2]

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


# Details How to deal with matrices? 
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

#Matrix and Tensor Operations in Adept2

- all arithmetic plus assignment but also boolean (?) 
- array slicing
- sum, mean, product, any or along dimension
- conditional
- >> Automatic differentiation of matrix multiplication is slow, and is missing for special matrices and for linear algebra routines.