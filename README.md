
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

## Matrix AD

The presentation "Efficient Automatic Differentiation of Matrix Functions" by Olsen of IBM 2012 presents very well the issues of matrix AD and its necessity also for matrix to scalar functions as they emerge in Machine Learning cases.


See also notes.md for some info.
