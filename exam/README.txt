/*-------------------------- DESCRIPTION     -----------------------------*/
This project is written as part of the examination of the course "Practical
Programming and Numerical Methods" taught at Aarhus University in the summer
2021.

The one-sided Jacobi algorithm for Singular Value Decomposition (SVD) of
a real matrix is implemented with success. It uses a sweep-method, such that
the iteration has converged when the diagonal elements are invariant under
further iterations. It resembles what has previously been done in the 4th
homework exercise -- where a two sided Jacobi method was implemented -- 
with the exception that this one sided method does not constrain it
self to a symmetric matrix; nor a square matrix.

For the case of a tall/thin matrix, a simple QR-decomposition, using a
Gram-Schmitt (GS) -method is utilised, to obtain

A = Q*R

for Q (orthogonal) and R (upper-triangular). When A is of size (n, m), then
Q is (n, m) and R is (m, m). In succession, then the R-matrix is subject for a one-sided
jacobi method, so

R = U'*D*VT,

and substitution gives

A = U*D*VT

for U'=Q*U (square and orthogonal), D (diagonal) and VT (orthogonal)

and then proceeds to do a one-sided jacobi method on the R-matrix; to obtain
A = Q*U*D*VT = U'*D*VT

A lot of tests are done;
- QR tests the properties of the Q and R matrices
- one-sided jacobi also tests for orthogonality as well the equality of the
  decomposition and given matrix.

I have found, that the U-matrix satisfy UT*U = IDENTITY, but U*UT does not.
This has been ignored, as the exercise (02: Linear Systems) only tests for
QT*Q = IDENTITY as well. I remember that back then, the Q*QT also struggled,
and therefore decided, that this was not an assigned task; just as the exercise
02.


/*-------------------------- SELF-EVALUATION -----------------------------*/
The algorithm has been successfully implemented in the programming language of
C, using the make-utility-tool as an intelligent builder. The algorithm has
been tested in its details; the decomposition of A in U*D*VT, and that the
factorisation products, U and V, truly are orthogonal matrices.
I think the system is a success.

10/10.



/*-------------------------- ASSIGNED PROBLEM -----------------------------*/

ONE-SIDED JACOBI ALGORITHM FOR SINGULAR VALUE DECOMPOSITION

INTRODUCTION:
The singular value decomposition (SVD) of a (real square, for simplicity)
matrix A is a representation of the matrix in the form
A = U D VT ,

where matrix D is diagonal with non-negative elements and matrices U and V
are orghogonal. The diagonal elements of matrix D can always be chosen
non-negative by multiplying the relevant columns of matrix U with (-1).
SVD can be used to solve a number of problems in linear algebra.

PROBLEM
Implement the one-sided Jacobi SVD algorithm.

ALGORITHM
In this method the elementary iteration is given as
A → A J(θ,p,q)

where the indices (p,q) are swept cyclicly (p=1..n, q=p+1..n) and where the
angle θ is chosen such that the columns number p and q of the matrix
AJ(θ,p,q) are orthogonal. One can show that the angle should be taken from the
following equation (you should use atan2 function),
tan(2θ)=2apTaq /(aqTaq - apTap)

where ai is the i-th column of matrix A (check this).
After the iterations converge and the matrix A'=AJ 
(where J is the accumulation of the individual rotations) 
has orthogonal columns, the SVD is simply given as

A = U*D*VT

where
V=J, Dii=||a'i||, ui=a'i/||a'i||,

where a'i is the i-th column of matrix A' and ui us the i-th column of matrix U.

