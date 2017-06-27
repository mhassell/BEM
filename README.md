# BEM

2d boundary element code in C++

Goals: Full steady state and time domain code for scalar problems (Laplace,
Helmholtz, Heat, and Wave).   

Some conventions I'm following are the use of std::vectors/std::arrays for pieces
that need user initialization (and where it was easier to copy from Matlab, like with 
quadrature rules).   Anything that needs linear algebra operations done on it or is 
part of another calculuation with matrices and vectors is done with the 
boost::numeric::ublas classes since they come with built in matrix operations.
This choice is because the boost classes don't have nice initializers like the 
STL does.  I don't like the mismatch, but that's what I have for now.  It would be nice to not have the type mismatch, but that's a future improvement.
