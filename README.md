# BEM

2d boundary element code in C++ for steady state and transient problems.

Goals: 

1. Full steady state and time domain code for scalar problems (Laplace, Helmholtz, Heat, and Wave).   

2. Eventually code FEM from scratch as well and work on BEM-FEM coupling, completing the conversion of my thesis to C++

3. Find a matrix library with linear algebra routines and nice constructors.

4. Eventually compile some of the libraries to be standalone and interface with Matlab - 
    see if there's any speed improvement.

	(Right now I have some of my doubts about the efficiency of this, since 
	Mex files have mxArray types for interfacing with Matlab, but I'd like to use 
	boost/STL or some other library for dealing with matrices.  I don't know
	if this back and forth copying would slow things down)

5. For time domain code, we're going to need an FFT library, I'll probably go with FFTW.

Some conventions I'm following are the use of std::vectors/std::arrays for pieces that need user initialization (and where it was easier to copy from Matlab, like with quadrature rules).   Anything that needs linear algebra operations done on it or is part of another calculuation with matrices and vectors is done with the boost::numeric::ublas classes since they come with built in matrix operations. This choice is because the boost classes don't have nice initializers like the STL does.  I'd like to improve this in the future.  I also don't like the difference in access operators between boost and STL.  It's just one more layer of syntax complexity between the user and program.
