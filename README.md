# BEM

2d boundary element code in C++ for steady state and transient problems.

Goals: 

1. Full steady state and time domain code for scalar problems (Laplace, Helmholtz, Heat, and Wave).   

2. Eventually code FEM from scratch as well and work on BEM-FEM coupling, completing the conversion of my thesis to C++

3. (DONE) Find a matrix library with linear algebra routines and nice constructors.

4. Eventually compile some of the libraries to be standalone and interface with Matlab - 
    see if there's any speed improvement.

	(Right now I have some of my doubts about the efficiency of this, since 
	Mex files have mxArray types for interfacing with Matlab, but I'd like to use 
	boost/STL or some other library for dealing with matrices.  I don't know
	if this back and forth copying would slow things down)

5. For time domain code, we're going to need an FFT library, I'll probably go with FFTW.

6. (DONE) Need to implement a Matlab like matrix library for my linear algebra needs.  This 
	would include a general solver (just with boost's LU routines), the kronecker product,
	and I would like a linspace that returns stl::vectors so I can generate some 
	bookeeping arrays a bit faster than now.  Some reshape routines would be helpful
	as well, with overloaded definitions for different array types. This would include
	something like a toRowVector() and toColumnVector() calls to do Matlab's equivalent of v(:).

7. I'm also hoping to eventually wrap all of the BEM code in a BEM:: namespce to keep it well separated from any future
	projects that may have overlaps in function names.   The matrix library described in #6 could then be pulled out
	and used across other projects in the future.	

(NOT ANYMORE) Some conventions I'm following are the use of std::vectors/std::arrays for pieces that need user initialization (and where it was easier to copy from Matlab, like with quadrature rules).   Anything that needs linear algebra operations done on it or is part of another calculuation with matrices and vectors is done with the boost::numeric::ublas classes since they come with built in matrix operations. This choice is because the boost classes don't have nice initializers like the STL does.  I'd like to improve this in the future.  I also don't like the difference in access operators between boost and STL.  It's just one more layer of syntax complexity between the user and program.

8/24/17 Update: I think I'm going to port a lot of the library to use eigen libraries.  These have nice constructors and tie in well with linear algebra routines.  Also, note to self: I'm putting the eigen header files in /usr/local/include

10/22/17 Update: Finally have Laplace code fully running.  I don't know if I'm ever going to do time domain code, since I made some decisions that would require reversing to implement that.
		 In lieu of function handles, it would be nice to use lambdas in the newer c++ standards, but I didn't need that for steady state code.
		 This, at the moment, is much slower than the equivalent Matlab code.  I'm not completely comparing apples to apples with some of my timing choices, but it is not competitive.  I'd like to work on this.

10/28/17 Update: I spent a few days looking for speed improvements using different profilers and was able to find some places to speed up the code.
		 Now the code is about 10x faster than Matlab!
