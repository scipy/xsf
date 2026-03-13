This is an implementation of the Mathieu fcns in C/C++.  The
implementation follows the prototype algos created in Matlab and
maintained on GitHub at
https://github.com/brorson/MathieuFcnsFourier.  This impl is a
header-only library for compatability with Scipy's xsf library.

The following Mathieu fcns are implemented:

*  Angular fcn ce(n,q,v)
*  Angular fcn se(n,q,v)
*  Radial (modified) fcn of first kind mc1(n,q,u)
*  Radial (modified) fcn of first kind ms1(n,q,u)
*  Radial (modified) fcn of second kind mc2(n,q,u)
*  Radial (modified) fcn of second kind ms2(n,q,u)

Here, n = fcn order, q = frequency (geometry) parmeter, v = angular
coord (radians), u = radial coord (au).

I also provide the following utility fcns:

*  Eigenvalue a_n(q)
*  Eigenvalue b_n(q)
*  Fourier coeffs A_n^k(q) for ce fcns
*  Fourier coeffs B_n^k(q) for se fcns

The goal is to provide a replacement of the Mathieu fcn suite used by
Scipy. 

These programs may be built the usual way on a Linux system using the
usual GNU build tools.  The main() function runs some simple sanity
checks on the functions.  In particular, it verifies some output
values against those computed by the Matlab programs.  I did a lot of
verification and accuracy testing on the Matlab implementations.
Therefore, tests run here just make sure the C implementation's
outputs match those from Matlab.  The code in main() also shows how to
invoke the various fcns.

Summer 2025, SDB
