#ifndef MATHIEU_FCNS_H
#define MATHIEU_FCNS_H

#include "../config.h"
#include "../error.h"
#include <math.h>
#include "matrix_utils.h"
#include "mathieu_coeffs.h"
#include "besseljyd.h"


/*
 *
 * This is part of the Mathieu function suite -- a reimplementation
 * of the Mathieu functions for Scipy.  This file holds the function
 * implementations themselves.  The prototype was written in Matlab
 * and validated.  This is a translation from Matlab to C.
 * 
 */

namespace xsf {
namespace mathieu {

  //==================================================================
  int mathieu_ce(int m, double q, double v, double *ce, double *ced) {
    // This computes the Mathieu fcn ce
    // Inputs:
    // m = Mathieu fcn order (scalar)
    // q = frequency parameter (scalar)
    // v = angle in radians (scalar)
    // Outputs:
    // ce = value of fcn for these inputs (scalar)
    // ced = value of fcn deriv w.r.t. v for these inputs (scalar)
    // Return code:  
    // Success = 0

    int retcode = SF_ERROR_OK;

    // Check input domain and flag any problems.
    if (m>500) {
      // Don't support absurdly larger orders for now.
      *ce = std::numeric_limits<double>::quiet_NaN();
      *ced = std::numeric_limits<double>::quiet_NaN();
      return SF_ERROR_DOMAIN;
    }
    // abs(q) > 1000 leads to low accuracy.
    if (abs(q)>1.0e3d) retcode = SF_ERROR_LOSS;
      
    
    // I find the peak Fourier coeff tracks m.  Therefore,
    // adjust the matrix size based on order m.  Later make this
    // a fcn of q also since the distribution of coeff mags allegedly
    // flattens out for large q.
    int N = m+25;  // N = size of recursion matrix to use.

    // Use different coeffs depending upon whether m is even or odd.
    if (m % 2 == 0) {
      // Even

      // Get coeff vector for even ce
      double *AA = (double *) calloc(N, sizeof(double));
      if (AA == NULL) return SF_ERROR_MEMORY;
      retcode = mathieu_coeffs_ee(N,q,m, AA);
      if (retcode != SF_ERROR_OK) {
	free(AA);
	return retcode;
      }
      
      // Local scope variables used in summing the Fourier series.
      double tt, td, cep, cem, cedp, cedm;
      cem = 0.0d; cep = 0.0d; cedm = 0.0d; cedp = 0.0d;
      
      // Sum from smallest to largest coeff.
      for (int k=(N-1); k>=0 ; k--) {
	tt = AA[k]*cos(2.0d*k*v); // Term for Mathieu ce
	if (tt<0) {
	  cem = cem + tt;  // Neg running sum
	} else {
	  cep = cep + tt;  // Pos running sum
	}

	td = -2.0d*k*AA[k]*sin(2.0d*k*v); // Term for deriv
	if (td<0) {
	  cedm = cedm + td;
	} else {
	  cedp = cedp + td;  
	}
      } // for (k=(N-1) ...

      // I should do a sort before doing the sum
      *ce = cep+cem;
      *ced = cedp+cedm;

      // Hack -- this makes sure the fcn has the right overall sign for q<0.
      // Someday combine this with the above sums into the same for loop.
      double s = 0.0d;
      for (int l = 0; l<N; l++) {
	s = s + AA[l];
      }
      *ce = SIGN(s)*(*ce);
      *ced = SIGN(s)*(*ced);
      free(AA);
      
    } else {
      // Odd

      // Get coeff vector for odd ce
      double *AA = (double *) calloc(N, sizeof(double));
      if (AA == NULL) return SF_ERROR_MEMORY;
      retcode = mathieu_coeffs_eo(N,q,m, AA);
      if (retcode != SF_ERROR_OK) {
	free(AA);
	return retcode;
      }
      
      // Local scope variables used in summing the Fourier series.
      double tt, td, cep, cem, cedp, cedm;
      cem = 0.0d; cep = 0.0d; cedm = 0.0d; cedp = 0.0d;
      
      // Perform Fourier sum on k = 0, 2, 4, ...
      for (int k=(N-1); k>=0 ; k--) {
	tt = AA[k]*cos((2.0d*k+1.0d)*v);  // Term for Mathieu ce
	if (tt<0) {
	  cem = cem + tt;  // Neg running sum
	} else {
	  cep = cep + tt;  // Pos running sum
	}

	td = -(2.0d*k+1.0d)*AA[k]*sin((2.0d*k+1.0d)*v); // Deriv.
	if (td<0) {
	  cedm = cedm + td;
	} else {
	  cedp = cedp + td;  
	}
      } // for (k=(N-1) ...

      // I should do a sort before doing the sum
      *ce = cep+cem;
      *ced = cedp+cedm;

      // Hack -- this makes sure the fcn has the right overall sign for q<0.
      // Someday combine this with the above sums into the same for loop.
      double s = 0.0d;
      for (int l = 0; l<N; l++) {
	s = s + AA[l];
      }
      *ce = SIGN(s)*(*ce);
      *ced = SIGN(s)*(*ced);

      free(AA);
    } // if (m % 2 == 0)

    return retcode;
  }


  //==================================================================
  int mathieu_se(int m, double q, double v, double *se, double *sed) {
    // This computes the Mathieu fcn se
    // Inputs:
    // m = Mathieu fcn order (scalar)
    // q = frequency parameter (scalar)
    // v = angle in radians (scalar)
    // Outputs:
    // se = value of fcn for these inputs (scalar)
    // sed = value of fcn deriv w.r.t. v for these inputs (scalar)
    // Return code:  
    // Success = 0

    int retcode = SF_ERROR_OK;

    // Check input domain and flag any problems.
    if (m>500) {
      // Don't support absurdly larger orders for now.
      *se = std::numeric_limits<double>::quiet_NaN();
      *sed = std::numeric_limits<double>::quiet_NaN();
      return SF_ERROR_DOMAIN;
    }
    // q>1000 leads to inaccuracy.
    if (abs(q)>1.0e3d) retcode = SF_ERROR_LOSS;
    
    // I find the peak Fourier coeff tracks m.  Therefore,
    // adjust the matrix size based on order m.  Later make this
    // a fcn of q also since the distribution of coeff mags allegedly
    // flattens out for large q.
    int N = m+25;  // N = size of recursion matrix to use.

    // Use different coeffs depending upon whether m is even or odd.
    if (m % 2 == 0) {
      // Even

      // Get coeff vector for even se
      double *BB = (double *) calloc(N, sizeof(double));
      if (BB == NULL) return SF_ERROR_MEMORY;
      retcode = mathieu_coeffs_oe(N,q,m, BB);
      if (retcode != SF_ERROR_OK) {
	free(BB);
	return retcode;
      }
      
      // Local scope variables used in summing the Fourier series.
      double tt, td, sep, sem, sedp, sedm;
      sem = 0.0d; sep = 0.0d; sedm = 0.0d; sedp = 0.0d;
      
      // Sum from smallest to largest coeff.
      for (int k=N; k>=1 ; k--) {
	tt = BB[k-1]*sin(2.0d*k*v); // Mathieu se term
	if (tt<0) {
	  sem = sem + tt;  // Neg running sum
	} else {
	  sep = sep + tt;  // Pos running sum
	}

	td = 2.0d*k*BB[k-1]*cos(2.0d*k*v); // Deriv term.
	if (td<0) {
	  sedm = sedm + td;
	} else {
	  sedp = sedp + td;  
	}
      } // for (k=(N-1) ...

      // I should do a sort before doing the sum
      *se = sep+sem;
      *sed = sedp+sedm;

      // Hack -- this makes sure the fcn has the right overall sign for q<0.
      // Someday combine this with the above sums into the same for loop.
      double s = 0.0d;
      for (int l = 0; l<N; l++) {
	s = s + BB[l];
      }
      *se = SIGN(s)*(*se);
      *sed = SIGN(s)*(*sed);
      
      free(BB);
    } else {
      // Odd

      // Get coeff vector for odd se
      double *BB = (double *) calloc(N, sizeof(double));
      if (BB == NULL) return SF_ERROR_MEMORY;
      retcode = mathieu_coeffs_oo(N,q,m, BB);
      if (retcode != SF_ERROR_OK) {
	free(BB);
	return retcode;
      }
      
      // Local scope variables used in summing the Fourier series.
      double tt, td, sep, sem, sedp, sedm;
      sem = 0.0d; sep = 0.0d; sedm = 0.0d; sedp = 0.0d;
      
      // Sum from smallest to largest coeff.
      for (int k=(N-1); k>=0 ; k--) {
	tt = BB[k]*sin((2.0d*k+1.0d)*v);  // Mathieu se term
	if (tt<0) {
	  sem = sem + tt;  // Neg running sum
	} else {
	  sep = sep + tt;  // Pos running sum
	}

	td = (2.0d*k+1.0d)*BB[k]*cos((2.0d*k+1.0d)*v);  // Deriv term.
	if (td<0) {
	  sedm = sedm + td;
	} else {
	  sedp = sedp + td;  
	}
      } // for (k=(N-1) ...

      // I should do a sort before doing the sum
      *se = sep+sem;
      *sed = sedp+sedm;

      // Hack -- this makes sure the fcn has the right overall sign for q<0.
      // Someday combine this with the above sums into the same for loop.
      double s = 0.0d;
      for (int l = 0; l<N; l++) {
	s = s + BB[l];
      }
      *se = SIGN(s)*(*se);
      *sed = SIGN(s)*(*sed);

      free(BB);
    }

    return retcode;
  } // int mathieu_se


  //==================================================================
  int mathieu_modmc1(int m, double q, double u, double *mc1, double *mc1d) {
    // This computes the Mathieu fcn modmc1
    // Inputs:
    // m = Mathieu fcn order (scalar)
    // q = frequency parameter (scalar)
    // u = radial coord (scalar)
    // Outputs:
    // mc1 = value of fcn for these inputs (scalar)
    // mc1d = value of fcn deriv w.r.t. u for these inputs (scalar)
    // Return code:  
    // Success = 0

    int retcode = SF_ERROR_OK;
    int c;        // Offset used in adaptive computation.

    // Check input domain and flag any problems
    if (m>500) {
      // Don't support absurdly larger orders for now.
      *mc1 = std::numeric_limits<double>::quiet_NaN();
      *mc1d = std::numeric_limits<double>::quiet_NaN();
      return SF_ERROR_DOMAIN;
    }
    if (q<0) return SF_ERROR_DOMAIN; // q<0 is unimplemented
    if (abs(q)>1.0e3d) retcode = SF_ERROR_LOSS;  // q>1000 is inaccurate
    if (m>15 && q>0.1d) retcode = SF_ERROR_LOSS;

    // I find the peak Fourier coeff tracks m.  Therefore,
    // adjust the matrix size based on order m.  Later make this
    // a fcn of q also since the distribution of coeff mags allegedly
    // flattens out for large q.
    int N = m+25;  // N = size of recursion matrix to use.

    // Utility vars.
    double sqq = sqrt(q);
    double exppu = exp(u);
    double expmu = exp(-u);
    double s = sqq*expmu;
    double t = sqq*exppu;
      
    // This is used for the adaptive computation.
    // I set the offset used in Bessel sum depending upon order m.
    // The offset depends upon exactly where in the [q,m] plane lives
    // the input args.
    // The idea comes from the book "Accurate Computation of Mathieu Functions",
    // Malcolm M. Bibby & Andrew F. Peterson.  Also used in the paper
    // "Accurate calculation of the modified Mathieu functions of
    // integer order", Van Buren & Boisvert.
    if ( (m>5 && q<.001) || 
	 (m>7 && q<.01) || 
	 (m>10 && q<.1) || 
	 (m>15 && q<1)  || 
	 (m>20 && q<10) || 
	 (m>30 && q<100) ) {
      c = m/2;
	} else {
      c = 0;
    }
    
    // Use different coeffs depending upon whether m is even or odd.
    if (m % 2 == 0) {
      // Even

      // Get coeff vector for even modmc1
      double *AA = (double *) calloc(N, sizeof(double));
      if (AA == NULL) return SF_ERROR_MEMORY;
      retcode = mathieu_coeffs_ee(N, q, m, AA);
      if (retcode != SF_ERROR_OK) {
	free(AA);
	return retcode;
      }
      
      // Local scope variables used in summing the Fourier series.
      // These are Float128 since some of the terms are near
      // equal amplitude, but different sign, and I want to
      // avoid catastrophic cancellation.
      _Float128 mc1p, mc1m, mc1dp, mc1dm;
      mc1p = 0.0; mc1m = 0.0; mc1dp = 0.0; mc1dm = 0.0;
      
      // Sum from smallest to largest coeff.
      for (int k=(N-1); k>=0 ; k--) {
	if (c==0) {
	  // Non-adaptive calc
	  double Jks = besselj(k,s);
	  double Jkt = besselj(k,t);
	  double Jdks = besseljd(k,s);
	  double Jdkt = besseljd(k,t);

	  _Float128 tt = AA[k]*(Jks*Jkt);
	  _Float128 ttd = AA[k]*(exppu*Jks*Jdkt - expmu*Jdks*Jkt);

	  // Even terms have + sign, odd terms have - sign
	  int sgn = (k%2 == 0) ? 1 : -1;
	  
	  // Do sum using separate sums for + and -
	  tt = sgn*tt;
	  if (tt<0) {
	    // Neg terms
	    mc1m = mc1m + tt;
	  } else {
	    // Pos terms
	    mc1p = mc1p + tt;
	  }
        
	  // Do sum using separate sums for + and -
	  ttd = sgn*ttd;
	  if (ttd<0) {
	    // Neg terms
	    mc1dm = mc1dm + ttd;
	  } else {
	    // Pos terms
	    mc1dp = mc1dp + ttd;
	  }

	} else {
	  // Adaptive calc
	  double Jkmcs = besselj(k-c,s);
	  double Jkpcs = besselj(k+c,s);
	  double Jkpct = besselj(k+c,t);
	  double Jkmct = besselj(k-c,t);

	  double Jdkmcs = besseljd(k-c,s);
	  double Jdkpcs = besseljd(k+c,s);
	  double Jdkpct = besseljd(k+c,t);
	  double Jdkmct = besseljd(k-c,t);

	  _Float128 tt = AA[k]*(Jkmcs*Jkpct + Jkpcs*Jkmct);
	  _Float128 ttd = AA[k]* 
	    (exppu*(Jkmcs*Jdkpct + Jkpcs*Jdkmct) -
	     expmu*(Jdkmcs*Jkpct + Jdkpcs*Jkmct)) ;

	  // Even terms have + sign, odd terms have - sign
	  int sgn = (k%2 == 0) ? 1 : -1;

	  // Do sum using separate sums for + and -
	  tt = sgn*tt;
	  if (tt<0) {
	    // Neg terms
	    mc1m = mc1m + tt;
	  } else {
	    // Pos terms
	    mc1p = mc1p + tt;
	  }
        
	  // Do sum using separate sums for + and -
	  ttd = sgn*ttd;
	  if (ttd<0) {
	    // Neg terms
	    mc1dm = mc1dm + ttd;
	  } else {
	    // Pos terms
	    mc1dp = mc1dp + ttd;
	  }
	  
	} // if (c==0)

      } // for (int k=(N-1) ...

      // Sum pos and neg terms to get final result
      *mc1 = static_cast<double>(mc1p+mc1m);
      *mc1d = static_cast<double>(mc1dp+mc1dm);

      // Do normalization.  Note normalization depends upon c.
      int sgn = m/2;
      if (sgn%2 == 0) {
	*mc1 = (*mc1)/AA[c];    
	*mc1d = sqq*(*mc1d)/AA[c];
      } else {
	*mc1 = -(*mc1)/AA[c];    
	*mc1d = -sqq*(*mc1d)/AA[c];
      }

      free(AA);
      
    } else {
      // Odd -- m = 1, 3, 5, 7 ...

      // Get coeff vector for even modmc1
      double *AA = (double *) calloc(N, sizeof(double));
      if (AA == NULL) return SF_ERROR_MEMORY;
      retcode = mathieu_coeffs_eo(N, q, m, AA);
      if (retcode != SF_ERROR_OK) {
	free(AA);
	return retcode;
      }
      
      // Variables used in summing the Fourier series.
      _Float128 mc1p, mc1m, mc1dp, mc1dm;
      mc1p = 0.0; mc1m = 0.0; mc1dp = 0.0; mc1dm = 0.0;
      
      // Sum from smallest to largest coeff.
      for (int k=(N-1); k>=0 ; k--) {
	if (c==0) {
	  // Non-adaptive calc
	  double Jks = besselj(k,s);
	  double Jkp1s = besselj(k+1,s);
	  double Jkt = besselj(k,t);
	  double Jkp1t = besselj(k+1,t);

	  double Jdks = besseljd(k,s);
	  double Jdkp1s = besseljd(k+1,s);
	  double Jdkt = besseljd(k,t);
	  double Jdkp1t = besseljd(k+1,t);
	
	  _Float128 tt = AA[k]*(Jks*Jkp1t + Jkp1s*Jkt);
	  _Float128 ttd = AA[k]*
	      (exppu*(Jks*Jdkp1t + Jkp1s*Jdkt) - 
	       expmu*(Jdks*Jkp1t + Jdkp1s*Jkt) ) ;
	  
	  int sgn = (k%2 == 0) ? 1 : -1;

	  // Do sum using separate sums for + and -
	  tt = sgn*tt;
	  if (tt<0) {
	    // Neg terms
	    mc1m = mc1m + tt;
	  } else {
	    // Pos terms
	    mc1p = mc1p + tt;
	  }
        
	  // Do sum using separate sums for + and -
	  ttd = sgn*ttd;
	  if (ttd<0) {
	    // Neg terms
	    mc1dm = mc1dm + ttd;
	  } else {
	    // Pos terms
	    mc1dp = mc1dp + ttd;
	  }
	  
	} else {
	  // Adaptive calc
	  double Jkmcs = besselj(k-c,s);
	  double Jkpcs = besselj(k+c+1,s);
	  double Jkpct = besselj(k+c+1,t);
	  double Jkmct = besselj(k-c,t);
	  
	  double Jdkmcs = besseljd(k-c,s);
	  double Jdkpcs = besseljd(k+c+1,s);
	  double Jdkpct = besseljd(k+c+1,t);
	  double Jdkmct = besseljd(k-c,t);

	  _Float128 tt = AA[k]*(Jkmcs*Jkpct + Jkpcs*Jkmct);
	  _Float128 ttd = AA[k]* 
	      (exppu*(Jkmcs*Jdkpct + Jkpcs*Jdkmct) -
	       expmu*(Jdkmcs*Jkpct + Jdkpcs*Jkmct)) ;

	  int sgn = (k%2 == 0) ? 1 : -1;

	  // Do sum using separate sums for + and -
	  tt = sgn*tt;
	  if (tt<0) {
	    // Neg terms
	    mc1m = mc1m + tt;
	  } else {
	    // Pos terms
	    mc1p = mc1p + tt;
	  }
        
	  // Do sum using separate sums for + and -
	  ttd = sgn*ttd;
	  if (ttd<0) {
	    // Neg terms
	    mc1dm = mc1dm + ttd;
	  } else {
	    // Pos terms
	    mc1dp = mc1dp + ttd;
	  }

	} // if (c==0)

      } // for (int k=(N-1) ...

      // Sum pos and neg terms to get final answer
      *mc1 = static_cast<double>(mc1p+mc1m);
      *mc1d = static_cast<double>(mc1dp+mc1dm);
      
      // Do normalization.  Note normalization depends upon c.
      int sgn = (m-1)/2;
      if (sgn%2 == 0) {
	*mc1 = (*mc1)/AA[c];    
	*mc1d = sqq*(*mc1d)/AA[c];
      } else {
	*mc1 = -(*mc1)/AA[c];    
	*mc1d = -sqq*(*mc1d)/AA[c];
      }

      free(AA);
    }

    return retcode;
  } // int mathieu_modmc1


  //==================================================================
  int mathieu_modms1(int m, double q, double u, double *ms1, double *ms1d) {
    // This computes the Mathieu fcn modms1
    // Inputs:
    // m = Mathieu fcn order (scalar)
    // q = frequency parameter (scalar)
    // u = radial coord (scalar)
    // Outputs:
    // ms1 = value of fcn for these inputs (scalar)
    // ms1d = value of fcn deriv w.r.t. u for these inputs (scalar)
    // Return code:  
    // Success = 0

    int retcode = SF_ERROR_OK;
    int c;        // Offset used in adaptive computation.
    
    // Check input domain and flag any problems
    if (m>500) {
      // Don't support absurdly larger orders for now.
      *ms1 = std::numeric_limits<double>::quiet_NaN();
      *ms1d = std::numeric_limits<double>::quiet_NaN();
      return SF_ERROR_DOMAIN;
    }
    if (q<0) return SF_ERROR_DOMAIN;  // q < 0 is currently unimplemented
    if (abs(q)>1.0e3d) retcode = SF_ERROR_LOSS;
    if (m>15 && q>0.1d) retcode = SF_ERROR_LOSS;

    // I find the peak Fourier coeff tracks m.  Therefore,
    // adjust the matrix size based on order m.  Later make this
    // a fcn of q also since the distribution of coeff mags allegedly
    // flattens out for large q.
    int N = m+25;  // N = size of recursion matrix to use.

    // Utility vars.
    double sqq = sqrt(q);
    double exppu = exp(u);
    double expmu = exp(-u);
    double s = sqq*expmu;
    double t = sqq*exppu;
      
    // This is used for the adaptive computation.
    // I set the offset used in Bessel fcn depending upon order m.
    // The idea comes from the book "Accurate Computation of Mathieu Functions",
    // Malcolm M. Bibby & Andrew F. Peterson.  Also used in the paper
    // "Accurate calculation of the modified Mathieu functions of
    // integer order", Van Buren & Boisvert.
    if ( (m>5 && q<.001) || 
	 (m>7 && q<.01) || 
	 (m>10 && q<.1) || 
	 (m>15 && q<1)  || 
	 (m>20 && q<10) || 
	 (m>30 && q<100) ) {
      c = m/2;
	} else {
      c = 0;
    }
    
    // Use different coeffs depending upon whether m is even or odd.
    if (m % 2 == 0) {
      // Even

      // Get coeff vector for even modms1
      double *BB = (double *) calloc(N, sizeof(double));
      if (BB == NULL) return SF_ERROR_MEMORY;
      retcode = mathieu_coeffs_oe(N, q, m, BB);
      if (retcode != SF_ERROR_OK) {
	free(BB);
	return retcode;
      }
      
      // Variables used in summing the Fourier series.
      // These are Float128 since some of the terms are near
      // equal amplitude, but different sign.
      _Float128 ms1p, ms1m, ms1dp, ms1dm;
      ms1p = 0.0; ms1m = 0.0; ms1dp = 0.0; ms1dm = 0.0;
      
      // Sum from smallest to largest coeff.
      for (int k=(N-1); k>=0 ; k--) {
	if (c==0) {
	  // Non-adaptive calc
	  double Jks = besselj(k,s);
	  double Jkp2t = besselj(k+2,t);
	  double Jkp2s = besselj(k+2,s);
	  double Jkt = besselj(k,t);
	  
	  double Jdks = besseljd(k,s);
	  double Jdkp2t = besseljd(k+2,t);
	  double Jdkp2s = besseljd(k+2,s);
	  double Jdkt = besseljd(k,t);

	  _Float128 tt = BB[k]*(Jks*Jkp2t - Jkp2s*Jkt);
	  _Float128 ttd = BB[k]*
              (exppu*(Jks*Jdkp2t - Jkp2s*Jdkt) -
               expmu*(Jdks*Jkp2t - Jdkp2s*Jkt));

	  // Even terms have + sign, odd terms have - sign
	  int sgn = (k%2 == 0) ? 1 : -1;
	  
	  // Do sum using separate sums for + and -
	  tt = sgn*tt;
	  if (tt<0) {
	    // Neg terms
	    ms1m = ms1m + tt;
	  } else {
	    // Pos terms
	    ms1p = ms1p + tt;
	  }
        
	  // Do sum using separate sums for + and -
	  ttd = sgn*ttd;
	  if (ttd<0) {
	    // Neg terms
	    ms1dm = ms1dm + ttd;
	  } else {
	    // Pos terms
	    ms1dp = ms1dp + ttd;
	  }

	} else {
	  // Adaptive calc
	  double Jkmcs = besselj(k-c,s);
	  double Jkpct = besselj(k+c+2,t);
	  double Jkpcs = besselj(k+c+2,s);
	  double Jkmct = besselj(k-c,t);
	  
	  double Jdkmcs = besseljd(k-c,s);
	  double Jdkpct = besseljd(k+c+2,t);
	  double Jdkpcs = besseljd(k+c+2,s);
	  double Jdkmct = besseljd(k-c,t);

	  _Float128 tt = BB[k]*(Jkmcs*Jkpct - Jkpcs*Jkmct);
	  _Float128 ttd = BB[k]*
            (exppu*(Jkmcs*Jdkpct - Jkpcs*Jdkmct) -
            expmu*(Jdkmcs*Jkpct - Jdkpcs*Jkmct));

	  // Even terms have + sign, odd terms have - sign
	  int sgn = (k%2 == 0) ? 1 : -1;

	  // Do sum using separate sums for + and -
	  tt = sgn*tt;
	  if (tt<0) {
	    // Neg terms
	    ms1m = ms1m + tt;
	  } else {
	    // Pos terms
	    ms1p = ms1p + tt;
	  }
        
	  // Do sum using separate sums for + and -
	  ttd = sgn*ttd;
	  if (ttd<0) {
	    // Neg terms
	    ms1dm = ms1dm + ttd;
	  } else {
	    // Pos terms
	    ms1dp = ms1dp + ttd;
	  }
	  
	} // if (c==0)

      } // for (int k=(N-1) ...

      // Sum pos and neg terms to get final result
      *ms1 = static_cast<double>(ms1p+ms1m);
      *ms1d = static_cast<double>(ms1dp+ms1dm);

      // Do normalization.  Note normalization depends upon c.
      int sgn = (m-2)/2;
      if (sgn%2 == 0) {
	*ms1 = (*ms1)/BB[c];    
	*ms1d = sqq*(*ms1d)/BB[c];
      } else {
	*ms1 = -(*ms1)/BB[c];    
	*ms1d = -sqq*(*ms1d)/BB[c];
      }

      free(BB);

    } else {
      // Odd -- m = 1, 3, 5, 7 ...

      // Get coeff vector for even modms1
      double *BB = (double *) calloc(N, sizeof(double));
      if (BB == NULL) return SF_ERROR_MEMORY;
      retcode = mathieu_coeffs_oo(N, q, m, BB);
      if (retcode != SF_ERROR_OK) {
	free(BB);
	return retcode;
      }
      
      // Variables used in summing the Fourier series.
      _Float128 ms1p, ms1m, ms1dp, ms1dm;
      ms1p = 0.0; ms1m = 0.0; ms1dp = 0.0; ms1dm = 0.0;
      
      // Sum from smallest to largest coeff.
      for (int k=(N-1); k>=0 ; k--) {
	if (c==0) {
	  // Non-adaptive calc
	  double Jks = besselj(k,s);
	  double Jkt = besselj(k,t);
	  double Jkp1s = besselj(k+1,s);
	  double Jkp1t = besselj(k+1,t);
	  
	  double Jdks = besseljd(k,s);
	  double Jdkt = besseljd(k,t);
	  double Jdkp1s = besseljd(k+1,s);
	  double Jdkp1t = besseljd(k+1,t);
	
	  _Float128 tt = BB[k]*(Jks*Jkp1t - Jkp1s*Jkt);
	  _Float128 ttd = BB[k]*
	    (exppu*(Jks*Jdkp1t - Jkp1s*Jdkt) - 
	     expmu*(Jdks*Jkp1t - Jdkp1s*Jkt));

	  int sgn = (k%2 == 0) ? 1 : -1;

	  // Do sum using separate sums for + and -
	  tt = sgn*tt;
	  if (tt<0) {
	    // Neg terms
	    ms1m = ms1m + tt;
	  } else {
	    // Pos terms
	    ms1p = ms1p + tt;
	  }
        
	  // Do sum using separate sums for + and -
	  ttd = sgn*ttd;
	  if (ttd<0) {
	    // Neg terms
	    ms1dm = ms1dm + ttd;
	  } else {
	    // Pos terms
	    ms1dp = ms1dp + ttd;
	  }
	  
	} else {
	  // Adaptive calc
	  double Jkmcs = besselj(k-c,s);
	  double Jkpcs = besselj(k+c+1,s);
	  double Jkpct = besselj(k+c+1,t);
	  double Jkmct = besselj(k-c,t);
	  
	  double Jdkmcs = besseljd(k-c,s);
	  double Jdkpcs = besseljd(k+c+1,s);
	  double Jdkpct = besseljd(k+c+1,t);
	  double Jdkmct = besseljd(k-c,t);

	  _Float128 tt = BB[k]*(Jkmcs*Jkpct - Jkpcs*Jkmct);
	  _Float128 ttd = BB[k]* 
	      (exppu*(Jkmcs*Jdkpct - Jkpcs*Jdkmct) -
	       expmu*(Jdkmcs*Jkpct - Jdkpcs*Jkmct)) ;

	  int sgn = (k%2 == 0) ? 1 : -1;

	  // Do sum using separate sums for + and -
	  tt = sgn*tt;
	  if (tt<0) {
	    // Neg terms
	    ms1m = ms1m + tt;
	  } else {
	    // Pos terms
	    ms1p = ms1p + tt;
	  }
        
	  // Do sum using separate sums for + and -
	  ttd = sgn*ttd;
	  if (ttd<0) {
	    // Neg terms
	    ms1dm = ms1dm + ttd;
	  } else {
	    // Pos terms
	    ms1dp = ms1dp + ttd;
	  }

	} // if (c==0)

      } // for (int k=(N-1) ...

      // Sum pos and neg terms to get final answer
      *ms1 = static_cast<double>(ms1p+ms1m);
      *ms1d = static_cast<double>(ms1dp+ms1dm);
      
      // Do normalization.  Note normalization depends upon c.
      int sgn = (m-1)/2;
      if (sgn%2 == 0) {
	*ms1 = (*ms1)/BB[c];    
	*ms1d = sqq*(*ms1d)/BB[c];
      } else {
	*ms1 = -(*ms1)/BB[c];    
	*ms1d = -sqq*(*ms1d)/BB[c];
      }
      
      free(BB);
    }

    return retcode;
  } // int mathieu_modms1


  //==================================================================
  int mathieu_modmc2(int m, double q, double u, double *mc2, double *mc2d) {
    // This computes the Mathieu fcn modmc2
    // Inputs:
    // m = Mathieu fcn order (scalar)
    // q = frequency parameter (scalar)
    // u = radial coord (scalar)
    // Outputs:
    // mc2 = value of fcn for these inputs (scalar)
    // mc2d = value of fcn deriv w.r.t. u for these inputs (scalar)
    // Return code:  
    // Success = 0

    int retcode = SF_ERROR_OK;
    int c;        // Offset used in adaptive computation.
    
    // Check input domain and flag any problems
    if (m>500) {
      // Don't support absurdly larger orders for now.
      *mc2 = std::numeric_limits<double>::quiet_NaN();
      *mc2d = std::numeric_limits<double>::quiet_NaN();
      return SF_ERROR_DOMAIN;
    }
    if (q<0) return SF_ERROR_DOMAIN;  // q<0 is currently unimplemented
    if (abs(q)>1.0e3d) retcode = SF_ERROR_LOSS;
    if (m>15 && q>0.1d) retcode = SF_ERROR_LOSS;

    // I find the peak Fourier coeff tracks m.  Therefore,
    // adjust the matrix size based on order m.  Later make this
    // a fcn of q also since the distribution of coeff mags allegedly
    // flattens out for large q.
    int N = m+25;  // N = size of recursion matrix to use.

    // Utility vars.
    double sqq = sqrt(q);
    double exppu = exp(u);
    double expmu = exp(-u);
    double s = sqq*expmu;
    double t = sqq*exppu;
      
    // This is used for the adaptive computation.
    // I set the offset used in Bessel fcn depending upon order m.
    // The idea comes from the book "Accurate Computation of Mathieu Functions",
    // Malcolm M. Bibby & Andrew F. Peterson.  Also used in the paper
    // "Accurate calculation of the modified Mathieu functions of
    // integer order", Van Buren & Boisvert.
    if ( (m>5 && q<.001) || 
	 (m>7 && q<.01) || 
	 (m>10 && q<.1) || 
	 (m>15 && q<1)  || 
	 (m>20 && q<10) || 
	 (m>30 && q<100) ) {
      c = m/2;
	} else {
      c = 0;
    }
    
    // Use different coeffs depending upon whether m is even or odd.
    if (m % 2 == 0) {
      // Even

      // Get coeff vector for even modmc2
      double *AA = (double *) calloc(N, sizeof(double));
      if (AA == NULL) return SF_ERROR_MEMORY;
      retcode = mathieu_coeffs_ee(N, q, m, AA);
      if (retcode != SF_ERROR_OK) {
	free(AA);
	return retcode;
      }
      
      // Variables used in summing the Fourier series.
      // These are Float128 since some of the terms are near
      // equal amplitude, but different sign.
      _Float128 mc2p, mc2m, mc2dp, mc2dm;
      
      // Sum from smallest to largest coeff.
      mc2p = 0.0; mc2m = 0.0; mc2dp = 0.0; mc2dm = 0.0;
      for (int k=(N-1); k>=0 ; k--) {
	if (c==0) {
	  // Non-adaptive calc
	  double Jks = besselj(k,s);
	  double Ykt = bessely(k,t);
	  double Jdks = besseljd(k,s);
	  double Ydkt = besselyd(k,t);

	  _Float128 tt = AA[k]*Jks*Ykt ;
	  _Float128 ttd = AA[k]*(exppu*Jks*Ydkt - expmu*Jdks*Ykt) ;

	  // Even terms have + sign, odd terms have - sign
	  int sgn = (k%2 == 0) ? 1 : -1;
	  
	  // Do sum using separate sums for + and -
	  tt = sgn*tt;
	  if (tt<0) {
	    // Neg terms
	    mc2m = mc2m + tt;
	  } else {
	    // Pos terms
	    mc2p = mc2p + tt;
	  }
        
	  // Do sum using separate sums for + and -
	  ttd = sgn*ttd;
	  if (ttd<0) {
	    // Neg terms
	    mc2dm = mc2dm + ttd;
	  } else {
	    // Pos terms
	    mc2dp = mc2dp + ttd;
	  }

	} else {
	  // Adaptive calc
	  double Jkmcs = besselj(k-c,s);
	  double Jkpcs = besselj(k+c,s);
	  double Ykpct = bessely(k+c,t);
	  double Ykmct = bessely(k-c,t);
	  
	  double Jdkmcs = besseljd(k-c,s);
	  double Jdkpcs = besseljd(k+c,s);
	  double Ydkpct = besselyd(k+c,t);
	  double Ydkmct = besselyd(k-c,t);

	  _Float128 tt = AA[k]*(Jkmcs*Ykpct + Jkpcs*Ykmct);
	  _Float128 ttd = AA[k]*
	    (exppu*(Jkmcs*Ydkpct + Jkpcs*Ydkmct) -
             expmu*(Jdkmcs*Ykpct + Jdkpcs*Ykmct)) ;

	  // Even terms have + sign, odd terms have - sign
	  int sgn = (k%2 == 0) ? 1 : -1;

	  // Do sum using separate sums for + and -
	  tt = sgn*tt;
	  if (tt<0) {
	    // Neg terms
	    mc2m = mc2m + tt;
	  } else {
	    // Pos terms
	    mc2p = mc2p + tt;
	  }
        
	  // Do sum using separate sums for + and -
	  ttd = sgn*ttd;
	  if (ttd<0) {
	    // Neg terms
	    mc2dm = mc2dm + ttd;
	  } else {
	    // Pos terms
	    mc2dp = mc2dp + ttd;
	  }
	  
	} // if (c==0)

      } // for (int k=(N-1) ...

      // Sum pos and neg terms to get final result
      *mc2 = static_cast<double>(mc2p+mc2m);
      *mc2d = static_cast<double>(mc2dp+mc2dm);

      // Do normalization.  Note normalization depends upon c.
      int sgn = m/2;
      if (sgn%2 == 0) {
	*mc2 = (*mc2)/AA[c];    
	*mc2d = sqq*(*mc2d)/AA[c];
      } else {
	*mc2 = -(*mc2)/AA[c];    
	*mc2d = -sqq*(*mc2d)/AA[c];
      }

      free(AA);
      
    } else {
      // Odd -- m = 1, 3, 5, 7 ...

      // Get coeff vector for odd mc2
      double *AA = (double *) calloc(N, sizeof(double));
      if (AA == NULL) return SF_ERROR_MEMORY;
      retcode = mathieu_coeffs_eo(N, q, m, AA);
      if (retcode != SF_ERROR_OK) {
	free(AA);
	return retcode;
      }
      
      // Variables used in summing the Fourier series.
      _Float128 mc2p, mc2m, mc2dp, mc2dm;
      mc2p = 0.0; mc2m = 0.0; mc2dp = 0.0; mc2dm = 0.0;
      
      // Sum from smallest to largest coeff.
      for (int k=(N-1); k>=0 ; k--) {
	if (c==0) {
	  // Non-adaptive calc
	  double Jks = besselj(k,s);
	  double Ykt = bessely(k,t);
	  double Jkp1s = besselj(k+1,s);
	  double Ykp1t = bessely(k+1,t);
	  
	  double Jdks = besseljd(k,s);
	  double Ydkt = besselyd(k,t);
	  double Jdkp1s = besseljd(k+1,s);
	  double Ydkp1t = besselyd(k+1,t);
	
	  _Float128 tt = AA[k]*(Jks*Ykp1t + Jkp1s*Ykt);
	  _Float128 ttd = AA[k]*
	    (exppu*(Jks*Ydkp1t + Jkp1s*Ydkt) -
	     expmu*(Jdks*Ykp1t + Jdkp1s*Ykt) ) ;

	  int sgn = (k%2 == 0) ? 1 : -1;

	  // Do sum using separate sums for + and -
	  tt = sgn*tt;
	  if (tt<0) {
	    // Neg terms
	    mc2m = mc2m + tt;
	  } else {
	    // Pos terms
	    mc2p = mc2p + tt;
	  }
        
	  // Do sum using separate sums for + and -
	  ttd = sgn*ttd;
	  if (ttd<0) {
	    // Neg terms
	    mc2dm = mc2dm + ttd;
	  } else {
	    // Pos terms
	    mc2dp = mc2dp + ttd;
	  }
	  
	} else {
	  // Adaptive calc
	  double Jkmcs = besselj(k-c,s);
	  double Jkpcs = besselj(k+c+1,s);
	  double Ykpct = bessely(k+c+1,t);
	  double Ykmct = bessely(k-c,t);
	  
	  double Jdkmcs = besseljd(k-c,s);
	  double Jdkpcs = besseljd(k+c+1,s);
	  double Ydkpct = besselyd(k+c+1,t);
	  double Ydkmct = besselyd(k-c,t);

	  _Float128 tt = AA[k]*(Jkmcs*Ykpct + Jkpcs*Ykmct);
	  _Float128 ttd = AA[k]*
	    (exppu*(Jkmcs*Ydkpct + Jkpcs*Ydkmct) -
	     expmu*(Jdkmcs*Ykpct + Jdkpcs*Ykmct) ) ;

	  int sgn = (k%2 == 0) ? 1 : -1;

	  // Do sum using separate sums for + and -
	  tt = sgn*tt;
	  if (tt<0) {
	    // Neg terms
	    mc2m = mc2m + tt;
	  } else {
	    // Pos terms
	    mc2p = mc2p + tt;
	  }
        
	  // Do sum using separate sums for + and -
	  ttd = sgn*ttd;
	  if (ttd<0) {
	    // Neg terms
	    mc2dm = mc2dm + ttd;
	  } else {
	    // Pos terms
	    mc2dp = mc2dp + ttd;
	  }

	} // if (c==0)

      } // for (int k=(N-1) ...

      // Sum pos and neg terms to get final answer
      *mc2 = static_cast<double>(mc2p+mc2m);
      *mc2d = static_cast<double>(mc2dp+mc2dm);
      
      // Do normalization.  Note normalization depends upon c.
      int sgn = (m-1)/2;
      if (sgn%2 == 0) {
	*mc2 = (*mc2)/AA[c];    
	*mc2d = sqq*(*mc2d)/AA[c];
      } else {
	*mc2 = -(*mc2)/AA[c];    
	*mc2d = -sqq*(*mc2d)/AA[c];
      }

      free(AA);
    }

    return retcode;
  } // int mathieu_modmc2


  //==================================================================
  int mathieu_modms2(int m, double q, double u, double *ms2, double *ms2d) {
    // This computes the Mathieu fcn modms2
    // Inputs:
    // m = Mathieu fcn order (scalar)
    // q = frequency parameter (scalar)
    // u = radial coord (scalar)
    // Outputs:
    // ms2 = value of fcn for these inputs (scalar)
    // ms2d = value of fcn deriv w.r.t. u for these inputs (scalar)
    // Return code:  
    // Success = 0

    int retcode = SF_ERROR_OK;
    int c;        // Offset used in adaptive computation.
    
    // Check input domain and flag any problems
    if (m>500) {
      // Don't support absurdly larger orders for now.
      *ms2 = std::numeric_limits<double>::quiet_NaN();
      *ms2d = std::numeric_limits<double>::quiet_NaN();
      return SF_ERROR_DOMAIN;
    }
    if (q<0) return SF_ERROR_DOMAIN;  // q<0 is currently unimplemented
    if (abs(q)>1.0e3d) retcode = SF_ERROR_LOSS;
    if (m>15 && q>0.1d) retcode = SF_ERROR_LOSS;

    // I find the peak Fourier coeff tracks m.  Therefore,
    // adjust the matrix size based on order m.  Later make this
    // a fcn of q also since the distribution of coeff mags allegedly
    // flattens out for large q.
    int N = m+25;  // N = size of recursion matrix to use.

    // Utility vars.
    double sqq = sqrt(q);
    double exppu = exp(u);
    double expmu = exp(-u);
    double s = sqq*expmu;
    double t = sqq*exppu;
      
    // This is used for the adaptive computation.
    // I set the offset used in Bessel fcn depending upon order m.
    // The idea comes from the book "Accurate Computation of Mathieu Functions",
    // Malcolm M. Bibby & Andrew F. Peterson.  Also used in the paper
    // "Accurate calculation of the modified Mathieu functions of
    // integer order", Van Buren & Boisvert.
    if ( (m>5 && q<.001) || 
	 (m>7 && q<.01) || 
	 (m>10 && q<.1) || 
	 (m>15 && q<1)  || 
	 (m>20 && q<10) || 
	 (m>30 && q<100) ) {
      c = m/2;
	} else {
      c = 0;
    }
    c = 0;
    
    // Use different coeffs depending upon whether m is even or odd.
    if (m % 2 == 0) {
      // Even

      // Get coeff vector for even modms2
      double *BB = (double *) calloc(N, sizeof(double));
      if (BB == NULL) return SF_ERROR_MEMORY;
      retcode = mathieu_coeffs_oe(N, q, m, BB);
      if (retcode != SF_ERROR_OK) {
	free(BB);
	return retcode;
      }
      
      // Variables used in summing the Fourier series.
      // These are Float128 since some of the terms are near
      // equal amplitude, but different sign.
      _Float128 ms2p, ms2m, ms2dp, ms2dm;
      ms2p = 0.0; ms2m = 0.0; ms2dp = 0.0; ms2dm = 0.0;
      
      // Sum from smallest to largest coeff.
      for (int k=(N-1); k>=0 ; k--) {
	if (c==0) {
	  // Non-adaptive calc
	  double Jks = besselj(k,s);
	  double Ykp2t = bessely(k+2,t);
	  double Jkp2s = besselj(k+2,s);
	  double Ykt = bessely(k,t);
	  
	  double Jdks = besseljd(k,s);
	  double Ydkp2t = besselyd(k+2,t);
	  double Jdkp2s = besseljd(k+2,s);
	  double Ydkt = besselyd(k,t);

	  _Float128 tt = BB[k]*(Jks*Ykp2t - Jkp2s*Ykt);
	  _Float128 ttd = BB[k]*
	    (exppu*(Jks*Ydkp2t - Jkp2s*Ydkt) - 
	     expmu*(Jdks*Ykp2t - Jdkp2s*Ykt));

	  // Even terms have + sign, odd terms have - sign
	  int sgn = (k%2 == 0) ? 1 : -1;
	  
	  // Do sum using separate sums for + and -
	  tt = sgn*tt;
	  if (tt<0) {
	    // Neg terms
	    ms2m = ms2m + tt;
	  } else {
	    // Pos terms
	    ms2p = ms2p + tt;
	  }
        
	  // Do sum using separate sums for + and -
	  ttd = sgn*ttd;
	  if (ttd<0) {
	    // Neg terms
	    ms2dm = ms2dm + ttd;
	  } else {
	    // Pos terms
	    ms2dp = ms2dp + ttd;
	  }

	} else {
	  // Adaptive calc
	  double Jkmcs = besselj(k-c,s);
	  double Ykpct = bessely(k+c+2,t);
	  double Jkpcs = besselj(k+c+2,s);
	  double Ykmct = bessely(k-c,t);
	  
	  double Jdkmcs = besseljd(k-c,s);
	  double Ydkpct = besselyd(k+c+2,t);
	  double Jdkpcs = besseljd(k+c+2,s);
	  double Ydkmct = besselyd(k-c,t);

	  _Float128 tt = BB[k]*(Jkmcs*Ykpct - Jkpcs*Ykmct);
	  _Float128 ttd = BB[k]*
	    (exppu*(Jkmcs*Ydkpct - Jkpcs*Ydkmct) - 
	     expmu*(Jdkmcs*Ykpct - Jdkpcs*Ykmct));

	  // Even terms have + sign, odd terms have - sign
	  int sgn = (k%2 == 0) ? 1 : -1;

	  // Do sum using separate sums for + and -
	  tt = sgn*tt;
	  if (tt<0) {
	    // Neg terms
	    ms2m = ms2m + tt;
	  } else {
	    // Pos terms
	    ms2p = ms2p + tt;
	  }
        
	  // Do sum using separate sums for + and -
	  ttd = sgn*ttd;
	  if (ttd<0) {
	    // Neg terms
	    ms2dm = ms2dm + ttd;
	  } else {
	    // Pos terms
	    ms2dp = ms2dp + ttd;
	  }
	  
	} // if (c==0)

      } // for (int k=(N-1) ...

      // Sum pos and neg terms to get final result
      *ms2 = static_cast<double>(ms2p+ms2m);
      *ms2d = static_cast<double>(ms2dp+ms2dm);

      // Do normalization.  Note normalization depends upon c.
      int sgn = (m-2)/2;
      if (sgn%2 == 0) {
	*ms2 = (*ms2)/BB[c];    
	*ms2d = sqq*(*ms2d)/BB[c];
      } else {
	*ms2 = -(*ms2)/BB[c];    
	*ms2d = -sqq*(*ms2d)/BB[c];
      }

      free(BB);

    } else {
      // Odd -- m = 1, 3, 5, 7 ...

      // Get coeff vector for even modms2
      double *BB = (double *) calloc(N, sizeof(double));
      if (BB == NULL) return SF_ERROR_MEMORY;
      retcode = mathieu_coeffs_oo(N, q, m, BB);
      if (retcode != SF_ERROR_OK) {
	free(BB);
	return retcode;
      }
      
      // Variables used in summing the Fourier series.
      _Float128 ms2p, ms2m, ms2dp, ms2dm;
      ms2p = 0.0; ms2m = 0.0; ms2dp = 0.0; ms2dm = 0.0;
      
      // Sum from smallest to largest coeff.
      for (int k=(N-1); k>=0 ; k--) {
	if (c==0) {
	  // Non-adaptive calc
	  double Jks = besselj(k,s);
	  double Ykt = bessely(k,t);
	  double Jkp1s = besselj(k+1,s);
	  double Ykp1t = bessely(k+1,t);
	  
	  double Jdks = besseljd(k,s);
	  double Ydkt = besselyd(k,t);
	  double Jdkp1s = besseljd(k+1,s);
	  double Ydkp1t = besselyd(k+1,t);
	
	  _Float128 tt = BB[k]*(Jks*Ykp1t - Jkp1s*Ykt);
	  _Float128 ttd = BB[k]*
	    (exppu*(Jks*Ydkp1t - Jkp1s*Ydkt) - 
	     expmu*(Jdks*Ykp1t - Jdkp1s*Ykt));

	  int sgn = (k%2 == 0) ? 1 : -1;

	  // Do sum using separate sums for + and -
	  tt = sgn*tt;
	  if (tt<0) {
	    // Neg terms
	    ms2m = ms2m + tt;
	  } else {
	    // Pos terms
	    ms2p = ms2p + tt;
	  }
        
	  // Do sum using separate sums for + and -
	  ttd = sgn*ttd;
	  if (ttd<0) {
	    // Neg terms
	    ms2dm = ms2dm + ttd;
	  } else {
	    // Pos terms
	    ms2dp = ms2dp + ttd;
	  }
	  
	} else {
	  // Adaptive calc
	  double Jkmcs = besselj(k-c,s);
	  double Jkpcs = besselj(k+c+1,s);
	  double Ykpct = bessely(k+c+1,t);
	  double Ykmct = bessely(k-c,t);
	  
	  double Jdkmcs = besseljd(k-c,s);
	  double Jdkpcs = besseljd(k+c+1,s);
	  double Ydkpct = besselyd(k+c+1,t);
	  double Ydkmct = besselyd(k-c,t);

	  _Float128 tt = BB[k]*(Jkmcs*Ykpct - Jkpcs*Ykmct);
	  _Float128 ttd = BB[k]*
	    (exppu*(Jkmcs*Ydkpct - Jkpcs*Ydkmct) - 
	     expmu*(Jdkmcs*Ykpct - Jdkpcs*Ykmct) ) ;

	  int sgn = (k%2 == 0) ? 1 : -1;

	  // Do sum using separate sums for + and -
	  tt = sgn*tt;
	  if (tt<0) {
	    // Neg terms
	    ms2m = ms2m + tt;
	  } else {
	    // Pos terms
	    ms2p = ms2p + tt;
	  }
        
	  // Do sum using separate sums for + and -
	  ttd = sgn*ttd;
	  if (ttd<0) {
	    // Neg terms
	    ms2dm = ms2dm + ttd;
	  } else {
	    // Pos terms
	    ms2dp = ms2dp + ttd;
	  }

	} // if (c==0)

      } // for (int k=(N-1) ...

      // Sum pos and neg terms to get final answer
      *ms2 = static_cast<double>(ms2p+ms2m);
      *ms2d = static_cast<double>(ms2dp+ms2dm);
      
      // Do normalization.  Note normalization depends upon c.
      int sgn = (m-1)/2;
      if (sgn%2 == 0) {
	*ms2 = (*ms2)/BB[c];    
	*ms2d = sqq*(*ms2d)/BB[c];
      } else {
	*ms2 = -(*ms2)/BB[c];    
	*ms2d = -sqq*(*ms2d)/BB[c];
      } 
      
      free(BB);
   }

    return retcode;
  } // int mathieu_modms2
 


} // namespace mathieu
} // namespace xsf

#endif  // #ifndef MATHIEU_FCNS_H

