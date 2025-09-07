#ifndef BESSELJYD_H
#define BESSELJYD_H

#include "../config.h"
#include "../bessel.h"


/*
 *
 * This is part of the Mathieu function suite -- a reimplementation
 * of the Mathieu functions for Scipy.  This file holds helpers
 * to the Bessel J and Y functions and also returns derivatives
 * of those fcns.
 * 
 */

namespace xsf {
namespace mathieu {

  //==================================================================
  double besselj(int k, double z) {
    // This is just a thin wrapper around the Bessel impl in the
    // std library.
    double v = (double) k;
    return xsf::cyl_bessel_j(v, z);
  }

  //==================================================================
  double bessely(int k, double z) {
    // This is just a thin wrapper around the Bessel impl in the
    // std library.
    double v = (double) k;
    return xsf::cyl_bessel_y(v, z);
  }

  //==================================================================
  double besseljd(int k, double z) {
    // This returns the derivative of besselj.  The deriv is
    // computed using common identities.
    double y;
    
    if (k == 0) {
      double v = 1.0;
      y = -besselj(v,z);
    } else {
      double kp1 = (double) (k+1);
      double km1 = (double) (k-1);      
      y = (besselj(km1,z)-besselj(kp1,z))/2.0;
    }

    // Must flip sign for negative k and odd k.
    if (k<0 && ((k % 2) != 0)) {
      y = -y;
    }

    return y;
  }

  //==================================================================
  double besselyd(int k, double z) {
    // This returns the derivative of besselj.  The deriv is
    // computed using common identities.
    double y;
    
    if (k == 0) {
      double v = 1.0;
      y = -bessely(v,z);
    } else {
      double kp1 = (double) (k+1);
      double km1 = (double) (k-1);      
      y = (bessely(km1,z)-bessely(kp1,z))/2.0;
    }

    // Must flip sign for negative k and odd k.
    if (k<0 && ((k % 2) != 0)) {
      y = -y;
    }
    
    return y;
  }
  
} // namespace xsf
} // namespace mathieu

#endif  // #ifndef BESSELJYD_H
