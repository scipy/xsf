#pragma once

#include "error.h"
#include "mathieu/matrix_utils.h"
#include "mathieu/make_matrix.h"
#include "mathieu/mathieu_coeffs.h"
#include "mathieu/mathieu_eigs.h"
#include "mathieu/besseljyd.h"
#include "mathieu/mathieu_fcns.h"

namespace xsf {


/* Characteristic values */
//-------------------------------------------------------------
template <typename T>
T cem_cva(T m, T q) {
  // This returns the even Mathieu characteristic value (eigenvalue) a.

  // Check for invalid Mathieu order.
  if ((m < 0) || (m != floor(m))) {
    set_error("mathieu_a", SF_ERROR_DOMAIN, NULL);
    return std::numeric_limits<T>::quiet_NaN();
  }
  
  // Must cast to correct types prior to fcn call.
  int im = static_cast<int>(m);
  double dq = static_cast<double>(q);
  double da;
    
  int retcode = xsf::mathieu::mathieu_a(im, dq, &da);
  if (retcode != SF_ERROR_OK) {
    set_error("mathieu_a", SF_ERROR_NO_RESULT, NULL);
    return std::numeric_limits<T>::quiet_NaN();
  }

  // Now cast back.
  T a = static_cast<T>(da);
  return a;
}

//-------------------------------------------------------------  
template <typename T>
T sem_cva(T m, T q) {
  // This returns the odd Mathieu characteristic value (eigenvalue) b.

  // Check for invalid Mathieu order.
  if ((m < 1) || (m != floor(m))) {
    set_error("mathieu_b", SF_ERROR_DOMAIN, NULL);
    return std::numeric_limits<T>::quiet_NaN();
  }

  // Must cast to correct types prior to fcn call.
  int im = static_cast<int>(m);
  double dq = static_cast<double>(q);
  double db;
    
  int retcode = xsf::mathieu::mathieu_b(im, dq, &db);
  if (retcode != SF_ERROR_OK) {
    set_error("mathieu_b", SF_ERROR_NO_RESULT, NULL);
    return std::numeric_limits<T>::quiet_NaN();
  }

  // Now cast back.
  T b = static_cast<T>(db);
  return b;
}


//---------------------------------------------------------------
/* Mathieu functions */
template <typename T>
void cem(T m, T q, T x, T &csf, T &csd) {

  if ((m < 0) || (m != floor(m))) {
    csf = std::numeric_limits<T>::quiet_NaN();
    csd = std::numeric_limits<T>::quiet_NaN();
    set_error("mathieu_ce", SF_ERROR_DOMAIN, NULL);
  } else {
    
    // Must cast to correct types prior to fcn call.
    int im = static_cast<int>(m);
    double dq = static_cast<double>(q);
    double dx = static_cast<double>(x);
    double dcsf;
    double dcsd;      
    
    // Call fcn and cast back.
    int retcode = xsf::mathieu::mathieu_ce(im, dq, dx, &dcsf, &dcsd);
    if (retcode != SF_ERROR_OK) {
      csf = std::numeric_limits<T>::quiet_NaN();
      csd = std::numeric_limits<T>::quiet_NaN();
      set_error("mathieu_ce", (sf_error_t) retcode, NULL);
    } else {
      csf = static_cast<T>(dcsf);
      csd = static_cast<T>(dcsd);
    }
  }
}
  

//---------------------------------------------------------------  
template <typename T>
void sem(T m, T q, T x, T &csf, T &csd) {

  if ((m < 1) || (m != floor(m))) {
    csf = std::numeric_limits<T>::quiet_NaN();
    csd = std::numeric_limits<T>::quiet_NaN();
    set_error("mathieu_sem", SF_ERROR_DOMAIN, NULL);
  } else {
    
    // Must cast to correct types prior to fcn call.
    int im = static_cast<int>(m);
    double dq = static_cast<double>(q);
    double dx = static_cast<double>(x);
    double dcsf;
    double dcsd;      
    
    // Call fcn and cast back.
    int retcode = xsf::mathieu::mathieu_se(im, dq, dx, &dcsf, &dcsd);
    if (retcode != SF_ERROR_OK) {
      csf = std::numeric_limits<T>::quiet_NaN();
      csd = std::numeric_limits<T>::quiet_NaN();
      set_error("mathieu_sem", (sf_error_t) retcode, NULL);
    } else {
      csf = static_cast<T>(dcsf);
      csd = static_cast<T>(dcsd);
    }
  }
}

//---------------------------------------------------------------  
template <typename T>
void mcm1(T m, T q, T x, T &f1r, T &d1r) {

  if ((m < 0) || (m != floor(m))) {
    f1r = std::numeric_limits<T>::quiet_NaN();
    d1r = std::numeric_limits<T>::quiet_NaN();
    set_error("mathieu_mcm1", SF_ERROR_DOMAIN, NULL);
  } else {
    
    // Must cast to correct types prior to fcn call.
    int im = static_cast<int>(m);
    double dq = static_cast<double>(q);
    double dx = static_cast<double>(x);
    double df1r;
    double dd1r;      
    
    // Call fcn and cast back.
    int retcode = xsf::mathieu::mathieu_modmc1(im, dq, dx, &df1r, &dd1r);
    if (retcode != SF_ERROR_OK) {
      f1r = std::numeric_limits<T>::quiet_NaN();
      d1r = std::numeric_limits<T>::quiet_NaN();
      set_error("mathieu_mcm1", (sf_error_t) retcode, NULL);
    } else {
      f1r = static_cast<T>(df1r);
      d1r = static_cast<T>(dd1r);
    }
  }
}

//---------------------------------------------------------------  
template <typename T>
void msm1(T m, T q, T x, T &f1r, T &d1r) {

  if ((m < 1) || (m != floor(m))) {
    f1r = std::numeric_limits<T>::quiet_NaN();
    d1r = std::numeric_limits<T>::quiet_NaN();
    set_error("mathieu_msm1", SF_ERROR_DOMAIN, NULL);
  } else {
    
    // Must cast to correct types prior to fcn call.
    int im = static_cast<int>(m);
    double dq = static_cast<double>(q);
    double dx = static_cast<double>(x);
    double df1r;
    double dd1r;      
    
    // Call fcn and cast back.
    int retcode = xsf::mathieu::mathieu_modms1(im, dq, dx, &df1r, &dd1r);
    if (retcode != SF_ERROR_OK) {
      f1r = std::numeric_limits<T>::quiet_NaN();
      d1r = std::numeric_limits<T>::quiet_NaN();
      set_error("mathieu_msm1", (sf_error_t) retcode, NULL);
    } else {
      f1r = static_cast<T>(df1r);
      d1r = static_cast<T>(dd1r);
    }
  }

}

//---------------------------------------------------------------  
template <typename T>
void mcm2(T m, T q, T x, T &f2r, T &d2r) {

  if ((m < 0) || (m != floor(m))) {
    f2r = std::numeric_limits<T>::quiet_NaN();
    d2r = std::numeric_limits<T>::quiet_NaN();
    set_error("mathieu_mcm2", SF_ERROR_DOMAIN, NULL);
  } else {
    
    // Must cast to correct types prior to fcn call.
    int im = static_cast<int>(m);
    double dq = static_cast<double>(q);
    double dx = static_cast<double>(x);
    double df2r;
    double dd2r;      
    
    // Call fcn and cast back.
    int retcode = xsf::mathieu::mathieu_modmc2(im, dq, dx, &df2r, &dd2r);
    if (retcode != SF_ERROR_OK) {
      f2r = std::numeric_limits<T>::quiet_NaN();
      d2r = std::numeric_limits<T>::quiet_NaN();
      set_error("mathieu_mcm2", (sf_error_t) retcode, NULL);
    } else {
      f2r = static_cast<T>(df2r);
      d2r = static_cast<T>(dd2r);
    }
  }

}

//---------------------------------------------------------------  
template <typename T>
void msm2(T m, T q, T x, T &f2r, T &d2r) {

  if ((m < 1) || (m != floor(m))) {
    f2r = std::numeric_limits<T>::quiet_NaN();
    d2r = std::numeric_limits<T>::quiet_NaN();
    set_error("mathieu_msm2", SF_ERROR_DOMAIN, NULL);
  } else {
    
    // Must cast to correct types prior to fcn call.
    int im = static_cast<int>(m);
    double dq = static_cast<double>(q);
    double dx = static_cast<double>(x);
    double df2r;
    double dd2r;      
    
    // Call fcn and cast back.
    int retcode = xsf::mathieu::mathieu_modms2(im, dq, dx, &df2r, &dd2r);
    if (retcode != SF_ERROR_OK) {
      f2r = std::numeric_limits<T>::quiet_NaN();
      d2r = std::numeric_limits<T>::quiet_NaN();
      set_error("mathieu_msm2", (sf_error_t) retcode, NULL);
    } else {
      f2r = static_cast<T>(df2r);
      d2r = static_cast<T>(dd2r);
    }
  }
}

  
} // namespace xsf
