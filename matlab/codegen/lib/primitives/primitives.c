/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * primitives.c
 *
 * Code generation for function 'primitives'
 *
 */

/* Include files */
#include "primitives.h"
#include "rt_nonfinite.h"
#include "rt_defines.h"
#include "rt_nonfinite.h"
#include <math.h>
#include <string.h>

/* Function Declarations */
static double absRelopProxies(double a, const creal_T b, double *y);

static double b_atan2(double y, double x);

static void b_log(creal_T *x);

static double b_xnrm2(int n, const double x[3]);

static void binary_expand_op(creal_T in1_data[], int in1_size[2], double in3,
                             double in4, const creal_T in7_data[],
                             const int *in7_size, double in8,
                             const creal_T in9_data[], const int *in9_size,
                             const creal_T in10_data[], const int *in10_size,
                             const creal_T in11_data[],
                             const creal_T in12_data[], const int *in12_size,
                             const creal_T in13_data[], const int *in13_size,
                             const creal_T in14_data[], const int *in14_size,
                             const creal_T in15_data[], const int *in15_size);

static void final_opt_time_pass(double v0, double a0, double sf, double vf,
                                creal_T final_opt_time_pass_var[4]);

static boolean_T ifWhileCond(const boolean_T x_data[], int x_size);

static int roots(const double c[5], creal_T r_data[]);

static double rt_powd_snf(double u0, double u1);

static int xdlahqr(int ilo, int ihi, double h_data[], const int h_size[2],
                   double wr_data[], int *wr_size, double wi_data[],
                   int *wi_size);

static double xdlanv2(double *a, double *b, double *c, double *d, double *rt1i,
                      double *rt2r, double *rt2i, double *cs, double *sn);

static double xnrm2(int n, const double x_data[], int ix0);

static int xzgebal(double A_data[], const int A_size[2], int *ihi,
                   double scale_data[], int *scale_size);

static void xzgehrd(double a_data[], const int a_size[2], int ilo, int ihi);

static double xzlarfg(int n, double *alpha1, double x[3]);

static void xzlascl(double cfrom, double cto, int m, double A_data[], int iA0);

/* Function Definitions */
/*
 *
 */
static double absRelopProxies(double a, const creal_T b, double *y)
{
  double b_a;
  double b_b;
  double x;
  boolean_T guard1;
  b_a = fabs(b.re);
  guard1 = false;
  if (b_a > 8.9884656743115785E+307) {
    guard1 = true;
  } else {
    b_b = fabs(b.im);
    if (b_b > 8.9884656743115785E+307) {
      guard1 = true;
    } else {
      x = fabs(a);
      if (b_a < b_b) {
        b_a /= b_b;
        *y = b_b * sqrt(b_a * b_a + 1.0);
      } else if (b_a > b_b) {
        b_b /= b_a;
        *y = b_a * sqrt(b_b * b_b + 1.0);
      } else if (rtIsNaN(b_b)) {
        *y = rtNaN;
      } else {
        *y = b_a * 1.4142135623730951;
      }
    }
  }
  if (guard1) {
    x = fabs(a) / 2.0;
    b_a = fabs(b.re / 2.0);
    b_b = fabs(b.im / 2.0);
    if (b_a < b_b) {
      b_a /= b_b;
      *y = b_b * sqrt(b_a * b_a + 1.0);
    } else if (b_a > b_b) {
      b_b /= b_a;
      *y = b_a * sqrt(b_b * b_b + 1.0);
    } else if (rtIsNaN(b_b)) {
      *y = rtNaN;
    } else {
      *y = b_a * 1.4142135623730951;
    }
  }
  if (x == *y) {
    x = b_atan2(0.0, a);
    *y = b_atan2(b.im, b.re);
    if (x == *y) {
      if (a != b.re) {
        if (x >= 0.0) {
          x = b.re;
          *y = a;
        } else {
          x = a;
          *y = b.re;
        }
      } else if (a < 0.0) {
        x = b.im;
        *y = 0.0;
      } else {
        x = 0.0;
        *y = b.im;
      }
      if (x == *y) {
        x = 0.0;
        *y = 0.0;
      }
    }
  }
  return x;
}

/*
 *
 */
static double b_atan2(double y, double x)
{
  double r;
  if (rtIsNaN(y) || rtIsNaN(x)) {
    r = rtNaN;
  } else if (rtIsInf(y) && rtIsInf(x)) {
    int i;
    int i1;
    if (y > 0.0) {
      i = 1;
    } else {
      i = -1;
    }
    if (x > 0.0) {
      i1 = 1;
    } else {
      i1 = -1;
    }
    r = atan2(i, i1);
  } else if (x == 0.0) {
    if (y > 0.0) {
      r = RT_PI / 2.0;
    } else if (y < 0.0) {
      r = -(RT_PI / 2.0);
    } else {
      r = 0.0;
    }
  } else {
    r = atan2(y, x);
  }
  return r;
}

/*
 *
 */
static void b_log(creal_T *x)
{
  if (x->im == 0.0) {
    if (x->re < 0.0) {
      x->re = log(fabs(x->re));
      x->im = 3.1415926535897931;
    } else {
      x->re = log(fabs(x->re));
      x->im = 0.0;
    }
  } else {
    double a;
    double b;
    boolean_T guard1;
    a = fabs(x->re);
    guard1 = false;
    if (a > 8.9884656743115785E+307) {
      guard1 = true;
    } else {
      b = fabs(x->im);
      if (b > 8.9884656743115785E+307) {
        guard1 = true;
      } else {
        if (a < b) {
          a /= b;
          a = b * sqrt(a * a + 1.0);
        } else if (a > b) {
          b /= a;
          a *= sqrt(b * b + 1.0);
        } else if (rtIsNaN(b)) {
          a = rtNaN;
        } else {
          a *= 1.4142135623730951;
        }
        b = x->re;
        x->re = log(a);
        x->im = b_atan2(x->im, b);
      }
    }
    if (guard1) {
      a = fabs(x->re / 2.0);
      b = fabs(x->im / 2.0);
      if (a < b) {
        a /= b;
        a = b * sqrt(a * a + 1.0);
      } else if (a > b) {
        b /= a;
        a *= sqrt(b * b + 1.0);
      } else if (rtIsNaN(b)) {
        a = rtNaN;
      } else {
        a *= 1.4142135623730951;
      }
      b = x->re;
      x->re = log(a) + 0.69314718055994529;
      x->im = b_atan2(x->im, b);
    }
  }
}

/*
 *
 */
static double b_xnrm2(int n, const double x[3])
{
  double y;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = fabs(x[1]);
    } else {
      double absxk;
      double scale;
      double t;
      scale = 3.3121686421112381E-170;
      absxk = fabs(x[1]);
      if (absxk > 3.3121686421112381E-170) {
        y = 1.0;
        scale = absxk;
      } else {
        t = absxk / 3.3121686421112381E-170;
        y = t * t;
      }
      absxk = fabs(x[2]);
      if (absxk > scale) {
        t = scale / absxk;
        y = y * t * t + 1.0;
        scale = absxk;
      } else {
        t = absxk / scale;
        y += t * t;
      }
      y = scale * sqrt(y);
    }
  }
  return y;
}

static void binary_expand_op(creal_T in1_data[], int in1_size[2], double in3,
                             double in4, const creal_T in7_data[],
                             const int *in7_size, double in8,
                             const creal_T in9_data[], const int *in9_size,
                             const creal_T in10_data[], const int *in10_size,
                             const creal_T in11_data[],
                             const creal_T in12_data[], const int *in12_size,
                             const creal_T in13_data[], const int *in13_size,
                             const creal_T in14_data[], const int *in14_size,
                             const creal_T in15_data[], const int *in15_size)
{
  creal_T in2[6];
  creal_T b_in5_data[4];
  creal_T c_in5_data[4];
  creal_T in5_data[4];
  creal_T varargin_1;
  double b_in8_re;
  double brm;
  double c_in8_re;
  double im;
  double in8_re;
  double r;
  double re;
  double varargout_1_tmp;
  int b_in8_re_tmp;
  int b_loop_ub;
  int i;
  int in8_re_tmp;
  int loop_ub;
  int stride_0_0;
  int stride_0_0_tmp;
  int stride_1_0;
  int stride_2_0_tmp;
  int stride_3_0;
  in8_re = in8 * -20.0;
  b_in8_re = in8 * -30.0;
  c_in8_re = in8 * -12.0;
  if (*in12_size == 1) {
    i = *in9_size;
  } else {
    i = *in12_size;
  }
  if (i == 1) {
    loop_ub = *in7_size;
  } else {
    loop_ub = i;
  }
  stride_0_0_tmp = (*in7_size != 1);
  stride_1_0 = (*in9_size != 1);
  stride_2_0_tmp = (*in12_size != 1);
  for (i = 0; i < loop_ub; i++) {
    varargin_1 = in7_data[i * stride_0_0_tmp];
    if ((varargin_1.im == 0.0) && (varargin_1.re >= 0.0)) {
      re = rt_powd_snf(varargin_1.re, 3.0);
      im = 0.0;
    } else if (varargin_1.re == 0.0) {
      re = 0.0;
      im = -rt_powd_snf(varargin_1.im, 3.0);
    } else {
      b_log(&varargin_1);
      r = 3.0 * varargin_1.re;
      varargout_1_tmp = 3.0 * varargin_1.im;
      if (r == 0.0) {
        re = cos(varargout_1_tmp);
        im = sin(varargout_1_tmp);
      } else if (varargout_1_tmp == 0.0) {
        re = exp(r);
        im = 0.0;
      } else if (rtIsInf(varargout_1_tmp) && rtIsInf(r) && (r < 0.0)) {
        re = 0.0;
        im = 0.0;
      } else {
        r = exp(r / 2.0);
        re = r * (r * cos(varargout_1_tmp));
        im = r * (r * sin(varargout_1_tmp));
      }
    }
    if (im == 0.0) {
      re = 1.0 / re;
      r = 0.0;
    } else if (re == 0.0) {
      re = 0.0;
      r = -(1.0 / im);
    } else {
      brm = fabs(re);
      r = fabs(im);
      if (brm > r) {
        varargout_1_tmp = im / re;
        r = re + varargout_1_tmp * im;
        re = (varargout_1_tmp * 0.0 + 1.0) / r;
        r = (0.0 - varargout_1_tmp) / r;
      } else if (r == brm) {
        if (re > 0.0) {
          varargout_1_tmp = 0.5;
        } else {
          varargout_1_tmp = -0.5;
        }
        if (im > 0.0) {
          r = 0.5;
        } else {
          r = -0.5;
        }
        re = (varargout_1_tmp + 0.0 * r) / brm;
        r = (0.0 * varargout_1_tmp - r) / brm;
      } else {
        varargout_1_tmp = re / im;
        r = im + varargout_1_tmp * re;
        re = varargout_1_tmp / r;
        r = (varargout_1_tmp * 0.0 - 1.0) / r;
      }
    }
    in8_re_tmp = i * stride_1_0;
    b_in8_re_tmp = i * stride_2_0_tmp;
    varargout_1_tmp =
        (((in8_re + in9_data[in8_re_tmp].re) + in10_data[in8_re_tmp].re) +
         in11_data[in8_re_tmp].re * 12.0) +
        in12_data[b_in8_re_tmp].re * 8.0;
    im = ((in9_data[in8_re_tmp].im + in10_data[in8_re_tmp].im) +
          in11_data[in8_re_tmp].im * 12.0) +
         in12_data[b_in8_re_tmp].im * 8.0;
    in5_data[i].re = -3.0 * (re * varargout_1_tmp - r * im);
    in5_data[i].im = -3.0 * (re * im + r * varargout_1_tmp);
  }
  if (*in12_size == 1) {
    i = *in14_size;
  } else {
    i = *in12_size;
  }
  if (i == 1) {
    b_loop_ub = *in13_size;
  } else {
    b_loop_ub = i;
  }
  stride_0_0 = (*in13_size != 1);
  stride_1_0 = (*in14_size != 1);
  for (i = 0; i < b_loop_ub; i++) {
    varargin_1 = in13_data[i * stride_0_0];
    re = varargin_1.re * varargin_1.re - varargin_1.im * varargin_1.im;
    r = varargin_1.re * varargin_1.im;
    im = r + r;
    if (im == 0.0) {
      re = 1.0 / re;
      r = 0.0;
    } else if (re == 0.0) {
      re = 0.0;
      r = -(1.0 / im);
    } else {
      brm = fabs(re);
      r = fabs(im);
      if (brm > r) {
        varargout_1_tmp = im / re;
        r = re + varargout_1_tmp * im;
        re = (varargout_1_tmp * 0.0 + 1.0) / r;
        r = (0.0 - varargout_1_tmp) / r;
      } else if (r == brm) {
        if (re > 0.0) {
          varargout_1_tmp = 0.5;
        } else {
          varargout_1_tmp = -0.5;
        }
        if (im > 0.0) {
          r = 0.5;
        } else {
          r = -0.5;
        }
        re = (varargout_1_tmp + 0.0 * r) / brm;
        r = (0.0 * varargout_1_tmp - r) / brm;
      } else {
        varargout_1_tmp = re / im;
        r = im + varargout_1_tmp * re;
        re = varargout_1_tmp / r;
        r = (varargout_1_tmp * 0.0 - 1.0) / r;
      }
    }
    in8_re_tmp = i * stride_1_0;
    b_in8_re_tmp = i * stride_2_0_tmp;
    in8_re = (((b_in8_re - in14_data[in8_re_tmp].re * 2.0) +
               in9_data[in8_re_tmp].re) +
              in11_data[in8_re_tmp].re * 16.0) +
             in12_data[b_in8_re_tmp].re * 14.0;
    im = (((0.0 - in14_data[in8_re_tmp].im * 2.0) + in9_data[in8_re_tmp].im) +
          in11_data[in8_re_tmp].im * 16.0) +
         in12_data[b_in8_re_tmp].im * 14.0;
    b_in5_data[i].re = 12.0 * (re * in8_re - r * im);
    b_in5_data[i].im = 12.0 * (re * im + r * in8_re);
  }
  if (*in15_size == 1) {
    if (*in12_size == 1) {
      i = *in10_size;
    } else {
      i = *in12_size;
    }
  } else {
    i = *in15_size;
  }
  if (i == 1) {
    stride_0_0 = *in7_size;
  } else {
    stride_0_0 = i;
  }
  stride_1_0 = (*in10_size != 1);
  stride_3_0 = (*in15_size != 1);
  for (i = 0; i < stride_0_0; i++) {
    int c_in8_re_tmp;
    varargin_1 = in7_data[i * stride_0_0_tmp];
    if ((varargin_1.im == 0.0) && (varargin_1.re >= 0.0)) {
      re = rt_powd_snf(varargin_1.re, 5.0);
      im = 0.0;
    } else if (varargin_1.re == 0.0) {
      re = 0.0;
      im = rt_powd_snf(varargin_1.im, 5.0);
    } else {
      b_log(&varargin_1);
      r = 5.0 * varargin_1.re;
      varargout_1_tmp = 5.0 * varargin_1.im;
      if (r == 0.0) {
        re = cos(varargout_1_tmp);
        im = sin(varargout_1_tmp);
      } else if (varargout_1_tmp == 0.0) {
        re = exp(r);
        im = 0.0;
      } else if (rtIsInf(varargout_1_tmp) && rtIsInf(r) && (r < 0.0)) {
        re = 0.0;
        im = 0.0;
      } else {
        r = exp(r / 2.0);
        re = r * (r * cos(varargout_1_tmp));
        im = r * (r * sin(varargout_1_tmp));
      }
    }
    if (im == 0.0) {
      re = 1.0 / re;
      r = 0.0;
    } else if (re == 0.0) {
      re = 0.0;
      r = -(1.0 / im);
    } else {
      brm = fabs(re);
      r = fabs(im);
      if (brm > r) {
        varargout_1_tmp = im / re;
        r = re + varargout_1_tmp * im;
        re = (varargout_1_tmp * 0.0 + 1.0) / r;
        r = (0.0 - varargout_1_tmp) / r;
      } else if (r == brm) {
        if (re > 0.0) {
          varargout_1_tmp = 0.5;
        } else {
          varargout_1_tmp = -0.5;
        }
        if (im > 0.0) {
          r = 0.5;
        } else {
          r = -0.5;
        }
        re = (varargout_1_tmp + 0.0 * r) / brm;
        r = (0.0 * varargout_1_tmp - r) / brm;
      } else {
        varargout_1_tmp = re / im;
        r = im + varargout_1_tmp * re;
        re = varargout_1_tmp / r;
        r = (varargout_1_tmp * 0.0 - 1.0) / r;
      }
    }
    in8_re_tmp = i * stride_1_0;
    b_in8_re_tmp = i * stride_2_0_tmp;
    c_in8_re_tmp = i * stride_3_0;
    in8_re = (((c_in8_re + in10_data[in8_re_tmp].re) +
               in11_data[in8_re_tmp].re * 6.0) +
              in12_data[b_in8_re_tmp].re * 6.0) +
             in15_data[c_in8_re_tmp].re;
    im = ((in10_data[in8_re_tmp].im + in11_data[in8_re_tmp].im * 6.0) +
          in12_data[b_in8_re_tmp].im * 6.0) +
         in15_data[c_in8_re_tmp].im;
    c_in5_data[i].re = -60.0 * (re * in8_re - r * im);
    c_in5_data[i].im = -60.0 * (re * im + r * in8_re);
  }
  in2[0].re = 0.0;
  in2[0].im = 0.0;
  in2[1].re = in3;
  in2[1].im = 0.0;
  in2[2].re = in4;
  in2[2].im = 0.0;
  if (loop_ub - 1 >= 0) {
    memcpy(&in2[3], &in5_data[0], (unsigned int)loop_ub * sizeof(creal_T));
  }
  if (b_loop_ub - 1 >= 0) {
    memcpy(&in2[4], &b_in5_data[0], (unsigned int)b_loop_ub * sizeof(creal_T));
  }
  if (stride_0_0 - 1 >= 0) {
    memcpy(&in2[5], &c_in5_data[0], (unsigned int)stride_0_0 * sizeof(creal_T));
  }
  in1_size[0] = 1;
  in1_size[1] = 6;
  for (i = 0; i < 6; i++) {
    in1_data[in1_size[0] * i] = in2[i];
  }
}

/*
 * function final_opt_time_pass_var = final_opt_time_pass(v0,a0,sf,vf)
 */
static void final_opt_time_pass(double v0, double a0, double sf, double vf,
                                creal_T final_opt_time_pass_var[4])
{
  creal_T b_tmp_data[4];
  creal_T c_tmp_data[4];
  creal_T tmp_data[4];
  double b_sf[5];
  double b_sf_tmp;
  double c_sf_tmp;
  double d_sf_tmp;
  double e_sf_tmp;
  double sf_tmp;
  /* FINAL_OPT_TIME_PASS */
  /*     FINAL_OPT_TIME_PASS_VAR = FINAL_OPT_TIME_PASS(V0,A0,SF,VF) */
  /*     This function was generated by the Symbolic Math Toolbox version 24.2.
   */
  /*     24-Oct-2025 11:35:05 */
  /* 'final_opt_time_pass:8' t0 =
   * roots([sf.^2.*1.2e+3,sf.*v0.*-9.6e+2-sf.*vf.*9.6e+2,a0.*sf.*-1.2e+2+v0.*vf.*3.36e+2+v0.^2.*1.92e+2+vf.^2.*1.92e+2,a0.*v0.*4.8e+1+a0.*vf.*3.2e+1,a0.^2.*3.0]);
   */
  sf_tmp = sf * sf * 1200.0;
  b_sf[0] = sf_tmp;
  b_sf_tmp = sf * v0 * -960.0 - sf * vf * 960.0;
  b_sf[1] = b_sf_tmp;
  c_sf_tmp = ((a0 * sf * -120.0 + v0 * vf * 336.0) + v0 * v0 * 192.0) +
             vf * vf * 192.0;
  b_sf[2] = c_sf_tmp;
  d_sf_tmp = a0 * v0 * 48.0 + a0 * vf * 32.0;
  b_sf[3] = d_sf_tmp;
  e_sf_tmp = a0 * a0 * 3.0;
  b_sf[4] = e_sf_tmp;
  /* 'final_opt_time_pass:9' t2 = t0(1); */
  roots(b_sf, tmp_data);
  /* 'final_opt_time_pass:10' t0 =
   * roots([sf.^2.*1.2e+3,sf.*v0.*-9.6e+2-sf.*vf.*9.6e+2,a0.*sf.*-1.2e+2+v0.*vf.*3.36e+2+v0.^2.*1.92e+2+vf.^2.*1.92e+2,a0.*v0.*4.8e+1+a0.*vf.*3.2e+1,a0.^2.*3.0]);
   */
  b_sf[0] = sf_tmp;
  b_sf[1] = b_sf_tmp;
  b_sf[2] = c_sf_tmp;
  b_sf[3] = d_sf_tmp;
  b_sf[4] = e_sf_tmp;
  /* 'final_opt_time_pass:11' t3 = t0(2); */
  roots(b_sf, b_tmp_data);
  /* 'final_opt_time_pass:12' t0 =
   * roots([sf.^2.*1.2e+3,sf.*v0.*-9.6e+2-sf.*vf.*9.6e+2,a0.*sf.*-1.2e+2+v0.*vf.*3.36e+2+v0.^2.*1.92e+2+vf.^2.*1.92e+2,a0.*v0.*4.8e+1+a0.*vf.*3.2e+1,a0.^2.*3.0]);
   */
  b_sf[0] = sf_tmp;
  b_sf[1] = b_sf_tmp;
  b_sf[2] = c_sf_tmp;
  b_sf[3] = d_sf_tmp;
  b_sf[4] = e_sf_tmp;
  /* 'final_opt_time_pass:13' t4 = t0(3); */
  roots(b_sf, c_tmp_data);
  /* 'final_opt_time_pass:14' t0 =
   * roots([sf.^2.*1.2e+3,sf.*v0.*-9.6e+2-sf.*vf.*9.6e+2,a0.*sf.*-1.2e+2+v0.*vf.*3.36e+2+v0.^2.*1.92e+2+vf.^2.*1.92e+2,a0.*v0.*4.8e+1+a0.*vf.*3.2e+1,a0.^2.*3.0]);
   */
  b_sf[0] = sf_tmp;
  b_sf[1] = b_sf_tmp;
  b_sf[2] = c_sf_tmp;
  b_sf[3] = d_sf_tmp;
  b_sf[4] = e_sf_tmp;
  /* 'final_opt_time_pass:15' t5 = t0(4); */
  /* 'final_opt_time_pass:16' final_opt_time_pass_var =
   * [1.0./t2;1.0./t3;1.0./t4;1.0./t5]; */
  if (tmp_data[0].im == 0.0) {
    final_opt_time_pass_var[0].re = 1.0 / tmp_data[0].re;
    final_opt_time_pass_var[0].im = 0.0;
  } else if (tmp_data[0].re == 0.0) {
    final_opt_time_pass_var[0].re = 0.0;
    final_opt_time_pass_var[0].im = -(1.0 / tmp_data[0].im);
  } else {
    c_sf_tmp = fabs(tmp_data[0].re);
    sf_tmp = fabs(tmp_data[0].im);
    if (c_sf_tmp > sf_tmp) {
      sf_tmp = tmp_data[0].im / tmp_data[0].re;
      b_sf_tmp = tmp_data[0].re + sf_tmp * tmp_data[0].im;
      final_opt_time_pass_var[0].re = (sf_tmp * 0.0 + 1.0) / b_sf_tmp;
      final_opt_time_pass_var[0].im = (0.0 - sf_tmp) / b_sf_tmp;
    } else if (sf_tmp == c_sf_tmp) {
      if (tmp_data[0].re > 0.0) {
        sf_tmp = 0.5;
      } else {
        sf_tmp = -0.5;
      }
      if (tmp_data[0].im > 0.0) {
        b_sf_tmp = 0.5;
      } else {
        b_sf_tmp = -0.5;
      }
      final_opt_time_pass_var[0].re = (sf_tmp + 0.0 * b_sf_tmp) / c_sf_tmp;
      final_opt_time_pass_var[0].im = (0.0 * sf_tmp - b_sf_tmp) / c_sf_tmp;
    } else {
      sf_tmp = tmp_data[0].re / tmp_data[0].im;
      b_sf_tmp = tmp_data[0].im + sf_tmp * tmp_data[0].re;
      final_opt_time_pass_var[0].re = sf_tmp / b_sf_tmp;
      final_opt_time_pass_var[0].im = (sf_tmp * 0.0 - 1.0) / b_sf_tmp;
    }
  }
  if (b_tmp_data[1].im == 0.0) {
    final_opt_time_pass_var[1].re = 1.0 / b_tmp_data[1].re;
    final_opt_time_pass_var[1].im = 0.0;
  } else if (b_tmp_data[1].re == 0.0) {
    final_opt_time_pass_var[1].re = 0.0;
    final_opt_time_pass_var[1].im = -(1.0 / b_tmp_data[1].im);
  } else {
    c_sf_tmp = fabs(b_tmp_data[1].re);
    sf_tmp = fabs(b_tmp_data[1].im);
    if (c_sf_tmp > sf_tmp) {
      sf_tmp = b_tmp_data[1].im / b_tmp_data[1].re;
      b_sf_tmp = b_tmp_data[1].re + sf_tmp * b_tmp_data[1].im;
      final_opt_time_pass_var[1].re = (sf_tmp * 0.0 + 1.0) / b_sf_tmp;
      final_opt_time_pass_var[1].im = (0.0 - sf_tmp) / b_sf_tmp;
    } else if (sf_tmp == c_sf_tmp) {
      if (b_tmp_data[1].re > 0.0) {
        sf_tmp = 0.5;
      } else {
        sf_tmp = -0.5;
      }
      if (b_tmp_data[1].im > 0.0) {
        b_sf_tmp = 0.5;
      } else {
        b_sf_tmp = -0.5;
      }
      final_opt_time_pass_var[1].re = (sf_tmp + 0.0 * b_sf_tmp) / c_sf_tmp;
      final_opt_time_pass_var[1].im = (0.0 * sf_tmp - b_sf_tmp) / c_sf_tmp;
    } else {
      sf_tmp = b_tmp_data[1].re / b_tmp_data[1].im;
      b_sf_tmp = b_tmp_data[1].im + sf_tmp * b_tmp_data[1].re;
      final_opt_time_pass_var[1].re = sf_tmp / b_sf_tmp;
      final_opt_time_pass_var[1].im = (sf_tmp * 0.0 - 1.0) / b_sf_tmp;
    }
  }
  if (c_tmp_data[2].im == 0.0) {
    final_opt_time_pass_var[2].re = 1.0 / c_tmp_data[2].re;
    final_opt_time_pass_var[2].im = 0.0;
  } else if (c_tmp_data[2].re == 0.0) {
    final_opt_time_pass_var[2].re = 0.0;
    final_opt_time_pass_var[2].im = -(1.0 / c_tmp_data[2].im);
  } else {
    c_sf_tmp = fabs(c_tmp_data[2].re);
    sf_tmp = fabs(c_tmp_data[2].im);
    if (c_sf_tmp > sf_tmp) {
      sf_tmp = c_tmp_data[2].im / c_tmp_data[2].re;
      b_sf_tmp = c_tmp_data[2].re + sf_tmp * c_tmp_data[2].im;
      final_opt_time_pass_var[2].re = (sf_tmp * 0.0 + 1.0) / b_sf_tmp;
      final_opt_time_pass_var[2].im = (0.0 - sf_tmp) / b_sf_tmp;
    } else if (sf_tmp == c_sf_tmp) {
      if (c_tmp_data[2].re > 0.0) {
        sf_tmp = 0.5;
      } else {
        sf_tmp = -0.5;
      }
      if (c_tmp_data[2].im > 0.0) {
        b_sf_tmp = 0.5;
      } else {
        b_sf_tmp = -0.5;
      }
      final_opt_time_pass_var[2].re = (sf_tmp + 0.0 * b_sf_tmp) / c_sf_tmp;
      final_opt_time_pass_var[2].im = (0.0 * sf_tmp - b_sf_tmp) / c_sf_tmp;
    } else {
      sf_tmp = c_tmp_data[2].re / c_tmp_data[2].im;
      b_sf_tmp = c_tmp_data[2].im + sf_tmp * c_tmp_data[2].re;
      final_opt_time_pass_var[2].re = sf_tmp / b_sf_tmp;
      final_opt_time_pass_var[2].im = (sf_tmp * 0.0 - 1.0) / b_sf_tmp;
    }
  }
  roots(b_sf, tmp_data);
  if (tmp_data[3].im == 0.0) {
    final_opt_time_pass_var[3].re = 1.0 / tmp_data[3].re;
    final_opt_time_pass_var[3].im = 0.0;
  } else if (tmp_data[3].re == 0.0) {
    final_opt_time_pass_var[3].re = 0.0;
    final_opt_time_pass_var[3].im = -(1.0 / tmp_data[3].im);
  } else {
    c_sf_tmp = fabs(tmp_data[3].re);
    sf_tmp = fabs(tmp_data[3].im);
    if (c_sf_tmp > sf_tmp) {
      sf_tmp = tmp_data[3].im / tmp_data[3].re;
      b_sf_tmp = tmp_data[3].re + sf_tmp * tmp_data[3].im;
      final_opt_time_pass_var[3].re = (sf_tmp * 0.0 + 1.0) / b_sf_tmp;
      final_opt_time_pass_var[3].im = (0.0 - sf_tmp) / b_sf_tmp;
    } else if (sf_tmp == c_sf_tmp) {
      if (tmp_data[3].re > 0.0) {
        sf_tmp = 0.5;
      } else {
        sf_tmp = -0.5;
      }
      if (tmp_data[3].im > 0.0) {
        b_sf_tmp = 0.5;
      } else {
        b_sf_tmp = -0.5;
      }
      final_opt_time_pass_var[3].re = (sf_tmp + 0.0 * b_sf_tmp) / c_sf_tmp;
      final_opt_time_pass_var[3].im = (0.0 * sf_tmp - b_sf_tmp) / c_sf_tmp;
    } else {
      sf_tmp = tmp_data[3].re / tmp_data[3].im;
      b_sf_tmp = tmp_data[3].im + sf_tmp * tmp_data[3].re;
      final_opt_time_pass_var[3].re = sf_tmp / b_sf_tmp;
      final_opt_time_pass_var[3].im = (sf_tmp * 0.0 - 1.0) / b_sf_tmp;
    }
  }
}

/*
 *
 */
static boolean_T ifWhileCond(const boolean_T x_data[], int x_size)
{
  boolean_T y;
  y = (x_size != 0);
  if (y) {
    int k;
    boolean_T exitg1;
    k = 0;
    exitg1 = false;
    while ((!exitg1) && (k <= x_size - 1)) {
      if (!x_data[k]) {
        y = false;
        exitg1 = true;
      } else {
        k++;
      }
    }
  }
  return y;
}

/*
 *
 */
static int roots(const double c[5], creal_T r_data[])
{
  double A_data[16];
  double a_data[16];
  double wi_data[4];
  double wr_data[4];
  int A_size[2];
  int b_i;
  int i;
  int j;
  int k1;
  int k2;
  int loop_ub_tmp;
  int nTrailingZeros_tmp;
  int r_size;
  memset(&r_data[0], 0, 4U * sizeof(creal_T));
  k1 = 1;
  while ((k1 <= 5) && (!(c[k1 - 1] != 0.0))) {
    k1++;
  }
  k2 = 5;
  while ((k2 >= k1) && (!(c[k2 - 1] != 0.0))) {
    k2--;
  }
  nTrailingZeros_tmp = 5 - k2;
  if (k1 < k2) {
    double ctmp[5];
    int companDim;
    boolean_T exitg1;
    companDim = k2 - k1;
    exitg1 = false;
    while ((!exitg1) && (companDim > 0)) {
      boolean_T exitg2;
      j = 0;
      exitg2 = false;
      while ((!exitg2) && (j + 1 <= companDim)) {
        ctmp[j] = c[k1 + j] / c[k1 - 1];
        if (rtIsInf(fabs(ctmp[j]))) {
          exitg2 = true;
        } else {
          j++;
        }
      }
      if (j + 1 > companDim) {
        exitg1 = true;
      } else {
        k1++;
        companDim--;
      }
    }
    if (companDim < 1) {
      r_size = 5 - k2;
    } else {
      creal_T eiga_data[4];
      loop_ub_tmp = companDim * companDim;
      memset(&a_data[0], 0, (unsigned int)loop_ub_tmp * sizeof(double));
      for (k1 = 0; k1 <= companDim - 2; k1++) {
        j = companDim * k1;
        a_data[j] = -ctmp[k1];
        a_data[(k1 + j) + 1] = 1.0;
      }
      a_data[companDim * (companDim - 1)] = -ctmp[companDim - 1];
      if (nTrailingZeros_tmp - 1 >= 0) {
        memset(&r_data[0], 0,
               (unsigned int)nTrailingZeros_tmp * sizeof(creal_T));
      }
      if (companDim == 1) {
        for (i = 0; i < companDim; i++) {
          eiga_data[i].re = a_data[i];
          eiga_data[i].im = 0.0;
        }
      } else {
        double absxk;
        double anrm;
        anrm = 0.0;
        k1 = 0;
        exitg1 = false;
        while ((!exitg1) && (k1 <= loop_ub_tmp - 1)) {
          absxk = fabs(a_data[k1]);
          if (rtIsNaN(absxk)) {
            anrm = rtNaN;
            exitg1 = true;
          } else {
            if (absxk > anrm) {
              anrm = absxk;
            }
            k1++;
          }
        }
        if (rtIsInf(anrm) || rtIsNaN(anrm)) {
          for (i = 0; i < companDim; i++) {
            eiga_data[i].re = rtNaN;
            eiga_data[i].im = 0.0;
          }
        } else {
          boolean_T guard1;
          boolean_T scalea;
          absxk = anrm;
          scalea = false;
          guard1 = false;
          if ((anrm > 0.0) && (anrm < 6.7178761075670888E-139)) {
            scalea = true;
            absxk = 6.7178761075670888E-139;
            guard1 = true;
          } else if (anrm > 1.4885657073574029E+138) {
            scalea = true;
            absxk = 1.4885657073574029E+138;
            guard1 = true;
          }
          if (guard1) {
            double cfromc;
            double ctoc;
            boolean_T notdone;
            cfromc = anrm;
            ctoc = absxk;
            notdone = true;
            while (notdone) {
              double cfrom1;
              double cto1;
              double mul;
              cfrom1 = cfromc * 2.0041683600089728E-292;
              cto1 = ctoc / 4.9896007738368E+291;
              if ((cfrom1 > ctoc) && (ctoc != 0.0)) {
                mul = 2.0041683600089728E-292;
                cfromc = cfrom1;
              } else if (cto1 > cfromc) {
                mul = 4.9896007738368E+291;
                ctoc = cto1;
              } else {
                mul = ctoc / cfromc;
                notdone = false;
              }
              for (j = 0; j < companDim; j++) {
                k1 = j * companDim - 1;
                for (b_i = 0; b_i < companDim; b_i++) {
                  i = (k1 + b_i) + 1;
                  a_data[i] *= mul;
                }
              }
            }
          }
          A_size[0] = companDim;
          A_size[1] = companDim;
          memcpy(&A_data[0], &a_data[0],
                 (unsigned int)loop_ub_tmp * sizeof(double));
          nTrailingZeros_tmp =
              xzgebal(A_data, A_size, &k1, wr_data, &loop_ub_tmp);
          xzgehrd(A_data, A_size, nTrailingZeros_tmp, k1);
          r_size = xdlahqr(nTrailingZeros_tmp, k1, A_data, A_size, wr_data,
                           &loop_ub_tmp, wi_data, &j);
          if (scalea) {
            i = companDim - r_size;
            xzlascl(absxk, anrm, i, wr_data, r_size + 1);
            xzlascl(absxk, anrm, i, wi_data, r_size + 1);
            if (r_size != 0) {
              xzlascl(absxk, anrm, nTrailingZeros_tmp - 1, wr_data, 1);
              xzlascl(absxk, anrm, nTrailingZeros_tmp - 1, wi_data, 1);
            }
          }
          if (r_size != 0) {
            for (b_i = nTrailingZeros_tmp; b_i <= r_size; b_i++) {
              wr_data[b_i - 1] = rtNaN;
              wi_data[b_i - 1] = 0.0;
            }
          }
          for (i = 0; i < loop_ub_tmp; i++) {
            eiga_data[i].re = wr_data[i];
            eiga_data[i].im = wi_data[i];
          }
        }
      }
      for (k1 = 0; k1 < companDim; k1++) {
        r_data[(k1 - k2) + 5] = eiga_data[k1];
      }
      r_size = (companDim - k2) + 5;
    }
  } else {
    r_size = 5 - k2;
  }
  return r_size;
}

static double rt_powd_snf(double u0, double u1)
{
  double y;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else {
    double d;
    y = fabs(u0);
    d = fabs(u1);
    if (rtIsInf(u1)) {
      if (y == 1.0) {
        y = 1.0;
      } else if (y > 1.0) {
        if (u1 > 0.0) {
          y = rtInf;
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = rtInf;
      }
    } else if (d == 0.0) {
      y = 1.0;
    } else if (d == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > floor(u1))) {
      y = rtNaN;
    } else {
      y = pow(u0, u1);
    }
  }
  return y;
}

/*
 *
 */
static int xdlahqr(int ilo, int ihi, double h_data[], const int h_size[2],
                   double wr_data[], int *wr_size, double wi_data[],
                   int *wi_size)
{
  double d;
  double h21;
  double h22;
  double rt1r;
  double rt2r;
  double s;
  double tst;
  int b_i;
  int b_k;
  int c_k;
  int i;
  int info;
  int j;
  *wr_size = h_size[0];
  *wi_size = *wr_size;
  info = 0;
  i = (unsigned char)(ilo - 1);
  for (b_i = 0; b_i < i; b_i++) {
    wr_data[b_i] = h_data[b_i + h_size[0] * b_i];
    wi_data[b_i] = 0.0;
  }
  i = ihi + 1;
  for (b_i = i; b_i <= *wr_size; b_i++) {
    wr_data[b_i - 1] = h_data[(b_i + h_size[0] * (b_i - 1)) - 1];
    wi_data[b_i - 1] = 0.0;
  }
  if (ilo == ihi) {
    wr_data[ilo - 1] = h_data[(ilo + h_size[0] * (ilo - 1)) - 1];
    wi_data[ilo - 1] = 0.0;
  } else {
    double smlnum;
    int kdefl;
    boolean_T exitg1;
    i = ihi - 3;
    for (j = ilo; j <= i; j++) {
      h_data[2] = 0.0;
      h_data[3] = 0.0;
    }
    if (ilo <= ihi - 2) {
      h_data[(ihi + h_size[0] * (ihi - 3)) - 1] = 0.0;
    }
    smlnum = 2.2250738585072014E-308 *
             ((double)((ihi - ilo) + 1) / 2.2204460492503131E-16);
    kdefl = 0;
    b_i = ihi - 1;
    exitg1 = false;
    while ((!exitg1) && (b_i + 1 >= ilo)) {
      int its;
      int l;
      int v_tmp;
      boolean_T converged;
      boolean_T exitg2;
      l = ilo;
      converged = false;
      its = 0;
      exitg2 = false;
      while ((!exitg2) && (its < 301)) {
        double aa;
        double h12;
        int k;
        int nr;
        boolean_T exitg3;
        k = b_i;
        exitg3 = false;
        while ((!exitg3) && (k + 1 > l)) {
          i = k + h_size[0] * (k - 1);
          d = fabs(h_data[i]);
          if (d <= smlnum) {
            exitg3 = true;
          } else {
            nr = k + h_size[0] * k;
            h21 = fabs(h_data[nr]);
            tst = fabs(h_data[i - 1]) + h21;
            if (tst == 0.0) {
              if (k - 1 >= ilo) {
                tst = fabs(h_data[(k + h_size[0] * (k - 2)) - 1]);
              }
              if (k + 2 <= ihi) {
                tst += fabs(h_data[nr + 1]);
              }
            }
            if (d <= 2.2204460492503131E-16 * tst) {
              h12 = fabs(h_data[nr - 1]);
              tst = fabs(h_data[i - 1] - h_data[nr]);
              aa = fmax(h21, tst);
              tst = fmin(h21, tst);
              s = aa + tst;
              if (fmin(d, h12) * (fmax(d, h12) / s) <=
                  fmax(smlnum, 2.2204460492503131E-16 * (tst * (aa / s)))) {
                exitg3 = true;
              } else {
                k--;
              }
            } else {
              k--;
            }
          }
        }
        l = k + 1;
        if (k + 1 > ilo) {
          h_data[k + h_size[0] * (k - 1)] = 0.0;
        }
        if (k + 1 >= b_i) {
          converged = true;
          exitg2 = true;
        } else {
          double v[3];
          int m;
          kdefl++;
          if (kdefl - kdefl / 20 * 20 == 0) {
            s = fabs(h_data[b_i + h_size[0] * (b_i - 1)]) +
                fabs(h_data[(b_i + h_size[0] * (b_i - 2)) - 1]);
            tst = 0.75 * s + h_data[b_i + h_size[0] * b_i];
            h12 = -0.4375 * s;
            h21 = s;
            h22 = tst;
          } else if (kdefl - kdefl / 10 * 10 == 0) {
            s = fabs(h_data[(k + h_size[0] * k) + 1]) +
                fabs(h_data[(k + h_size[0] * (k + 1)) + 2]);
            tst = 0.75 * s + h_data[k + h_size[0] * k];
            h12 = -0.4375 * s;
            h21 = s;
            h22 = tst;
          } else {
            nr = b_i + h_size[0] * (b_i - 1);
            tst = h_data[nr - 1];
            h21 = h_data[nr];
            h12 = h_data[(b_i + h_size[0] * b_i) - 1];
            h22 = h_data[b_i + h_size[0] * b_i];
          }
          s = ((fabs(tst) + fabs(h12)) + fabs(h21)) + fabs(h22);
          if (s == 0.0) {
            rt1r = 0.0;
            h21 = 0.0;
            rt2r = 0.0;
            h12 = 0.0;
          } else {
            tst /= s;
            h21 /= s;
            h12 /= s;
            h22 /= s;
            aa = (tst + h22) / 2.0;
            tst = (tst - aa) * (h22 - aa) - h12 * h21;
            h21 = sqrt(fabs(tst));
            if (tst >= 0.0) {
              rt1r = aa * s;
              rt2r = rt1r;
              h21 *= s;
              h12 = -h21;
            } else {
              rt1r = aa + h21;
              rt2r = aa - h21;
              if (fabs(rt1r - h22) <= fabs(rt2r - h22)) {
                rt1r *= s;
                rt2r = rt1r;
              } else {
                rt2r *= s;
                rt1r = rt2r;
              }
              h21 = 0.0;
              h12 = 0.0;
            }
          }
          m = b_i - 1;
          exitg3 = false;
          while ((!exitg3) && (m >= k + 1)) {
            nr = m + h_size[0] * (m - 1);
            aa = h_data[nr - 1];
            h22 = aa - rt2r;
            tst = h_data[nr];
            s = (fabs(h22) + fabs(h12)) + fabs(tst);
            tst /= s;
            v_tmp = m + h_size[0] * m;
            v[0] =
                (tst * h_data[v_tmp - 1] + h22 * (h22 / s)) - h21 * (h12 / s);
            v[1] = tst * (((aa + h_data[v_tmp]) - rt1r) - rt2r);
            v[2] = tst * h_data[v_tmp + 1];
            s = (fabs(v[0]) + fabs(v[1])) + fabs(v[2]);
            v[0] /= s;
            v[1] /= s;
            v[2] /= s;
            if ((m == k + 1) ||
                (fabs(h_data[m - 1]) * (fabs(v[1]) + fabs(v[2])) <=
                 2.2204460492503131E-16 * fabs(v[0]) *
                     ((fabs(h_data[0]) + fabs(h_data[nr - 1])) +
                      fabs(h_data[v_tmp])))) {
              exitg3 = true;
            } else {
              m--;
            }
          }
          for (b_k = m; b_k <= b_i; b_k++) {
            nr = (b_i - b_k) + 2;
            if (nr >= 3) {
              nr = 3;
            }
            if (b_k > m) {
              v_tmp = ((b_k - 2) * *wr_size + b_k) - 1;
              for (c_k = 0; c_k < nr; c_k++) {
                v[c_k] = h_data[v_tmp + c_k];
              }
            }
            tst = v[0];
            aa = xzlarfg(nr, &tst, v);
            if (b_k > m) {
              h_data[(b_k + h_size[0] * (b_k - 2)) - 1] = tst;
              h_data[b_k + h_size[0] * (b_k - 2)] = 0.0;
              if (b_k < b_i) {
                h_data[b_k + 1] = 0.0;
              }
            } else if (m > k + 1) {
              h_data[b_k - 1] *= 1.0 - aa;
            }
            d = v[1];
            tst = aa * v[1];
            if (nr == 3) {
              h22 = v[2];
              h12 = aa * v[2];
              for (j = b_k; j <= b_i + 1; j++) {
                i = b_k + h_size[0] * (j - 1);
                rt2r = h_data[i - 1];
                rt1r = h_data[i];
                s = h_data[i + 1];
                h21 = (rt2r + d * rt1r) + h22 * s;
                rt2r -= h21 * aa;
                h_data[i - 1] = rt2r;
                rt1r -= h21 * tst;
                h_data[i] = rt1r;
                s -= h21 * h12;
                h_data[i + 1] = s;
              }
              if (b_k + 3 <= b_i + 1) {
                i = b_k;
              } else {
                i = b_i - 2;
              }
              for (j = k + 1; j <= i + 3; j++) {
                v_tmp = (j + h_size[0] * (b_k - 1)) - 1;
                rt2r = h_data[v_tmp];
                nr = (j + h_size[0] * b_k) - 1;
                rt1r = h_data[nr];
                c_k = (j + h_size[0] * (b_k + 1)) - 1;
                s = h_data[c_k];
                h21 = (rt2r + d * rt1r) + h22 * s;
                rt2r -= h21 * aa;
                h_data[v_tmp] = rt2r;
                rt1r -= h21 * tst;
                h_data[nr] = rt1r;
                s -= h21 * h12;
                h_data[c_k] = s;
              }
            } else if (nr == 2) {
              for (j = b_k; j <= b_i + 1; j++) {
                i = b_k + h_size[0] * (j - 1);
                h22 = h_data[i - 1];
                rt2r = h_data[i];
                h21 = h22 + d * rt2r;
                h22 -= h21 * aa;
                h_data[i - 1] = h22;
                rt2r -= h21 * tst;
                h_data[i] = rt2r;
              }
              for (j = k + 1; j <= b_i + 1; j++) {
                i = (j + h_size[0] * (b_k - 1)) - 1;
                h22 = h_data[i];
                v_tmp = (j + h_size[0] * b_k) - 1;
                rt2r = h_data[v_tmp];
                h21 = h22 + d * rt2r;
                h22 -= h21 * aa;
                h_data[i] = h22;
                rt2r -= h21 * tst;
                h_data[v_tmp] = rt2r;
              }
            }
          }
          its++;
        }
      }
      if (!converged) {
        info = b_i + 1;
        exitg1 = true;
      } else {
        if (l == b_i + 1) {
          wr_data[b_i] = h_data[b_i + h_size[0] * b_i];
          wi_data[b_i] = 0.0;
        } else if (l == b_i) {
          i = b_i + h_size[0] * b_i;
          d = h_data[i - 1];
          v_tmp = b_i + h_size[0] * (b_i - 1);
          h22 = h_data[v_tmp];
          rt2r = h_data[i];
          wr_data[b_i - 1] =
              xdlanv2(&h_data[(b_i + h_size[0] * (b_i - 1)) - 1], &d, &h22,
                      &rt2r, &wi_data[b_i - 1], &rt1r, &s, &tst, &h21);
          wr_data[b_i] = rt1r;
          wi_data[b_i] = s;
          h_data[i - 1] = d;
          h_data[v_tmp] = h22;
          h_data[i] = rt2r;
        }
        kdefl = 0;
        b_i = l - 2;
      }
    }
    if ((info != 0) && (*wr_size > 2)) {
      for (j = 3; j <= *wr_size; j++) {
        for (b_i = j; b_i <= *wr_size; b_i++) {
          h_data[(b_i + h_size[0] * (j - 3)) - 1] = 0.0;
        }
      }
    }
  }
  return info;
}

/*
 *
 */
static double xdlanv2(double *a, double *b, double *c, double *d, double *rt1i,
                      double *rt2r, double *rt2i, double *cs, double *sn)
{
  double rt1r;
  if (*c == 0.0) {
    *cs = 1.0;
    *sn = 0.0;
  } else if (*b == 0.0) {
    *cs = 0.0;
    *sn = 1.0;
    rt1r = *d;
    *d = *a;
    *a = rt1r;
    *b = -*c;
    *c = 0.0;
  } else {
    rt1r = *a - *d;
    if ((rt1r == 0.0) && ((*b < 0.0) != (*c < 0.0))) {
      *cs = 1.0;
      *sn = 0.0;
    } else {
      double bcmax;
      double bcmis;
      double p;
      double scale;
      double sigma;
      double z;
      int count;
      int i;
      p = 0.5 * rt1r;
      bcmis = fabs(*b);
      sigma = fabs(*c);
      bcmax = fmax(bcmis, sigma);
      if (!(*b < 0.0)) {
        count = 1;
      } else {
        count = -1;
      }
      if (!(*c < 0.0)) {
        i = 1;
      } else {
        i = -1;
      }
      bcmis = fmin(bcmis, sigma) * (double)count * (double)i;
      scale = fmax(fabs(p), bcmax);
      z = p / scale * p + bcmax / scale * bcmis;
      if (z >= 8.8817841970012523E-16) {
        *a = sqrt(scale) * sqrt(z);
        if (!(p < 0.0)) {
          rt1r = *a;
        } else {
          rt1r = -*a;
        }
        z = p + rt1r;
        *a = *d + z;
        *d -= bcmax / z * bcmis;
        scale = fabs(z);
        if (sigma < scale) {
          bcmax = sigma / scale;
          scale *= sqrt(bcmax * bcmax + 1.0);
        } else if (sigma > scale) {
          scale /= sigma;
          scale = sigma * sqrt(scale * scale + 1.0);
        } else if (rtIsNaN(scale)) {
          scale = rtNaN;
        } else {
          scale = sigma * 1.4142135623730951;
        }
        *cs = z / scale;
        *sn = *c / scale;
        *b -= *c;
        *c = 0.0;
      } else {
        sigma = *b + *c;
        scale = fmax(fabs(rt1r), fabs(sigma));
        count = 0;
        while ((scale >= 7.4428285367870146E+137) && (count <= 20)) {
          sigma *= 1.3435752215134178E-138;
          rt1r *= 1.3435752215134178E-138;
          scale = fmax(fabs(rt1r), fabs(sigma));
          count++;
        }
        while ((scale <= 1.3435752215134178E-138) && (count <= 20)) {
          sigma *= 7.4428285367870146E+137;
          rt1r *= 7.4428285367870146E+137;
          scale = fmax(fabs(rt1r), fabs(sigma));
          count++;
        }
        bcmis = fabs(sigma);
        scale = fabs(rt1r);
        if (bcmis < scale) {
          bcmax = bcmis / scale;
          scale *= sqrt(bcmax * bcmax + 1.0);
        } else if (bcmis > scale) {
          scale /= bcmis;
          scale = bcmis * sqrt(scale * scale + 1.0);
        } else if (rtIsNaN(scale)) {
          scale = rtNaN;
        } else {
          scale = bcmis * 1.4142135623730951;
        }
        *cs = sqrt(0.5 * (bcmis / scale + 1.0));
        if (!(sigma < 0.0)) {
          count = 1;
        } else {
          count = -1;
        }
        *sn = -(0.5 * rt1r / (scale * *cs)) * (double)count;
        sigma = *a * *cs + *b * *sn;
        bcmax = -*a * *sn + *b * *cs;
        scale = *c * *cs + *d * *sn;
        bcmis = -*c * *sn + *d * *cs;
        *b = bcmax * *cs + bcmis * *sn;
        *c = -sigma * *sn + scale * *cs;
        rt1r =
            0.5 * ((sigma * *cs + scale * *sn) + (-bcmax * *sn + bcmis * *cs));
        *a = rt1r;
        *d = rt1r;
        if (*c != 0.0) {
          if (*b != 0.0) {
            if ((*b < 0.0) == (*c < 0.0)) {
              bcmax = sqrt(fabs(*b));
              bcmis = sqrt(fabs(*c));
              *a = bcmax * bcmis;
              if (!(*c < 0.0)) {
                p = *a;
              } else {
                p = -*a;
              }
              scale = 1.0 / sqrt(fabs(*b + *c));
              *a = rt1r + p;
              *d = rt1r - p;
              *b -= *c;
              *c = 0.0;
              bcmax *= scale;
              bcmis *= scale;
              rt1r = *cs * bcmax - *sn * bcmis;
              *sn = *cs * bcmis + *sn * bcmax;
              *cs = rt1r;
            }
          } else {
            *b = -*c;
            *c = 0.0;
            rt1r = *cs;
            *cs = -*sn;
            *sn = rt1r;
          }
        }
      }
    }
  }
  rt1r = *a;
  *rt2r = *d;
  if (*c == 0.0) {
    *rt1i = 0.0;
    *rt2i = 0.0;
  } else {
    *rt1i = sqrt(fabs(*b)) * sqrt(fabs(*c));
    *rt2i = -*rt1i;
  }
  return rt1r;
}

/*
 *
 */
static double xnrm2(int n, const double x_data[], int ix0)
{
  double y;
  int k;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = fabs(x_data[ix0 - 1]);
    } else {
      double scale;
      int kend;
      scale = 3.3121686421112381E-170;
      kend = ix0 + n;
      for (k = ix0; k < kend; k++) {
        double absxk;
        absxk = fabs(x_data[k - 1]);
        if (absxk > scale) {
          double t;
          t = scale / absxk;
          y = y * t * t + 1.0;
          scale = absxk;
        } else {
          double t;
          t = absxk / scale;
          y += t * t;
        }
      }
      y = scale * sqrt(y);
    }
  }
  return y;
}

/*
 *
 */
static int xzgebal(double A_data[], const int A_size[2], int *ihi,
                   double scale_data[], int *scale_size)
{
  double scale;
  int b_k;
  int c_tmp;
  int exitg5;
  int i;
  int ilo;
  int ix;
  int iy;
  int k;
  int kend;
  int n_tmp;
  boolean_T converged;
  boolean_T notdone;
  *scale_size = A_size[0];
  for (i = 0; i < *scale_size; i++) {
    scale_data[i] = 1.0;
  }
  k = 0;
  *ihi = *scale_size;
  notdone = true;
  do {
    exitg5 = 0;
    if (notdone) {
      int exitg4;
      notdone = false;
      c_tmp = *ihi;
      do {
        exitg4 = 0;
        if (c_tmp > 0) {
          boolean_T exitg6;
          converged = false;
          ix = 0;
          exitg6 = false;
          while ((!exitg6) && (ix <= (unsigned char)*ihi - 1)) {
            if ((ix + 1 == c_tmp) ||
                (!(A_data[(c_tmp + A_size[0] * ix) - 1] != 0.0))) {
              ix++;
            } else {
              converged = true;
              exitg6 = true;
            }
          }
          if (converged) {
            c_tmp--;
          } else {
            scale_data[*ihi - 1] = c_tmp;
            if (c_tmp != *ihi) {
              ix = (c_tmp - 1) * *scale_size;
              iy = (*ihi - 1) * *scale_size;
              i = (unsigned char)*ihi;
              for (b_k = 0; b_k < i; b_k++) {
                ilo = ix + b_k;
                scale = A_data[ilo];
                n_tmp = iy + b_k;
                A_data[ilo] = A_data[n_tmp];
                A_data[n_tmp] = scale;
              }
              for (b_k = 0; b_k < *scale_size; b_k++) {
                ilo = b_k * *scale_size;
                n_tmp = (c_tmp + ilo) - 1;
                scale = A_data[n_tmp];
                i = (*ihi + ilo) - 1;
                A_data[n_tmp] = A_data[i];
                A_data[i] = scale;
              }
            }
            exitg4 = 1;
          }
        } else {
          exitg4 = 2;
        }
      } while (exitg4 == 0);
      if (exitg4 == 1) {
        if (*ihi == 1) {
          ilo = 1;
          *ihi = 1;
          exitg5 = 1;
        } else {
          (*ihi)--;
          notdone = true;
        }
      }
    } else {
      notdone = true;
      while (notdone) {
        boolean_T exitg6;
        notdone = false;
        c_tmp = k;
        exitg6 = false;
        while ((!exitg6) && (c_tmp + 1 <= *ihi)) {
          boolean_T exitg7;
          converged = false;
          ix = k;
          exitg7 = false;
          while ((!exitg7) && (ix + 1 <= *ihi)) {
            if ((ix + 1 == c_tmp + 1) ||
                (!(A_data[ix + A_size[0] * c_tmp] != 0.0))) {
              ix++;
            } else {
              converged = true;
              exitg7 = true;
            }
          }
          if (converged) {
            c_tmp++;
          } else {
            scale_data[k] = c_tmp + 1;
            if (c_tmp + 1 != k + 1) {
              ix = c_tmp * *scale_size;
              kend = k * *scale_size;
              i = (unsigned char)*ihi;
              for (b_k = 0; b_k < i; b_k++) {
                ilo = ix + b_k;
                scale = A_data[ilo];
                n_tmp = kend + b_k;
                A_data[ilo] = A_data[n_tmp];
                A_data[n_tmp] = scale;
              }
              ix = kend + c_tmp;
              iy = kend + k;
              kend = *scale_size - k;
              for (b_k = 0; b_k < kend; b_k++) {
                ilo = b_k * *scale_size;
                n_tmp = ix + ilo;
                scale = A_data[n_tmp];
                i = iy + ilo;
                A_data[n_tmp] = A_data[i];
                A_data[i] = scale;
              }
            }
            k++;
            notdone = true;
            exitg6 = true;
          }
        }
      }
      ilo = k + 1;
      converged = false;
      exitg5 = 2;
    }
  } while (exitg5 == 0);
  if (exitg5 != 1) {
    boolean_T exitg3;
    exitg3 = false;
    while ((!exitg3) && (!converged)) {
      int exitg2;
      converged = true;
      ix = k;
      do {
        exitg2 = 0;
        if (ix + 1 <= *ihi) {
          double absxk;
          double c;
          double ca;
          double r;
          double t;
          kend = *ihi - k;
          c_tmp = ix * *scale_size;
          c = xnrm2(kend, A_data, (c_tmp + k) + 1);
          iy = k * *scale_size + ix;
          r = 0.0;
          if (kend >= 1) {
            if (kend == 1) {
              r = fabs(A_data[iy]);
            } else {
              scale = 3.3121686421112381E-170;
              kend = (iy + (kend - 1) * *scale_size) + 1;
              for (b_k = iy + 1; *scale_size < 0 ? b_k >= kend : b_k <= kend;
                   b_k += *scale_size) {
                absxk = fabs(A_data[b_k - 1]);
                if (absxk > scale) {
                  t = scale / absxk;
                  r = r * t * t + 1.0;
                  scale = absxk;
                } else {
                  t = absxk / scale;
                  r += t * t;
                }
              }
              r = scale * sqrt(r);
            }
          }
          kend = 1;
          if (*ihi > 1) {
            scale = fabs(A_data[c_tmp]);
            for (b_k = 2; b_k <= *ihi; b_k++) {
              t = fabs(A_data[(c_tmp + b_k) - 1]);
              if (t > scale) {
                kend = b_k;
                scale = t;
              }
            }
          }
          ca = fabs(A_data[(kend + A_size[0] * ix) - 1]);
          n_tmp = *scale_size - k;
          if (n_tmp < 1) {
            kend = 0;
          } else {
            kend = 1;
            if (n_tmp > 1) {
              scale = fabs(A_data[iy]);
              for (b_k = 2; b_k <= n_tmp; b_k++) {
                t = fabs(A_data[iy + (b_k - 1) * *scale_size]);
                if (t > scale) {
                  kend = b_k;
                  scale = t;
                }
              }
            }
          }
          scale = fabs(A_data[ix + A_size[0] * ((kend + k) - 1)]);
          if ((c == 0.0) || (r == 0.0)) {
            ix++;
          } else {
            double f;
            int exitg1;
            absxk = r / 2.0;
            f = 1.0;
            t = c + r;
            do {
              exitg1 = 0;
              if ((c < absxk) &&
                  (fmax(f, fmax(c, ca)) < 4.9896007738368E+291) &&
                  (fmin(r, fmin(absxk, scale)) > 2.0041683600089728E-292)) {
                if (rtIsNaN(((((c + f) + ca) + r) + absxk) + scale)) {
                  exitg1 = 1;
                } else {
                  f *= 2.0;
                  c *= 2.0;
                  ca *= 2.0;
                  r /= 2.0;
                  absxk /= 2.0;
                  scale /= 2.0;
                }
              } else {
                absxk = c / 2.0;
                while ((absxk >= r) &&
                       (fmax(r, scale) < 4.9896007738368E+291) &&
                       (fmin(fmin(f, c), fmin(absxk, ca)) >
                        2.0041683600089728E-292)) {
                  f /= 2.0;
                  c /= 2.0;
                  absxk /= 2.0;
                  ca /= 2.0;
                  r *= 2.0;
                  scale *= 2.0;
                }
                if ((!(c + r >= 0.95 * t)) &&
                    ((!(f < 1.0)) || (!(scale_data[ix] < 1.0)) ||
                     (!(f * scale_data[ix] <= 1.0020841800044864E-292))) &&
                    ((!(f > 1.0)) || (!(scale_data[ix] > 1.0)) ||
                     (!(scale_data[ix] >= 9.9792015476736E+291 / f)))) {
                  scale = 1.0 / f;
                  scale_data[ix] *= f;
                  kend = iy + 1;
                  i = (iy + *scale_size * (n_tmp - 1)) + 1;
                  for (b_k = kend; *scale_size < 0 ? b_k >= i : b_k <= i;
                       b_k += *scale_size) {
                    A_data[b_k - 1] *= scale;
                  }
                  i = c_tmp + *ihi;
                  for (b_k = c_tmp + 1; b_k <= i; b_k++) {
                    A_data[b_k - 1] *= f;
                  }
                  converged = false;
                }
                exitg1 = 2;
              }
            } while (exitg1 == 0);
            if (exitg1 == 1) {
              exitg2 = 2;
            } else {
              ix++;
            }
          }
        } else {
          exitg2 = 1;
        }
      } while (exitg2 == 0);
      if (exitg2 != 1) {
        exitg3 = true;
      }
    }
  }
  return ilo;
}

/*
 *
 */
static void xzgehrd(double a_data[], const int a_size[2], int ilo, int ihi)
{
  double work_data[4];
  double tau_data[3];
  int b_i;
  int c_i;
  int ia;
  int im1n_tmp;
  int k;
  int n_tmp;
  n_tmp = a_size[0];
  if ((ihi - ilo) + 1 > 1) {
    int i;
    i = (unsigned char)(ilo - 1);
    if (i - 1 >= 0) {
      memset(&tau_data[0], 0, (unsigned int)i * sizeof(double));
    }
    if (ihi <= n_tmp - 1) {
      memset(&tau_data[ihi + -1], 0,
             (unsigned int)(n_tmp - ihi) * sizeof(double));
    }
    if (n_tmp - 1 >= 0) {
      memset(&work_data[0], 0, (unsigned int)n_tmp * sizeof(double));
    }
    for (b_i = ilo; b_i < ihi; b_i++) {
      double a;
      double alpha1;
      double beta1;
      int alpha1_tmp;
      int exitg1;
      int i1;
      int in;
      int iv0;
      int jA;
      int knt;
      int lastc;
      int lastv;
      int n;
      boolean_T exitg2;
      im1n_tmp = (b_i - 1) * n_tmp;
      in = b_i * n_tmp;
      alpha1_tmp = b_i + a_size[0] * (b_i - 1);
      alpha1 = a_data[alpha1_tmp];
      c_i = b_i + 2;
      if (c_i > n_tmp) {
        c_i = n_tmp;
      }
      c_i += im1n_tmp;
      n = (ihi - b_i) - 1;
      tau_data[b_i - 1] = 0.0;
      if (n + 1 > 0) {
        beta1 = xnrm2(n, a_data, c_i);
        if (beta1 != 0.0) {
          a = fabs(alpha1);
          beta1 = fabs(beta1);
          if (a < beta1) {
            a /= beta1;
            beta1 *= sqrt(a * a + 1.0);
          } else if (a > beta1) {
            beta1 /= a;
            beta1 = a * sqrt(beta1 * beta1 + 1.0);
          } else if (rtIsNaN(beta1)) {
            beta1 = rtNaN;
          } else {
            beta1 = a * 1.4142135623730951;
          }
          if (alpha1 >= 0.0) {
            beta1 = -beta1;
          }
          if (fabs(beta1) < 1.0020841800044864E-292) {
            knt = 0;
            i = c_i + n;
            do {
              knt++;
              for (k = c_i; k < i; k++) {
                a_data[k - 1] *= 9.9792015476736E+291;
              }
              beta1 *= 9.9792015476736E+291;
              alpha1 *= 9.9792015476736E+291;
            } while ((fabs(beta1) < 1.0020841800044864E-292) && (knt < 20));
            a = fabs(alpha1);
            beta1 = fabs(xnrm2(n, a_data, c_i));
            if (a < beta1) {
              a /= beta1;
              beta1 *= sqrt(a * a + 1.0);
            } else if (a > beta1) {
              beta1 /= a;
              beta1 = a * sqrt(beta1 * beta1 + 1.0);
            } else if (rtIsNaN(beta1)) {
              beta1 = rtNaN;
            } else {
              beta1 = a * 1.4142135623730951;
            }
            if (alpha1 >= 0.0) {
              beta1 = -beta1;
            }
            tau_data[b_i - 1] = (beta1 - alpha1) / beta1;
            a = 1.0 / (alpha1 - beta1);
            for (k = c_i; k < i; k++) {
              a_data[k - 1] *= a;
            }
            for (k = 0; k < knt; k++) {
              beta1 *= 1.0020841800044864E-292;
            }
            alpha1 = beta1;
          } else {
            tau_data[b_i - 1] = (beta1 - alpha1) / beta1;
            a = 1.0 / (alpha1 - beta1);
            i = c_i + n;
            for (k = c_i; k < i; k++) {
              a_data[k - 1] *= a;
            }
            alpha1 = beta1;
          }
        }
      }
      a_data[alpha1_tmp] = 1.0;
      iv0 = b_i + im1n_tmp;
      im1n_tmp = in + 1;
      a = tau_data[b_i - 1];
      if (a != 0.0) {
        lastv = n;
        c_i = iv0 + n;
        while ((lastv + 1 > 0) && (a_data[c_i] == 0.0)) {
          lastv--;
          c_i--;
        }
        lastc = ihi;
        exitg2 = false;
        while ((!exitg2) && (lastc > 0)) {
          knt = in + lastc;
          ia = knt;
          do {
            exitg1 = 0;
            if ((n_tmp > 0) && (ia <= knt + lastv * n_tmp)) {
              if (a_data[ia - 1] != 0.0) {
                exitg1 = 1;
              } else {
                ia += n_tmp;
              }
            } else {
              lastc--;
              exitg1 = 2;
            }
          } while (exitg1 == 0);
          if (exitg1 == 1) {
            exitg2 = true;
          }
        }
      } else {
        lastv = -1;
        lastc = 0;
      }
      if (lastv + 1 > 0) {
        if (lastc != 0) {
          i = (unsigned char)lastc;
          memset(&work_data[0], 0, (unsigned int)i * sizeof(double));
          knt = iv0;
          i = (in + n_tmp * lastv) + 1;
          for (k = im1n_tmp; n_tmp < 0 ? k >= i : k <= i; k += n_tmp) {
            i1 = k + lastc;
            for (ia = k; ia < i1; ia++) {
              c_i = ia - k;
              work_data[c_i] += a_data[ia - 1] * a_data[knt];
            }
            knt++;
          }
        }
        if (!(-a == 0.0)) {
          jA = in;
          i = (unsigned char)(lastv + 1);
          for (c_i = 0; c_i < i; c_i++) {
            beta1 = a_data[iv0 + c_i];
            if (beta1 != 0.0) {
              beta1 *= -a;
              i1 = jA + 1;
              knt = lastc + jA;
              for (im1n_tmp = i1; im1n_tmp <= knt; im1n_tmp++) {
                a_data[im1n_tmp - 1] += work_data[(im1n_tmp - jA) - 1] * beta1;
              }
            }
            jA += n_tmp;
          }
        }
      }
      jA = (b_i + in) + 1;
      if (a != 0.0) {
        lastv = n + 1;
        c_i = iv0 + n;
        while ((lastv > 0) && (a_data[c_i] == 0.0)) {
          lastv--;
          c_i--;
        }
        lastc = (n_tmp - b_i) - 1;
        exitg2 = false;
        while ((!exitg2) && (lastc + 1 > 0)) {
          c_i = jA + lastc * n_tmp;
          ia = c_i;
          do {
            exitg1 = 0;
            if (ia <= (c_i + lastv) - 1) {
              if (a_data[ia - 1] != 0.0) {
                exitg1 = 1;
              } else {
                ia++;
              }
            } else {
              lastc--;
              exitg1 = 2;
            }
          } while (exitg1 == 0);
          if (exitg1 == 1) {
            exitg2 = true;
          }
        }
      } else {
        lastv = 0;
        lastc = -1;
      }
      if (lastv > 0) {
        if (lastc + 1 != 0) {
          if (lastc >= 0) {
            memset(&work_data[0], 0,
                   (unsigned int)(lastc + 1) * sizeof(double));
          }
          knt = 0;
          i = jA + n_tmp * lastc;
          for (k = jA; n_tmp < 0 ? k >= i : k <= i; k += n_tmp) {
            beta1 = 0.0;
            i1 = k + lastv;
            for (ia = k; ia < i1; ia++) {
              beta1 += a_data[ia - 1] * a_data[(iv0 + ia) - k];
            }
            work_data[knt] += beta1;
            knt++;
          }
        }
        if (!(-a == 0.0)) {
          for (c_i = 0; c_i <= lastc; c_i++) {
            beta1 = work_data[c_i];
            if (beta1 != 0.0) {
              beta1 *= -a;
              i = lastv + jA;
              for (im1n_tmp = jA; im1n_tmp < i; im1n_tmp++) {
                a_data[im1n_tmp - 1] += a_data[(iv0 + im1n_tmp) - jA] * beta1;
              }
            }
            jA += n_tmp;
          }
        }
      }
      a_data[alpha1_tmp] = alpha1;
    }
  }
}

/*
 *
 */
static double xzlarfg(int n, double *alpha1, double x[3])
{
  double tau;
  int k;
  tau = 0.0;
  if (n > 0) {
    double xnorm;
    xnorm = b_xnrm2(n - 1, x);
    if (xnorm != 0.0) {
      double a_tmp;
      a_tmp = fabs(*alpha1);
      xnorm = fabs(xnorm);
      if (a_tmp < xnorm) {
        a_tmp /= xnorm;
        xnorm *= sqrt(a_tmp * a_tmp + 1.0);
      } else if (a_tmp > xnorm) {
        xnorm /= a_tmp;
        xnorm = a_tmp * sqrt(xnorm * xnorm + 1.0);
      } else if (rtIsNaN(xnorm)) {
        xnorm = rtNaN;
      } else {
        xnorm = a_tmp * 1.4142135623730951;
      }
      if (*alpha1 >= 0.0) {
        xnorm = -xnorm;
      }
      if (fabs(xnorm) < 1.0020841800044864E-292) {
        int knt;
        knt = 0;
        do {
          knt++;
          for (k = 2; k <= n; k++) {
            x[k - 1] *= 9.9792015476736E+291;
          }
          xnorm *= 9.9792015476736E+291;
          *alpha1 *= 9.9792015476736E+291;
        } while ((fabs(xnorm) < 1.0020841800044864E-292) && (knt < 20));
        a_tmp = fabs(*alpha1);
        xnorm = fabs(b_xnrm2(n - 1, x));
        if (a_tmp < xnorm) {
          a_tmp /= xnorm;
          xnorm *= sqrt(a_tmp * a_tmp + 1.0);
        } else if (a_tmp > xnorm) {
          xnorm /= a_tmp;
          xnorm = a_tmp * sqrt(xnorm * xnorm + 1.0);
        } else if (rtIsNaN(xnorm)) {
          xnorm = rtNaN;
        } else {
          xnorm = a_tmp * 1.4142135623730951;
        }
        if (*alpha1 >= 0.0) {
          xnorm = -xnorm;
        }
        tau = (xnorm - *alpha1) / xnorm;
        a_tmp = 1.0 / (*alpha1 - xnorm);
        for (k = 2; k <= n; k++) {
          x[k - 1] *= a_tmp;
        }
        for (k = 0; k < knt; k++) {
          xnorm *= 1.0020841800044864E-292;
        }
        *alpha1 = xnorm;
      } else {
        tau = (xnorm - *alpha1) / xnorm;
        a_tmp = 1.0 / (*alpha1 - xnorm);
        for (k = 2; k <= n; k++) {
          x[k - 1] *= a_tmp;
        }
        *alpha1 = xnorm;
      }
    }
  }
  return tau;
}

/*
 *
 */
static void xzlascl(double cfrom, double cto, int m, double A_data[], int iA0)
{
  double cfromc;
  double ctoc;
  int i;
  boolean_T notdone;
  cfromc = cfrom;
  ctoc = cto;
  notdone = true;
  while (notdone) {
    double cfrom1;
    double cto1;
    double mul;
    cfrom1 = cfromc * 2.0041683600089728E-292;
    cto1 = ctoc / 4.9896007738368E+291;
    if ((fabs(cfrom1) > fabs(ctoc)) && (ctoc != 0.0)) {
      mul = 2.0041683600089728E-292;
      cfromc = cfrom1;
    } else if (fabs(cto1) > fabs(cfromc)) {
      mul = 4.9896007738368E+291;
      ctoc = cto1;
    } else {
      mul = ctoc / cfromc;
      notdone = false;
    }
    for (i = 0; i < m; i++) {
      int b_i;
      b_i = (iA0 + i) - 1;
      A_data[b_i] *= mul;
    }
  }
}

void primitives_initialize(void)
{
}

void primitives_terminate(void)
{
}

/*
 * function [coeffsT2, v2, T2, coeffsT1, v1, T1] = student_pass_primitive(v0,
 * a0, sf, vfmin, vfmax, Tmin, Tmax)
 */
void student_pass_primitive(double v0, double a0, double sf, double vfmin,
                            double vfmax, double Tmin, double Tmax,
                            creal_T coeffsT2_data[], int coeffsT2_size[2],
                            creal_T v2_data[], int v2_size[1],
                            creal_T T2_data[], int T2_size[1],
                            creal_T coeffsT1_data[], int coeffsT1_size[2],
                            creal_T v1_data[], int v1_size[1],
                            creal_T T1_data[], int T1_size[1])
{
  creal_T dcv[6];
  creal_T Tvmax_data[4];
  creal_T Tvmin_data[4];
  creal_T varargin_2;
  double t7;
  double t7_tmp;
  double v_star;
  int Tvmax_size;
  int Tvmin_size;
  int i;
  boolean_T b_Tvmax_data[4];
  boolean_T guard1;
  /*  */
  /*            Agent Logic */
  /*          Pass Primitive */
  /*               2025 */
  /*  */
  /* 'student_pass_primitive:11' if a0 >= 0 */
  if (a0 >= 0.0) {
    /* 'student_pass_primitive:12' Tvmin = final_opt_time_pass(v0, a0, sf,
     * vfmin); */
    final_opt_time_pass(v0, a0, sf, vfmin, Tvmax_data);
    Tvmin_size = 4;
    memcpy(&Tvmin_data[0], &Tvmax_data[0], 4U * sizeof(creal_T));
    /* 'student_pass_primitive:13' Tvmax = final_opt_time_pass(v0, a0, sf,
     * vfmax); */
    final_opt_time_pass(v0, a0, sf, vfmax, Tvmax_data);
    Tvmax_size = 4;
  } else {
    /* 'student_pass_primitive:14' else */
    /* 'student_pass_primitive:15' T_star = time_min_vel(a0, sf); */
    /* TIME_MIN_VEL */
    /*     TIME_MIN_VEL_VAR = TIME_MIN_VEL(A0,SF) */
    /*     This function was generated by the Symbolic Math Toolbox
     * version 24.2. */
    /*     24-Oct-2025 11:35:05 */
    /* 'time_min_vel:8' t2 = a0.*sf; */
    /* 'time_min_vel:9' t3 = 1.0./a0; */
    /* 'time_min_vel:10' t5 = sqrt(1.5e+1); */
    /* 'time_min_vel:11' t4 = -t2; */
    /* 'time_min_vel:12' t6 = sqrt(t4); */
    /* 'time_min_vel:13' t7 = t3.*t5.*t6; */
    v_star = a0 * sf;
    t7_tmp = sqrt(-v_star);
    t7 = 1.0 / a0 * 3.872983346207417 * t7_tmp;
    /* 'time_min_vel:14' time_min_vel_var = [t7;-t7]; */
    /* 'student_pass_primitive:16' v_star = min_vel(v0, a0, sf); */
    /* MIN_VEL */
    /*     MIN_VEL_VAR = MIN_VEL(V0,A0,SF) */
    /*     This function was generated by the Symbolic Math Toolbox
     * version 24.2. */
    /*     24-Oct-2025 11:35:05 */
    /* 'min_vel:8' t2 = a0.*sf; */
    /* 'min_vel:9' t3 = -t2; */
    /* 'min_vel:10' min_vel_var
     * = 1.0./sqrt(t3).*(sqrt(1.5e+1).*t2.*2.0+sqrt(t3).*v0.*7.0).*(-1.0./8.0);
     */
    v_star = 1.0 / t7_tmp *
             (3.872983346207417 * v_star * 2.0 + t7_tmp * v0 * 7.0) * -0.125;
    /* 'student_pass_primitive:17' if v_star < vfmin && vfmin < vfmax */
    if ((v_star < vfmin) && (vfmin < vfmax)) {
      /* 'student_pass_primitive:18' Tvmin = final_opt_time_pass(v0, a0, sf,
       * vfmin); */
      final_opt_time_pass(v0, a0, sf, vfmin, Tvmax_data);
      Tvmin_size = 4;
      memcpy(&Tvmin_data[0], &Tvmax_data[0], 4U * sizeof(creal_T));
      /* 'student_pass_primitive:19' Tvmax = final_opt_time_pass(v0, a0, sf,
       * vfmax); */
      final_opt_time_pass(v0, a0, sf, vfmax, Tvmax_data);
      Tvmax_size = 4;
    } else if ((vfmin < v_star) && (v_star < vfmax)) {
      creal_T b_t7[2];
      /* 'student_pass_primitive:20' elseif vfmin < v_star && v_star < vfmax */
      /* 'student_pass_primitive:21' Tvmin = T_star; */
      b_t7[0].re = t7;
      b_t7[0].im = 0.0;
      b_t7[1].re = -t7;
      b_t7[1].im = 0.0;
      Tvmin_size = 2;
      Tvmin_data[0] = b_t7[0];
      Tvmin_data[1] = b_t7[1];
      /* 'student_pass_primitive:22' Tvmax = final_opt_time_pass(v0, a0, sf,
       * vfmax); */
      final_opt_time_pass(v0, a0, sf, vfmax, Tvmax_data);
      Tvmax_size = 4;
    } else {
      /* 'student_pass_primitive:23' else */
      /* 'student_pass_primitive:24' Tvmin = 0.; */
      Tvmin_size = 1;
      Tvmin_data[0].re = 0.0;
      Tvmin_data[0].im = 0.0;
      /* 'student_pass_primitive:25' Tvmax = 0.; */
      Tvmax_size = 1;
      Tvmax_data[0].re = 0.0;
      Tvmax_data[0].im = 0.0;
    }
  }
  /*  FREE-FLOW PRIMITIVE */
  /* 'student_pass_primitive:30' if Tmin == 0. && Tmax == 0. */
  if ((Tmin == 0.0) && (Tmax == 0.0)) {
    /* 'student_pass_primitive:31' T1 = Tvmin; */
    T1_size[0] = Tvmin_size;
    memcpy(&T1_data[0], &Tvmin_data[0],
           (unsigned int)Tvmin_size * sizeof(creal_T));
    /* 'student_pass_primitive:32' T2 = Tvmax; */
    T2_size[0] = Tvmax_size;
    memcpy(&T2_data[0], &Tvmax_data[0],
           (unsigned int)Tvmax_size * sizeof(creal_T));
    /*  FREE-FLOW PRIMITIVE */
  } else {
    boolean_T p;
    /* 'student_pass_primitive:34' else */
    /* 'student_pass_primitive:35' T1 = max(Tmin, Tvmax); */
    T1_size[0] = Tvmax_size;
    for (i = 0; i < Tvmax_size; i++) {
      varargin_2 = Tvmax_data[i];
      if (rtIsNaN(varargin_2.re) || rtIsNaN(varargin_2.im)) {
        p = false;
      } else if (rtIsNaN(Tmin)) {
        p = true;
      } else {
        v_star = absRelopProxies(Tmin, varargin_2, &t7_tmp);
        p = (v_star < t7_tmp);
      }
      if (p) {
        T1_data[i] = varargin_2;
      } else {
        T1_data[i].re = Tmin;
        T1_data[i].im = 0.0;
      }
    }
    /* 'student_pass_primitive:36' T2 = min(Tmax, Tvmin); */
    T2_size[0] = Tvmin_size;
    for (i = 0; i < Tvmin_size; i++) {
      varargin_2 = Tvmin_data[i];
      if (rtIsNaN(varargin_2.re) || rtIsNaN(varargin_2.im)) {
        p = false;
      } else if (rtIsNaN(Tmax)) {
        p = true;
      } else {
        v_star = absRelopProxies(Tmax, varargin_2, &t7_tmp);
        p = (v_star > t7_tmp);
      }
      if (p) {
        T2_data[i] = varargin_2;
      } else {
        T2_data[i].re = Tmax;
        T2_data[i].im = 0.0;
      }
    }
  }
  /* 'student_pass_primitive:39' if Tvmax ~= 0. & Tmax & Tvmax <= Tvmin & T1 >
   * 0. & T1 <= T2 */
  for (i = 0; i < Tvmax_size; i++) {
    b_Tvmax_data[i] = ((Tvmax_data[i].re != 0.0) || (Tvmax_data[i].im != 0.0));
  }
  guard1 = false;
  if (ifWhileCond(b_Tvmax_data, Tvmax_size) && (Tmax != 0.0)) {
    int loop_ub;
    int stride_0_0;
    if (Tvmin_size == 1) {
      loop_ub = Tvmax_size;
    } else {
      loop_ub = Tvmin_size;
    }
    stride_0_0 = (Tvmax_size != 1);
    Tvmax_size = (Tvmin_size != 1);
    for (i = 0; i < loop_ub; i++) {
      b_Tvmax_data[i] =
          (Tvmax_data[i * stride_0_0].re <= Tvmin_data[i * Tvmax_size].re);
    }
    if (ifWhileCond(b_Tvmax_data, loop_ub)) {
      loop_ub = T1_size[0];
      for (i = 0; i < loop_ub; i++) {
        b_Tvmax_data[i] = (T1_data[i].re > 0.0);
      }
      if (ifWhileCond(b_Tvmax_data, T1_size[0])) {
        if (T2_size[0] == 1) {
          Tvmin_size = T1_size[0];
        } else {
          Tvmin_size = T2_size[0];
        }
        stride_0_0 = (T1_size[0] != 1);
        Tvmax_size = (T2_size[0] != 1);
        for (i = 0; i < Tvmin_size; i++) {
          b_Tvmax_data[i] =
              (T1_data[i * stride_0_0].re <= T2_data[i * Tvmax_size].re);
        }
        if (ifWhileCond(b_Tvmax_data, Tvmin_size)) {
          creal_T b_tmp_data[4];
          creal_T t4_data[4];
          creal_T t4_tmp_data[4];
          creal_T t5_data[4];
          creal_T tmp_data[4];
          double b_sf_re;
          double brm;
          double c_sf_re;
          double re;
          double s;
          double sf_re;
          double varargin_1_im;
          double varargin_1_re;
          /* 'student_pass_primitive:40' v1 = final_opt_vel_pass(v0, a0, sf,
           * T1); */
          /* FINAL_OPT_VEL_PASS */
          /*     FINAL_OPT_VEL_PASS_VAR = FINAL_OPT_VEL_PASS(V0,A0,SF,T) */
          /*     This function was generated by the Symbolic Math Toolbox
           * version 24.2. */
          /*     24-Oct-2025 11:35:04 */
          /* 'final_opt_vel_pass:8' final_opt_vel_pass_var =
           * ((sf.*-1.5e+1+T.*v0.*7.0+T.^2.*a0).*(-1.0./8.0))./T; */
          t7 = sf * -15.0;
          v1_size[0] = T1_size[0];
          for (i = 0; i < loop_ub; i++) {
            varargin_1_re = T1_data[i].re;
            varargin_1_im = T1_data[i].im;
            v_star = varargin_1_re * varargin_1_im;
            sf_re = ((t7 + v0 * varargin_1_re * 7.0) +
                     (varargin_1_re * varargin_1_re -
                      varargin_1_im * varargin_1_im) *
                         a0) *
                    -0.125;
            t7_tmp =
                (v0 * varargin_1_im * 7.0 + (v_star + v_star) * a0) * -0.125;
            if (varargin_1_im == 0.0) {
              if (t7_tmp == 0.0) {
                v1_data[i].re = sf_re / varargin_1_re;
                v1_data[i].im = 0.0;
              } else if (sf_re == 0.0) {
                v1_data[i].re = 0.0;
                v1_data[i].im = t7_tmp / varargin_1_re;
              } else {
                v1_data[i].re = sf_re / varargin_1_re;
                v1_data[i].im = t7_tmp / varargin_1_re;
              }
            } else if (varargin_1_re == 0.0) {
              if (sf_re == 0.0) {
                v1_data[i].re = t7_tmp / varargin_1_im;
                v1_data[i].im = 0.0;
              } else if (t7_tmp == 0.0) {
                v1_data[i].re = 0.0;
                v1_data[i].im = -(sf_re / varargin_1_im);
              } else {
                v1_data[i].re = t7_tmp / varargin_1_im;
                v1_data[i].im = -(sf_re / varargin_1_im);
              }
            } else {
              brm = fabs(varargin_1_re);
              v_star = fabs(varargin_1_im);
              if (brm > v_star) {
                s = varargin_1_im / varargin_1_re;
                v_star = varargin_1_re + s * varargin_1_im;
                v1_data[i].re = (sf_re + s * t7_tmp) / v_star;
                v1_data[i].im = (t7_tmp - s * sf_re) / v_star;
              } else if (v_star == brm) {
                if (varargin_1_re > 0.0) {
                  s = 0.5;
                } else {
                  s = -0.5;
                }
                if (varargin_1_im > 0.0) {
                  v_star = 0.5;
                } else {
                  v_star = -0.5;
                }
                v1_data[i].re = (sf_re * s + t7_tmp * v_star) / brm;
                v1_data[i].im = (t7_tmp * s - sf_re * v_star) / brm;
              } else {
                s = varargin_1_re / varargin_1_im;
                v_star = varargin_1_im + s * varargin_1_re;
                v1_data[i].re = (s * sf_re + t7_tmp) / v_star;
                v1_data[i].im = (s * t7_tmp - sf_re) / v_star;
              }
            }
          }
          /* 'student_pass_primitive:41' v2 = final_opt_vel_pass(v0, a0, sf,
           * T2); */
          /* FINAL_OPT_VEL_PASS */
          /*     FINAL_OPT_VEL_PASS_VAR = FINAL_OPT_VEL_PASS(V0,A0,SF,T) */
          /*     This function was generated by the Symbolic Math Toolbox
           * version 24.2. */
          /*     24-Oct-2025 11:35:04 */
          /* 'final_opt_vel_pass:8' final_opt_vel_pass_var =
           * ((sf.*-1.5e+1+T.*v0.*7.0+T.^2.*a0).*(-1.0./8.0))./T; */
          Tvmin_size = T2_size[0];
          v2_size[0] = T2_size[0];
          for (i = 0; i < Tvmin_size; i++) {
            varargin_1_re = T2_data[i].re;
            varargin_1_im = T2_data[i].im;
            v_star = varargin_1_re * varargin_1_im;
            sf_re = ((t7 + v0 * varargin_1_re * 7.0) +
                     (varargin_1_re * varargin_1_re -
                      varargin_1_im * varargin_1_im) *
                         a0) *
                    -0.125;
            t7_tmp =
                (v0 * varargin_1_im * 7.0 + (v_star + v_star) * a0) * -0.125;
            if (varargin_1_im == 0.0) {
              if (t7_tmp == 0.0) {
                v2_data[i].re = sf_re / varargin_1_re;
                v2_data[i].im = 0.0;
              } else if (sf_re == 0.0) {
                v2_data[i].re = 0.0;
                v2_data[i].im = t7_tmp / varargin_1_re;
              } else {
                v2_data[i].re = sf_re / varargin_1_re;
                v2_data[i].im = t7_tmp / varargin_1_re;
              }
            } else if (varargin_1_re == 0.0) {
              if (sf_re == 0.0) {
                v2_data[i].re = t7_tmp / varargin_1_im;
                v2_data[i].im = 0.0;
              } else if (t7_tmp == 0.0) {
                v2_data[i].re = 0.0;
                v2_data[i].im = -(sf_re / varargin_1_im);
              } else {
                v2_data[i].re = t7_tmp / varargin_1_im;
                v2_data[i].im = -(sf_re / varargin_1_im);
              }
            } else {
              brm = fabs(varargin_1_re);
              v_star = fabs(varargin_1_im);
              if (brm > v_star) {
                s = varargin_1_im / varargin_1_re;
                v_star = varargin_1_re + s * varargin_1_im;
                v2_data[i].re = (sf_re + s * t7_tmp) / v_star;
                v2_data[i].im = (t7_tmp - s * sf_re) / v_star;
              } else if (v_star == brm) {
                if (varargin_1_re > 0.0) {
                  s = 0.5;
                } else {
                  s = -0.5;
                }
                if (varargin_1_im > 0.0) {
                  v_star = 0.5;
                } else {
                  v_star = -0.5;
                }
                v2_data[i].re = (sf_re * s + t7_tmp * v_star) / brm;
                v2_data[i].im = (t7_tmp * s - sf_re * v_star) / brm;
              } else {
                s = varargin_1_re / varargin_1_im;
                v_star = varargin_1_im + s * varargin_1_re;
                v2_data[i].re = (s * sf_re + t7_tmp) / v_star;
                v2_data[i].im = (s * t7_tmp - sf_re) / v_star;
              }
            }
          }
          /* 'student_pass_primitive:42' coeffsT1 = coef_list_fun(v0, a0, sf,
           * v1, 0., T1); */
          /* COEF_LIST_FUN */
          /*     COEF_LIST_VAR = COEF_LIST_FUN(V0,A0,SF,VF,AF,T) */
          /*     This function was generated by the Symbolic Math Toolbox
           * version 24.2. */
          /*     24-Oct-2025 11:34:51 */
          /* 'coef_list_fun:8' t2 = T.^2; */
          /* 'coef_list_fun:9' t3 = af.*t2; */
          /* 'coef_list_fun:10' t4 = a0.*t2.*3.0; */
          /* 'coef_list_fun:11' t5 = -t3; */
          /* 'coef_list_fun:12' coef_list_var =
           * [0.0,v0,a0,1.0./T.^3.*(sf.*-2.0e+1+t4+t5+T.*v0.*1.2e+1+T.*vf.*8.0).*-3.0,1.0./t2.^2.*(sf.*-3.0e+1-t3.*2.0+t4+T.*v0.*1.6e+1+T.*vf.*1.4e+1).*1.2e+1,1.0./T.^5.*(sf.*-1.2e+1+t5+T.*v0.*6.0+T.*vf.*6.0+a0.*t2).*-6.0e+1];
           */
          for (i = 0; i < loop_ub; i++) {
            varargin_1_re = T1_data[i].re;
            varargin_1_im = T1_data[i].im;
            sf_re =
                varargin_1_re * varargin_1_re - varargin_1_im * varargin_1_im;
            Tvmin_data[i].re = sf_re;
            v_star = varargin_1_re * varargin_1_im;
            v_star += v_star;
            Tvmin_data[i].im = v_star;
            t7_tmp = 0.0 * sf_re;
            Tvmax_data[i].re = t7_tmp;
            t7 = 0.0 * v_star;
            Tvmax_data[i].im = t7;
            sf_re *= a0;
            t4_tmp_data[i].re = sf_re;
            v_star *= a0;
            t4_tmp_data[i].im = v_star;
            t4_data[i].re = 3.0 * sf_re;
            t4_data[i].im = 3.0 * v_star;
            t5_data[i].re = -t7_tmp;
            t5_data[i].im = -t7;
            tmp_data[i].re = v0 * varargin_1_re;
            tmp_data[i].im = v0 * varargin_1_im;
          }
          Tvmax_size = T1_size[0];
          for (i = 0; i < loop_ub; i++) {
            sf_re = T1_data[i].re;
            v_star = v1_data[i].im;
            t7_tmp = T1_data[i].im;
            t7 = v1_data[i].re;
            b_tmp_data[i].re = sf_re * t7 - t7_tmp * v_star;
            b_tmp_data[i].im = sf_re * v_star + t7_tmp * t7;
          }
          if (T1_size[0] == 1) {
            i = Tvmax_size;
          } else {
            i = T1_size[0];
          }
          if ((T1_size[0] == Tvmax_size) && (T1_size[0] == i) &&
              (T1_size[0] == Tvmax_size) && (T1_size[0] == i) &&
              (T1_size[0] == Tvmax_size) && (i == T1_size[0]) &&
              (T1_size[0] == i)) {
            sf_re = sf * -20.0;
            b_sf_re = sf * -30.0;
            c_sf_re = sf * -12.0;
            dcv[0].re = 0.0;
            dcv[0].im = 0.0;
            dcv[1].re = v0;
            dcv[1].im = 0.0;
            dcv[2].re = a0;
            dcv[2].im = 0.0;
            for (i = 0; i < Tvmax_size; i++) {
              varargin_2 = T1_data[i];
              if ((varargin_2.im == 0.0) && (varargin_2.re >= 0.0)) {
                re = rt_powd_snf(varargin_2.re, 3.0);
                t7 = 0.0;
              } else if (varargin_2.re == 0.0) {
                re = 0.0;
                t7 = -rt_powd_snf(varargin_2.im, 3.0);
              } else {
                b_log(&varargin_2);
                v_star = 3.0 * varargin_2.re;
                t7_tmp = 3.0 * varargin_2.im;
                if (v_star == 0.0) {
                  re = cos(t7_tmp);
                  t7 = sin(t7_tmp);
                } else if (t7_tmp == 0.0) {
                  re = exp(v_star);
                  t7 = 0.0;
                } else if (rtIsInf(t7_tmp) && rtIsInf(v_star) &&
                           (v_star < 0.0)) {
                  re = 0.0;
                  t7 = 0.0;
                } else {
                  v_star = exp(v_star / 2.0);
                  re = v_star * (v_star * cos(t7_tmp));
                  t7 = v_star * (v_star * sin(t7_tmp));
                }
              }
              if (t7 == 0.0) {
                re = 1.0 / re;
                t7 = 0.0;
              } else if (re == 0.0) {
                re = 0.0;
                t7 = -(1.0 / t7);
              } else {
                brm = fabs(re);
                v_star = fabs(t7);
                if (brm > v_star) {
                  s = t7 / re;
                  v_star = re + s * t7;
                  re = (s * 0.0 + 1.0) / v_star;
                  t7 = (0.0 - s) / v_star;
                } else if (v_star == brm) {
                  if (re > 0.0) {
                    s = 0.5;
                  } else {
                    s = -0.5;
                  }
                  if (t7 > 0.0) {
                    v_star = 0.5;
                  } else {
                    v_star = -0.5;
                  }
                  re = (s + 0.0 * v_star) / brm;
                  t7 = (0.0 * s - v_star) / brm;
                } else {
                  s = re / t7;
                  v_star = t7 + s * re;
                  re = s / v_star;
                  t7 = (s * 0.0 - 1.0) / v_star;
                }
              }
              v_star = (((sf_re + t4_data[i].re) + t5_data[i].re) +
                        tmp_data[i].re * 12.0) +
                       b_tmp_data[i].re * 8.0;
              t7_tmp =
                  ((t4_data[i].im + t5_data[i].im) + tmp_data[i].im * 12.0) +
                  b_tmp_data[i].im * 8.0;
              dcv[i + 3].re = -3.0 * (re * v_star - t7 * t7_tmp);
              dcv[i + 3].im = -3.0 * (re * t7_tmp + t7 * v_star);
            }
            for (i = 0; i < loop_ub; i++) {
              varargin_1_re = Tvmin_data[i].re;
              varargin_1_im = Tvmin_data[i].im;
              re =
                  varargin_1_re * varargin_1_re - varargin_1_im * varargin_1_im;
              v_star = varargin_1_re * varargin_1_im;
              t7 = v_star + v_star;
              if (t7 == 0.0) {
                re = 1.0 / re;
                t7 = 0.0;
              } else if (re == 0.0) {
                re = 0.0;
                t7 = -(1.0 / t7);
              } else {
                brm = fabs(re);
                v_star = fabs(t7);
                if (brm > v_star) {
                  s = t7 / re;
                  v_star = re + s * t7;
                  re = (s * 0.0 + 1.0) / v_star;
                  t7 = (0.0 - s) / v_star;
                } else if (v_star == brm) {
                  if (re > 0.0) {
                    s = 0.5;
                  } else {
                    s = -0.5;
                  }
                  if (t7 > 0.0) {
                    v_star = 0.5;
                  } else {
                    v_star = -0.5;
                  }
                  re = (s + 0.0 * v_star) / brm;
                  t7 = (0.0 * s - v_star) / brm;
                } else {
                  s = re / t7;
                  v_star = t7 + s * re;
                  re = s / v_star;
                  t7 = (s * 0.0 - 1.0) / v_star;
                }
              }
              sf_re = (((b_sf_re - Tvmax_data[i].re * 2.0) + t4_data[i].re) +
                       tmp_data[i].re * 16.0) +
                      b_tmp_data[i].re * 14.0;
              t7_tmp = (((0.0 - Tvmax_data[i].im * 2.0) + t4_data[i].im) +
                        tmp_data[i].im * 16.0) +
                       b_tmp_data[i].im * 14.0;
              dcv[i + 4].re = 12.0 * (re * sf_re - t7 * t7_tmp);
              dcv[i + 4].im = 12.0 * (re * t7_tmp + t7 * sf_re);
            }
            for (i = 0; i < Tvmax_size; i++) {
              varargin_2 = T1_data[i];
              if ((varargin_2.im == 0.0) && (varargin_2.re >= 0.0)) {
                re = rt_powd_snf(varargin_2.re, 5.0);
                t7 = 0.0;
              } else if (varargin_2.re == 0.0) {
                re = 0.0;
                t7 = rt_powd_snf(varargin_2.im, 5.0);
              } else {
                b_log(&varargin_2);
                v_star = 5.0 * varargin_2.re;
                t7_tmp = 5.0 * varargin_2.im;
                if (v_star == 0.0) {
                  re = cos(t7_tmp);
                  t7 = sin(t7_tmp);
                } else if (t7_tmp == 0.0) {
                  re = exp(v_star);
                  t7 = 0.0;
                } else if (rtIsInf(t7_tmp) && rtIsInf(v_star) &&
                           (v_star < 0.0)) {
                  re = 0.0;
                  t7 = 0.0;
                } else {
                  v_star = exp(v_star / 2.0);
                  re = v_star * (v_star * cos(t7_tmp));
                  t7 = v_star * (v_star * sin(t7_tmp));
                }
              }
              if (t7 == 0.0) {
                re = 1.0 / re;
                t7 = 0.0;
              } else if (re == 0.0) {
                re = 0.0;
                t7 = -(1.0 / t7);
              } else {
                brm = fabs(re);
                v_star = fabs(t7);
                if (brm > v_star) {
                  s = t7 / re;
                  v_star = re + s * t7;
                  re = (s * 0.0 + 1.0) / v_star;
                  t7 = (0.0 - s) / v_star;
                } else if (v_star == brm) {
                  if (re > 0.0) {
                    s = 0.5;
                  } else {
                    s = -0.5;
                  }
                  if (t7 > 0.0) {
                    v_star = 0.5;
                  } else {
                    v_star = -0.5;
                  }
                  re = (s + 0.0 * v_star) / brm;
                  t7 = (0.0 * s - v_star) / brm;
                } else {
                  s = re / t7;
                  v_star = t7 + s * re;
                  re = s / v_star;
                  t7 = (s * 0.0 - 1.0) / v_star;
                }
              }
              sf_re = (((c_sf_re + t5_data[i].re) + tmp_data[i].re * 6.0) +
                       b_tmp_data[i].re * 6.0) +
                      t4_tmp_data[i].re;
              t7_tmp = ((t5_data[i].im + tmp_data[i].im * 6.0) +
                        b_tmp_data[i].im * 6.0) +
                       t4_tmp_data[i].im;
              dcv[i + 5].re = -60.0 * (re * sf_re - t7 * t7_tmp);
              dcv[i + 5].im = -60.0 * (re * t7_tmp + t7 * sf_re);
            }
            coeffsT1_size[0] = 1;
            coeffsT1_size[1] = 6;
            memcpy(&coeffsT1_data[0], &dcv[0], 6U * sizeof(creal_T));
          } else {
            binary_expand_op(coeffsT1_data, coeffsT1_size, v0, a0, T1_data,
                             &T1_size[0], sf, t4_data, &T1_size[0], t5_data,
                             &T1_size[0], tmp_data, b_tmp_data, &Tvmax_size,
                             Tvmin_data, &T1_size[0], Tvmax_data, &T1_size[0],
                             t4_tmp_data, &T1_size[0]);
          }
          /* 'student_pass_primitive:43' coeffsT2 = coef_list_fun(v0, a0, sf,
           * v2, 0., T2); */
          /* COEF_LIST_FUN */
          /*     COEF_LIST_VAR = COEF_LIST_FUN(V0,A0,SF,VF,AF,T) */
          /*     This function was generated by the Symbolic Math Toolbox
           * version 24.2. */
          /*     24-Oct-2025 11:34:51 */
          /* 'coef_list_fun:8' t2 = T.^2; */
          /* 'coef_list_fun:9' t3 = af.*t2; */
          /* 'coef_list_fun:10' t4 = a0.*t2.*3.0; */
          /* 'coef_list_fun:11' t5 = -t3; */
          /* 'coef_list_fun:12' coef_list_var =
           * [0.0,v0,a0,1.0./T.^3.*(sf.*-2.0e+1+t4+t5+T.*v0.*1.2e+1+T.*vf.*8.0).*-3.0,1.0./t2.^2.*(sf.*-3.0e+1-t3.*2.0+t4+T.*v0.*1.6e+1+T.*vf.*1.4e+1).*1.2e+1,1.0./T.^5.*(sf.*-1.2e+1+t5+T.*v0.*6.0+T.*vf.*6.0+a0.*t2).*-6.0e+1];
           */
          for (i = 0; i < Tvmin_size; i++) {
            varargin_1_re = T2_data[i].re;
            varargin_1_im = T2_data[i].im;
            sf_re =
                varargin_1_re * varargin_1_re - varargin_1_im * varargin_1_im;
            Tvmin_data[i].re = sf_re;
            v_star = varargin_1_re * varargin_1_im;
            v_star += v_star;
            Tvmin_data[i].im = v_star;
            t7_tmp = 0.0 * sf_re;
            Tvmax_data[i].re = t7_tmp;
            t7 = 0.0 * v_star;
            Tvmax_data[i].im = t7;
            sf_re *= a0;
            t4_tmp_data[i].re = sf_re;
            v_star *= a0;
            t4_tmp_data[i].im = v_star;
            t4_data[i].re = 3.0 * sf_re;
            t4_data[i].im = 3.0 * v_star;
            t5_data[i].re = -t7_tmp;
            t5_data[i].im = -t7;
            tmp_data[i].re = v0 * varargin_1_re;
            tmp_data[i].im = v0 * varargin_1_im;
          }
          Tvmax_size = T2_size[0];
          for (i = 0; i < Tvmin_size; i++) {
            sf_re = T2_data[i].re;
            v_star = v2_data[i].im;
            t7_tmp = T2_data[i].im;
            t7 = v2_data[i].re;
            b_tmp_data[i].re = sf_re * t7 - t7_tmp * v_star;
            b_tmp_data[i].im = sf_re * v_star + t7_tmp * t7;
          }
          if (T2_size[0] == 1) {
            i = Tvmax_size;
          } else {
            i = T2_size[0];
          }
          if ((T2_size[0] == Tvmax_size) && (T2_size[0] == i) &&
              (T2_size[0] == Tvmax_size) && (T2_size[0] == i) &&
              (T2_size[0] == Tvmax_size) && (i == T2_size[0]) &&
              (T2_size[0] == i)) {
            sf_re = sf * -20.0;
            b_sf_re = sf * -30.0;
            c_sf_re = sf * -12.0;
            dcv[0].re = 0.0;
            dcv[0].im = 0.0;
            dcv[1].re = v0;
            dcv[1].im = 0.0;
            dcv[2].re = a0;
            dcv[2].im = 0.0;
            for (i = 0; i < Tvmax_size; i++) {
              varargin_2 = T2_data[i];
              if ((varargin_2.im == 0.0) && (varargin_2.re >= 0.0)) {
                re = rt_powd_snf(varargin_2.re, 3.0);
                t7 = 0.0;
              } else if (varargin_2.re == 0.0) {
                re = 0.0;
                t7 = -rt_powd_snf(varargin_2.im, 3.0);
              } else {
                b_log(&varargin_2);
                v_star = 3.0 * varargin_2.re;
                t7_tmp = 3.0 * varargin_2.im;
                if (v_star == 0.0) {
                  re = cos(t7_tmp);
                  t7 = sin(t7_tmp);
                } else if (t7_tmp == 0.0) {
                  re = exp(v_star);
                  t7 = 0.0;
                } else if (rtIsInf(t7_tmp) && rtIsInf(v_star) &&
                           (v_star < 0.0)) {
                  re = 0.0;
                  t7 = 0.0;
                } else {
                  v_star = exp(v_star / 2.0);
                  re = v_star * (v_star * cos(t7_tmp));
                  t7 = v_star * (v_star * sin(t7_tmp));
                }
              }
              if (t7 == 0.0) {
                re = 1.0 / re;
                t7 = 0.0;
              } else if (re == 0.0) {
                re = 0.0;
                t7 = -(1.0 / t7);
              } else {
                brm = fabs(re);
                v_star = fabs(t7);
                if (brm > v_star) {
                  s = t7 / re;
                  v_star = re + s * t7;
                  re = (s * 0.0 + 1.0) / v_star;
                  t7 = (0.0 - s) / v_star;
                } else if (v_star == brm) {
                  if (re > 0.0) {
                    s = 0.5;
                  } else {
                    s = -0.5;
                  }
                  if (t7 > 0.0) {
                    v_star = 0.5;
                  } else {
                    v_star = -0.5;
                  }
                  re = (s + 0.0 * v_star) / brm;
                  t7 = (0.0 * s - v_star) / brm;
                } else {
                  s = re / t7;
                  v_star = t7 + s * re;
                  re = s / v_star;
                  t7 = (s * 0.0 - 1.0) / v_star;
                }
              }
              v_star = (((sf_re + t4_data[i].re) + t5_data[i].re) +
                        tmp_data[i].re * 12.0) +
                       b_tmp_data[i].re * 8.0;
              t7_tmp =
                  ((t4_data[i].im + t5_data[i].im) + tmp_data[i].im * 12.0) +
                  b_tmp_data[i].im * 8.0;
              dcv[i + 3].re = -3.0 * (re * v_star - t7 * t7_tmp);
              dcv[i + 3].im = -3.0 * (re * t7_tmp + t7 * v_star);
            }
            for (i = 0; i < Tvmin_size; i++) {
              varargin_1_re = Tvmin_data[i].re;
              varargin_1_im = Tvmin_data[i].im;
              re =
                  varargin_1_re * varargin_1_re - varargin_1_im * varargin_1_im;
              v_star = varargin_1_re * varargin_1_im;
              t7 = v_star + v_star;
              if (t7 == 0.0) {
                re = 1.0 / re;
                t7 = 0.0;
              } else if (re == 0.0) {
                re = 0.0;
                t7 = -(1.0 / t7);
              } else {
                brm = fabs(re);
                v_star = fabs(t7);
                if (brm > v_star) {
                  s = t7 / re;
                  v_star = re + s * t7;
                  re = (s * 0.0 + 1.0) / v_star;
                  t7 = (0.0 - s) / v_star;
                } else if (v_star == brm) {
                  if (re > 0.0) {
                    s = 0.5;
                  } else {
                    s = -0.5;
                  }
                  if (t7 > 0.0) {
                    v_star = 0.5;
                  } else {
                    v_star = -0.5;
                  }
                  re = (s + 0.0 * v_star) / brm;
                  t7 = (0.0 * s - v_star) / brm;
                } else {
                  s = re / t7;
                  v_star = t7 + s * re;
                  re = s / v_star;
                  t7 = (s * 0.0 - 1.0) / v_star;
                }
              }
              sf_re = (((b_sf_re - Tvmax_data[i].re * 2.0) + t4_data[i].re) +
                       tmp_data[i].re * 16.0) +
                      b_tmp_data[i].re * 14.0;
              t7_tmp = (((0.0 - Tvmax_data[i].im * 2.0) + t4_data[i].im) +
                        tmp_data[i].im * 16.0) +
                       b_tmp_data[i].im * 14.0;
              dcv[i + 4].re = 12.0 * (re * sf_re - t7 * t7_tmp);
              dcv[i + 4].im = 12.0 * (re * t7_tmp + t7 * sf_re);
            }
            for (i = 0; i < Tvmax_size; i++) {
              varargin_2 = T2_data[i];
              if ((varargin_2.im == 0.0) && (varargin_2.re >= 0.0)) {
                re = rt_powd_snf(varargin_2.re, 5.0);
                t7 = 0.0;
              } else if (varargin_2.re == 0.0) {
                re = 0.0;
                t7 = rt_powd_snf(varargin_2.im, 5.0);
              } else {
                b_log(&varargin_2);
                v_star = 5.0 * varargin_2.re;
                t7_tmp = 5.0 * varargin_2.im;
                if (v_star == 0.0) {
                  re = cos(t7_tmp);
                  t7 = sin(t7_tmp);
                } else if (t7_tmp == 0.0) {
                  re = exp(v_star);
                  t7 = 0.0;
                } else if (rtIsInf(t7_tmp) && rtIsInf(v_star) &&
                           (v_star < 0.0)) {
                  re = 0.0;
                  t7 = 0.0;
                } else {
                  v_star = exp(v_star / 2.0);
                  re = v_star * (v_star * cos(t7_tmp));
                  t7 = v_star * (v_star * sin(t7_tmp));
                }
              }
              if (t7 == 0.0) {
                re = 1.0 / re;
                t7 = 0.0;
              } else if (re == 0.0) {
                re = 0.0;
                t7 = -(1.0 / t7);
              } else {
                brm = fabs(re);
                v_star = fabs(t7);
                if (brm > v_star) {
                  s = t7 / re;
                  v_star = re + s * t7;
                  re = (s * 0.0 + 1.0) / v_star;
                  t7 = (0.0 - s) / v_star;
                } else if (v_star == brm) {
                  if (re > 0.0) {
                    s = 0.5;
                  } else {
                    s = -0.5;
                  }
                  if (t7 > 0.0) {
                    v_star = 0.5;
                  } else {
                    v_star = -0.5;
                  }
                  re = (s + 0.0 * v_star) / brm;
                  t7 = (0.0 * s - v_star) / brm;
                } else {
                  s = re / t7;
                  v_star = t7 + s * re;
                  re = s / v_star;
                  t7 = (s * 0.0 - 1.0) / v_star;
                }
              }
              sf_re = (((c_sf_re + t5_data[i].re) + tmp_data[i].re * 6.0) +
                       b_tmp_data[i].re * 6.0) +
                      t4_tmp_data[i].re;
              t7_tmp = ((t5_data[i].im + tmp_data[i].im * 6.0) +
                        b_tmp_data[i].im * 6.0) +
                       t4_tmp_data[i].im;
              dcv[i + 5].re = -60.0 * (re * sf_re - t7 * t7_tmp);
              dcv[i + 5].im = -60.0 * (re * t7_tmp + t7 * sf_re);
            }
            coeffsT2_size[0] = 1;
            coeffsT2_size[1] = 6;
            memcpy(&coeffsT2_data[0], &dcv[0], 6U * sizeof(creal_T));
          } else {
            binary_expand_op(coeffsT2_data, coeffsT2_size, v0, a0, T2_data,
                             &T2_size[0], sf, t4_data, &T2_size[0], t5_data,
                             &T2_size[0], tmp_data, b_tmp_data, &Tvmax_size,
                             Tvmin_data, &T2_size[0], Tvmax_data, &T2_size[0],
                             t4_tmp_data, &T2_size[0]);
          }
        } else {
          guard1 = true;
        }
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }
  } else {
    guard1 = true;
  }
  if (guard1) {
    /* 'student_pass_primitive:44' else */
    /* 'student_pass_primitive:45' coeffsT1 = [0.;0.;0.;0.;0.;0.]; */
    coeffsT1_size[0] = 6;
    coeffsT1_size[1] = 1;
    /* 'student_pass_primitive:46' coeffsT2 = [0.;0.;0.;0.;0.;0.]; */
    coeffsT2_size[0] = 6;
    coeffsT2_size[1] = 1;
    memset(&coeffsT1_data[0], 0, 6U * sizeof(creal_T));
    memset(&coeffsT2_data[0], 0, 6U * sizeof(creal_T));
    /* 'student_pass_primitive:47' Tvmin = 0.; */
    /* 'student_pass_primitive:48' Tvmax = 0.; */
    /* 'student_pass_primitive:49' T_star = 0.; */
    /* 'student_pass_primitive:50' v_star = 0.; */
    /* 'student_pass_primitive:51' v1 = 0.; */
    v1_size[0] = 1;
    v1_data[0].re = 0.0;
    v1_data[0].im = 0.0;
    /* 'student_pass_primitive:52' v2 = 0.; */
    v2_size[0] = 1;
    v2_data[0].re = 0.0;
    v2_data[0].im = 0.0;
  }
}

/*
 * function [coefs,maxsf,tf] = student_stop_primitive(v0,a0,sf)
 */
void student_stop_primitive(double v0, double a0, double sf, double coefs[6],
                            double *maxsf, double *tf)
{
  int i;
  /*  */
  /*           Agent Logic */
  /*          Stop Primitive */
  /*               2025 */
  /*  */
  /*  Impossible cases: negative velocity and null final position --> */
  /*  return 0 */
  /* 'student_stop_primitive:11' if (v0 <= 0) || (sf == 0) */
  if ((v0 <= 0.0) || (sf == 0.0)) {
    /* 'student_stop_primitive:12' coefs = zeros(1,6); */
    for (i = 0; i < 6; i++) {
      coefs[i] = 0.0;
    }
    /*  Returns a vector of 6 zeros for each coefficient */
    /* 'student_stop_primitive:13' tf = 0.; */
    *tf = 0.0;
    /* 'student_stop_primitive:14' maxsf = 0.; */
    *maxsf = 0.0;
    /*  Also set other returns to 0 */
  } else {
    double coefs_tmp;
    double t2;
    double t3;
    double t4;
    double t4_tmp;
    /* 'student_stop_primitive:15' else */
    /*  Can't reach a specified sf using given inputs (v0, a0) */
    /* 'student_stop_primitive:17' if 4*v0^2 + 5*a0*sf < 0 */
    t2 = 4.0 * (v0 * v0);
    if (t2 + 5.0 * a0 * sf < 0.0) {
      /* 'student_stop_primitive:18' maxsf = - (4*v0^2)/(5*a0); */
      *maxsf = -t2 / (5.0 * a0);
      /* 'student_stop_primitive:19' tf = (10*maxsf)/(2*v0); */
      *tf = 10.0 * *maxsf / (2.0 * v0);
    } else {
      /* 'student_stop_primitive:20' else */
      /*  Valid case: compute final time with defined function */
      /* 'student_stop_primitive:22' maxsf = sf; */
      *maxsf = sf;
      /* 'student_stop_primitive:23' tf = final_opt_time_stop(v0, a0, maxsf); */
      /* FINAL_OPT_TIME_STOP */
      /*     OUT1 = FINAL_OPT_TIME_STOP(V0,A0,SF) */
      /*     This function was generated by the Symbolic Math Toolbox
       * version 24.2. */
      /*     24-Oct-2025 11:34:57 */
      /* 'final_opt_time_stop:8' out1 =
       * (sf.*1.0e+1)./(v0.*2.0+sqrt(a0.*sf.*5.0+v0.^2.*4.0)); */
      *tf = sf * 10.0 / (v0 * 2.0 + sqrt(a0 * sf * 5.0 + t2));
    }
    /*  Evaluate coefficients (consider vf = af = 0) */
    /* 'student_stop_primitive:27' coefs = coef_list_fun(v0, a0, 0., 0., 0.,
     * tf); */
    /* COEF_LIST_FUN */
    /*     COEF_LIST_VAR = COEF_LIST_FUN(V0,A0,SF,VF,AF,T) */
    /*     This function was generated by the Symbolic Math Toolbox
     * version 24.2. */
    /*     24-Oct-2025 11:34:51 */
    /* 'coef_list_fun:8' t2 = T.^2; */
    t2 = *tf * *tf;
    /* 'coef_list_fun:9' t3 = af.*t2; */
    t3 = 0.0 * t2;
    /* 'coef_list_fun:10' t4 = a0.*t2.*3.0; */
    t4_tmp = a0 * t2;
    t4 = t4_tmp * 3.0;
    /* 'coef_list_fun:11' t5 = -t3; */
    /* 'coef_list_fun:12' coef_list_var =
     * [0.0,v0,a0,1.0./T.^3.*(sf.*-2.0e+1+t4+t5+T.*v0.*1.2e+1+T.*vf.*8.0).*-3.0,1.0./t2.^2.*(sf.*-3.0e+1-t3.*2.0+t4+T.*v0.*1.6e+1+T.*vf.*1.4e+1).*1.2e+1,1.0./T.^5.*(sf.*-1.2e+1+t5+T.*v0.*6.0+T.*vf.*6.0+a0.*t2).*-6.0e+1];
     */
    coefs[0] = 0.0;
    coefs[1] = v0;
    coefs[2] = a0;
    coefs_tmp = *tf * v0;
    coefs[3] = 1.0 / rt_powd_snf(*tf, 3.0) *
               (((t4 - t3) + coefs_tmp * 12.0) + *tf * 0.0 * 8.0) * -3.0;
    coefs[4] =
        1.0 / (t2 * t2) *
        ((((-0.0 - t3 * 2.0) + t4) + coefs_tmp * 16.0) + *tf * 0.0 * 14.0) *
        12.0;
    coefs[5] = 1.0 / rt_powd_snf(*tf, 5.0) *
               (((-t3 + coefs_tmp * 6.0) + *tf * 0.0 * 6.0) + t4_tmp) * -60.0;
  }
}

/* End of code generation (primitives.c) */
