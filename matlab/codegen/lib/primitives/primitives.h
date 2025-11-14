/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * primitives.h
 *
 * Code generation for function 'primitives'
 *
 */

#ifndef PRIMITIVES_H
#define PRIMITIVES_H

/* Include files */
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
extern double a_opt(double t, double v0, double a0, double sf, double vf,
                    double af, double T);

extern void coef_list_fun(double v0, double a0, double sf, double vf, double af,
                          double T, double coef_list_var[6]);

extern void primitives_initialize(void);

extern void primitives_terminate(void);

extern void student_pass_primitive(double v0, double a0, double sf,
                                   double vfmin, double vfmax, double Tmin,
                                   double Tmax, double coeffsT2_data[],
                                   int coeffsT2_size[2], double *v2, double *T2,
                                   double coeffsT1_data[], int coeffsT1_size[2],
                                   double *v1, double *T1);

extern void student_stop_primitive(double v0, double a0, double sf,
                                   double coefs[6], double *maxsf, double *tf);

extern double v_opt(double t, double v0, double a0, double sf, double vf,
                    double af, double T);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (primitives.h) */
