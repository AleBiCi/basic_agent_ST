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
extern void primitives_initialize(void);

extern void primitives_terminate(void);

extern void student_pass_primitive(
    double v0, double a0, double sf, double vfmin, double vfmax, double Tmin,
    double Tmax, creal_T coeffsT2_data[], int coeffsT2_size[2],
    creal_T v2_data[], int v2_size[1], creal_T T2_data[], int T2_size[1],
    creal_T coeffsT1_data[], int coeffsT1_size[2], creal_T v1_data[],
    int v1_size[1], creal_T T1_data[], int T1_size[1]);

extern void student_stop_primitive(double v0, double a0, double sf,
                                   double coefs[6], double *maxsf, double *tf);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (primitives.h) */
