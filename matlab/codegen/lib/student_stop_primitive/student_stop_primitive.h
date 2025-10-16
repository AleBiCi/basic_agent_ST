/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * student_stop_primitive.h
 *
 * Code generation for function 'student_stop_primitive'
 *
 */

#ifndef STUDENT_STOP_PRIMITIVE_H
#define STUDENT_STOP_PRIMITIVE_H

/* Include files */
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
extern void student_stop_primitive(double v0, double a0, double sf,
                                   double coefs[6], double *maxsf, double *tf);

extern void student_stop_primitive_initialize(void);

extern void student_stop_primitive_terminate(void);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (student_stop_primitive.h) */
