/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_student_stop_primitive_api.h
 *
 * Code generation for function 'student_stop_primitive'
 *
 */

#ifndef _CODER_STUDENT_STOP_PRIMITIVE_API_H
#define _CODER_STUDENT_STOP_PRIMITIVE_API_H

/* Include files */
#include "emlrt.h"
#include "mex.h"
#include "tmwtypes.h"
#include <string.h>

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void student_stop_primitive(real_T v0, real_T a0, real_T sf, real_T coefs[6],
                            real_T *maxsf, real_T *tf);

void student_stop_primitive_api(const mxArray *const prhs[3], int32_T nlhs,
                                const mxArray *plhs[3]);

void student_stop_primitive_atexit(void);

void student_stop_primitive_initialize(void);

void student_stop_primitive_terminate(void);

void student_stop_primitive_xil_shutdown(void);

void student_stop_primitive_xil_terminate(void);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (_coder_student_stop_primitive_api.h) */
