/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_primitives_api.h
 *
 * Code generation for function 'student_pass_primitive'
 *
 */

#ifndef _CODER_PRIMITIVES_API_H
#define _CODER_PRIMITIVES_API_H

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
void primitives_atexit(void);

void primitives_initialize(void);

void primitives_terminate(void);

void primitives_xil_shutdown(void);

void primitives_xil_terminate(void);

void student_pass_primitive(real_T v0, real_T a0, real_T sf, real_T vfmin,
                            real_T vfmax, real_T Tmin, real_T Tmax,
                            creal_T coeffsT2_data[], int32_T coeffsT2_size[2],
                            creal_T v2_data[], int32_T v2_size[1],
                            creal_T T2_data[], int32_T T2_size[1],
                            creal_T coeffsT1_data[], int32_T coeffsT1_size[2],
                            creal_T v1_data[], int32_T v1_size[1],
                            creal_T T1_data[], int32_T T1_size[1]);

void student_pass_primitive_api(const mxArray *const prhs[7], int32_T nlhs,
                                const mxArray *plhs[6]);

void student_stop_primitive(real_T v0, real_T a0, real_T sf, real_T coefs[6],
                            real_T *maxsf, real_T *tf);

void student_stop_primitive_api(const mxArray *const prhs[3], int32_T nlhs,
                                const mxArray *plhs[3]);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (_coder_primitives_api.h) */
