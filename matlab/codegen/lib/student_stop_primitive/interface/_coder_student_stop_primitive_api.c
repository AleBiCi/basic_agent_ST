/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_student_stop_primitive_api.c
 *
 * Code generation for function 'student_stop_primitive'
 *
 */

/* Include files */
#include "_coder_student_stop_primitive_api.h"
#include "_coder_student_stop_primitive_mex.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;

emlrtContext emlrtContextGlobal = {
    true,                                                 /* bFirstTime */
    false,                                                /* bInitialized */
    131659U,                                              /* fVersionInfo */
    NULL,                                                 /* fErrorFunction */
    "student_stop_primitive",                             /* fFunctionName */
    NULL,                                                 /* fRTCallStack */
    false,                                                /* bDebugMode */
    {2045744189U, 2170104910U, 2743257031U, 4284093946U}, /* fSigWrd */
    NULL                                                  /* fSigMem */
};

/* Function Declarations */
static real_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                 const emlrtMsgIdentifier *parentId);

static const mxArray *b_emlrt_marshallOut(const real_T u);

static real_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                 const emlrtMsgIdentifier *msgId);

static void emlrtExitTimeCleanupDtorFcn(const void *r);

static real_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                               const char_T *identifier);

static const mxArray *emlrt_marshallOut(real_T u[6]);

/* Function Definitions */
static real_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                 const emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = c_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static const mxArray *b_emlrt_marshallOut(const real_T u)
{
  const mxArray *m;
  const mxArray *y;
  y = NULL;
  m = emlrtCreateDoubleScalar(u);
  emlrtAssign(&y, m);
  return y;
}

static real_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                 const emlrtMsgIdentifier *msgId)
{
  static const int32_T dims = 0;
  real_T ret;
  emlrtCheckBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 0U,
                          (const void *)&dims);
  ret = *(real_T *)emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static void emlrtExitTimeCleanupDtorFcn(const void *r)
{
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

static real_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                               const char_T *identifier)
{
  emlrtMsgIdentifier thisId;
  real_T y;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

static const mxArray *emlrt_marshallOut(real_T u[6])
{
  static const int32_T iv[2] = {0, 0};
  static const int32_T iv1[2] = {1, 6};
  const mxArray *m;
  const mxArray *y;
  y = NULL;
  m = emlrtCreateNumericArray(2, (const void *)&iv[0], mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m, &u[0]);
  emlrtSetDimensions((mxArray *)m, &iv1[0], 2);
  emlrtAssign(&y, m);
  return y;
}

void student_stop_primitive_api(const mxArray *const prhs[3], int32_T nlhs,
                                const mxArray *plhs[3])
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  real_T(*coefs)[6];
  real_T a0;
  real_T maxsf;
  real_T sf;
  real_T tf;
  real_T v0;
  st.tls = emlrtRootTLSGlobal;
  coefs = (real_T(*)[6])mxMalloc(sizeof(real_T[6]));
  /* Marshall function inputs */
  v0 = emlrt_marshallIn(&st, emlrtAliasP(prhs[0]), "v0");
  a0 = emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "a0");
  sf = emlrt_marshallIn(&st, emlrtAliasP(prhs[2]), "sf");
  /* Invoke the target function */
  student_stop_primitive(v0, a0, sf, *coefs, &maxsf, &tf);
  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(*coefs);
  if (nlhs > 1) {
    plhs[1] = b_emlrt_marshallOut(maxsf);
  }
  if (nlhs > 2) {
    plhs[2] = b_emlrt_marshallOut(tf);
  }
}

void student_stop_primitive_atexit(void)
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtPushHeapReferenceStackR2021a(
      &st, false, NULL, (void *)&emlrtExitTimeCleanupDtorFcn, NULL, NULL, NULL);
  emlrtEnterRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  student_stop_primitive_xil_terminate();
  student_stop_primitive_xil_shutdown();
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

void student_stop_primitive_initialize(void)
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, NULL);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

void student_stop_primitive_terminate(void)
{
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (_coder_student_stop_primitive_api.c) */
