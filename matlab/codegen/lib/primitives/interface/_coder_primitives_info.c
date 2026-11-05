/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_primitives_info.c
 *
 * Code generation for function 'coef_list_fun'
 *
 */

/* Include files */
#include "_coder_primitives_info.h"
#include "emlrt.h"
#include "tmwtypes.h"

/* Function Declarations */
static const mxArray *c_emlrtMexFcnResolvedFunctionsI(void);

/* Function Definitions */
static const mxArray *c_emlrtMexFcnResolvedFunctionsI(void)
{
  const mxArray *nameCaptureInfo;
  const char_T *data[10] = {
      "789ced9add72da4614c7d76dfa71d10fa699767a9167c8246e3038770a888f0884413802"
      "a28e2260c1027d21090cb49d66dadef7017ad717685e218fd0bbde74"
      "c6579d642657ed0b54489641aa77c0111696ba7bb35efdb57bce9e957e5ece0aec15cb7b"
      "00804f8053be2c3bf5c7e7edc479fd0ef016bfbee7bbbee7bd1dbc07",
      "6e5daaff7c5e7754c58453d36928820c2f7a7655595404c5accf34087468a8d204766da5"
      "274ab02eca90596dd08b969c5b912e1a0b69f177e6047686cc5806fa"
      "89b1f4505a6d5cc4e30c31df5b1bc6e339221e099ffe94fc3af3903b36a06e7084040d43"
      "50babaca65a13134558d3ba68b759a2b13f92253af11259263c84c85",
      "ce567882a62bdc912eca2acf40992b3e21b25c0d6aaac1b50543ecf0421f2a26cfd43959"
      "3025a1cd19e6b8bbb86258a3f29ad54f34c509bc2bbbf39d069cefed"
      "35f375f59eb5a612af6a266f5a2b63bb635f77fd7816d08ff7917e384a571db725b8bd75"
      "fe0e69cfab87b5ce97c4f7ae0cd6c7f7d30de7ebaf97f77fe8d4f29f",
      "b61496bdbff60bbf8769cf2dbbb217f43dfd02612fe1d31f173bc7a74c523f22d583d683"
      "5ebb3e48a5d8ecd28fa33576d6f90110edb0c63f43f48f2bdfb580f3"
      "45fd7f4ff8f48e0a7bbc241a26df1b2b60775cff23a0bd53a43daf1ed6fa7ae26a113dae"
      "3c7ffd4f12f31c6c9fe7a559919a1798a9208a05e9b095cbd70ab354",
      "3e3e3c7f85e8bf691c7f448c9ff0e9bbe5393fb867efe6c27f7e2ed9575ace4477dffe1c"
      "69cfabef70df6e2f765c39ff779fc59c07db7f4fe9f4fd99381ba606"
      "add2e323cd54eff7dbba51c09c8f3ae783eedf51714ef8ea25872c9f2d37a2caf739d29e"
      "570f9fef4e5ccf9737b6fbf85fdf5430dfc1f6f9aeb22a4b120d4627",
      "5bad119b6db2b983513983f9fe7fe73bcecf78efc3f999ebb187f3334ec1f999ab8d7f86"
      "e81fd5dfeb2ed735c1306ed079eac21dfb7a54f7ed37fc3c7511df58"
      "9fa796eff431dfc1f6f9fe95d09006eda644d42ae94c4e6852ad475506f33db27c0fba5f"
      "ff68cd7c5ddda68e2c2afc044a8b7654f7eb13a43daf1ed6faaec675",
      "01f4b8f2fc87df24cc73b07d9e8fcb658d18d3d464aa2655b694aacd93b54731cac39c21"
      "fac795e7cf02cef703e47c1dc5457854f93d42daf3ea61ade70abaed"
      "12577ebf204dcc6f700ddf33e44b5435d3d0ea6cb5caa429ea98d89f1ee630bfa3caefa0"
      "cfcd676be6ebeacb7c80459f65ba25b2f9966f91f6bc7af8f91637be",
      "0ee2e3caf7a7bff098ef60f3e7f57384bd844fcf768bb374aa911fe9aabccfce0f92fd71"
      "960698ef51e57bd07c0b3e1ff5de87cf47afc71e3e1f750a3e1fbdda"
      "f871fdeec5cbf31bf3fdbaed16fe7efd5acf4963fdfd7aeedd9f30e7c1f6dfd339dd1bc9"
      "b34a63383e11e7994a931d28d5748cf2ea98f36f37ef4d7ff7fd377f",
      "60633eb29cff1e69cfabef2e3fe3ac755c39fff2ce3798f360fb9c3f9d3f38ac0c4f86a3"
      "5993ea50cd7b64b7d8a348ccf9a8731ee769ae660fe769709e260c7b"
      "384f136cfc7f017dc67535",
      ""};
  nameCaptureInfo = NULL;
  emlrtNameCaptureMxArrayR2016a(&data[0], 19760U, &nameCaptureInfo);
  return nameCaptureInfo;
}

mxArray *emlrtMexFcnProperties(void)
{
  mxArray *xEntryPoints;
  mxArray *xInputs;
  mxArray *xResult;
  const char_T *propFieldName[9] = {"Version",
                                    "ResolvedFunctions",
                                    "Checksum",
                                    "EntryPoints",
                                    "CoverageInfo",
                                    "IsPolymorphic",
                                    "PropertyList",
                                    "UUID",
                                    "ClassEntryPointIsHandle"};
  const char_T *epFieldName[8] = {
      "QualifiedName",    "NumberOfInputs", "NumberOfOutputs", "ConstantInputs",
      "ResolvedFilePath", "TimeStamp",      "Constructor",     "Visible"};
  xEntryPoints =
      emlrtCreateStructMatrix(1, 9, 8, (const char_T **)&epFieldName[0]);
  xInputs = emlrtCreateLogicalMatrix(1, 6);
  emlrtSetField(xEntryPoints, 0, "QualifiedName",
                emlrtMxCreateString("coef_list_fun"));
  emlrtSetField(xEntryPoints, 0, "NumberOfInputs",
                emlrtMxCreateDoubleScalar(6.0));
  emlrtSetField(xEntryPoints, 0, "NumberOfOutputs",
                emlrtMxCreateDoubleScalar(1.0));
  emlrtSetField(xEntryPoints, 0, "ConstantInputs", xInputs);
  emlrtSetField(
      xEntryPoints, 0, "ResolvedFilePath",
      emlrtMxCreateString(
          "C:\\Users\\Alessandro\\Desktop\\UNITN\\MAGISTRALE\\SECONDO_"
          "ANNO\\Primo_Sem\\IVAD\\Repos\\basic_agent_ST\\matlab\\coef_list_fun."
          "m"));
  emlrtSetField(xEntryPoints, 0, "TimeStamp",
                emlrtMxCreateDoubleScalar(740291.46527777775));
  emlrtSetField(xEntryPoints, 0, "Constructor",
                emlrtMxCreateLogicalScalar(false));
  emlrtSetField(xEntryPoints, 0, "Visible", emlrtMxCreateLogicalScalar(true));
  xInputs = emlrtCreateLogicalMatrix(1, 2);
  emlrtSetField(xEntryPoints, 1, "QualifiedName",
                emlrtMxCreateString("coeffs_s_opt"));
  emlrtSetField(xEntryPoints, 1, "NumberOfInputs",
                emlrtMxCreateDoubleScalar(2.0));
  emlrtSetField(xEntryPoints, 1, "NumberOfOutputs",
                emlrtMxCreateDoubleScalar(1.0));
  emlrtSetField(xEntryPoints, 1, "ConstantInputs", xInputs);
  emlrtSetField(
      xEntryPoints, 1, "ResolvedFilePath",
      emlrtMxCreateString(
          "C:\\Users\\Alessandro\\Desktop\\UNITN\\MAGISTRALE\\SECONDO_"
          "ANNO\\Primo_Sem\\IVAD\\Repos\\basic_agent_ST\\matlab\\coeffs_s_opt."
          "m"));
  emlrtSetField(xEntryPoints, 1, "TimeStamp",
                emlrtMxCreateDoubleScalar(740291.46527777775));
  emlrtSetField(xEntryPoints, 1, "Constructor",
                emlrtMxCreateLogicalScalar(false));
  emlrtSetField(xEntryPoints, 1, "Visible", emlrtMxCreateLogicalScalar(true));
  xInputs = emlrtCreateLogicalMatrix(1, 2);
  emlrtSetField(xEntryPoints, 2, "QualifiedName",
                emlrtMxCreateString("coeffs_v_opt"));
  emlrtSetField(xEntryPoints, 2, "NumberOfInputs",
                emlrtMxCreateDoubleScalar(2.0));
  emlrtSetField(xEntryPoints, 2, "NumberOfOutputs",
                emlrtMxCreateDoubleScalar(1.0));
  emlrtSetField(xEntryPoints, 2, "ConstantInputs", xInputs);
  emlrtSetField(
      xEntryPoints, 2, "ResolvedFilePath",
      emlrtMxCreateString(
          "C:\\Users\\Alessandro\\Desktop\\UNITN\\MAGISTRALE\\SECONDO_"
          "ANNO\\Primo_Sem\\IVAD\\Repos\\basic_agent_ST\\matlab\\coeffs_v_opt."
          "m"));
  emlrtSetField(xEntryPoints, 2, "TimeStamp",
                emlrtMxCreateDoubleScalar(740291.46527777775));
  emlrtSetField(xEntryPoints, 2, "Constructor",
                emlrtMxCreateLogicalScalar(false));
  emlrtSetField(xEntryPoints, 2, "Visible", emlrtMxCreateLogicalScalar(true));
  xInputs = emlrtCreateLogicalMatrix(1, 2);
  emlrtSetField(xEntryPoints, 3, "QualifiedName",
                emlrtMxCreateString("coeffs_a_opt"));
  emlrtSetField(xEntryPoints, 3, "NumberOfInputs",
                emlrtMxCreateDoubleScalar(2.0));
  emlrtSetField(xEntryPoints, 3, "NumberOfOutputs",
                emlrtMxCreateDoubleScalar(1.0));
  emlrtSetField(xEntryPoints, 3, "ConstantInputs", xInputs);
  emlrtSetField(
      xEntryPoints, 3, "ResolvedFilePath",
      emlrtMxCreateString(
          "C:\\Users\\Alessandro\\Desktop\\UNITN\\MAGISTRALE\\SECONDO_"
          "ANNO\\Primo_Sem\\IVAD\\Repos\\basic_agent_ST\\matlab\\coeffs_a_opt."
          "m"));
  emlrtSetField(xEntryPoints, 3, "TimeStamp",
                emlrtMxCreateDoubleScalar(740291.46527777775));
  emlrtSetField(xEntryPoints, 3, "Constructor",
                emlrtMxCreateLogicalScalar(false));
  emlrtSetField(xEntryPoints, 3, "Visible", emlrtMxCreateLogicalScalar(true));
  xInputs = emlrtCreateLogicalMatrix(1, 2);
  emlrtSetField(xEntryPoints, 4, "QualifiedName",
                emlrtMxCreateString("coeffs_j_opt"));
  emlrtSetField(xEntryPoints, 4, "NumberOfInputs",
                emlrtMxCreateDoubleScalar(2.0));
  emlrtSetField(xEntryPoints, 4, "NumberOfOutputs",
                emlrtMxCreateDoubleScalar(1.0));
  emlrtSetField(xEntryPoints, 4, "ConstantInputs", xInputs);
  emlrtSetField(
      xEntryPoints, 4, "ResolvedFilePath",
      emlrtMxCreateString(
          "C:\\Users\\Alessandro\\Desktop\\UNITN\\MAGISTRALE\\SECONDO_"
          "ANNO\\Primo_Sem\\IVAD\\Repos\\basic_agent_ST\\matlab\\coeffs_j_opt."
          "m"));
  emlrtSetField(xEntryPoints, 4, "TimeStamp",
                emlrtMxCreateDoubleScalar(740291.46527777775));
  emlrtSetField(xEntryPoints, 4, "Constructor",
                emlrtMxCreateLogicalScalar(false));
  emlrtSetField(xEntryPoints, 4, "Visible", emlrtMxCreateLogicalScalar(true));
  xInputs = emlrtCreateLogicalMatrix(1, 3);
  emlrtSetField(xEntryPoints, 5, "QualifiedName",
                emlrtMxCreateString("student_stop_primitive"));
  emlrtSetField(xEntryPoints, 5, "NumberOfInputs",
                emlrtMxCreateDoubleScalar(3.0));
  emlrtSetField(xEntryPoints, 5, "NumberOfOutputs",
                emlrtMxCreateDoubleScalar(3.0));
  emlrtSetField(xEntryPoints, 5, "ConstantInputs", xInputs);
  emlrtSetField(
      xEntryPoints, 5, "ResolvedFilePath",
      emlrtMxCreateString(
          "C:\\Users\\Alessandro\\Desktop\\UNITN\\MAGISTRALE\\SECONDO_"
          "ANNO\\Primo_Sem\\IVAD\\Repos\\basic_agent_ST\\matlab\\student_stop_"
          "pr"
          "imitive.m"));
  emlrtSetField(xEntryPoints, 5, "TimeStamp",
                emlrtMxCreateDoubleScalar(740290.78478009254));
  emlrtSetField(xEntryPoints, 5, "Constructor",
                emlrtMxCreateLogicalScalar(false));
  emlrtSetField(xEntryPoints, 5, "Visible", emlrtMxCreateLogicalScalar(true));
  xInputs = emlrtCreateLogicalMatrix(1, 2);
  emlrtSetField(xEntryPoints, 6, "QualifiedName",
                emlrtMxCreateString("student_stop_primitive_j0"));
  emlrtSetField(xEntryPoints, 6, "NumberOfInputs",
                emlrtMxCreateDoubleScalar(2.0));
  emlrtSetField(xEntryPoints, 6, "NumberOfOutputs",
                emlrtMxCreateDoubleScalar(3.0));
  emlrtSetField(xEntryPoints, 6, "ConstantInputs", xInputs);
  emlrtSetField(
      xEntryPoints, 6, "ResolvedFilePath",
      emlrtMxCreateString(
          "C:\\Users\\Alessandro\\Desktop\\UNITN\\MAGISTRALE\\SECONDO_"
          "ANNO\\Primo_Sem\\IVAD\\Repos\\basic_agent_ST\\matlab\\student_stop_"
          "pr"
          "imitive_j0.m"));
  emlrtSetField(xEntryPoints, 6, "TimeStamp",
                emlrtMxCreateDoubleScalar(739966.76707175921));
  emlrtSetField(xEntryPoints, 6, "Constructor",
                emlrtMxCreateLogicalScalar(false));
  emlrtSetField(xEntryPoints, 6, "Visible", emlrtMxCreateLogicalScalar(true));
  xInputs = emlrtCreateLogicalMatrix(1, 7);
  emlrtSetField(xEntryPoints, 7, "QualifiedName",
                emlrtMxCreateString("student_pass_primitive"));
  emlrtSetField(xEntryPoints, 7, "NumberOfInputs",
                emlrtMxCreateDoubleScalar(7.0));
  emlrtSetField(xEntryPoints, 7, "NumberOfOutputs",
                emlrtMxCreateDoubleScalar(6.0));
  emlrtSetField(xEntryPoints, 7, "ConstantInputs", xInputs);
  emlrtSetField(
      xEntryPoints, 7, "ResolvedFilePath",
      emlrtMxCreateString(
          "C:\\Users\\Alessandro\\Desktop\\UNITN\\MAGISTRALE\\SECONDO_"
          "ANNO\\Primo_Sem\\IVAD\\Repos\\basic_agent_ST\\matlab\\student_pass_"
          "pr"
          "imitive.m"));
  emlrtSetField(xEntryPoints, 7, "TimeStamp",
                emlrtMxCreateDoubleScalar(740290.7849884259));
  emlrtSetField(xEntryPoints, 7, "Constructor",
                emlrtMxCreateLogicalScalar(false));
  emlrtSetField(xEntryPoints, 7, "Visible", emlrtMxCreateLogicalScalar(true));
  xInputs = emlrtCreateLogicalMatrix(1, 5);
  emlrtSetField(xEntryPoints, 8, "QualifiedName",
                emlrtMxCreateString("student_pass_primitive_j0"));
  emlrtSetField(xEntryPoints, 8, "NumberOfInputs",
                emlrtMxCreateDoubleScalar(5.0));
  emlrtSetField(xEntryPoints, 8, "NumberOfOutputs",
                emlrtMxCreateDoubleScalar(3.0));
  emlrtSetField(xEntryPoints, 8, "ConstantInputs", xInputs);
  emlrtSetField(
      xEntryPoints, 8, "ResolvedFilePath",
      emlrtMxCreateString(
          "C:\\Users\\Alessandro\\Desktop\\UNITN\\MAGISTRALE\\SECONDO_"
          "ANNO\\Primo_Sem\\IVAD\\Repos\\basic_agent_ST\\matlab\\student_pass_"
          "pr"
          "imitive_j0.m"));
  emlrtSetField(xEntryPoints, 8, "TimeStamp",
                emlrtMxCreateDoubleScalar(740290.7849884259));
  emlrtSetField(xEntryPoints, 8, "Constructor",
                emlrtMxCreateLogicalScalar(false));
  emlrtSetField(xEntryPoints, 8, "Visible", emlrtMxCreateLogicalScalar(true));
  xResult =
      emlrtCreateStructMatrix(1, 1, 9, (const char_T **)&propFieldName[0]);
  emlrtSetField(xResult, 0, "Version",
                emlrtMxCreateString("24.2.0.2773142 (R2024b) Update 2"));
  emlrtSetField(xResult, 0, "ResolvedFunctions",
                (mxArray *)c_emlrtMexFcnResolvedFunctionsI());
  emlrtSetField(xResult, 0, "Checksum",
                emlrtMxCreateString("nld6RTtCjr480tBsr8L64B"));
  emlrtSetField(xResult, 0, "EntryPoints", xEntryPoints);
  return xResult;
}

/* End of code generation (_coder_primitives_info.c) */
