/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_student_stop_primitive_info.c
 *
 * Code generation for function 'student_stop_primitive'
 *
 */

/* Include files */
#include "_coder_student_stop_primitive_info.h"
#include "emlrt.h"
#include "tmwtypes.h"

/* Function Declarations */
static const mxArray *c_emlrtMexFcnResolvedFunctionsI(void);

/* Function Definitions */
static const mxArray *c_emlrtMexFcnResolvedFunctionsI(void)
{
  const mxArray *nameCaptureInfo;
  const char_T *data[6] = {
      "789ced56416fd3301476d040bb009590b8f003b8cdb7dd431b68046d47938c31848c9338"
      "ab496c47b153755738234e3b7142da3fe28f4cda1f206d9a2e891ab5"
      "52513610eff2fcfcd9f99e3fbfe40568e64003003c02b91d3dccfdd281ced2df0355abe3"
      "5a6d9d565d0eee83bdcabe02ffbef49ee08acc541e70ccc86aa72f18",
      "e5982bfb3c262021524453e22f908046c4a68c58e560388fd8cb12b40ae6d07cdc9d102f"
      "b452069289bcc9302a072b3df6b5f5e7dddb528fcb063d3a35fc83f1"
      "113a922412e2884889b99f089762ee4da847926c40618fc85089185a5e42a7d9044683a3"
      "17507f63a0637d6c1aa8371adaf6c8e9f6a13334ed211ce8af4ccb1e",
      "670ba045327d7d8174ce05348ff51e1c935848e862493d84cf0857c8b221c32ac22e942a"
      "f5e73332634371421955744a0e58599759c3b9b7d5e5c9065d0a3cc8"
      "ee3e422256486537b8480994f3f8b4631e0f1af3c8115fa46e44fe5c3dfc6ce4abe277a9"
      "1ed6dc415e0c60f33d3cde5297babf59bfbff0c6afab499b7cd75fbf",
      "74dae42becb6f8767d9f9f36f0756ab8b2de1be3a81b064ef079649cccdc941ebe354afd"
      "6703cfa63c4043dcd6f3fff78bf5fd226e38f7b6bad4ff1beaba14b8"
      "274880222a150a520e6eaf4f801debe047235f15bf4b7550d17ed921dafa9e3d6bb93f7c"
      "bbbc78de265f61ff7a7f905eff5d7876ca4ece5f0bec394168d8a7fd",
      "eedfdf1f7e03883e4b9c",
      ""};
  nameCaptureInfo = NULL;
  emlrtNameCaptureMxArrayR2016a(&data[0], 3432U, &nameCaptureInfo);
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
      emlrtCreateStructMatrix(1, 1, 8, (const char_T **)&epFieldName[0]);
  xInputs = emlrtCreateLogicalMatrix(1, 3);
  emlrtSetField(xEntryPoints, 0, "QualifiedName",
                emlrtMxCreateString("student_stop_primitive"));
  emlrtSetField(xEntryPoints, 0, "NumberOfInputs",
                emlrtMxCreateDoubleScalar(3.0));
  emlrtSetField(xEntryPoints, 0, "NumberOfOutputs",
                emlrtMxCreateDoubleScalar(3.0));
  emlrtSetField(xEntryPoints, 0, "ConstantInputs", xInputs);
  emlrtSetField(
      xEntryPoints, 0, "ResolvedFilePath",
      emlrtMxCreateString("/Users/alessandrobianchiceriani/Desktop/"
                          "Scrivania_MPB/ALE_VARIE_DONTTOUCH/UNITN/MAGISTRALE/"
                          "Secondo_Anno/IVAD/Repos/basic_agent_S"
                          "T/matlab/student_stop_primitive.m"));
  emlrtSetField(xEntryPoints, 0, "TimeStamp",
                emlrtMxCreateDoubleScalar(739906.546412037));
  emlrtSetField(xEntryPoints, 0, "Constructor",
                emlrtMxCreateLogicalScalar(false));
  emlrtSetField(xEntryPoints, 0, "Visible", emlrtMxCreateLogicalScalar(true));
  xResult =
      emlrtCreateStructMatrix(1, 1, 9, (const char_T **)&propFieldName[0]);
  emlrtSetField(xResult, 0, "Version",
                emlrtMxCreateString("24.2.0.2712019 (R2024b)"));
  emlrtSetField(xResult, 0, "ResolvedFunctions",
                (mxArray *)c_emlrtMexFcnResolvedFunctionsI());
  emlrtSetField(xResult, 0, "Checksum",
                emlrtMxCreateString("GHx47UXU03F3AHNANNvmwF"));
  emlrtSetField(xResult, 0, "EntryPoints", xEntryPoints);
  return xResult;
}

/* End of code generation (_coder_student_stop_primitive_info.c) */
