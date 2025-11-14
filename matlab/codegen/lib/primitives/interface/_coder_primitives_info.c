/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_primitives_info.c
 *
 * Code generation for function 'a_opt'
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
  const char_T *data[8] = {
      "789ced59cb6ed340149da0162af10a5081f80644ad86beb24172e3a4b594a46dec246d01"
      "398e33c1d3f8917ac651fa110809a91bd85009898f80357fc13fb044"
      "6283f37063875a4965e4a623dfcde4fad83977eedcf8e4ce80045f480000ee8181fd7a36"
      "18ef0efde470bc01fc368e2786e3ad31dfb57930e77bcec5df0f47c5",
      "3408ec928163c83a3c7fb261eac8900d229eb421b02036b50e6cf49126d2a0887428789d"
      "62cfd3731ee8dce941bdcf19152a2dc1d681a5e251849ad739cfc742"
      "e2e2f9ce4d998f2f01f9488ee1afb26f9832861666640d622c1b0dcbac23d95054a440cb"
      "f980180ee21631db8ca058a8e35c90a5c2ee26c3e6b352852df15989",
      "db298ae24e39b3cd948bbc58640aec162f8825e70646804e7e1ba6c41a86c9f01596634a"
      "b06d62a62e63a448f25b68104910195d269a5c6730b11bbd2b6d1963"
      "a96d211d11d4814bba372fdd80794f9b974713f2e2e24d67ed35c96c1389382bd80f0978"
      "e3a8858ce366601c03a461da750dfebf7a380be4f3e3b3540f17acc1",
      "a018c0e475b83f655ec6c7d1fd0bfdf1fbcb27284abef98dfc8328f95cbb2abeb0bfe7c7"
      "017cc9311cab58e5b2cb95ee7a71953db08e116c95b39ba3387627f0"
      "4c8a0304f8517d7fac1717eb453b60ded3e6e5ce84bcb878ff0de5fc61903a50f3f2d742"
      "f25f562740c83af818c8e7c767a90ebcb9770582567df8b174b61825",
      "9f6bb4ebc3724a663954dcdf4b57bbfbea1aaacaa452cac4fa40bb3ed402e63d6d5e82fa"
      "cee4107125e1aaf4e07748bed3403e3f3e4bebee95025a75e0e7bbd7"
      "4fa3e4738d761d7891db1156ec75d5e6f7d59200855596cd6fe4621da05d07c2d6d7c309"
      "7971f1d19e86f3861a6d2b5ddb7ee173209f1f9fa57af8770d96745a",
      "75e2c3d7c5db51f2b946bb4e348e05cbe0ccf2a19159b335bc5e5e3b24f9b85fa05e27c2"
      "ee278d9f438de7c5c5151336250d6122356d035c5f7df814c8e7c767"
      "a90e7cb91f6e2845f53efbb610ad3e9cf2563a4a3ed768d7079cd3f1b65225bb1c6a6e75"
      "5b5ab5a091c256ac0f34e80376d866ec7cba1712f0c6510b19477c3e",
      "7de9f3e9de1a447e3ecd6eff51a3e47bbeba928c92cf35daf5c2c8a74f8ed2b983bd5626"
      "9be28fb8d4614a8dcfa7e9d78bb89fb81c5fdc4ff8e71ff71303a35d"
      "1f68ed27fe023df11545",
      ""};
  nameCaptureInfo = NULL;
  emlrtNameCaptureMxArrayR2016a(&data[0], 11272U, &nameCaptureInfo);
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
      emlrtCreateStructMatrix(1, 5, 8, (const char_T **)&epFieldName[0]);
  xInputs = emlrtCreateLogicalMatrix(1, 7);
  emlrtSetField(xEntryPoints, 0, "QualifiedName", emlrtMxCreateString("a_opt"));
  emlrtSetField(xEntryPoints, 0, "NumberOfInputs",
                emlrtMxCreateDoubleScalar(7.0));
  emlrtSetField(xEntryPoints, 0, "NumberOfOutputs",
                emlrtMxCreateDoubleScalar(1.0));
  emlrtSetField(xEntryPoints, 0, "ConstantInputs", xInputs);
  emlrtSetField(
      xEntryPoints, 0, "ResolvedFilePath",
      emlrtMxCreateString("/Users/alessandrobianchiceriani/Desktop/"
                          "Scrivania_MPB/ALE_VARIE_DONTTOUCH/UNITN/MAGISTRALE/"
                          "Secondo_Anno/IVAD/Repos/basic_agent_S"
                          "T/matlab/a_opt.m"));
  emlrtSetField(xEntryPoints, 0, "TimeStamp",
                emlrtMxCreateDoubleScalar(739935.491238426));
  emlrtSetField(xEntryPoints, 0, "Constructor",
                emlrtMxCreateLogicalScalar(false));
  emlrtSetField(xEntryPoints, 0, "Visible", emlrtMxCreateLogicalScalar(true));
  xInputs = emlrtCreateLogicalMatrix(1, 6);
  emlrtSetField(xEntryPoints, 1, "QualifiedName",
                emlrtMxCreateString("coef_list_fun"));
  emlrtSetField(xEntryPoints, 1, "NumberOfInputs",
                emlrtMxCreateDoubleScalar(6.0));
  emlrtSetField(xEntryPoints, 1, "NumberOfOutputs",
                emlrtMxCreateDoubleScalar(1.0));
  emlrtSetField(xEntryPoints, 1, "ConstantInputs", xInputs);
  emlrtSetField(
      xEntryPoints, 1, "ResolvedFilePath",
      emlrtMxCreateString("/Users/alessandrobianchiceriani/Desktop/"
                          "Scrivania_MPB/ALE_VARIE_DONTTOUCH/UNITN/MAGISTRALE/"
                          "Secondo_Anno/IVAD/Repos/basic_agent_S"
                          "T/matlab/coef_list_fun.m"));
  emlrtSetField(xEntryPoints, 1, "TimeStamp",
                emlrtMxCreateDoubleScalar(739935.491238426));
  emlrtSetField(xEntryPoints, 1, "Constructor",
                emlrtMxCreateLogicalScalar(false));
  emlrtSetField(xEntryPoints, 1, "Visible", emlrtMxCreateLogicalScalar(true));
  xInputs = emlrtCreateLogicalMatrix(1, 7);
  emlrtSetField(xEntryPoints, 2, "QualifiedName",
                emlrtMxCreateString("student_pass_primitive"));
  emlrtSetField(xEntryPoints, 2, "NumberOfInputs",
                emlrtMxCreateDoubleScalar(7.0));
  emlrtSetField(xEntryPoints, 2, "NumberOfOutputs",
                emlrtMxCreateDoubleScalar(6.0));
  emlrtSetField(xEntryPoints, 2, "ConstantInputs", xInputs);
  emlrtSetField(
      xEntryPoints, 2, "ResolvedFilePath",
      emlrtMxCreateString("/Users/alessandrobianchiceriani/Desktop/"
                          "Scrivania_MPB/ALE_VARIE_DONTTOUCH/UNITN/MAGISTRALE/"
                          "Secondo_Anno/IVAD/Repos/basic_agent_S"
                          "T/matlab/student_pass_primitive.m"));
  emlrtSetField(xEntryPoints, 2, "TimeStamp",
                emlrtMxCreateDoubleScalar(739928.64553240745));
  emlrtSetField(xEntryPoints, 2, "Constructor",
                emlrtMxCreateLogicalScalar(false));
  emlrtSetField(xEntryPoints, 2, "Visible", emlrtMxCreateLogicalScalar(true));
  xInputs = emlrtCreateLogicalMatrix(1, 3);
  emlrtSetField(xEntryPoints, 3, "QualifiedName",
                emlrtMxCreateString("student_stop_primitive"));
  emlrtSetField(xEntryPoints, 3, "NumberOfInputs",
                emlrtMxCreateDoubleScalar(3.0));
  emlrtSetField(xEntryPoints, 3, "NumberOfOutputs",
                emlrtMxCreateDoubleScalar(3.0));
  emlrtSetField(xEntryPoints, 3, "ConstantInputs", xInputs);
  emlrtSetField(
      xEntryPoints, 3, "ResolvedFilePath",
      emlrtMxCreateString("/Users/alessandrobianchiceriani/Desktop/"
                          "Scrivania_MPB/ALE_VARIE_DONTTOUCH/UNITN/MAGISTRALE/"
                          "Secondo_Anno/IVAD/Repos/basic_agent_S"
                          "T/matlab/student_stop_primitive.m"));
  emlrtSetField(xEntryPoints, 3, "TimeStamp",
                emlrtMxCreateDoubleScalar(739912.47503472224));
  emlrtSetField(xEntryPoints, 3, "Constructor",
                emlrtMxCreateLogicalScalar(false));
  emlrtSetField(xEntryPoints, 3, "Visible", emlrtMxCreateLogicalScalar(true));
  xInputs = emlrtCreateLogicalMatrix(1, 7);
  emlrtSetField(xEntryPoints, 4, "QualifiedName", emlrtMxCreateString("v_opt"));
  emlrtSetField(xEntryPoints, 4, "NumberOfInputs",
                emlrtMxCreateDoubleScalar(7.0));
  emlrtSetField(xEntryPoints, 4, "NumberOfOutputs",
                emlrtMxCreateDoubleScalar(1.0));
  emlrtSetField(xEntryPoints, 4, "ConstantInputs", xInputs);
  emlrtSetField(
      xEntryPoints, 4, "ResolvedFilePath",
      emlrtMxCreateString("/Users/alessandrobianchiceriani/Desktop/"
                          "Scrivania_MPB/ALE_VARIE_DONTTOUCH/UNITN/MAGISTRALE/"
                          "Secondo_Anno/IVAD/Repos/basic_agent_S"
                          "T/matlab/v_opt.m"));
  emlrtSetField(xEntryPoints, 4, "TimeStamp",
                emlrtMxCreateDoubleScalar(739935.491238426));
  emlrtSetField(xEntryPoints, 4, "Constructor",
                emlrtMxCreateLogicalScalar(false));
  emlrtSetField(xEntryPoints, 4, "Visible", emlrtMxCreateLogicalScalar(true));
  xResult =
      emlrtCreateStructMatrix(1, 1, 9, (const char_T **)&propFieldName[0]);
  emlrtSetField(xResult, 0, "Version",
                emlrtMxCreateString("24.2.0.2712019 (R2024b)"));
  emlrtSetField(xResult, 0, "ResolvedFunctions",
                (mxArray *)c_emlrtMexFcnResolvedFunctionsI());
  emlrtSetField(xResult, 0, "Checksum",
                emlrtMxCreateString("CjjAt2mck50voG0hyrKjpC"));
  emlrtSetField(xResult, 0, "EntryPoints", xEntryPoints);
  return xResult;
}

/* End of code generation (_coder_primitives_info.c) */
