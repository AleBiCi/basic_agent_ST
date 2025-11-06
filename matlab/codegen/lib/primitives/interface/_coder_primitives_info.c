/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_primitives_info.c
 *
 * Code generation for function 'student_pass_primitive'
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
      "789ced59cb6ed340149da0029578455474c11a56881a425b14b1caabb1db246d63a76d8a"
      "90e3c7b876e2573c8e49f909242436b00189351bfe8435ff004b2436"
      "380f37b6a99554466e62f96e26d7c7ceb973e7664eae07a4886a0a00701b8cecd7a3d178"
      "6beca7c7e315e0353f9e1a8fd77dbe6357c192e739077f3b1e394d35",
      "61df1c392aa3c0b327794d91544635a9531d0203224db6203f4404498694a440d2edd406"
      "9eb2e582ce9c0134f85c1021d7217b0a3044348950763b67f9584e9d"
      "3fdfa519f3f125201f691ffeb2f40a6b2068208c9121428cca1b1a2b312a274a1c34ec0f"
      "125684a8636a3a46728664d91718baba97c77295127d90ab1325bab8",
      "5ba3a8dd4601c71a3582aa61d55c9920a9ba7d0346423bbfbc46e75455c388835c11ab43"
      "5d4318cb2089a39913a89a3449610a63ca0c8b21b3c70faee80c42b4"
      "6e488a644a165c53dc79e907cc7bd6bcac4cc98b830bf6dacbb4a69bb469afe03024e08e"
      "a315328e6b81718c105eebb132fc7ff5f03990cf8bcf533d9cb306a3",
      "6200d3d7e1ce8c79f18f93fb97872381ff11a3e4fbf173f746947c8e5d165fd8dff36a00"
      "5fda873feb1e2a790dcf56e4c34cb9b2bf5eeed676347c12c7de149e"
      "697180003faaef4ff4e27cbdd003e63d6b5e6e4ec98b830f7728fb0f036d41d9cddf0ac9"
      "7f519d0021ebe043209f179fa73a70e7de1188b8ea43edebc39528f9",
      "1c8bbb3ef48ebb9d0cde3db2f6baf8b6d25345b6ca09f9441fe2ae0fad8079cf9a97a0be"
      "333d461c49b82c3df81d92ef7d209f179fa775774b415c75e08dbe7a"
      "3f4a3ec7165507ee05f0a57d38434aeb4f2c719bed6f72acce57ccf649f708243a10771d"
      "085b5f77a7e4c5c127ef34ec1d6af25a6961fb854f817c5e7c9eeae1",
      "df355853a2dad7f0887562e7dbf7e751f239b6a83a316bbfa0bd6ef3c7cc0155dc10b6b3"
      "a7e56cb351579e26fd42ec7522ecfb24ff39943f2f0ece6950a06509"
      "99b4d053c1e2eac3c7403e2f3e4f75e0c9fdf8855254fbd98b88f5e11dd87c10259f6371"
      "d7070675ac0c4fe0fdada65ae5f84a7e47682ac5441fe2a00fc8669b",
      "b3f3e94148c01d472b641cc9f9f485cfa7076b10f9f9742e62bd78bcb19e8e92cfb1b8eb"
      "855ac99eb6b35bcdfd4ea19421dac5cc71462c25fd44ecf522e9272e"
      "c697f413def927fdc4c8e2ae0f71ed27fe0293c01862",
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
      emlrtCreateStructMatrix(1, 2, 8, (const char_T **)&epFieldName[0]);
  xInputs = emlrtCreateLogicalMatrix(1, 7);
  emlrtSetField(xEntryPoints, 0, "QualifiedName",
                emlrtMxCreateString("student_pass_primitive"));
  emlrtSetField(xEntryPoints, 0, "NumberOfInputs",
                emlrtMxCreateDoubleScalar(7.0));
  emlrtSetField(xEntryPoints, 0, "NumberOfOutputs",
                emlrtMxCreateDoubleScalar(6.0));
  emlrtSetField(xEntryPoints, 0, "ConstantInputs", xInputs);
  emlrtSetField(
      xEntryPoints, 0, "ResolvedFilePath",
      emlrtMxCreateString("/Users/alessandrobianchiceriani/Desktop/"
                          "Scrivania_MPB/ALE_VARIE_DONTTOUCH/UNITN/MAGISTRALE/"
                          "Secondo_Anno/IVAD/Repos/basic_agent_S"
                          "T/matlab/student_pass_primitive.m"));
  emlrtSetField(xEntryPoints, 0, "TimeStamp",
                emlrtMxCreateDoubleScalar(739914.51457175927));
  emlrtSetField(xEntryPoints, 0, "Constructor",
                emlrtMxCreateLogicalScalar(false));
  emlrtSetField(xEntryPoints, 0, "Visible", emlrtMxCreateLogicalScalar(true));
  xInputs = emlrtCreateLogicalMatrix(1, 3);
  emlrtSetField(xEntryPoints, 1, "QualifiedName",
                emlrtMxCreateString("student_stop_primitive"));
  emlrtSetField(xEntryPoints, 1, "NumberOfInputs",
                emlrtMxCreateDoubleScalar(3.0));
  emlrtSetField(xEntryPoints, 1, "NumberOfOutputs",
                emlrtMxCreateDoubleScalar(3.0));
  emlrtSetField(xEntryPoints, 1, "ConstantInputs", xInputs);
  emlrtSetField(
      xEntryPoints, 1, "ResolvedFilePath",
      emlrtMxCreateString("/Users/alessandrobianchiceriani/Desktop/"
                          "Scrivania_MPB/ALE_VARIE_DONTTOUCH/UNITN/MAGISTRALE/"
                          "Secondo_Anno/IVAD/Repos/basic_agent_S"
                          "T/matlab/student_stop_primitive.m"));
  emlrtSetField(xEntryPoints, 1, "TimeStamp",
                emlrtMxCreateDoubleScalar(739912.47503472224));
  emlrtSetField(xEntryPoints, 1, "Constructor",
                emlrtMxCreateLogicalScalar(false));
  emlrtSetField(xEntryPoints, 1, "Visible", emlrtMxCreateLogicalScalar(true));
  xResult =
      emlrtCreateStructMatrix(1, 1, 9, (const char_T **)&propFieldName[0]);
  emlrtSetField(xResult, 0, "Version",
                emlrtMxCreateString("24.2.0.2712019 (R2024b)"));
  emlrtSetField(xResult, 0, "ResolvedFunctions",
                (mxArray *)c_emlrtMexFcnResolvedFunctionsI());
  emlrtSetField(xResult, 0, "Checksum",
                emlrtMxCreateString("c8uHyDgCIjkuXNutMsmNnH"));
  emlrtSetField(xResult, 0, "EntryPoints", xEntryPoints);
  return xResult;
}

/* End of code generation (_coder_primitives_info.c) */
