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
  const char_T *data[10] = {
      "789ced9b4d73d3461880d70c1f3db4e02903ffa1171a27b1a9b929fe50446419fc417010"
      "08d95ed94a56966c492e61a6bfa13fa3a75e3a3d73603a1d66b8c081"
      "81e9f4e003c30fe8ad9d1e90e4288e84776c46b21289ddcb66f56af7fd58f9f1eb771590"
      "622a2900c065306d3fb0d3fe9ba371faa83f07bccd2f4ff9aea7bcb7",
      "830be0fc5cf9cf477d471d18f089311d0c44051ecfecaa8a3c100746e35083600475158d"
      "61d7914832820d5981f59303ce1e29e513a2e3812db2ff2ef461e7a0"
      "6e2a60d4d76716a29383e378bcc3f87b7ec9788c31f148fbe40f4a0f0bb7f8a60e473a4f"
      "21a8ebe2a03b52f922d40f0c55e39b1cd3e0f80a4533f5468d624b7c",
      "bd54a872c5aa40715c9567ee5145be063555e7dba22e7704b1070786506ff08a6820b1cd"
      "eb86d9b5af68a2ae0bda485664431ec31bcacccf2701fdbcbac04f57"
      "2e597b8904553304c3da11c720e7ba6bc7e380765cc4da31957455b38d6078fbabfbc6c0"
      "779f2b5fedfece89aabdb960715caf2ce9a7bf9fddff95d37f50fe4e",
      "45a9ef8df2f2f728f5b9edb4f405fd7c5ec7e84bfbe43777e94cbf9e330f1969adbc878a"
      "798a4385ed991d7716e8596407c08ca35aff4be1b916d0cfaf17f8e9"
      "ca1dde585fd0c218227b7c5a1c7f15501fc2eaf3ca57bbaf27a339057852f9fd9df49af0"
      "1b84cfef5a3197dde8655b92cc2095d172951a2b154a84df71e3f7e3",
      "807e5ec2fa3995b8c83e2d5ebf08a8af87d5e795af761f67a84e2aa79fadbd259c06e173"
      "7ac79073a859dadaccf5a43c5da7b26663ab4df2ecd8713ae8f3f2ed"
      "023f5df9ec17bec59b59d924b6759311569f571e55ddc48dea51d924b13c4fd1cf09cf41"
      "f83cd7e8b6f274f3705081ddad3da6b6759b2b73770b84e771e379d0",
      "ba09eefc23ed937754280948d60d41320720be751305abcf2b5fedbe7aa2e9203ca9fcce"
      "bcfffeb728f5b92de9fc2e8bf0805daf0c0d98cdb05256cd8bda53b6"
      "9c1c7e4f30f3978de32166fdb44f7e1afc16f6d7dcb42df2e766cec99b654e7cf3f2b3f1"
      "3d3d3faacbd4c7e2caf56717feff334a7d6e4b3ad7b73798dcbec8e6",
      "e51fb3cdaad832bb2cd2e804e5e513cc7cc2f5f9fe5e5be0af2bffb432e0603db65c37b1"
      "fabcf2a8eb2dee1e2795ebbffcfb0fe13a089febbb3b3976871e9932"
      "b5aeb4606b5cca64d87282eae713ccfca4729dd45d3e4f1fa9bb90bacb2af591ba4b38eb"
      "27e3f7b8cb71dd5aeb4cbd3f6e1be45c8f6b5e7e26df1fb7a31af1fb",
      "e3ef23e6f923314bf272b0029e8ffb43b596dbc87746d47e71af092b1d8d49505efea5f0"
      "9ce4e39fa78fe4e3241f5fa53e928f87b3fe04333f9e75152fbfcfd0"
      "39a8631839075d415e1ef139e85f1173fd8fd67fbf46a9cf6d71e5fab2e7599b55eb511c"
      "56ef518d9ea4c10cab37d6e93d40b81e57ae07cdcf71f14dfbfa1981",
      "2cab2d43e29a9f6bbe31f0dde7caa3e2f9349a3794a472fca775c271bb85fe7f43c32eea"
      "f7abbbed22baad6cdfbf39bc5fced209aaaf4c30f309c7e7fb4bea2c"
      "defb489d255c7da4ce326da4ceb2dcfa1f01e9e11332",
      ""};
  nameCaptureInfo = NULL;
  emlrtNameCaptureMxArrayR2016a(&data[0], 19536U, &nameCaptureInfo);
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
      emlrtCreateStructMatrix(1, 8, 8, (const char_T **)&epFieldName[0]);
  xInputs = emlrtCreateLogicalMatrix(1, 7);
  emlrtSetField(xEntryPoints, 0, "QualifiedName", emlrtMxCreateString("a_opt"));
  emlrtSetField(xEntryPoints, 0, "NumberOfInputs",
                emlrtMxCreateDoubleScalar(7.0));
  emlrtSetField(xEntryPoints, 0, "NumberOfOutputs",
                emlrtMxCreateDoubleScalar(1.0));
  emlrtSetField(xEntryPoints, 0, "ConstantInputs", xInputs);
  emlrtSetField(
      xEntryPoints, 0, "ResolvedFilePath",
      emlrtMxCreateString(
          "C:\\Users\\Alessandro\\Desktop\\UNITN\\MAGISTRALE\\SECONDO_"
          "ANNO\\IVAD\\Repos\\basic_agent_ST\\matlab\\a_opt.m"));
  emlrtSetField(xEntryPoints, 0, "TimeStamp",
                emlrtMxCreateDoubleScalar(739966.76216435188));
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
      emlrtMxCreateString(
          "C:\\Users\\Alessandro\\Desktop\\UNITN\\MAGISTRALE\\SECONDO_"
          "ANNO\\IVAD\\Repos\\basic_agent_ST\\matlab\\coef_list_fun.m"));
  emlrtSetField(xEntryPoints, 1, "TimeStamp",
                emlrtMxCreateDoubleScalar(739966.76216435188));
  emlrtSetField(xEntryPoints, 1, "Constructor",
                emlrtMxCreateLogicalScalar(false));
  emlrtSetField(xEntryPoints, 1, "Visible", emlrtMxCreateLogicalScalar(true));
  xInputs = emlrtCreateLogicalMatrix(1, 7);
  emlrtSetField(xEntryPoints, 2, "QualifiedName", emlrtMxCreateString("s_opt"));
  emlrtSetField(xEntryPoints, 2, "NumberOfInputs",
                emlrtMxCreateDoubleScalar(7.0));
  emlrtSetField(xEntryPoints, 2, "NumberOfOutputs",
                emlrtMxCreateDoubleScalar(1.0));
  emlrtSetField(xEntryPoints, 2, "ConstantInputs", xInputs);
  emlrtSetField(
      xEntryPoints, 2, "ResolvedFilePath",
      emlrtMxCreateString(
          "C:\\Users\\Alessandro\\Desktop\\UNITN\\MAGISTRALE\\SECONDO_"
          "ANNO\\IVAD\\Repos\\basic_agent_ST\\matlab\\s_opt.m"));
  emlrtSetField(xEntryPoints, 2, "TimeStamp",
                emlrtMxCreateDoubleScalar(739966.76215277775));
  emlrtSetField(xEntryPoints, 2, "Constructor",
                emlrtMxCreateLogicalScalar(false));
  emlrtSetField(xEntryPoints, 2, "Visible", emlrtMxCreateLogicalScalar(true));
  xInputs = emlrtCreateLogicalMatrix(1, 7);
  emlrtSetField(xEntryPoints, 3, "QualifiedName",
                emlrtMxCreateString("student_pass_primitive"));
  emlrtSetField(xEntryPoints, 3, "NumberOfInputs",
                emlrtMxCreateDoubleScalar(7.0));
  emlrtSetField(xEntryPoints, 3, "NumberOfOutputs",
                emlrtMxCreateDoubleScalar(6.0));
  emlrtSetField(xEntryPoints, 3, "ConstantInputs", xInputs);
  emlrtSetField(
      xEntryPoints, 3, "ResolvedFilePath",
      emlrtMxCreateString(
          "C:\\Users\\Alessandro\\Desktop\\UNITN\\MAGISTRALE\\SECONDO_"
          "ANNO\\IVAD\\Repos\\basic_agent_ST\\matlab\\student_pass_primitive."
          "m"));
  emlrtSetField(xEntryPoints, 3, "TimeStamp",
                emlrtMxCreateDoubleScalar(739934.75704861106));
  emlrtSetField(xEntryPoints, 3, "Constructor",
                emlrtMxCreateLogicalScalar(false));
  emlrtSetField(xEntryPoints, 3, "Visible", emlrtMxCreateLogicalScalar(true));
  xInputs = emlrtCreateLogicalMatrix(1, 5);
  emlrtSetField(xEntryPoints, 4, "QualifiedName",
                emlrtMxCreateString("student_pass_primitive_j0"));
  emlrtSetField(xEntryPoints, 4, "NumberOfInputs",
                emlrtMxCreateDoubleScalar(5.0));
  emlrtSetField(xEntryPoints, 4, "NumberOfOutputs",
                emlrtMxCreateDoubleScalar(3.0));
  emlrtSetField(xEntryPoints, 4, "ConstantInputs", xInputs);
  emlrtSetField(
      xEntryPoints, 4, "ResolvedFilePath",
      emlrtMxCreateString(
          "C:\\Users\\Alessandro\\Desktop\\UNITN\\MAGISTRALE\\SECONDO_"
          "ANNO\\IVAD\\Repos\\basic_agent_ST\\matlab\\student_pass_primitive_"
          "j0."
          "m"));
  emlrtSetField(xEntryPoints, 4, "TimeStamp",
                emlrtMxCreateDoubleScalar(739966.76650462963));
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
          "ANNO\\IVAD\\Repos\\basic_agent_ST\\matlab\\student_stop_primitive."
          "m"));
  emlrtSetField(xEntryPoints, 5, "TimeStamp",
                emlrtMxCreateDoubleScalar(739966.732650463));
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
          "ANNO\\IVAD\\Repos\\basic_agent_ST\\matlab\\student_stop_primitive_"
          "j0."
          "m"));
  emlrtSetField(xEntryPoints, 6, "TimeStamp",
                emlrtMxCreateDoubleScalar(739966.76707175921));
  emlrtSetField(xEntryPoints, 6, "Constructor",
                emlrtMxCreateLogicalScalar(false));
  emlrtSetField(xEntryPoints, 6, "Visible", emlrtMxCreateLogicalScalar(true));
  xInputs = emlrtCreateLogicalMatrix(1, 7);
  emlrtSetField(xEntryPoints, 7, "QualifiedName", emlrtMxCreateString("v_opt"));
  emlrtSetField(xEntryPoints, 7, "NumberOfInputs",
                emlrtMxCreateDoubleScalar(7.0));
  emlrtSetField(xEntryPoints, 7, "NumberOfOutputs",
                emlrtMxCreateDoubleScalar(1.0));
  emlrtSetField(xEntryPoints, 7, "ConstantInputs", xInputs);
  emlrtSetField(
      xEntryPoints, 7, "ResolvedFilePath",
      emlrtMxCreateString(
          "C:\\Users\\Alessandro\\Desktop\\UNITN\\MAGISTRALE\\SECONDO_"
          "ANNO\\IVAD\\Repos\\basic_agent_ST\\matlab\\v_opt.m"));
  emlrtSetField(xEntryPoints, 7, "TimeStamp",
                emlrtMxCreateDoubleScalar(739966.76216435188));
  emlrtSetField(xEntryPoints, 7, "Constructor",
                emlrtMxCreateLogicalScalar(false));
  emlrtSetField(xEntryPoints, 7, "Visible", emlrtMxCreateLogicalScalar(true));
  xResult =
      emlrtCreateStructMatrix(1, 1, 9, (const char_T **)&propFieldName[0]);
  emlrtSetField(xResult, 0, "Version",
                emlrtMxCreateString("24.2.0.2773142 (R2024b) Update 2"));
  emlrtSetField(xResult, 0, "ResolvedFunctions",
                (mxArray *)c_emlrtMexFcnResolvedFunctionsI());
  emlrtSetField(xResult, 0, "Checksum",
                emlrtMxCreateString("L4UIQb0D0FTZqJRYckSnCC"));
  emlrtSetField(xResult, 0, "EntryPoints", xEntryPoints);
  return xResult;
}

/* End of code generation (_coder_primitives_info.c) */
