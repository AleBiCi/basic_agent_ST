/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * main.c
 *
 * Code generation for function 'main'
 *
 */

/*************************************************************************/
/* This automatically generated example C main file shows how to call    */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/

/* Include files */
#include "main.h"
#include "primitives.h"
#include "rt_nonfinite.h"

/* Function Declarations */
static double argInit_real_T(void);

/* Function Definitions */
static double argInit_real_T(void)
{
  return 0.0;
}

int main(int argc, char **argv)
{
  (void)argc;
  (void)argv;
  /* The initialize function is being called automatically from your entry-point
   * function. So, a call to initialize is not included here. */
  /* Invoke the entry-point functions.
You can call entry-point functions multiple times. */
  main_a_opt();
  main_coef_list_fun();
  main_s_opt();
  main_student_pass_primitive();
  main_student_pass_primitive_j0();
  main_student_stop_primitive();
  main_student_stop_primitive_j0();
  main_v_opt();
  /* Terminate the application.
You do not need to do this more than one time. */
  primitives_terminate();
  return 0;
}

void main_a_opt(void)
{
  double t_tmp;
  /* Initialize function 'a_opt' input arguments. */
  t_tmp = argInit_real_T();
  /* Call the entry-point 'a_opt'. */
  t_tmp = a_opt(t_tmp, t_tmp, t_tmp, t_tmp, t_tmp, t_tmp, t_tmp);
}

void main_coef_list_fun(void)
{
  double coef_list_var[6];
  double v0_tmp;
  /* Initialize function 'coef_list_fun' input arguments. */
  v0_tmp = argInit_real_T();
  /* Call the entry-point 'coef_list_fun'. */
  coef_list_fun(v0_tmp, v0_tmp, v0_tmp, v0_tmp, v0_tmp, v0_tmp, coef_list_var);
}

void main_s_opt(void)
{
  double t_tmp;
  /* Initialize function 's_opt' input arguments. */
  t_tmp = argInit_real_T();
  /* Call the entry-point 's_opt'. */
  t_tmp = s_opt(t_tmp, t_tmp, t_tmp, t_tmp, t_tmp, t_tmp, t_tmp);
}

void main_student_pass_primitive(void)
{
  double coeffsT1_data[36];
  double coeffsT2_data[36];
  double T1;
  double T2;
  double v0_tmp;
  double v1;
  double v2;
  int coeffsT1_size[2];
  int coeffsT2_size[2];
  /* Initialize function 'student_pass_primitive' input arguments. */
  v0_tmp = argInit_real_T();
  /* Call the entry-point 'student_pass_primitive'. */
  student_pass_primitive(v0_tmp, v0_tmp, v0_tmp, v0_tmp, v0_tmp, v0_tmp, v0_tmp,
                         coeffsT2_data, coeffsT2_size, &v2, &T2, coeffsT1_data,
                         coeffsT1_size, &v1, &T1);
}

void main_student_pass_primitive_j0(void)
{
  double coefsj0[6];
  double tfj0;
  double vfj0;
  /* Initialize function 'student_pass_primitive_j0' input arguments. */
  vfj0 = argInit_real_T();
  /* Call the entry-point 'student_pass_primitive_j0'. */
  student_pass_primitive_j0(vfj0, vfj0, vfj0, vfj0, vfj0, coefsj0, &vfj0,
                            &tfj0);
}

void main_student_stop_primitive(void)
{
  double coefs[6];
  double maxsf;
  double tf;
  double v0_tmp;
  /* Initialize function 'student_stop_primitive' input arguments. */
  v0_tmp = argInit_real_T();
  /* Call the entry-point 'student_stop_primitive'. */
  student_stop_primitive(v0_tmp, v0_tmp, v0_tmp, coefs, &maxsf, &tf);
}

void main_student_stop_primitive_j0(void)
{
  double coefsj0[6];
  double sfj0;
  double tfj0;
  double v0_tmp;
  /* Initialize function 'student_stop_primitive_j0' input arguments. */
  v0_tmp = argInit_real_T();
  /* Call the entry-point 'student_stop_primitive_j0'. */
  student_stop_primitive_j0(v0_tmp, v0_tmp, coefsj0, &sfj0, &tfj0);
}

void main_v_opt(void)
{
  double t_tmp;
  /* Initialize function 'v_opt' input arguments. */
  t_tmp = argInit_real_T();
  /* Call the entry-point 'v_opt'. */
  t_tmp = v_opt(t_tmp, t_tmp, t_tmp, t_tmp, t_tmp, t_tmp, t_tmp);
}

/* End of code generation (main.c) */
