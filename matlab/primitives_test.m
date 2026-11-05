% 
% Test script for primitives codegen
% 

t = 0.5; v0val = 3.0; a0val = 5.9; sfval = 100.0; vfval = 10.0; vfminval = 2.0; vfmaxval = 15.0; afval = 3.0; Tminval = 100.0; Tmaxval = 300.0;

% s_opt
s_optimal = s_opt(t, v0val, a0val, sfval, vfval, afval, Tmaxval);

% a_opt
a_optimal = a_opt(t, v0val, a0val, sfval, vfval, afval, Tmaxval);

% v_opt
v_optimal = v_opt(t, v0val, a0val, sfval, vfval, afval, Tmaxval);

% coef_list_fun
coeffs = coef_list_fun(v0val, a0val, sfval, vfval, afval, Tmaxval);

% Pass primitive
[coeffsT2, v2, T2, coeffsT1, v1, T1] = student_pass_primitive(v0val, a0val, sfval, vfminval, vfmaxval, Tminval, Tmaxval);

% Stop primitive
[coefs,maxsf,tf] = student_stop_primitive(v0val, -a0val, sfval);

% Pass primitive w/ j0 = 0
[coefsj0, vj0, Tj0] = student_pass_primitive_j0(v0val, a0val, sfval, vfminval, vfmaxval);

% Stop primitive w/ j0 = 0
[coefsj0, xfj0, Tj0] = student_stop_primitive_j0(v0val, -a0val);