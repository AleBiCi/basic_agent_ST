% 
% Test script for primitives codegen
% 

v0 = 3.0; a0 = 5.9; sf = 100.0; vfmin = 2.0; vfmax = 15.0; Tmin = 100.0; Tmax = 300.0;

% Pass primitive
% [coeffsT2, v2, T2, coeffsT1, v1, T1] = student_pass_primitive(v0, a0, sf, vfmin, vfmax, Tmin, Tmax);
[coeffsT2, v2, T2, coeffsT1, v1, T1] = student_pass_primitive(3.0, 5.9, 100.0, 2.0, 15.0, 100.0, 300.0);

% Stop primitive
% [coefs,maxsf,tf] = student_stop_primitive(v0,a0,sf);
[coefs,maxsf,tf] = student_stop_primitive(3.0, 5.0, 100.0);