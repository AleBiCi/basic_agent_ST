%%
%   Rebuild the primitives library
%
%   1) runs the derivation scripts, which rewrite every matlabFunction output
%      (s/v/a/j_opt, coeffs_*_opt, coef_list_fun, total_cost, final_opt_*, time_min_vel, min_vel)
%   2) generates the C library in codegen/lib/primitives (primitives.c/.h, single file)
%
%   Replaces the GUI project primitives.prj, whose file paths still point to the old
%   Mac and Windows folders. Run from the matlab/ folder.
%%

student_def_optimal_control;    % clears the workspace first
student_def_stop_primitive;
student_def_stop_primitive_j0;
student_def_pass_primitive;
student_def_pass_primitive_j0;
rehash;                         % pick up the .m files just written

cfg = coder.config('lib');
cfg.TargetLang          = 'C';
cfg.FilePartitionMethod = 'SingleFile';
cfg.GenerateReport      = true;
cfg.GenCodeOnly         = true;     % C sources only: CMake compiles them (no MATLAB toolchain needed)

d  = 0;             % scalar double
c6 = zeros(1,6);    % coefficient vector [c0 = 0, c1 = v0, c2 = a0, c3 = j0, c4, c5]

codegen('-config', cfg, '-d', 'codegen/lib/primitives', '-o', 'primitives', ...
    'coef_list_fun',             '-args', {d, d, d, d, d, d}, ...
    'coeffs_s_opt',              '-args', {d, c6}, ...
    'coeffs_v_opt',              '-args', {d, c6}, ...
    'coeffs_a_opt',              '-args', {d, c6}, ...
    'coeffs_j_opt',              '-args', {d, c6}, ...
    'student_stop_primitive',    '-args', {d, d, d}, ...
    'student_stop_primitive_j0', '-args', {d, d}, ...
    'student_pass_primitive',    '-args', {d, d, d, d, d, d, d}, ...
    'student_pass_primitive_j0', '-args', {d, d, d, d, d});

%% Quick numeric checks against the slides
% The derivation leaves symbolic variables named like the generated files: drop them so the
% names below resolve to the .m functions
clear coeffs_s_opt coeffs_v_opt coeffs_a_opt coeffs_j_opt coef_list_fun
% coeffs_a_opt / coeffs_j_opt must match a_opt / j_opt (slide 35)
m = coef_list_fun(10, 1, 90, 25, 0, 4);
assert(abs(coeffs_a_opt(1.3, m) - a_opt(1.3, 10, 1, 90, 25, 0, 4)) < 1e-9);
assert(abs(coeffs_j_opt(0, m) - m(4)) < 1e-12);
% Stop primitive ends at sf with v = a = 0 (slide 50)
[m, smax, tf] = student_stop_primitive(10, -1, 40);
assert(abs(coeffs_s_opt(tf, m) - smax) < 1e-6 && abs(coeffs_v_opt(tf, m)) < 1e-6);
% T* > 0 and v* from slide 58 for a0 < 0
assert(time_min_vel(-1, 40) > 0);
assert(abs(min_vel(10, -1, 40) - (2*sqrt(15)*sqrt(40) - 70)/8) < 1e-9);
disp('generate_primitives: OK');
