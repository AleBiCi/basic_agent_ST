%%
%           Agent Logic
%         Pass Primitive
%              2025
%
%
%% Determine the optimal final velocity to pass --> CHECKED
eq_vel = diff(subs(total_cost_var, af, 0.), vf) == 0;
final_opt_vel_pass_var = solve(eq_vel, vf);
final_opt_vel_pass_fun = matlabFunction(final_opt_vel_pass_var,'Vars',[v0,a0,sf,T],'File','final_opt_vel_pass.m');
% - Use functions 'solve', 'diff', 'subs', and 'matlabFunction'

%% Determine the optimal time to reach to final velocity --> CHECKED
syms z; % We create this variable for more stable solution of the equations

eq_time_z = subs(total_cost_var, [af, T], [0., 1/z]);
eq_time_z_prime = simplify(diff(eq_time_z, vf)) == 0;
final_opt_time_pass_var = 1./solve(eq_time_z_prime, z); % select fourth solution for this equation
final_opt_time_pass_fun = matlabFunction(final_opt_time_pass_var(4),'Vars',[v0,a0,sf,vf],'File','final_opt_time_pass.m');
% - Use the function solve 
%   Use 'subs' function to solve the system using the frequency instead of the time

%% Determine the time to reach the minimum velocity --> CHECKED
time_min_vel_roots = solve(diff(final_opt_vel_pass_var, T) == 0, T);
% Keep the positive root T* = sqrt(15 sf/(-a0)) (slide 58): test the sign at a0 = -1, sf = 1
time_min_vel_var = time_min_vel_roots(double(subs(time_min_vel_roots, [a0, sf], [-1, 1])) > 0);
time_min_vel_fun = matlabFunction(time_min_vel_var,'Vars',[a0,sf],'File','time_min_vel.m');
min_vel_var = simplify(subs(final_opt_vel_pass_var, T, time_min_vel_var));
min_vel_fun = matlabFunction(min_vel_var,'Vars',[v0,a0,sf],'File','min_vel.m');
% - Use functions 'solve', 'diff', 'subs', and 'matlabFunction'
