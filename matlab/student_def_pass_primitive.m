%%
%           Agent Logic
%         Pass Primitive
%              2025
%
%
%% Determine the optimal final velocity to pass
final_opt_vel_pass_var = solve(diff(subs(total_cost_var, af, 0.), vf) == 0, vf);
final_opt_vel_pass_fun = matlabFunction(final_opt_vel_pass_var,'Vars',[v0,a0,sf,T],'File','final_opt_vel_pass.m');
% - Use functions 'solve', 'diff', 'subs', and 'matlabFunction'

%% Determine the optimal time to reach to final velocity
syms z; % We create this variable for more stable solution of the equations

final_opt_time_pass_var = 1./solve(diff(subs(total_cost_var, [af, T], [0., 1/z]), z) == 0, z);
final_opt_time_pass_fun = matlabFunction(final_opt_time_pass_var,'Vars',[v0,a0,sf,vf],'File','final_opt_time_pass.m');
% - Use the function solve 
%   Use 'subs' function to solve the system using the frequency instead of the time

%% Determine the time to reach the minimum velocity
time_min_vel_var = solve(diff(final_opt_vel_pass_var, T) == 0, T);
time_min_vel_fun = matlabFunction(time_min_vel_var,'Vars',[v0,a0,sf,T],'File','time_min_vel.m');
min_vel_var = simplify(subs(final_opt_vel_pass_var, T, time_min_vel_var(2)));
min_vel_fun = matlabFunction(min_vel_var,'Vars',[v0,a0,sf,T],'File','min_vel.m');
% - Use functions 'solve', 'diff', 'subs', and 'matlabFunction'
