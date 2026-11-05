%%
%           Agent Logic
%          Pass Primitive J0
%               2025
%
%
%% Pass primitive with j0 = 0 %%
syms j0;
syms z;

% Determine the optimal final velocity to pass function of j0 = 0 %%
final_opt_vel_pass_j0_var = solve(subs(sol_opt.j,t,0)==j0,vf);
final_opt_vel_pass_j0_fun = matlabFunction(final_opt_vel_pass_j0_var,'Vars',[v0,a0,j0,sf,af,T],'File','final_opt_vel_pass_j0.m');

% Determine the optimal time to reach to final velocity with j0 = 0 %%
final_opt_time_pass_j0_var = 1./solve(final_opt_time_pass_fun(v0,a0,sf,subs(final_opt_vel_pass_j0_var,[T,af,j0],[1/z,0,0]))==1/z,z);
final_opt_time_pass_j0_fun = matlabFunction(final_opt_time_pass_j0_var,'Vars',[v0,a0,sf],'File','final_opt_time_pass_j0.m');
