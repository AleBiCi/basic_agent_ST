%%
%           Agent Logic
%         Stop Primitive J0
%              2025
%
%
%% Stop primitive with j0 = 0 %%

%% Determine the optimal final position imposing j0 %%
% Solve for t = 0 since we're computing the optimal position with j0 = 0 at
% the current time instant
opt_pos_eq = subs(sol_opt.j, [t, vf, af], [0., 0., 0.]) == 0;
final_opt_pos_j0_var = simplify(solve(opt_pos_eq, sf));
final_opt_pos_j0_fun = matlabFunction(final_opt_pos_j0_var, 'Vars', [v0, a0, T], 'File', 'final_opt_pos_j0.m');
% - Use functions 'solve', 'diff', 'subs', and 'matlabFunction'

%% Determine the optimal time to stop with j0 = 0 %%
opt_time_stop_eq = final_opt_time_stop(v0, a0, final_opt_pos_j0_fun(v0, a0, T)) == T;
final_opt_time_stop_j0_var = simplify(solve(opt_time_stop_eq, T));
final_opt_time_stop_j0_fun = matlabFunction(final_opt_time_stop_j0_var(2), 'Vars', [v0, a0], 'File', 'final_opt_time_stop_j0.m');
% - Use the function solve  
%   Use the function final_opt_time_stop_fun 
%   Use 'subs' function to solve the system using the frequency instead of the time

