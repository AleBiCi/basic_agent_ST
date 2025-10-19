%%
%           Agent Logic
%         Stop Primitive
%              2025
%
%
%% Determine the optimal time to stop %%
syms z; % We create this variable for more stable solution of the equations

%% Simbolic definition of the optimal time to stop
final_opt_time_stop_var = 1./solve(diff(subs(total_cost_var, [vf, af, T], [0., 0., 1/z]), z) == 0, z);
% - Use the function solve 
%   Use 'subs' function to solve the system using the frequency instead of the time
% subs T (tf) w/ 1/z, derive wrt z and solve ode == 0 wrt z, then invert
% the result to get final result (Txf)

%% Export the solution in a matlab function
final_opt_time_stop_fun = matlabFunction(final_opt_time_stop_var(1), 'Vars', [v0, a0, sf], 'File', 'final_opt_time_stop.m');