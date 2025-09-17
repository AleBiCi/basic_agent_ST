%%
%             Agent Logic
%       Optimal Control Solution
%                2025
%
%%
clear; clc; close all

%% Define a symbolic variables use "syms" (e.g. syms s(t)) and the cost function.
% Symbolic variables
syms s(t) v(t) a(t) u(t) l1(t) l2(t) l3(t);

%% Define the system model by mean the relation between the variables. 
% Plant equations
% ode1 = ...
% ode2 = ...
% ode3 = ...
% - Use "diff()" function to derive the function (e.g. diff(s) == v)

%% Define the Lagrangian and the Hamiltonian.
% Lagrangian
% L = 

% Hamiltonian
% H = L + ...
% - Use "rhs()" function for getting the right hand side 

%% Solve the optimal control problem (Solving the Hamiltonian). 
% Derivative of Hamiltonian w.r.t. the control
% du = ...
syms opt_u;
% Optimal solution
% opt_u = ...
% - use "functionalDerivative()" to derive a function with respect to another function (e.g. du = functionalDerivative(H,u)) 
% - Use "solve" and "subs()" to replace variable inside the equation

%% Write the second optimality condition.
% Second Optimality Condition
% Dl1 = ...
% Dl2 = ...
% Dl3 = ...
% - Use "diff(,t)" to derive w.r.t. the time 
%   and "functionalDerivative()" to detive w.r.t. a variable 

%% Substitute the optimal solution opt_u to state equations
% New system equation
% ode3s = ...
% - Replace the ode3 equation with the solution of the optimal control problem

%% Define the Boundary condition on initial state and Boundary condition on final state
% Boundary condition on initial and final states  
% ICs = ... 
% FCs = ... 
% - Write the condition as a string

%% Find the solution of the OCP imposing the boundary condition
% Solution of the optimal control problem
% sol_opt = ...
% - Use the function "dsolve([], , )" to obtain a solution of the optimal
%   control problem

disp('Optimal polynomial longitudinal position:');
% pretty(sol_opt.s)

disp('Optimal polynomial velocity:');
% pretty(sol_opt.v)

disp('Optimal polynomial acceleration:');
% pretty(sol_opt.a)

% Assign to functions the solutions found

%% Get the optimal control solution
% Obtain optimal control solution
% sol_opt.j = ...
% - Use the "subs" function on the opt_u with the value of l3

%% Export the solution in a matlab function
syms t v0 a0 sf vf af T  
% Create the matlab function for the solutions
% s_opt_fun = matlabFunction(sol_opt.s,'Vars',[t,v0,a0,sf,vf,af,T],'File','s_opt.m');
% v_opt_fun = ...
% a_opt_fun = ...
% j_opt_fun = ...
% - Use the matlabFunction function to generate a matlab function using a
%   symbolic function

%% Create the matlab function for the solution using the coeffs
c = sym('c',[1 6]);
% coeffs_s_opt = ...
% coeffs_v_opt = ...
% coeffs_a_opt = ...
% coeffs_j_opt = ...
% coeffs_s_opt_fun = ...
% coeffs_v_opt_fun = ...
% coeffs_a_opt_fun = ...
% coeffs_j_opt_fun = ...

%% Export the coefficent list in a matlab function
% the coeffs are moltiplied by [1,2,6,24,120] to obtain the value of c1, c2, c3, c4, c5
% coef_list_var = [0,coeffs(sol_opt.s,t) .* [1,2,6,24,120]];
% coef_list_fun = ...
% - Use the matlabFunction function to generate a matlab function using a
%   symbolic function 

%% Export the total cost in a matlab function 
% total_cost_var = simplify(int(sol_opt.j^2,t,0,T));
% total_cost_fun = ...
% - Use the matlabFunction function to generate a matlab function using a
%   symbolic function

%% Inspect solutions 

Tmax  = 4.;
t_list = linspace(0,Tmax,100);

v0val = 10;
a0val = 1;
xfval = 90;
vfval = 25;
afval = 0.;

%% The functions work both for vector (t_list) or number (t)
% s_list = s_opt_fun(t_list, v0val, a0val, xfval, vfval, afval, Tmax);
% v_list = ...
% a_list = ...
% j_list = ...

%% The functions that used the coeflist
coeffs = coef_list_fun(v0val, a0val, xfval, vfval, afval, Tmax);
% coeffs_s_list = coeffs_s_opt_fun(t_list, coeffs);
% coeffs_v_list = ...
% coeffs_a_list = ...
% coeffs_j_list = ...

figure(1)
%% Position
% subplot(2,4,1)
% plot(t_list, s_list,'b');hold on;
% plot(t_list, coeffs_s_list,'r--');
% grid on
% xlabel('Time (s)','Interpreter','latex');
% ylabel('Position $(m)$','Interpreter','latex');
    
%% Velocity
% subplot(2,4,2)
% ...
    
%% Acceleration
% subplot(2,4,3)
% ...

%% Control
% subplot(2,4,4)
% ...
    
%% Velocity on position
% subplot(2,4,6)
% ...

%% Acceleration on position
% subplot(2,4,7)
% ...

%% Control on position
% subplot(2,4,8)
% ...