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
ode1 = diff(s) == v;
ode2 = diff(v) == a;
ode3 = diff(a) == u;
% - Use "diff()" function to derive the function (e.g. diff(s) == v)

%% Define the Lagrangian and the Hamiltonian.
% Lagrangian
L = u^2;

% Hamiltonian
H = L + l1*rhs(ode1) + l2*rhs(ode2)+ l3*rhs(ode3);
% - Use "rhs()" function for getting the right hand side 

%% Solve the optimal control problem (Solving the Hamiltonian). 
% Derivative of Hamiltonian w.r.t. the control
du = functionalDerivative(H,u);
syms opt_u;
% Optimal solution
opt_u = solve(subs(du,u(t),opt_u) == 0,opt_u);  % use subs because trying to solve using function of time results in warning and void solution
% - use "functionalDerivative()" to derive a function with respect to another function (e.g. du = functionalDerivative(H,u)) 
% - Use "solve" and "subs()" to replace variable inside the equation

%% Write the second optimality condition.
% Second Optimality Condition
Dl1 = diff(l1,t) == -functionalDerivative(H,s);  % can omit t also,  since it's the only variable in l1,l2,l3
Dl2 = diff(l2,t) == -functionalDerivative(H,v);
Dl3 = diff(l3,t) == -functionalDerivative(H,a);
% - Use "diff(,t)" to derive w.r.t. the time
%   and "functionalDerivative()" to detive w.r.t. a variable (for the H)

%% Substitute the optimal solution opt_u to state equations
% New system equation
ode3s = diff(a) == subs(rhs(ode3), u, opt_u);  % inside the ode3 replace u with opt_u
% - Replace the ode3 equation with the solution of the optimal control problem

%% Define the Boundary condition on initial state and Boundary condition on final state
% Boundary condition on initial and final states, written as strings
ICs = 's(0)=0, v(0)=v0, a(0)=a0';  % automatically creates a new symoblic variable (v0, a0) for each starting condition
FCs = 's(T)=sf, v(T)=vf, a(T)=af';
% - Write the condition as a string

%% Find the solution of the OCP imposing the boundary condition
% Solution of the optimal control problem
sol_opt = dsolve([ode1, ode2, ode3s, Dl1, Dl2, Dl3], ICs, FCs);
% - Use the function "dsolve([], , )" to obtain a solution of the optimal
%   control problem

disp('Optimal polynomial longitudinal position:');
pretty(sol_opt.s)

disp('Optimal polynomial velocity:');
pretty(sol_opt.v)

disp('Optimal polynomial acceleration:');
pretty(sol_opt.a)

% Assign to functions the solutions found

%% Get the optimal control solution
% Obtain optimal control solution
sol_opt.j = subs(opt_u, l3(t), sol_opt.l3);
% - Use the "subs" function on the opt_u with the value of l3

%% Export the solution in a matlab function
syms t v0 a0 sf vf af T  
% Create the matlab function for the solutions
% Can also call simplify on sol_opt.* to keep them shorter before writing
% in a file
s_opt_fun = matlabFunction(sol_opt.s,'Vars',[t,v0,a0,sf,vf,af,T],'File','s_opt.m');
v_opt_fun = matlabFunction(sol_opt.v,'Vars',[t,v0,a0,sf,vf,af,T],'File','v_opt.m');
a_opt_fun = matlabFunction(sol_opt.a,'Vars',[t,v0,a0,sf,vf,af,T],'File','a_opt.m');
j_opt_fun = matlabFunction(sol_opt.j,'Vars',[t,v0,a0,sf,vf,af,T],'File','j_opt.m');
% - Use matlabFunction to generate a matlab function using a
%   symbolic function

%% Create the matlab function for the solution using the coeffs
c = sym('c',[1 6]);  % called 'm' in the slides
% idea: write optimal solution equations with explicit coefficients
% c1,...,c5
% Note: c0 = 0

coeffs_s_opt = c(1)*t + 1/2*c(2)*t^2 + 1/6*c(3)*t^3 + 1/24*c(4)*t^4 + 1/120*c(5)*t^5 + c(6)*t^5;
% Get further equations by differentiating s
coeffs_v_opt = diff(coeffs_s_opt, t);
coeffs_a_opt = diff(coeffs_v_opt, t);
coeffs_j_opt = diff(coeffs_a_opt, t);

coeffs_s_opt_fun = matlabFunction(coeffs_s_opt,'Vars',{t,c},'File','coeffs_s_opt.m');
coeffs_v_opt_fun = matlabFunction(coeffs_v_opt,'Vars',{t,c},'File','coeffs_v_opt.m');
coeffs_a_opt_fun = matlabFunction(coeffs_a_opt,'Vars',{t,c},'File','coeffs_a_opt.m');
coeffs_j_opt_fun = matlabFunction(coeffs_j_opt,'Vars',{t,c},'File','coeffs_j_opt.m');

%% Export the coefficent list in a matlab function
% the coeffs are moltiplied by [1,2,6,24,120] to obtain the value of c1, c2, c3, c4, c5
coef_list_var = [0,coeffs(sol_opt.s,t) .* [1,2,6,24,120]]; % also includes c0=0 at the beginning
coef_list_fun = matlabFunction(coef_list_var,'Vars',[v0,a0,sf,vf,af,T],'File','coef_list_fun.m');  % !! doesn't depend on t, only on T
% - Use the matlabFunction function to generate a matlab function using a
%   symbolic function 

%% Export the total cost in a matlab function 
total_cost_var = simplify(int(sol_opt.j^2,t,0,T));
total_cost_fun = matlabFunction(total_cost_var,'Vars',[v0,a0,sf,vf,af,T],'File','total_cost.m');
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
s_list = s_opt_fun(t_list, v0val, a0val, xfval, vfval, afval, Tmax);
v_list = v_opt_fun(t_list, v0val, a0val, xfval, vfval, afval, Tmax);
a_list = a_opt_fun(t_list, v0val, a0val, xfval, vfval, afval, Tmax);
j_list = j_opt_fun(t_list, v0val, a0val, xfval, vfval, afval, Tmax);

%% The functions that used the coeflist
coeffs = coef_list_fun(v0val, a0val, xfval, vfval, afval, Tmax);
coeffs_s_list = coeffs_s_opt_fun(t_list, coeffs);
coeffs_v_list = coeffs_v_opt_fun(t_list, coeffs);
coeffs_a_list = coeffs_v_opt_fun(t_list, coeffs);
coeffs_j_list = coeffs_a_opt_fun(t_list, coeffs);

% figure(1)
% %% Position
% subplot(2,4,1)
% plot(t_list, s_list,'b');hold on;
% plot(t_list, coeffs_s_list,'r--');
% grid on
% xlabel('Time (s)','Interpreter','latex');
% ylabel('Position $(m)$','Interpreter','latex');
% 
% %% Velocity
% subplot(2,4,2)
% plot(t_list, v_list,'b','LineWidth',1.5); hold on;
% plot(t_list, coeffs_v_list,'r--','LineWidth',1.2);
% grid on
% xlabel('Time (s)','Interpreter','latex');
% ylabel('Velocity $(m/s)$','Interpreter','latex');
% 
% %% Acceleration
% subplot(2,4,3)
% plot(t_list, a_list,'b','LineWidth',1.5); hold on;
% plot(t_list, coeffs_a_list,'r--','LineWidth',1.2);
% grid on
% xlabel('Time (s)','Interpreter','latex');
% ylabel('Acceleration $(m/s^2)$','Interpreter','latex');
% 
% %% Control (jerk)
% subplot(2,4,4)
% plot(t_list, j_list,'b','LineWidth',1.5); hold on;
% plot(t_list, coeffs_j_list,'r--','LineWidth',1.2);
% grid on
% xlabel('Time (s)','Interpreter','latex');
% ylabel('Jerk $(m/s^3)$','Interpreter','latex');
% 
% %% Velocity vs. Position
% subplot(2,4,6)
% plot(s_list, v_list,'b','LineWidth',1.5);
% grid on
% xlabel('Position $(m)$','Interpreter','latex');
% ylabel('Velocity $(m/s)$','Interpreter','latex');
% 
% %% Acceleration vs. Position
% subplot(2,4,7)
% plot(s_list, a_list,'b','LineWidth',1.5);
% grid on
% xlabel('Position $(m)$','Interpreter','latex');
% ylabel('Acceleration $(m/s^2)$','Interpreter','latex');
% 
% %% Control vs. Position
% subplot(2,4,8)
% plot(s_list, j_list,'b','LineWidth',1.5);
% grid on
% xlabel('Position $(m)$','Interpreter','latex');
% ylabel('Jerk $(m/s^3)$','Interpreter','latex');