clc
close all
clear

%% Variable Definition

J = [66.66 0 0; 0 66.66 0;
    0 0 66.66]; %spacecraft moments of inertia
mu = 398600; %[km^3/s^2]
h = 700; %[km] initial orbit height (wrt Earth's surface) 
Re = 6371; %[km]
a = Re + h; %[km] Semi-major axis of the orbit
n = sqrt(mu/a^3);%[rad/s] angular rate of the spacecraft around the Earth

Td_prem = [1e-4; 1e-4; 1e-4]; %[N] preliminary simplified disturbance torque

%% Constraint values
T_max = 1;  % [Nm]
pointing_accuracy = deg2rad(2);  % [rad]
settling_time = 90;  % [s]

%% Optimisation with Sequential Quadratic Programming

x0 = [0.6, 10];
lb = [0.1, 4];
ub = [1.4, 14];

A = [];
b = [];
Aeq = [];
beq = [];

nonlcon = @(x)simp_constraints(x, J, n, Td_prem, T_max, pointing_accuracy, ...
    settling_time);

fun = @(x)simp_del_ang_mom(x, J, n, Td_prem, T_max, pointing_accuracy, ...
    settling_time);

options = optimoptions('fmincon','Display','iter');

% x_opt = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon, options)

%% Optimisation with Genetic Algorightm

nvars = 2;
options = optimoptions('ga','Display','iter', 'PopulationSize', 50);
[x,fval,exitflag,output] = ga(fun,nvars,A,b,Aeq,beq,lb,ub,nonlcon, options)
