clc
close all
clear

%% Variable Definition

J = [66.66 0 0; 0 66.66 0;
    0 0 66.66];  % spacecraft moments of inertia
mu = 398600;  % [km^3/s^2]
h = 700;  % [km] initial orbit height (wrt Earth's surface) 
Re = 6371;  % [km]
a = Re + h;  % [km] Semi-major axis of the orbit
n = sqrt(mu/a^3);  % [rad/s] angular rate of the spacecraft around the Earth
isp = 285;  % [s]

Td_prem = [1e-4; 1e-4; 1e-4]; %[N] preliminary simplified disturbance torque

%% Constraint values
T_max = 1;  % [Nm]
pointing_accuracy = deg2rad(2);  % [rad]
settling_time = 90;  % [rad]
prop_mass = 0.005;  % [kg]
energy_stored = (8.64e6 - 50*60*710)/10;  % [J]

%% Optimisation with Interior Point Algorithm

x0 = [1.4, 4, 1e-6];  % [0.4145, 12.7255, 4.9622e-04];
lb = [0.1, 4, 1e-6];
ub = [1.4, 14, 2e-3];

A = [];
b = [];
Aeq = [];
beq = [];

nonlcon = @(x)constraints(x, J, n, Td_prem, T_max, pointing_accuracy, ...
    settling_time);

fun = @(x)del_ang_mom(x, J, n, Td_prem);

options = optimoptions('fmincon','Display','iter');

%x_opt = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon, options)

%% Optimisation with Sequential Quadratic Programming

options = optimoptions('fmincon','Display','iter', 'Algorithm', 'sqp');

x_opt = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon, options)

%% Optimisation with Generalised Pattern Search

options = optimoptions('patternsearch','Display','iter');

%x_opt = patternsearch(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon, options)

%% Optimisation with Genetic Algorightm

nvars = 3;
options = optimoptions('ga', 'Display', 'iter', 'PopulationSize', 100);
%[x,fval,exitflag,output] = ga(fun, nvars, A, b, Aeq, beq, lb, ub, nonlcon, options)

