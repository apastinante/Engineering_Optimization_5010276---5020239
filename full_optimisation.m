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
pointing_accuracy = deg2rad(1);  % [rad]
settling_time = 90;  % [rad]
prop_mass = 0.005;  % [kg]
energy_stored = (8.64e6 - 50*60*710)/10;  % [J]

%% Optimisation with Sequential Quadratic Programming

x0 = [0.6, 10, 1e-4];
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

x_opt = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon, options)
