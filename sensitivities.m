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

%% Design Variable Bounds
bounds = [0.3, 1.5; 4, 20; 1e-5, 5e-3];

%% Simulation

h_array = [1e-12, 5e-12, 1e-11, 5e-11, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1e-1, 5e-1];
kp = 0.6;
kd = 10;
ki = 1e-4;
assessment = ones(length(h_array), 3);


for i = 1:length(h_array)
    display(i)
    tic
    assessment(i, :) = central_difference(h_array(i), kp, kd, J, n, Td_prem, T_max, pointing_accuracy, ...
    settling_time, ki);
    toc
end

%% Plotting

figure(1)
semilogx(h_array, assessment(:, 1))
xlabel('h Step [-]')
ylabel('$\partial y/\partial k_p [Nms]$', 'interpreter', 'latex', 'FontSize', 14)
grid on
figure(2)
semilogx(h_array, assessment(:, 2))
xlabel('h Step [-]')
ylabel('$\partial y/ \partial k_d [Nms]$', 'interpreter', 'latex', 'FontSize', 14)
grid on
figure(3)
semilogx(h_array(1:9), assessment(1:9, 3))
xlabel('h Step [-]')
ylabel('$\partial y/ \partial k_i [Nms]$', 'interpreter', 'latex', 'FontSize', 14)
grid on


function [grad] = central_difference(h, kp, kd, J, n, Td_prem, ...
    T_max, pointing_accuracy, settling_time, ki)
top = final_objective_function(kp+h, kd, J, n, Td_prem, T_max, ...
    pointing_accuracy, settling_time, ki);
bottom = final_objective_function( ...
    kp-h, kd, J, n, Td_prem, T_max, pointing_accuracy, settling_time, ...
    ki);
dkp = (top(1) - bottom(1))/2/h;
top = final_objective_function(kp, kd+h, J, n, Td_prem, T_max, ...
    pointing_accuracy, settling_time, ki);
bottom = final_objective_function( ...
    kp, kd-h, J, n, Td_prem, T_max, pointing_accuracy, settling_time, ...
    ki);
dkd = (top(1) - bottom(1))/2/h;
if h < 2e-2
    top = final_objective_function(kp, kd, J, n, Td_prem, T_max, ...
        pointing_accuracy, settling_time, ki+h);
    bottom = final_objective_function( ...
        kp, kd, J, n, Td_prem, T_max, pointing_accuracy, settling_time, ...
        ki-h);
    dki = (top(1) - bottom(1))/2/h;
else
    dki = 0;
end
grad = [dkp, dkd, dki];
end



