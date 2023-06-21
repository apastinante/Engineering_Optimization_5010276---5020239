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
bounds = [0.1, 2; 1, 15; 1e-5, 5e-3];

%% Simulation
n_sim = 50;
pool = ones(n_sim, 3);
assessment = ones(n_sim, 6);

for i = 1:n_sim
    kp = (bounds(1, 2) - bounds(1, 1))*rand() + bounds(1, 1);
    kd = (bounds(2, 2) - bounds(2, 1))*rand() + bounds(2, 1);
    ki = (bounds(3, 2) - bounds(3, 1))*rand() + bounds(3, 1);
    pool(i, :) = [kp, kd, ki];
end

for i = 1:n_sim
    display(i)
    tic
    assessment(i, :) = full_objective_function(pool(i, 1), pool(i, 2), J, n, Td_prem, T_max, pointing_accuracy, ...
    settling_time, pool(i, 3), isp, prop_mass, energy_stored);
    toc
end

%% Plotting

assessment

bool = assessment(:, end) > 0;
bool2 = assessment(:, end) == 0;

figure(1)
scatter(pool(bool, 1), assessment(bool, 1), 'Filled', 'b')
xlabel('$K_p$ [-]', 'interpreter', 'latex')
ylabel('Delivered Angular Momentum [Nms]')
hold on
scatter(pool(bool2, 1), assessment(bool2, 1), 'Filled', 'r')
figure(2)
scatter(pool(bool, 2), assessment(bool, 1), 'Filled', 'b')
xlabel('$K_d$ [-]', 'interpreter', 'latex')
ylabel('Delivered Angular Momentum [Nms]')
hold on
scatter(pool(bool2, 2), assessment(bool2, 1), 'Filled', 'r')
figure(3)
scatter(pool(bool, 3), assessment(bool, 1), 'Filled', 'b')
hold on
scatter(pool(bool2, 3), assessment(bool2, 1), 'Filled', 'r')
xlabel('$K_i$ [-]', 'interpreter', 'latex')
ylabel('Delivered Angular Momentum [Nms]')


x = 1:4;
y = zeros(4,1);
for i=1:length(assessment(:, 1))
    if assessment(i, 2) > 0
        y(1) = y(1) + 1;
    elseif assessment(i, 3) > 0
        y(2) = y(2) + 1;
    elseif assessment(i, 4) > 0
        y(3) = y(3) + 1;
    elseif assessment(i, 5) > 0
        y(4) = y(4) + 1;
    end
end

y = y./sum(y);
figure(4)
bar(x, y)
xlabel('Contraint Number [-]')
ylabel('Percentage of violations [-]')





