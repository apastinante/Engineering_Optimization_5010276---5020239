clc
clear
close all

%%                               Variable Definition
%% spacecraft related
J = [66.66 0 0; 0 66.66 0;
    0 0 66.66]; %spacecraft moments of inertia
                               
%% orbit related 
mu = 398600; %[km^3/s^2]
h = 700; %[km] initial orbit height (wrt Earth's surface) 
Re = 6371; %[km]
a = Re + h; %[km] Semi-major axis of the orbit
n = sqrt(mu/a^3);%[rad/s] angular rate of the spacecraft around the Earth
T = 2*pi*sqrt(a^3/mu); %[s] orbital period
v = sqrt(mu/a)*1000; %[m/s] orbital velocity

%% disturbance torque related
Td_prem = [1e-4; 1e-4; 1e-4]; %[N] preliminary simplified disturbance torque
u_v = [1; 0; 0]; %unit vector in velocity direction

%% simulation/ feedback related
% nominal -> del_ang_mom = 13.7373 Nms
% kp = 0.6;
% kd = 10;
% ki = 1e-4;

% optimal simplified problem -> del_ang_mom = 8.843966 Nms
% kp = 0.2492;
% kd = 6.0768;
% ki = 0;

% fitness_simple = objective_function(kp, kd, J, n, Td_prem, ...
%     T_max, pointing_accuracy, settling_time)

% optimal full problem -> del_ang_mom = 16.76 Nms (full) 
%                     and del_ang_mom = 9.9847 Nms (simplified)
kp = 0.2669; %0.349384116049841;
kd = 7.4259; %8.550631994521687;
ki = 1.877867588671635e-04; %1.026123024369114e-04;

% optimal simplified problem ga -> del_ang_mom = 9.9853 NMS
kp = 0.3977; % 0.3593
kd = 10.1746; % 9.1939
ki = 0;

% optimal full problem ga -> del_ang_mom = 18.2793 Nms (full)
%                              and del_ang_mom = (simplified)
kp = 0.4145;
kd = 12.7255;
ki = 4.9622e-04;

%% Constraint values
T_max = 1;  % [Nm]
pointing_accuracy = deg2rad(2);  % [rad]
settling_time = 90;  % [s]


%%                             Initial Conditions
%% 
% u0 = deg2rad(0.0209); %[rad] initial argument of latitude
% RAAN0 = deg2rad(0); %[rad] initial right ascension of ascending node

r0 = deg2rad(80); %initial roll
p0 = deg2rad(80); %initial pitch
y0 = deg2rad(80); %initial yaw

q0 = [r0 p0 y0]; %initial quaternions

w0 = [0 0 0]; %initial angular rates
y0 = [w0 q0];
sampling_time = 0.3;
sim_length = 120;
tspan = 0:sampling_time:sim_length;
opts = odeset('RelTol', 1e-3, 'AbsTol', 1e-3);

theta_ref = [0, 0, 0];

tau = 0;
taus = zeros(length(tspan), 3);
ts = zeros(1, length(tspan));
ys = zeros(length(tspan), 6);
ys(1, :) = y0;

%%                                Main Program
%% 

% fitness_simple = objective_function(kp, kd, J, n, Td_prem, ...
%     T_max, pointing_accuracy, settling_time)

fitness_full = final_objective_function(kp, kd, J, n, Td_prem, T_max, ...
    pointing_accuracy, settling_time, ki)

commands_ref = 0*ones(length(tspan), 1)';
tic

for i = 2:length(tspan)
    [t, y] = ode23(@(t, y) odefunc(t, y, J, n, Td_prem, tau), ...
        [tspan(i-1) tspan(i)], y0, opts);
    y0 = y(end, :);
    ts(i) = t(end);
    ys(i, :) = y(end, :);
    w = y0(1:3)';
    theta = y0(4:6)';
    theta_dot = euler_dot(theta, w, n);
    bool = tspan <= t(end);
    theta_error = ys(bool, 4:6) - theta_ref;
    t_pid = ts(bool);
    tau = pid_controller(theta, kp, kd, ki, theta_dot, theta_ref, ...
        t_pid, theta_error);
    taus(i, :) = tau;
end
toc

del_ang_mom = trapz(tspan, vecnorm(taus, 2, 2))

ws = ys(:, 1:3)'; %obtain the angular rates vector from y at every instance of time
w_arr_x = ws(1, :);
w_arr_y = ws(2, :);
w_arr_z = ws(3, :);
qs = ys(:, 4:6)'; %obtain the quaternion vector from y at every instance of time
roll_arr = rad2deg(qs(1, :));
pitch_arr = rad2deg(qs(2, :));
yaw_arr = rad2deg(qs(3, :));

%%                                 Plotting
%% 
t_orb = tspan;

figure(1)
%
plot(t_orb, roll_arr, t_orb, pitch_arr, t_orb, yaw_arr)
hold on
plot(tspan, commands_ref, '--')
hold on
yline(2, 'k')
hold on
yline(-2, 'k')
hold on
xline(90, 'k')
xlabel('Time [s]');
ylabel('Angle [deg]');
legend('Roll', 'Pitch', 'Yaw', 'Reference');
grid minor 
%}

%%
figure(2)
%angular rates
%wx
subplot(3,1,1)
plot(t_orb, w_arr_x)
xlabel('Time [s]')
ylabel('wx [rad/s]')
%ylim([-0.006, 0.006])
grid minor
%wy
subplot(3,1,2)
plot(t_orb, w_arr_y)
xlabel('Time [s]')
ylabel('wy [rad/s]')
%ylim([-0.006, 0.006])
grid minor 
%wz
subplot(3,1,3)
plot(t_orb, w_arr_z)
xlabel('Time [s]')
ylabel('wz [rad/s]')
%ylim([-0.006, 0.006])
grid minor 

length(t_orb)
length(tau(:,1))

% Torques
figure(3)
subplot(3,1,1)
area(t_orb, taus(:, 1), 'EdgeColor', 'b', 'FaceColor', 'b');
xlabel('Time [s]');
ylabel('Tx');
grid minor
subplot(3,1,2)
area(t_orb, taus(:, 2), 'EdgeColor', 'b', 'FaceColor', 'b');
xlabel('Time [s]');
ylabel('Ty');
grid minor
subplot(3,1,3)
area(t_orb, taus(:, 3), 'EdgeColor', 'b', 'FaceColor', 'b');
xlabel('Time [s]');
ylabel('Tz');
grid minor
%}

%%                            Function Definition
%% Mathematical operations

function [C] = C1(angle)
%Computes a right handed euler rotation matrix around the x axis
%angle must be given in radians
C = [1 0 0; 0 cos(angle) sin(angle); 0 -sin(angle) cos(angle)];
end

function [C] = C2(angle)
%Computes a right handed euler rotation matrix around the y axis
%angle must be given in radians
C = [cos(angle) 0 -sin(angle); 0 1 0; sin(angle) 0 cos(angle)];
end

function [C] = C3(angle)
%Computes a right handed euler rotation matrix around the z axis
%angle must be given in radians
C = [cos(angle) sin(angle) 0; -sin(angle) cos(angle) 0; 0 0 1];
end

function [p_tilda] = skew_sym(p)
%computes the skew_symmetric matrix of a 3x1 vector (4th quaternion must be
%removed)
p_tilda = [0 -p(3) p(2); p(3) 0 -p(1); -p(2) p(1) 0];
end

%% Spacecraft EOM related
function [C] = rot_mat(theta)
%for a given set of quaternions, computes the rotation matrix between the
%body reference frame and the LVLH reference frame. 
C = C1(theta(1))*C2(theta(2))*C3(theta(3));
end 

function [wr] = omega_r(theta, w, n)
%computes the angular rate vector with respect to the geocentric inertial 
%frame which is required for the kinematic differential equation in terms
%of quaternions
C = rot_mat(theta);
cy = C(:, 2);
wr = w + n*cy;
end

function [theta_dot] = euler_dot(theta, w, n)
%provides the derivative with respect to time of a given vector of
%quaternions, given also the angular rate of the s/c and the angular rate
%vector
wr = omega_r(theta, w, n);
A = 1/cos(theta(2))*[cos(theta(2)) sin(theta(1))*sin(theta(2)) cos(theta(1))*sin(theta(2));
                     0 cos(theta(1))*cos(theta(2)) -sin(theta(1))*cos(theta(2));
                     0 sin(theta(1)) cos(theta(1))];
theta_dot = A*wr;
end

function [w_dot] = omega_dot(w, J, tau, Td)
%uses the dynamic eqn to compute the angular acceleration vector given as
%inputs, respectively: angular rate, principal moments of inertia matrix,
%momentum wheel moment ofinertia matrix, momentum wheekl angular velocity,
%input torque, and disturbance torque. 
w_dot = J\(cross(-w, J*w) + tau + Td);
end

%% Disturbance related
function [T_gg] = grav_grad(theta, n, J)
%Computes the gravity gradient disturbance torque 
C = rot_mat(theta);
c_z = C(:, 3);
a = J*c_z;
%T_gg = cross(3*n^2*c_z, a);
T_gg = 3*n^2*skew_sym(c_z)*a;
end

%% Feedback Laws

function [tau] = pd_controller(theta, kp, kd, theta_dot, theta_ref)
tau = -(kp*(theta - theta_ref') + kd*(theta_dot));
end

function [tau] = pid_controller(theta, kp, kd, ki, theta_dot, theta_ref, ...
    t_pid, theta_error)
error_r = trapz(t_pid, theta_error(:, 1));
error_p = trapz(t_pid, theta_error(:, 2));
error_y = trapz(t_pid, theta_error(:, 3));
error_integral = [error_r, error_p, error_y]';
tau = -(kp*(theta - theta_ref') + kd*(theta_dot) + ki*error_integral);
end

%% ODE45
function [dydt] = odefunc(t, y, J, n, Td_prem, tau)
%displays both differential equations (angular rates and quaternions)  in 
%an array as a function of time

%w -> y(1)-y(3)
%q -> y(4)-y(7)
w = y(1:3);
theta = y(4:6);

T_gg = grav_grad(theta, n, J);
Td_prem = Td_prem + T_gg;

%System of differential eqns
dydt = zeros(6, 1);

%differential eqns in vector form
w_dot = omega_dot(w, J, tau, Td_prem);
theta_dot = euler_dot(theta, w, n);

%differntial eqns in components
dydt(1) = w_dot(1);
dydt(2) = w_dot(2);
dydt(3) = w_dot(3);

dydt(4) = theta_dot(1);
dydt(5) = theta_dot(2);
dydt(6) = theta_dot(3);
end