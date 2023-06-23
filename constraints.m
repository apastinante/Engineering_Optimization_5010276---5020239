function [c, ceq] = constraints(x, J, n, Td_prem, ...
    T_max, pointing_accuracy, settling_time)

%%                             Initial Conditions
%% 

kp = x(1);
kd = x(2);
ki = x(3);

tot_ang_mom = 0;
thrust_constraints = zeros(3, 1);
pointing_constraints = zeros(3, 1);

initial_devs = [30, 50, 80];

for j=1:3
  
dev = initial_devs(j);
r0 = deg2rad(dev); %initial roll
p0 = deg2rad(dev); %initial pitch
y0 = deg2rad(dev); %initial yaw

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

del_ang_mom = trapz(tspan, vecnorm(taus, 2, 2));
C_t_i = torque_constraint(taus, T_max);
C_acc_i = accuracy_constraint(ts, ys(:, 4:6), pointing_accuracy, settling_time);

tot_ang_mom = tot_ang_mom + del_ang_mom;
thrust_constraints(j) = C_t_i;
pointing_constraints(j) = C_acc_i;
end

C_t = max(thrust_constraints);
C_acc = max(pointing_constraints);

c = [C_t, C_acc];
ceq = [];

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

function [C_t] = torque_constraint(taus, T_max)
tau_req = max(taus, [], 1);
C_t = (max(tau_req) - T_max)/T_max;
end

function [C_acc] = accuracy_constraint(ts, theta, pointing_accuracy, settling_time)
bool = ts>settling_time;
aux = max(abs(theta(bool, :)), [], 2);
accuracy = max(aux);
C_acc = (accuracy - pointing_accuracy)/pointing_accuracy;
end

end