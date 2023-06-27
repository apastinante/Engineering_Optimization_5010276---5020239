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

T_max = 1;  % [Nm]
pointing_accuracy = deg2rad(2);  % [rad]
settling_time = 90;  % [s]

%% Optimisation with Homemade SLP

x0 = [0.2492    6.0766]; 
lb = [0.1, 4];
ub = [1.4, 14];

A = [];
b = [];
Aeq = [];
beq = [];


nonlcon = @(x)simp_constraints(x, J, n, Td_prem, T_max, pointing_accuracy, ...
    settling_time);

lcon = @(x)dgslp(x, J, n, Td_prem, ...
    T_max, pointing_accuracy, settling_time)

dfun = @(x)dfslp(x, J, n, Td_prem, ...
    T_max, pointing_accuracy, settling_time)

fun = @(x)simp_del_ang_mom(x, J, n, Td_prem, T_max, pointing_accuracy, ...
    settling_time);

%options = optimoptions('linprog','Algorithm','dual-simplex')

%x_opt = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon, options)



% Initialize variables


% Define convergence criteria
max_iterations = 300;   % Maximum number of iterations
tolerance = 1e-6;       % Convergence tolerance for the objective function

% Initialize iteration counter and objective function value
iteration = 0;
prev_obj_value = inf;
x = x0;

while iteration < max_iterations
    % Solve the linear programming subproblem
    % Use appropriate linear programming formulation based on your problem
    
    % Define objective coefficients for linear programming
    f = dfun(x);  
    
    % Define inequality constraints for linear programming
    A = lcon(x);  % Placeholder values, replace with actual constraint matrix
    b = [0 ; 0];  % Placeholder values, replace with actual constraint vector
    
    % Define equality constraints for linear programming
    Aeq = [];  % Placeholder values, replace with actual constraint matrix
    beq = [];  % Placeholder values, replace with actual constraint vector
    
    
    % Solve the linear programming problem
    [x, obj_value, exitflag, output] = linprog(f, A, b, Aeq, beq, lb, ub);
    
    % Check convergence criteria
    if abs(obj_value - prev_obj_value) < tolerance
        break;  % Convergence achieved, exit the loop
    end
    
    % Update variables
    Kp = x(1);
    Kd = x(2);
    %Ki = x(3);
    
    % Update iteration counter and previous objective function value
    iteration = iteration + 1;
    prev_obj_value = obj_value;
end

% Retrieve optimal solution
optimal_Kp = x(1);
optimal_Kd = x(2);
%optimal_Ki = Ki;

% Display results
fprintf('Optimal Solution:\n');
fprintf('Kp: %.4f\n', optimal_Kp);
fprintf('Kd: %.4f\n', optimal_Kd);
fprintf('Ki: %.4f\n', optimal_Ki);