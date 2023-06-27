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

x_opt = [0.1776 4.00]
fun_opt = fun(x_opt)

fprintf('Opt_func: %.4f\n', fun_opt);
%options = optimoptions('linprog','Algorithm','dual-simplex')

%x_opt = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon, options)



% Initialize variables


% Define convergence criteria
max_iterations = 50;   % Maximum number of iterations
tolerance = 1e-6;       % Convergence tolerance for the objective function

% Initialize iteration counter and objective function value
iteration = 0;
prev_obj_value = inf;
x = x0;

% Make a list to store the objective function and the x values at each iteration
obj_values = [];
x_values = [];

% Make a list to store the constraint values at each iteration
constraint_values = [];


while iteration < max_iterations
    % Solve the linear programming subproblem
    % Use appropriate linear programming formulation based on your problem
    
    % Define objective coefficients for linear programming
    f = dfun(x);  
    
    % Define inequality constraints for linear programming
    A = lcon(x);  
    b = -nonlcon(x);  
    
    % Define equality constraints for linear programming
    Aeq = [];  
    beq = [];  
    
    
    % Solve the linear programming problem
    [x, obj_value, exitflag, output] = linprog(f, A, b, Aeq, beq, lb, ub);
    
    % Check convergence criteria
    if abs(obj_value - prev_obj_value) < tolerance
        break;  % Convergence achieved, exit the loop
    end
    
    % Save x and objective function value for plotting
    x_values = [x_values , x];
    obj_values = [obj_values , obj_value];
    
    % Save constraint values for plotting
    constraint_values = [constraint_values; nonlcon(x)];

    % Update iteration counter and previous objective function value
    iteration = iteration + 1;
    prev_obj_value = obj_value;
end

% Retrieve optimal solution
optimal_Kp = x(1);
optimal_Kd = x(2);

% Display results
fprintf('Optimal Solution:\n');
fprintf('Kp: %.4f\n', optimal_Kp);
fprintf('Kd: %.4f\n', optimal_Kd);

% Make a meshgrid with the lb and ub values for each variable
[xq,yq] = meshgrid(lb(1):0.01:ub(1), lb(2):0.01:ub(2));

% Plot the constraint values per constraint (2D array) as a function of the x (2D array) values, 
% so it has to be a 2D colormap for each constraint
c_tau = constraint_values(:, 1);
c_acc = constraint_values(:, 2);
c_tau_surf = griddata(x_values(1, :), x_values(2, :), c_tau, xq, yq);
c_acc_surf = griddata(x_values(1, :), x_values(2, :), c_acc, xq, yq);

figure(1);
% Color map for constraint 1 at each combination of x(1) and x(2) value
plot3(xq, yq, c_tau_surf);
xlabel('Kp');
ylabel('Kd');
zlabel('Constraint 1');
title('Constraint 1 vs. Kp and Kd');
colorbar;
figure(2);
% Color map for constraint 2 at each combination of x(1) and x(2) value
plot3(xq, yq, c_acc_surf);
xlabel('Kp');
ylabel('Kd');
zlabel('Constraint 2');
title('Constraint 2 vs. Kp and Kd');
colorbar;

% Show the objective function value per iteration
figure(3);
plot(obj_values);
xlabel('Iteration');
ylabel('Objective Function Value');
title('Objective Function Value vs. Iteration');






