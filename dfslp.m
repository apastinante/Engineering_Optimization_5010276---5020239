function dF = dfslp(x, J, n, Td_prem, ...
    T_max, pointing_accuracy, settling_time)
% Derivatives of objective function of valve spring design  - Exercise 1

% Input:
% x  : design point "[D, d]" for which derivatives are computed.

% Output:
% dF  : [1x2] gradient (row)vector "[dfdx1 dfdx2]" of objective function.

% Note: Constant parameter values are read within the function springobj2. 
% Forward finite diffence gradients of objective function and constraints

% Finite diffence step
hx = 1.0e-6;

% Gradient of objective function
fx = simp_del_ang_mom(x, J, n, Td_prem, ...
    T_max, pointing_accuracy, settling_time);
fx1plush = simp_del_ang_mom([x(1)+hx, x(2)],J, n, Td_prem, ...
    T_max, pointing_accuracy, settling_time);
fx2plush = simp_del_ang_mom([x(1), x(2)+hx], J, n, Td_prem, ...
    T_max, pointing_accuracy, settling_time);
fx1minush = simp_del_ang_mom([x(1)-hx, x(2)],J, n, Td_prem, ...
    T_max, pointing_accuracy, settling_time);
fx2minush = simp_del_ang_mom([x(1), x(2)-hx], J, n, Td_prem, ...
    T_max, pointing_accuracy, settling_time);

dfdx1 = (fx1plush - fx )/(hx);
dfdx2 = (fx2plush - fx )/(hx);

% dfdx1 = (fx1plush - 2*fx + fx1minush)/(2*hx);
% dfdx2 = (fx2plush - 2*fx + fx2minush)/(2*hx);
dF = [dfdx1 dfdx2];

end 