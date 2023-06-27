function flin = flinx0(x, x0);
% Return linearized objective function 
% Input:
%   x  : ([1x2] row) design variables (Kp and Kd)
%   x0 : ([1x2] row) design point where objective is linearized.
% Output:
%   flin  : [1x1] scalar of linearized objective function value



% Analysis of valve spring in linearization point x0:
[del_ang_mom_x0] = simp_del_ang_mom(x0, J, n, Td_prem, ...
    T_max, pointing_accuracy, settling_time);
% Objective function in linearization point:
fx0 = del_ang_mom_x0;
 
% Gradient of objective function in linearization point:
dF = dfslp(x0, J, n, Td_prem, ...
    T_max, pointing_accuracy, settling_time);
 
% Linearized objective function:
flin = fx0 + dF*(x - x0)';
    
end 