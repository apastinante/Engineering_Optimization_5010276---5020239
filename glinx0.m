function glin = glinx0(x, x0);
% Return linearized objective function 
% Input:
%   x  : ([1x2] row) design variables (Kp and Kd)
%   x0 : ([1x2] row) design point where objective is linearized.
% Output:
%   flin  : [1x1] scalar of linearized objective function value



% Analysis of valve spring in linearization point x0:
cx0 = simp_constraints(x0, , J, n, Td_prem, ...
    T_max, pointing_accuracy, settling_time);
% Objective function in linearization point:
gx0 = cx0;
 
% Gradient of objective function in linearization point:
dG = dgslp(x0, J, n, Td_prem, ...
    T_max, pointing_accuracy, settling_time);
 
% Linearized objective function:
glin = gx0 + dG*(x - x0)';
    
end 