function dG = dgslp(x, J, n, Td_prem, ...
    T_max, pointing_accuracy, settling_time)

% Input:
% x  : gains "[Kp, Kd]" for which derivatives are computed.

% Output:
% dG  : [5X2] matrix with gradients of 5 constraints:
%        "[dg1dx1  dg1dx2
%          ....     ....
%          dg5dx1  dg5dx2]" 

% Note: Constant parameter values are read within the function  springcon3.

% Forward finite diffence gradients of objective function and constraints

% Finite diffence step
hx = 1.0e-6;

% Constraint gradients 
gx = simp_constraints(x, J, n, Td_prem, ...
    T_max, pointing_accuracy, settling_time);  % Note: gx is [1X2] row vector
gx1plush = simp_constraints([x(1)+hx, x(2)], J, n, Td_prem, ...
    T_max, pointing_accuracy, settling_time);
gx2plush = simp_constraints([x(1), x(2)+hx], J, n, Td_prem, ...
    T_max, pointing_accuracy, settling_time);
dgdx1 = (gx1plush - gx)./hx;
dgdx2 = (gx2plush - gx)./hx;
dG = [dgdx1' dgdx2'];

end 