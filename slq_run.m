
lb = [0.26 ; 7];
ub = [0.32 ; 8];

% Iterate over the design space 
err_dom = [];

corr_dom = [];
obj_values = [];
iterations = [],
opt_answers = [];
for kp = lb(1):0.02:ub(1)
    for kd = lb(2):0.1:ub(2)
        
        % Try/except to simulate the SLP optimization
        try 
            answers = sqlp_opt_simp(kp,kd);
            if answers(1) == -2
                disp('No point satisfies constraints at kp = ' + string(kp) + ' kd = ' + string(kd))
                err_dom = [err_dom ; [kp kd]];
            elseif answers(2) == 0
                disp('No point satisfies constraints at kp = ' + string(kp) + ' kd = ' + string(kd))
                err_dom = [err_dom ; [kp kd]];
            else 
                % If the optimization is successful, add the point to the correct domain
                corr_dom = [corr_dom ; [kp kd]];

                % Add the objective value to the list of objective values
                obj_values = [obj_values ; answers(6)];

                % Add the number of iterations to the list of iterations
                iterations = [iterations ; answers(3)];

                % Add the optimal answers to the list of optimal answers
                opt_answers = [opt_answers ; answers(4:5)];
            end 
        catch
            disp('Error at kp = ' + string(kp) + ' kd = ' + string(kd))
            err_dom = [err_dom ; [kp kd]];
        end


    end
end

figure(1)
% Plot the initial points at which the optimization failed
scatter(err_dom(:,1),err_dom(:,2),'r')
hold on

% Plot the points at which the optimization was successful
scatter(corr_dom(:,1),corr_dom(:,2),'b')
xlabel('Kp')
ylabel('Kd')
title('Optimization Successfull Domains')
legend('Optimization Failed','Optimization Successful')

figure(2)
% Plot the objective values
scatter3(corr_dom(:,1),corr_dom(:,2),obj_values,'b')
xlabel('Kp')
ylabel('Kd')
zlabel('Objective Value')
title('Objective Values')

figure(3)
% Plot the number of iterations
scatter3(corr_dom(:,1),corr_dom(:,2),iterations,'b')
xlabel('Kp')
ylabel('Kd')
zlabel('Number of Iterations')
title('Number of Iterations')

figure(4)
% Plot the optimal answers and project them onto the Kp-optimal Kd plane and the Kd-optimal Kd plane
scatter3(corr_dom(:,1),corr_dom(:,2),opt_answers(:,1),'b')
hold on
scatter3(corr_dom(:,1),ub(2)*ones(size(corr_dom(:,2))),opt_answers(:,1),'k')
scatter3(ub(1)*ones(size(corr_dom(:,1))),corr_dom(:,2),opt_answers(:,1),'k')
xlabel('Kp')
ylabel('Kd')
xlim([lb(1) ub(1)])
ylim([lb(2) ub(2)])
zlabel('Optimal Kd')
title('Optimal Kd')
hold off

figure(5)
% Plot the optimal answers and project them onto the Kp-optimal Kp plane and the Kd-optimal Kp plane
scatter3(corr_dom(:,1),corr_dom(:,2),opt_answers(:,2),'b')
hold on
scatter3(corr_dom(:,1),ub(2)*ones(size(corr_dom(:,2))),opt_answers(:,2),'k')
scatter3(ub(1)*ones(size(corr_dom(:,1))),corr_dom(:,2),opt_answers(:,2),'k')
xlim([lb(1) ub(1)])
ylim([lb(2) ub(2)])
xlabel('Kp')
ylabel('Kd')
zlabel('Optimal Kp')
title('Optimal Kp')

