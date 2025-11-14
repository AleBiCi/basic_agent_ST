%% IMPORT THE TEST DATA AS TABLES AND PLOT THEM
clc, clear opt_contr_data;
opt_contr_data = importfile_vel_acc_coefs("../log_internal/acc_vel_test.csv");


%% OPTIMAL CONTROL VELOCITY AND ACCELERATION
% Save acc and vel data on another array
acc_capped = zeros(height(opt_contr_data.acc));
req_acc_capped = zeros(height(opt_contr_data.req_acc));

vel_capped = zeros(height(opt_contr_data.vel));
req_vel_capped = zeros(height(opt_contr_data.req_vel));


% Cap the acceleration and velocity values to a maximum
for i=1:height(opt_contr_data.acc)
    acc_capped(i) = opt_contr_data.acc(i);
    req_acc_capped(i) = opt_contr_data.req_acc(i);

    vel_capped(i) = opt_contr_data.vel(i);
    req_vel_capped(i) = opt_contr_data.req_vel(i);
    
    if abs(acc_capped(i)) > 1e2
        if sign(acc_capped(i)) == -1
            acc_capped(i) = -1e2;
        else
            acc_capped(i) = 1e2;
        end
    end

    if abs(req_acc_capped(i)) > 2e2
        if sign(req_acc_capped(i)) == -1
            req_acc_capped(i) = -2e2;
        else
            req_acc_capped(i) = 2e2;
        end
    end

    if abs(vel_capped(i)) > 1e2
        if sign(vel_capped(i)) == -1
            vel_capped(i) = -1e2;
        else
            vel_capped(i) = 1e2;
        end
    end

    if abs(req_vel_capped(i)) > 2e2
        if sign(req_vel_capped(i)) == -1
            req_vel_capped(i) = -2e2;
        else
            req_vel_capped(i) = 2e2;
        end
    end
end

%% Optimal control Plots

% Acceleration plot
figure(1), clf, hold on, grid on, grid minor;
plot(opt_contr_data.time, acc_capped, 'r');
plot(opt_contr_data.time, req_acc_capped, 'b--');
legend('acc', 'req_acc', 'Location', 'best');

% Velocity plot
figure(2), clf, hold on, grid on, grid minor;
plot(opt_contr_data.time, vel_capped, 'r');
plot(opt_contr_data.time, req_vel_capped, 'b--');
legend('vel', 'req_vel', 'Location', 'best');
