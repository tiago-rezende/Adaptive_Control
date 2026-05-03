%%
% Experiment 9 –
clear
close all

%% Common parameters
a = 1;
b = 5;
s = tf('s');
l_star = 3.3/b;
k_star = (6.7 - a)/b;

model_tf = 3.3/(s+6.7);


%% Case 1: Constant reference
time = linspace(0,50,1000);
reference = 15 * ones(size(time));

[model_output, real_output, u_opt] = simulate_system(time, reference, a, b, k_star, l_star, model_tf);

[output_d, error_d, k_d, l_d, u_d] = mrac_direct(time, reference, model_output, a, b, 0.001, 0.002);
[output_i, error_i, a_i, b_i, u_i] = mrac_indirect(time, reference, model_output, a, b, 0.05, 0.02);

plot_results(time, model_output, real_output, output_d, output_i, error_d, error_i, u_opt, u_d, u_i, k_star, k_d, a_i, l_star, l_d, b_i, 'Case 1');

%% Case 2: Sinusoidal reference
time = linspace(0,10,1000);
reference = 2*sin(10*time) + 5*sin(3*time);

[model_output, real_output, u_opt] = simulate_system(time, reference, a, b, k_star, l_star, model_tf);

[output_d, error_d, k_d, l_d, u_d] = mrac_direct(time, reference, model_output, a, b, 0.1, 0.2);
[output_i, error_i, a_i, b_i, u_i] = mrac_indirect(time, reference, model_output, a, b, 5, 2);

plot_results(time, model_output, real_output, output_d, output_i, error_d, error_i, u_opt, u_d, u_i, k_star, k_d, a_i, l_star, l_d, b_i, 'Case 2');


%% Function definitions

function [model_output, real_output, u_opt] = simulate_system(time, reference, a, b, k_star, l_star, model_tf)
    s = tf('s');
    model_output = lsim(model_tf, reference, time);
    system_tf = l_star*b/(s+a+b*k_star);
    real_output = lsim(system_tf, reference, time);
    u_opt = -k_star * real_output' + l_star * reference;
end

function [output, error, k_est, l_est, u_control] = mrac_direct(time, reference, model_output, a, b, gamma1, gamma2)
    n = length(time);
    output = zeros(1,n);
    error = zeros(1,n);
    k_est = zeros(1,n);
    l_est = zeros(1,n);
    u_control = zeros(1,n);

    dt = time(end)/n;
    delta_output = 0;

    for i = 2:n
        u_control(i) = -k_est(i-1)*output(i-1) + l_est(i-1)*reference(i);
        output(i) = output(i-1) + dt*(-delta_output/a + (b/a)*u_control(i));
        delta_output = output(i) - output(i-1);

        error(i) = output(i) - model_output(i);
        k_est(i) = k_est(i-1) + gamma1 * error(i) * output(i);
        l_est(i) = l_est(i-1) - gamma2 * error(i) * reference(i);
    end
end

function [output, error, a_est, b_est, u_control] = mrac_indirect(time, reference, model_output, a, b, gamma1, gamma2)
    n = length(time);
    output = zeros(1,n);
    error = zeros(1,n);
    a_est = zeros(1,n);
    b_est = ones(1,n);
    u_control = zeros(1,n);

    param_a_model = 6.7;
    param_b_model = 3.3;
    b_min = 1;

    dt = time(end)/n;
    delta_output = 0;

    for i = 2:n
        u_control(i) = -output(i-1)*(param_a_model + a_est(i-1))/b_est(i-1) + param_b_model*reference(i-1)/b_est(i-1);

        output(i) = output(i-1) + dt*(-(output(i-1) - delta_output)/a + (b/a)*u_control(i));
        delta_output = output(i) - output(i-1);

        error(i) = output(i) - model_output(i);

        a_est(i) = a_est(i-1) + gamma1 * error(i) * output(i);

        if abs(b_est(i-1)) > b_min || (b_est(i-1) == b_min && error(i)*u_control(i) >= 0)
            b_est(i) = b_est(i-1) + gamma2 * error(i) * u_control(i);
        else
            b_est(i) = b_est(i-1);
        end
    end
end

function plot_results(time, model_output, real_output, output_direct, output_indirect, error_direct, error_indirect, u_opt, u_direct, u_indirect, k_star, k_direct, k_indirect, l_star, l_direct, l_indirect, label)

    figure
    plot(time, model_output, 'r', time, real_output, 'b--', time, real_output - model_output, 'LineWidth', 1.5);
    legend('Reference model', 'Optimal control output','Error');
    title([label ' - Optimal control']);
    grid on;

    figure
    plot(time, model_output, 'r', time, output_direct, 'b', time, error_direct, 'LineWidth', 1.5);
    legend('Reference model', 'Direct MRAC','Error');
    title([label ' - Direct MRAC']);
    grid on;

    figure
    plot(time, model_output, 'r', time, output_indirect, 'b', time, error_indirect, 'LineWidth', 1.5);
    legend('Reference model', 'Indirect MRAC','Error');
    title([label ' - Indirect MRAC']);
    grid on;

    figure
    plot(time, u_opt, 'r', time, u_direct, 'b', time, u_indirect, 'g--', 'LineWidth', 1.5);
    legend('Optimal','Direct MRAC','Indirect MRAC');
    title([label ' - Control signals']);
    grid on;

    k_ind = (6.7 + k_indirect)./l_indirect;
    k_star_vec = k_star * ones(size(time));

    figure
    plot(time, k_star_vec, 'r', time, k_direct, 'b', time, k_ind, 'g', 'LineWidth', 1.5);
    legend('Optimal','Direct MRAC','Indirect MRAC');
    title([label ' - K(t)']);
    grid on;

    l_ind = 3.3 ./ l_indirect;
    l_star_vec = l_star * ones(size(time));

    figure
    plot(time, l_star_vec, 'r', time, l_direct, 'b', time, l_ind, 'g', 'LineWidth', 1.5);
    legend('Optimal','Direct MRAC','Indirect MRAC');
    title([label ' - L(t)']);
    grid on;
end
