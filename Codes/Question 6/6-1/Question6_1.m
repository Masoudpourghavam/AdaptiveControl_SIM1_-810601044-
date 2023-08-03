% Adaptive Control - Simulation 1
% Masoud Pourghavam
% Student Number: 810601044
% Question 6-1

%% --------------------------------------------- %%
clear all;
close all;
clc;

%% Load data from Excel file
data = readtable('Battery.xlsx');
y = data.Current;
u = data.Voltage;

%% Determine maximum order of the model to be tested
max_order = 10;

%% Initialize variables to store cost functions
J_white = zeros(max_order, 1);
J_color = zeros(max_order, 1);

%% Loop through different model orders
for n = 1:max_order
    % Estimate ARMAX model with white noise
    model_white = armax(iddata(y, u), [n, n, n, 0]);
    yhat_white = sim(model_white, iddata([], u));
    J_white(n) = sum((y - yhat_white.y).^2);

    % Estimate ARMAX model with color noise
    model_color = armax(iddata(y, u), [n, n, n, n]);
    yhat_color = sim(model_color, iddata([], u));
    J_color(n) = sum((y - yhat_color.y).^2);
end

%% Plot cost functions for different model orders
subplot(2,1,1)
plot(1:max_order, J_white)
title('Cost function for ARMAX model with white noise')
xlabel('Model order')
ylabel('Cost')
subplot(2,1,2)
plot(1:max_order, J_color)
title('Cost function for ARMAX model with color noise')
xlabel('Model order')
ylabel('Cost')

%% Define number of regressors
ny = 1; % output regressor order
nu = 1; % input regressor order
ne = 1; % noise regressor order

%% Set RLS algorithm parameters
lambda = 0.99; % forgetting factor
delta = 1e-6; % regularization factor

%% Initialize RLS algorithm variables
theta = zeros(ny + nu + ne, 1); % parameter vector
P = delta * eye(ny + nu + ne); % covariance matrix
y_prev = 0; % previous output value
u_prev = 0; % previous input value

%% Loop through data samples and update RLS algorithm variables
for k = 1:length(y)
    % Construct regressor vector
    phi_k = [-y_prev; u_prev; 1];
    
    % Compute RLS algorithm gain
    K_k = P * phi_k / (lambda + phi_k' * P * phi_k);
    
    % Compute prediction error
    e_k = y(k) - theta' * phi_k;
    
    % Update parameter vector and covariance matrix
    theta = theta + K_k * e_k;
    P = (P - K_k * phi_k' * P) / lambda;
    
    % Update previous output and input values
    y_prev = y(k);
    u_prev = u(k);
end

%% Define true model parameters
a_1_true = -1.5;
b_0_true = 0.6;
b_1_true = 0.3;

%% Print estimated model parameters
fprintf('Estimated ARMAX parameters:\n');
fprintf('a_1 = %.4f\n', -theta(1));
fprintf('b_0 = %.4f\n', theta(2));
fprintf('b_1 = %.4f\n', theta(3));

%% Plot estimated and true model parameters
figure;
bar([a_1_true b_0_true b_1_true; -theta(1) theta(2) theta(3)]);
xlabel('Parameter');
ylabel('Value');
legend('True', 'Estimated');
