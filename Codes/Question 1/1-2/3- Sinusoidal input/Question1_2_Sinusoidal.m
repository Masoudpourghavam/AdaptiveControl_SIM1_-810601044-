% Adaptive Control - Simulation 1
% Masoud Pourghavam
% Student Number: 810601044
% Question 1-2 Sinusoidal input

%% --------------------------------------------- %%
clear all;
close all;
clc;

%% Define the parameters of system
k_1 = 0.6;
k_2 = 0.26;
k_3 = 0.4;
b_1 = 0.26;
b_2 = 0.6;
m_1 = 1;
m_2 = 1;
Samples = 1000;
noise_variance = 0.058;
q = 8;
s = tf('s');

%% Define Transfer Function:
TF = (k_3)/(-(k_3^2)+(m_2*(s^2)+b_2*s+(k_2+k_3))*(m_1*(s^2)+b_1*s+(k_1+k_3)));

%% Define Band width:
bw = bandwidth(TF);

%% Define Sample Time:
SampleTime = 6.28/(5*bw);

%% Define Discrete Transfer Function:
discrete_TF = c2d(TF , SampleTime ,'zoh');

%% Define nominators and denominators:
[num,den] = tfdata(discrete_TF);
num = cell2mat(num);
den = cell2mat(den);

%% Define a1...a4 and b1...b4 and theta matrix of coefficients:
a_1 = den(2); a_2 = den(3); a_3 = den(4); a_4 = den(5); b_1 = ...
    num(2); b_2 = num(3); b_3 = num(4); b_4 = num(5);
theta = [a_1 ; a_2 ; a_3 ; a_4 ; b_1 ; b_2 ; b_3 ; b_4];

%% Define a white noise:
white_noise = sqrt(noise_variance)*randn(Samples,1); white_noise = white_noise - mean(white_noise);

%% Define an input:
u = zeros(Samples,1);
omega = 0.032;
 for t = 1:Samples
   u(t,1) = sin(omega*t*SampleTime);
end

y_with_noise = zeros(length(u),1); y_with_noise(1,1) = 0;
y_with_noise(2,1) = [-y_with_noise(1,1) 0 0 0 u(1,1) ...
    0 0 0]*theta+white_noise(2,1);
y_with_noise(3,1) = [-y_with_noise(2,1) -y_with_noise(1,1) ...
    0 0 u(2,1) u(1,1) 0 0]*theta+white_noise(3,1);
y_with_noise(4,1) = [-y_with_noise(3,1) -y_with_noise(2,1) -y_with_noise(1,1) ...
    0 u(3,1) u(2,1) u(1,1) 0]*theta+white_noise(4,1);
for t = 5:length(y_with_noise)
    y_with_noise(t,1) = [-y_with_noise(t-1,1) -y_with_noise(t-2,1) ...
        -y_with_noise(t-3,1) -y_with_noise(t-4,1) u(t-1,1) u(t-2,1) u(t-3,1) ...
        u(t-4,1)]*theta+white_noise(t,1);
end

Phi = zeros(length(y_with_noise),q);
Phi(2,:) = [-y_with_noise(1,1) 0 0 0 u(1,1) 0 0 0];
Phi(3,:) = [-y_with_noise(2,1) -y_with_noise(1,1) 0 0 u(2,1) u(1,1) 0 0];
Phi(4,:) = [-y_with_noise(3,1) -y_with_noise(2,1) -y_with_noise(1,1) ...
    0 u(3,1) u(2,1) u(1,1) 0];
for t = 5:length(y_with_noise)
    Phi(t,:) = [-y_with_noise(t-1,1) -y_with_noise(t-2,1) ...
        -y_with_noise(t-3,1) -y_with_noise(t-4,1) u(t-1,1) u(t-2,1) u(t-3,1) u(t-4,1)];
end

theta_ls = (inv((transpose(Phi))*(Phi)))*(transpose(Phi))*y_with_noise;

y_without_noise = zeros(length(u),1); y_without_noise(1,1) = 0;
y_without_noise(2,1) = [-y_without_noise(1,1) 0 0 0 u(1,1) 0 0 0]*theta;
y_without_noise(3,1) = [-y_without_noise(2,1) -y_without_noise(1,1) ...
    0 0 u(2,1) u(1,1) 0 0]*theta;
y_without_noise(4,1) = [-y_without_noise(3,1) -y_without_noise(2,1) ...
    -y_without_noise(1,1) 0 u(3,1) u(2,1) u(1,1) 0]*theta;
for t = 5:length(y_without_noise)
    y_without_noise(t,1) = [-y_without_noise(t-1,1) -y_without_noise(t-2,1) ...
        -y_without_noise(t-3,1) -y_without_noise(t-4,1) u(t-1,1) ...
        u(t-2,1) u(t-3,1) u(t-4,1)]*theta;
end

z = tf('z');
TF_LS = (theta_ls(5)*z^3+theta_ls(6)*z^2+theta_ls(7)* ... 
    (z)+theta_ls(8))/(z^4+theta_ls(1)*z^3+theta_ls(2)*z^2+theta_ls(3)*z+theta_ls(4));


sample_number = zeros(Samples,1);
for t = 1:Samples
    sample_number(t,1) = t;
end

y_hat = Phi*theta_ls;
error = y_with_noise - y_hat;
Loss_Function = sum((y_with_noise - y_hat).^2);
sigma_hat = ((error')*error)/(Samples-q);
cov_thetahat_ls = sigma_hat*((inv((transpose(Phi))*(Phi))));

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure()
bode(discrete_TF, "black")
title('Discerete system')

figure()
bode(TF_LS, "black")
title('Peredicted system')

figure()
bode(discrete_TF)
hold on
bode(TF_LS)
legend('Discerete real system','Peredicted system')

figure()
plot(sample_number,error,"black")
title('Error of estimation')

figure()
plot(sample_number,y_with_noise,"black")
title('Real output')

figure()
plot(sample_number,y_hat,"green")
title('Predicted output')

figure()
plot(sample_number,y_with_noise,"black")
hold on
plot(sample_number,y_hat,"green")
xlabel('Samples')
ylabel('Outputs')
title('Real and predicted systems output')
legend('Real output','Prediction Output')