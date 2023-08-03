% Adaptive Control - Simulation 1
% Masoud Pourghavam
% Student Number: 810601044
% Question 1-3

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
a_1 = den(2); a_2 = den(3); a_3 = den(4); a_4 = den(5); b_1  ...
    = num(2); b_2 = num(3); b_3 = num(4); b_4 = num(5);
theta = [a_1 ; a_2 ; a_3 ; a_4 ; b_1 ; b_2 ; b_3 ; b_4];

%% Define a white noise:
white_noise = sqrt(noise_variance)*randn(Samples,1);
white_noise = white_noise - mean(white_noise);

%% Define an input:
u = sqrt(noise_variance)*randn(Samples,1);

y_with_noise=zeros(Samples,1);
y_with_noise(1,1)=0;
y_with_noise(2,1)=[-y_with_noise(1,1) 0 0 0 u(1,1) ...
    0 0 0]*theta+white_noise(2,1);
y_with_noise(3,1)=[-y_with_noise(2,1) -y_with_noise(1,1) ...
    0 0 u(2,1) u(1,1) 0 0]*theta+white_noise(3,1);
y_with_noise(4,1)=[-y_with_noise(3,1) -y_with_noise(2,1) ...
    -y_with_noise(1,1) 0 u(3,1) u(2,1) u(1,1) 0]*theta+white_noise(4,1);
for t = 5:length(y_with_noise)
    y_with_noise(t,1)=[-y_with_noise(t-1,1) -y_with_noise(t-2,1) ...
        -y_with_noise(t-3,1) -y_with_noise(t-4,1) u(t-1,1) u(t-2,1) u(t-3,1) ...
        u(t-4,1)]*theta+white_noise(t,1);
end
Sample_number = zeros(Samples,1);
for t = 1:Samples
    Sample_number(t,1)=t;
end

phi = zeros(1,q); theta_hat = zeros(q,1); Phi_o = zeros(Samples,q);

for L = 2:Samples
    for i=1:4
          if L-i<=0
            y_i = 0;
            u_i = 0;
          else 
            y_i = -y_with_noise(L-i,1);
            u_i = u(L-i,1);
          end
        phi(1,i) = y_i;
        phi(1,i+4) = u_i;
    end

    Phi_o(L,:) = phi;
end

theta_hat = (((Phi_o')*(Phi_o))^-1)*(Phi_o')*y_with_noise;
y_hat = zeros(Samples,1);
y_hat = Phi_o*theta_hat;
error = zeros(Samples,1);
error = y_with_noise - y_hat;
LossFunction = error'*error;



q_overparameterized = 9;
phi = zeros(1,q_overparameterized);
theta_hat_overparametrized = zeros(q_overparameterized,1);
Phi_overparametrized = zeros(Samples,q_overparameterized);

%%

for L = 2:Samples
    for i=2:5
          if L-i<=0
            y_i = 0;
            u_i = 0;
          else 
            y_i = -y_with_noise(L-i,1);
            u_i = u(L-i,1);
          end
        phi(1,i) = y_i;
        phi(1,i+4) = u_i;
    end
    phi(1,1) = -y_with_noise(L-1,1);
    Phi_overparametrized(L,:) = phi;
end

theta_hat_overparametrized = (((Phi_overparametrized')* ... 
    (Phi_overparametrized))^-1)*(Phi_overparametrized')*y_with_noise;

y_hat_overparametrized = zeros(Samples,1);
y_hat_overparametrized = Phi_overparametrized*theta_hat_overparametrized;



q_underparameterized = 7;
phi = zeros(1,q_underparameterized);
theta_hat_underparametrized = zeros(q_underparameterized,1);
Phi_underparametrized = zeros(Samples,q_underparameterized);

for L = 2:Samples
    for i=2:4
          if L-i<=0
            y_i = 0;
            u_i = 0;
          else 
            y_i = -y_with_noise(L-i,1);
            u_i = u(L-i,1);
          end
        phi(1,i) = y_i;
        phi(1,i+3) = u_i;
    end
    phi(1,1) = -y_with_noise(L-1,1);
    Phi_underparametrized(L,:) = phi;
end
theta_hat_underparametrized = (((Phi_underparametrized')* ... 
    (Phi_underparametrized))^-1)*(Phi_underparametrized')*y_with_noise;
y_hat_underparameterized = zeros(Samples,1);
y_hat_underparameterized = Phi_underparametrized*theta_hat_underparametrized;



V_theta_hat_overparametrized = 0;
for t = 1:Samples
    V_theta_hat_overparametrized = V_theta_hat_overparametrized  ... 
        + (y_hat_overparametrized(t,1)-y_with_noise(t,1))^2;
end
error_overparametrized = zeros(Samples,1);
error_overparametrized = y_with_noise-y_hat_overparametrized;

V_theta_hat_underparametrized = 0;
for t = 1:Samples
    V_theta_hat_underparametrized = V_theta_hat_underparametrized+ ... 
        (y_hat_underparameterized(t,1)-y_with_noise(t,1))^2;
end
error_underparameterized = zeros(Samples,1);
error_underparameterized = y_with_noise-y_hat_underparameterized;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure()
plot(Sample_number , y_with_noise,"black")
title('Real output with noise')

figure()
plot(Sample_number,y_hat,"green")
title('Predicted output of system with correct parameter')

figure()
plot(Sample_number,y_with_noise,"black")
hold on 
plot(Sample_number,y_hat,"green")
xlabel('Samples')
ylabel('Outputs')
title('Real and correct parameter outputs')
legend('Real output with noise','Predicted output of system with correct parameter')

figure()
plot(Sample_number,y_with_noise,"black")
title('Real output')

figure()
plot(Sample_number,y_hat_overparametrized,"green")
title('Overparameter output')

figure()
plot(Sample_number,y_with_noise,"black")
hold on 
plot(Sample_number,y_hat_overparametrized,"green")
xlabel('Samples')
ylabel('Outputs')
title('Real and over parameter systems output')
legend('Real output','Over parameter predicted output')

figure()
plot(Sample_number,y_with_noise,"black")
title('Real output')

figure()
plot(Sample_number,y_hat_underparameterized,"green")
title('Under parameter output')

figure()
plot(Sample_number,y_with_noise,"black")
hold on 
plot(Sample_number,y_hat_underparameterized,"green")
xlabel('Samples')
ylabel('Outputs')
title('Real and under parameter systems output')
legend('Real output of system','Predicted output of under parameterized system')

figure()
plot(Sample_number,y_with_noise,"black")
title('Real output')

figure()
plot(Sample_number,y_hat_underparameterized,"green")
title('Under parameter output')

figure()
plot(Sample_number,y_with_noise,"black")
hold on 
plot(Sample_number,y_hat_underparameterized,"green")
xlabel('Samples')
ylabel('Outputs')
title('Real and under parameter systems output')
legend('Real output of system','Predicted output of under parameterized system')

figure() 
plot(Sample_number,error,"black")
title('Correct parameter error')

figure() 
plot(Sample_number,error_overparametrized,"black")
title('Over parameter error')

figure() 
plot(Sample_number,error_underparameterized,"black")
title('Under parameter error')

figure() 
plot(Sample_number,error,"black--")
hold on 
plot(Sample_number,error_overparametrized,"green--")
hold on 
plot(Sample_number,error_underparameterized,'m--')
legend('Correct parameter error','Over parameter error','Under parameter error')