% Adaptive Control - Simulation 1
% Masoud Pourghavam
% Student Number: 810601044
% Question 4-1 Kalman Filter

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
q = 8;
noise_variance = 0.058;
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
[num , den] = tfdata(discrete_TF);
num = cell2mat(num);
den = cell2mat(den);

a_1 = den(2); a_2 = den(3); a_3 = den(4); a_4 = den(5); b_1 = num(2); b_2 = num(3); b_3 = num(4); b_4 = num(5);

theta = [a_1 ; a_2 ; a_3 ; a_4 ; b_1 ; b_2 ; b_3 ; b_4];
white_noise = sqrt(noise_variance)*randn(Samples,1);
white_noise = white_noise - mean(white_noise);
u_in = sqrt(noise_variance)*randn(Samples,1);

y_out = zeros(Samples,1);
y_out(1,1) = 0;
y_out(2,1) = [-y_out(1,1) 0 0 0 u_in(1,1) 0 0 0]*theta+white_noise(2,1);
y_out(3,1) = [-y_out(2,1) -y_out(1,1) 0 0 u_in(2,1) u_in(1,1) 0 0]*theta+white_noise(3,1);
y_out(4,1) = [-y_out(3,1) -y_out(2,1) -y_out(1,1) 0 u_in(3,1) ...
    u_in(2,1) u_in(1,1) 0]*theta+white_noise(4,1);
for t = 5:length(y_out)
    y_out(t,1) = [-y_out(t-1,1) -y_out(t-2,1) -y_out(t-3,1) ...
        -y_out(t-4,1) u_in(t-1,1) u_in(t-2,1) u_in(t-3,1) u_in(t-4,1)]*theta+white_noise(t,1);
end

%%

init_value = 2000*eye(q); phi_t = zeros(q,1); H = zeros(1,q); L = 1; e = 10; LF = 0;

a_1_hat = zeros(Samples,1); a_2_hat = zeros(Samples,1);
a_3_hat = zeros(Samples,1); a_4_hat = zeros(Samples,1);
b_1_hat = zeros(Samples,1); b_2_hat = zeros(Samples,1);
b_3_hat = zeros(Samples,1); b_4_hat = zeros(Samples,1);

sample_number = zeros(Samples,1);
for t = 1:Samples
    sample_number(t,1) = t;
end

U = eye(q); PHI = zeros(Samples,q);

%%

 for L = 2:Samples
    for i = 1:4
        if L-i<=0
            y = 0;
            u = 0;
        else 
            y = -y_out(L-i,1);
            u = u_in(L-i,1);
        end
        H(1,i) = y;
        H(1,i+4) = u;
    end
   
    init_value_0 = U*init_value*(U'); phi_t_0 = U*phi_t;
    y_tilda = y_out(L,1)-(H*phi_t_0);
    N = (H*init_value*(H')) + noise_variance;
    k_N = init_value_0*(H')/N;
    init_value = (eye(q)-k_N*H)*init_value_0;
    phi_t = phi_t_0 + (k_N*y_tilda);
    PHI(L,:) = phi_t;
    
    a_1_hat(L,1) = phi_t(1); a_2_hat(L,1) = phi_t(2);
    a_3_hat(L,1) = phi_t(3); a_4_hat(L,1) = phi_t(4);
    b_1_hat(L,1) = phi_t(5); b_2_hat(L,1) = phi_t(6);
    b_3_hat(L,1) = phi_t(7); b_4_hat(L,1) = phi_t(8);
    
    LF = LF + 0.5*((y_out(L,1)-(H*phi_t))^2);
    e = LF;
    
 end
 

%%

a_1_real = zeros(Samples,1); a_2_real = zeros(Samples,1);
a_3_real = zeros(Samples,1); a_4_real = zeros(Samples,1);
b_1_real = zeros(Samples,1); b_2_real = zeros(Samples,1);
b_3_real = zeros(Samples,1); b_4_real = zeros(Samples,1);

a_1_real(:,1) = a_1; a_2_real(:,1) = a_2;
a_3_real(:,1) = a_3; a_4_real(:,1) = a_4;
b_1_real(:,1) = b_1; b_2_real(:,1) = b_2;
b_3_real(:,1) = b_3; b_4_real(:,1) = b_4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure()
plot(sample_number,a_1_real,"black")
hold on
plot(sample_number,a_1_hat,"green")
xlim ([0 1000])
legend('a_1 real','a_1 hat')
xlabel('Samples')
ylabel('Output')

figure()
plot(sample_number,a_2_real,"black")
hold on
plot(sample_number,a_2_hat,"green")
xlim ([0 1000])
legend('a_2 real','a_2 hat')
xlabel('Samples')
ylabel('Output')

figure()
plot(sample_number,a_3_real,"black")
hold on
plot(sample_number,a_3_hat,"green")
xlim ([0 1000])
legend('a_3 real','a_3 hat')
xlabel('Samples')
ylabel('Output')

figure()
plot(sample_number,a_4_real,"black")
hold on
plot(sample_number,a_4_hat,"green")
xlim ([0 1000])
legend('a_4 real','a_4 hat')
xlabel('Samples')
ylabel('Output') 

figure()
plot(sample_number,b_1_real,"black")
hold on
plot(sample_number,b_1_hat,"green")
xlim ([0 1000])
legend('b_1 real','b_1 hat')
xlabel('Samples')
ylabel('Output')

figure()
plot(sample_number,b_2_real,"black")
hold on
plot(sample_number,b_2_hat,"green")
xlim ([0 1000])
legend('b_2 real','b_2 hat')
xlabel('Samples')
ylabel('Output')

figure()
plot(sample_number,b_3_real,"black")
hold on
plot(sample_number,b_3_hat,"green")
xlim ([0 1000])
legend('b_3 real','b_3 hat')
xlabel('Samples')
ylabel('Output')

figure()
plot(sample_number,b_4_real,"black")
hold on
plot(sample_number,b_4_hat,"green")
xlim ([0 1000])
legend('b_4 real','b_4 hat')
xlabel('Samples')
ylabel('Output')