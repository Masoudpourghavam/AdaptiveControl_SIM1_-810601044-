% Adaptive Control - Simulation 1
% Masoud Pourghavam
% Student Number: 810601044
% Question 3-1

%% --------------------------------------------- %%
clear all;
close all;
clc;

%% Define the parameters of system
Samples = 1000;
noise_variance = 0.058;
q = 8;

%% Define Transfer Function:
z = tf('z');
TF = 0.59899*(z+1.253)*(z+0.221)/((z+0.4178)*(z+0.6611)*(z+0.9687));

%% Define Band width:
bw = bandwidth(TF);

%% Define Sample Time:
SampleTime = 6.28/(5*bw);

%% Define nominators and denominators:
[num,den] = tfdata(TF);
num = cell2mat(num);
den = cell2mat(den);

%% Define a1...a4 and b1...b4 and theta matrix of coefficients:
a_1 = den(2); a_2 = den(3); a_3 = den(4); ...
    b_1 = num(2); b_2 = num(3); b_3 = num(4);
theta = [a_1 ; a_2 ; a_3 ; b_1 ; b_2 ; b_3];

%% Define a white noise:
white_noise = sqrt(noise_variance)*randn(Samples,1);
white_noise = white_noise - mean(white_noise);
NOISE = white_noise;

u_in = sqrt(noise_variance)*randn(Samples,1);

y_out = zeros(length(u_in),1);
y_out(1,1) = 0;
y_out(2,1) = [-y_out(1,1) 0 0  u_in(1,1) 0 0 ]*theta+NOISE(2,1);
y_out(3,1) = [-y_out(2,1) -y_out(1,1) 0  u_in(2,1) u_in(1,1) 0 ]*theta+NOISE(3,1);
y_out(4,1) = [-y_out(3,1) -y_out(2,1) ...
    -y_out(1,1)  u_in(3,1) u_in(2,1) u_in(1,1) ]*theta+NOISE(4,1);
for t = 5:length(y_out)
    y_out(t,1) = [-y_out(t-1,1) -y_out(t-2,1) -y_out(t-3,1) ...
        u_in(t-1,1) u_in(t-2,1) u_in(t-3,1)]*theta+NOISE(t,1);
end

%% Feedback:
K_c = 0.1;
u_in = -(K_c)*(y_out)

theta_hat = zeros(q,1); init_value = 1000*eye(q); phi_t = zeros(1,q); L = 1; e = 10; LF = 0;

a_1_hat = zeros(Samples,1); a_2_hat = zeros(Samples,1);
a_3_hat = zeros(Samples,1); 
b_1_hat = zeros(Samples,1); b_2_hat = zeros(Samples,1);
b_3_hat = zeros(Samples,1); 

sample_number = zeros(Samples,1);
for t=1:Samples
    sample_number(t,1)=t;
end

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
        phi_t(1,i) = y;
        phi_t(1,i+4) = u;
    end
    init_value = init_value-((init_value*(phi_t')*phi_t*init_value)/(1+phi_t*init_value*(phi_t')));
    gain = init_value*(phi_t');
    theta_hat = theta_hat+gain*(y_out(L,1)-(phi_t*theta_hat));
    
a_1_hat(L,1) = theta_hat(1); a_2_hat(L,1) = theta_hat(2);
a_3_hat(L,1) = theta_hat(3); 
b_1_hat(L,1) = theta_hat(5); b_2_hat(L,1) = theta_hat(6);
b_3_hat(L,1) = theta_hat(7); 
    
    LF = LF+0.5*((y_out(L,1)-(phi_t*theta_hat))^2);
    e=LF;
end


a_1_real = zeros(Samples,1); a_2_real = zeros(Samples,1);
a_3_real = zeros(Samples,1); 
b_1_real = zeros(Samples,1); b_2_real = zeros(Samples,1);
b_3_real = zeros(Samples,1); 

a_1_real(:,1) = a_1; a_2_real(:,1) = a_2;
a_3_real(:,1) = a_3; 
b_1_real(:,1) = b_1; b_2_real(:,1) = b_2;
b_3_real(:,1) = b_3; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure()
plot(sample_number,a_1_real,"black")
hold on
plot(sample_number,a_1_hat,"green")
legend('a_1 real','a_1 hat')
xlabel('Samples')

figure()
plot(sample_number,a_2_real,"black")
hold on
plot(sample_number,a_2_hat,"green")
legend('a_2 real','a_2 hat')
xlabel('Samples')

figure()
plot(sample_number,a_3_real,"black")
hold on
plot(sample_number,a_3_hat,"green")
legend('a_3 real','a_3 hat')
xlabel('Samples')

figure()
plot(sample_number,b_1_real,"black")
hold on
plot(sample_number,b_1_hat,"green")
legend('b_1 real','b_1 hat')
xlabel('Samples')

figure()
plot(sample_number,b_2_real,"black")
hold on
plot(sample_number,b_2_hat,"green")
legend('b_2 real','b_2 hat')
xlabel('Samples')

figure()
plot(sample_number,b_3_real,"black")
hold on
plot(sample_number,b_3_hat,"green")
legend('b_3 real','b_3 hat')
xlabel('Samples')
