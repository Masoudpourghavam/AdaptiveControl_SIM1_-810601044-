% Adaptive Control - Simulation 1
% Masoud Pourghavam
% Student Number: 810601044
% Question 2-6

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
Samples = 5000;
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
a_1 = den(2); a_2 = den(3); a_3 = den(4); a_4 = den(5); ...
    b_1 = num(2); b_2 = num(3); b_3 = num(4); b_4 = num(5);
theta = [a_1 ; a_2 ; a_3 ; a_4 ; b_1 ; b_2 ; b_3 ; b_4];

%% Define a white and color noise:
white_noise = sqrt(noise_variance)*randn(Samples,1);
white_noise = white_noise - mean(white_noise);

color_noise = zeros(Samples,1);
color_noise(1,1) = white_noise(1,1);
color_coeff = 0.32;

for t = 2:Samples
color_noise(t,1) = white_noise(t,1)+color_coeff*white_noise(t-1,1);
end
NOISE = color_noise;

u_in = sqrt(noise_variance)*randn(Samples,1);

y_out = zeros(length(u_in),1);
y_out(1,1) = 0;
y_out(2,1) = [-y_out(1,1) 0 0 0 u_in(1,1) 0 0 0]*theta+NOISE(2,1);
y_out(3,1) = [-y_out(2,1) -y_out(1,1) 0 0 u_in(2,1) u_in(1,1) 0 0]*theta+NOISE(3,1);
y_out(4,1) = [-y_out(3,1) -y_out(2,1) ...
    -y_out(1,1) 0 u_in(3,1) u_in(2,1) u_in(1,1) 0]*theta+NOISE(4,1);
for t = 5:length(y_out)
    y_out(t,1) = [-y_out(t-1,1) -y_out(t-2,1) -y_out(t-3,1) -y_out(t-4,1) ...
        u_in(t-1,1) u_in(t-2,1) u_in(t-3,1) u_in(t-4,1)]*theta+NOISE(t,1);
end

theta_hat = zeros(q,1); init_value = 1000*eye(q); phi_t = zeros(1,q); L = 1; e = 10; LF = 0;

a_1_hat = zeros(Samples,1); a_2_hat = zeros(Samples,1);
a_3_hat = zeros(Samples,1); a_4_hat = zeros(Samples,1);
b_1_hat = zeros(Samples,1); b_2_hat = zeros(Samples,1);
b_3_hat = zeros(Samples,1); b_4_hat = zeros(Samples,1);

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
a_3_hat(L,1) = theta_hat(3); a_4_hat(L,1) = theta_hat(4);
b_1_hat(L,1) = theta_hat(5); b_2_hat(L,1) = theta_hat(6);
b_3_hat(L,1) = theta_hat(7); b_4_hat(L,1) = theta_hat(8);
    
    LF = LF+0.5*((y_out(L,1)-(phi_t*theta_hat))^2);
    e=LF;
end


a_1_real = zeros(Samples,1); a_2_real = zeros(Samples,1);
a_3_real = zeros(Samples,1); a_4_real = zeros(Samples,1);
b_1_real = zeros(Samples,1); b_2_real = zeros(Samples,1);
b_3_real = zeros(Samples,1); b_4_real = zeros(Samples,1);

a_1_real(:,1) = a_1; a_2_real(:,1) = a_2;
a_3_real(:,1) = a_3; a_4_real(:,1) = a_4;
b_1_real(:,1) = b_1; b_2_real(:,1) = b_2;
b_3_real(:,1) = b_3; b_4_real(:,1) = b_4;

%%

ep = zeros(Samples,1); q_els = q+1;
theta_els = zeros(q_els,1); theta_els = [theta;color_coeff];
theta_hat_els = zeros(q_els,1); init_value_els = 1000*eye(q_els);
phi_t_els = zeros(1,q_els);
i=1; e_els = 10; LF_els = 0;
    
a_1_hat_els = zeros(Samples,1); a_2_hat_els = zeros(Samples,1);
a_3_hat_els = zeros(Samples,1); a_4_hat_els = zeros(Samples,1);
b_1_hat_els = zeros(Samples,1); b_2_hat_els = zeros(Samples,1);
b_3_hat_els = zeros(Samples,1); b_4_hat_els = zeros(Samples,1);
color_coeff_hat = zeros(Samples,1);
    
    while(e_els > 10^-12)
         i=i+1;
         
         for q=1:4
             if i-q<=0
                 y=0;
                 u=0;
             else 
            y = -y_out(i-q,1);
            u = u_in(i-q,1);
             end
             phi_t_els(1,q) = y;
             phi_t_els(1,q+4) = u;
         end
          phi_t_els(1,end) = ep(i-1,1);
          ep(i,1) = y_out(i,1)-phi_t_els*theta_hat_els;
          init_value_els = init_value_els-((init_value_els*(phi_t_els')*phi_t_els*init_value_els)/(1+phi_t_els*init_value_els*(phi_t_els')));
          theta_hat_els = theta_hat_els+(init_value_els*(phi_t_els')*ep(i,1));
          
a_1_hat_els(i,1) = theta_hat_els(1); a_2_hat_els(i,1) = theta_hat_els(2);
a_3_hat_els(i,1) = theta_hat_els(3); a_4_hat_els(i,1) = theta_hat_els(4);
b_1_hat_els(i,1) = theta_hat_els(5); b_2_hat_els(i,1) = theta_hat_els(6);
b_3_hat_els(i,1) = theta_hat_els(7); b_4_hat_els(i,1) = theta_hat_els(8);
color_coeff_hat(i,1) = theta_hat_els(9);
          
          LF_els = LF_els+0.5*((y_out(i,1)-(phi_t_els*theta_hat_els))^2);
          e_els = LF_els;
          if i>=Samples
              break
          end   
    end
color_coeff_real = zeros(Samples,1);
color_coeff_real(:,1) = color_coeff;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure()
plot(sample_number,a_1_real,"black")
hold on
plot(sample_number,a_1_hat_els,"green")
legend('a_1 real','a_1 hat els')

figure()
plot(sample_number,a_2_real,"black")
hold on
plot(sample_number,a_2_hat_els,"green")
legend('a_2 real','a_2 hat els')

figure()
plot(sample_number,a_3_real,"black")
hold on
plot(sample_number,a_3_hat_els,"green")
legend('a_3 real','a_3 hat els')

figure()
plot(sample_number,a_4_real,"black")
hold on
plot(sample_number,a_4_hat_els,"green")
legend('a_4 real','a_4 hat els')
 
figure()
plot(sample_number,b_1_real,"black")
hold on
plot(sample_number,b_1_hat_els,"green")
legend('b_1 real','b_1 hat els')

figure()
plot(sample_number,b_2_real,"black")
hold on
plot(sample_number,b_2_hat_els,"green")
legend('b_2 real','b_2 hat els')

figure()
plot(sample_number,b_3_real,"black")
hold on
plot(sample_number,b_3_hat_els,"green")
legend('b_3 real','b_3 hat els')

figure()
plot(sample_number,b_4_real,"black")
hold on
plot(sample_number,b_4_hat_els,"green")
legend('b_4 real','b_4 hat els')

figure
plot(sample_number,color_coeff_real, "black")
hold on
plot(sample_number,color_coeff_hat,"green")
legend('color coefficinet_{real}','color coefficient_{hat}')