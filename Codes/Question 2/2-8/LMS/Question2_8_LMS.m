% Adaptive Control - Simulation 1
% Masoud Pourghavam
% Student Number: 810601044
% Question 2-8 LMS

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
[num,den] = tfdata(discrete_TF);
num = cell2mat(num);
den = cell2mat(den);

a_1 = den(2); a_2 = den(3);
a_3 = den(4); a_4 = den(5);
b_1 = num(2); b_2 = num(3);
b_3 = num(4); b_4 = num(5);

theta = [a_1 ; a_2 ; a_3 ; a_4 ; b_1 ; b_2 ; b_3 ; b_4];
white_noise = sqrt(noise_variance)*randn(Samples,1);
white_noise = white_noise - mean(white_noise);
NOISE = white_noise;
u_in = sqrt(noise_variance)*randn(Samples,1);

y_out = zeros(length(u_in),1);
y_out(1,1)=0;
y_out(2,1)=[-y_out(1,1) 0 0 0 u_in(1,1) 0 0 0]*theta+NOISE(2,1);
y_out(3,1)=[-y_out(2,1) -y_out(1,1) 0 0 u_in(2,1) u_in(1,1) 0 0]*theta+NOISE(3,1);
y_out(4,1)=[-y_out(3,1) -y_out(2,1) -y_out(1,1) 0 ...
    u_in(3,1) u_in(2,1) u_in(1,1) 0]*theta+NOISE(4,1);
for t=5:length(y_out)
    y_out(t,1)=[-y_out(t-1,1) -y_out(t-2,1) -y_out(t-3,1) ...
        -y_out(t-4,1) u_in(t-1,1) u_in(t-2,1) u_in(t-3,1) u_in(t-4,1)]*theta+NOISE(t,1);
end

theta_hat = zeros(q,1); phi_t = zeros(1,q); L = 1; E = 10; LF = 0; e = 0; e_param = 0;

a_1_hat = zeros(Samples,1); a_2_hat = zeros(Samples,1); 
a_3_hat = zeros(Samples,1); a_4_hat = zeros(Samples,1);
b_1_hat = zeros(Samples,1); b_2_hat = zeros(Samples,1);
b_3_hat = zeros(Samples,1); b_4_hat = zeros(Samples,1);

sample_number = zeros(Samples,1);
for t=1:Samples
    sample_number(t,1) = t;
end
%%

for L=2:Samples
    for i=1:4
        if L-i<=0
            y=0;
            u=0;
        else 
            y=-y_out(L-i,1);
            u=u_in(L-i,1);
        end
        phi_t(1,i)=y;
        phi_t(1,i+4)=u;
    end
   
    gamma=0.01;

    gain=(gamma*(phi_t'));
    theta_hat=theta_hat+gain*(y_out(L,1)-(phi_t*theta_hat));
    
    a_1_hat(L,1)=theta_hat(1); a_2_hat(L,1)=theta_hat(2);
    a_3_hat(L,1)=theta_hat(3); a_4_hat(L,1)=theta_hat(4);
    b_1_hat(L,1)=theta_hat(5); b_2_hat(L,1)=theta_hat(6);
    b_3_hat(L,1)=theta_hat(7); b_4_hat(L,1)=theta_hat(8);
    
    LF=(abs(y_out(L,1)-(phi_t*theta_hat)));
    e=e+(LF^2);
    e_param = e_param+(a_1_hat(L,1)-a_1)^2+(a_2_hat(L,1)-a_2)^2+ ...
        (a_3_hat(L,1)-a_3)^2+(a_4_hat(L,1)-a_4)^2+...
        (b_1_hat(L,1)-b_1)^2+(b_2_hat(L,1)-b_2)^2+(b_3_hat(L,1)-b_3)^2+(b_4_hat(L,1)-b_4)^2;
    E=LF;
  
end

%%

a_1_real = zeros(Samples,1);
for t=1:Samples
    a_1_real(t,1) = a_1;
end
a_2_real=zeros(Samples,1);
for t=1:Samples
    a_2_real(t,1) = a_2;
end

a_3_real=zeros(Samples,1);
for t=1:Samples
    a_3_real(t,1) = a_3;
end

a_4_real=zeros(Samples,1);
for t=1:Samples
    a_4_real(t,1) = a_4;
end

b_1_real = zeros(Samples,1);
for t=1:Samples
    b_1_real(t,1) = b_1;
end
b_2_real=zeros(Samples,1);
for t=1:Samples
    b_2_real(t,1) = b_2;
end 

b_3_real=zeros(Samples,1);
for t=1:Samples
    b_3_real(t,1) = b_3;
end

b_4_real=zeros(Samples,1);
for t=1:Samples
    b_4_real(t,1) = b_4;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure()
plot(sample_number,a_1_real,"black")
hold on
plot(sample_number,a_1_hat,"green")
ylim ([-5 5])
legend('a_1 real','a_1 hat')
xlabel('Samples')
ylabel('Output')

figure()
plot(sample_number,a_2_real,"black")
hold on
plot(sample_number,a_2_hat,"green")
ylim ([-5 5])
legend('a_2 real','a_2 hat')
xlabel('Samples')
ylabel('Output')

figure()
plot(sample_number,a_3_real,"black")
hold on
plot(sample_number,a_3_hat,"green")
ylim ([-5 5])
legend('a_3 real','a_3 hat')
xlabel('Samples')
ylabel('Output')

figure()
plot(sample_number,a_4_real,"black")
hold on
plot(sample_number,a_4_hat,"green")
ylim ([-5 5])
legend('a_4 real','a_4 hat')
xlabel('Samples')
ylabel('Output') 

figure()
plot(sample_number,b_1_real,"black")
hold on
plot(sample_number,b_1_hat,"green")
ylim ([-2 2])
legend('b_1 real','b_1 hat')
xlabel('Samples')
ylabel('Output')

figure()
plot(sample_number,b_2_real,"black")
hold on
plot(sample_number,b_2_hat,"green")
ylim ([-2 2])
legend('b_2 real','b_2 hat')
xlabel('Samples')
ylabel('Output')

figure()
plot(sample_number,b_3_real,"black")
hold on
plot(sample_number,b_3_hat,"green")
ylim ([-2 2])
legend('b_3 real','b_3 hat')
xlabel('Samples')
ylabel('Output')

figure()
plot(sample_number,b_4_real,"black")
hold on
plot(sample_number,b_4_hat,"green")
ylim ([-2 2])
legend('b_4 real','b_4 hat')
xlabel('Samples')
ylabel('Output')