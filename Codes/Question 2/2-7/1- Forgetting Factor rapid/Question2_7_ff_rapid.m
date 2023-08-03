% Adaptive Control - Simulation 1
% Masoud Pourghavam
% Student Number: 810601044
% Question 2-7 forgetting factor rapid

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
Samples = 6000;
change = 0.04;
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

a_1=zeros(Samples,1);
for t=1:Samples
    if t<=Samples/4
        a_1(t,1) = den(2);
    elseif t>Samples/4 && t<=2*Samples/4
        a_1(t,1) = den(2)+change*den(2);
    elseif t>2*Samples/4 && t<=3*Samples/4
        a_1(t,1) = den(2);
    elseif  t>3*Samples/4 && t<=4*Samples/4
        a_1(t,1) = den(2)-change*den(2);
    elseif t>4*Samples/4 && t<=Samples
        a_1(t,1) = den(2);
    end
end

a_2=zeros(Samples,1);
for t=1:Samples
    if t<=Samples/4
        a_2(t,1) = den(3);
    elseif t>1*Samples/4 && t<=2*Samples/4
        a_2(t,1) = den(3)+change*den(3);
    elseif t>2*Samples/4 && t<=3*Samples/4
        a_2(t,1) = den(3);
    elseif  t>3*Samples/4 && t<=4*Samples/4
        a_2(t,1) = den(3)-change*den(3);
    elseif t>4*Samples/4 && t<=Samples
        a_2(t,1) = den(3);
    end
end

a_3=zeros(Samples,1);
for t=1:Samples
    if t<=Samples/4
        a_3(t,1) = den(4);
    elseif t>1*Samples/4 && t<=2*Samples/4
        a_3(t,1) = den(4)+change*den(4);
    elseif t>2*Samples/4 && t<=3*Samples/4
        a_3(t,1) = den(4);
    elseif  t>3*Samples/4 && t<=4*Samples/4
        a_3(t,1) = den(4)-change*den(4);
    elseif t>4*Samples/4 && t<=Samples
        a_3(t,1) = den(4);
    end
end

a_4=zeros(Samples,1);
for t=1:Samples
    if t<=Samples/4
        a_4(t,1) = den(5);
    elseif t>1*Samples/4 && t<=2*Samples/4
        a_4(t,1) = den(5)+change*den(5);
    elseif t>2*Samples/4 && t<=3*Samples/4
        a_4(t,1) = den(5);
    elseif  t>3*Samples/4 && t<=4*Samples/4
        a_4(t,1) = den(5)-change*den(5);
    elseif t>4*Samples/4 && t<=Samples
        a_4(t,1) = den(5);
    end
end

b_1=zeros(Samples,1);
for t=1:Samples
    if t<=Samples/5
        b_1(t,1) = num(2);
    elseif t>Samples/4 && t<=2*Samples/4
        b_1(t,1) = num(2)+change*num(2);
    elseif t>2*Samples/4 && t<=3*Samples/4
        b_1(t,1) = num(2);
    elseif  t>3*Samples/4 && t<=4*Samples/4
        b_1(t,1) = num(2)-change*num(2);
    elseif t>4*Samples/4 && t<=Samples
        b_1(t,1) = num(2);
    end
end

b_2=zeros(Samples,1);
for t=1:Samples
    if t<=Samples/4
        b_2(t,1) = num(3);
    elseif t>Samples/4 && t<=2*Samples/4
        b_2(t,1) = num(3)+change*num(3);
    elseif t>2*Samples/4 && t<=3*Samples/4
        b_2(t,1) = num(3);
    elseif  t>3*Samples/4 && t<=4*Samples/4
        b_2(t,1) = num(3)-change*num(3);
    elseif t>4*Samples/4 && t<=Samples
        b_2(t,1) = num(3);
    end
end

b_3=zeros(Samples,1);
for t=1:Samples
    if t<=Samples/4
        b_3(t,1) = num(4);
    elseif t>Samples/4 && t<=2*Samples/4
        b_3(t,1) = num(4)+change*num(4);
    elseif t>2*Samples/4 && t<=3*Samples/4
        b_3(t,1) = num(4);
    elseif  t>3*Samples/4 && t<=4*Samples/4
        b_3(t,1) = num(4)-change*num(4);
    elseif t>4*Samples/4 && t<=5*Samples/4
        b_3(t,1) = num(4);
    end
end

b_4=zeros(Samples,1);
for t=1:Samples
    if t<=Samples/4
        b_4(t,1) = num(5);
    elseif t>Samples/4 && t<=2*Samples/4
        b_4(t,1) = num(5)+change*num(5);
    elseif t>2*Samples/4 && t<=3*Samples/4
        b_4(t,1) = num(5);
    elseif  t>3*Samples/4 && t<=4*Samples/4
        b_4(t,1) = num(5)-change*num(5);
    elseif t>4*Samples/4 && t<=5*Samples/4
        b_4(t,1) = num(5);
    end
end

theta = [a_1 , a_2 , a_3 , a_4 , b_1 , b_2 , b_3 , b_4];

noise_variance = 0.01;
white_noise = sqrt(noise_variance)*randn(Samples,1);
white_noise = white_noise - mean(white_noise);
NOISE = white_noise;

u_in=sqrt(noise_variance)*randn(Samples,1);

%%

y_out = zeros(length(u_in),1);
y_out(1,1) = 0;
y_out(2,1) = [-y_out(1,1) 0 0 0 u_in(1,1) 0 0 0]*(theta(2,:)')+NOISE(2,1);
y_out(3,1) = [-y_out(2,1) -y_out(1,1) 0 0 u_in(2,1) u_in(1,1) 0 0]*(theta(3,:)')+NOISE(3,1);
y_out(4,1) = [-y_out(3,1) -y_out(2,1) -y_out(1,1) 0 u_in(3,1) ... 
    u_in(2,1) u_in(1,1) 0]*(theta(4,:)')+NOISE(4,1);
for t=5:length(y_out)
    y_out(t,1) = [-y_out(t-1,1) -y_out(t-2,1) -y_out(t-3,1) ...
        -y_out(t-4,1) u_in(t-1,1) u_in(t-2,1) u_in(t-3,1) u_in(t-4,1)]*(theta(t,:)')+NOISE(t,1);
end

%%

theta_hat = zeros(q,1); init_value=2000*eye(q); phi_t=zeros(1,q); LF=0;

a_1_hat = zeros(Samples,1); a_2_hat=zeros(Samples,1);
a_3_hat = zeros(Samples,1); a_4_hat=zeros(Samples,1);
b_1_hat = zeros(Samples,1); b_2_hat=zeros(Samples,1);
b_3_hat = zeros(Samples,1); b_4_hat=zeros(Samples,1);

sample_number = zeros(Samples,1);
for t=1:Samples
    sample_number(t,1) = t;
end

%%

 for L = 1:Samples
   
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
    
    Forgetting_Coeff = 0.999;

    gain = (init_value*(phi_t'))/(Forgetting_Coeff+phi_t*init_value*(phi_t'));
    init_value = (init_value-((gain)*phi_t*init_value))/Forgetting_Coeff;
    theta_hat=theta_hat+gain*(y_out(L,1)-(phi_t*theta_hat));
    
    a_1_hat(L,1) = theta_hat(1); a_2_hat(L,1) = theta_hat(2);
    a_3_hat(L,1) = theta_hat(3); a_4_hat(L,1) = theta_hat(4);
    b_1_hat(L,1) = theta_hat(5); b_2_hat(L,1) = theta_hat(6);
    b_3_hat(L,1) = theta_hat(7); b_4_hat(L,1) = theta_hat(8);
    LF = (abs(y_out(L,1)-(phi_t*theta_hat)));
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure()
plot(sample_number,theta(:,1),"black")
hold on
plot(sample_number,a_1_hat,"green")
legend('a_1 real','a_1 hat')

figure()
plot(sample_number,theta(:,2),"black")
hold on
plot(sample_number,a_2_hat,"green")
legend('a_2 real','a_2 hat')

figure()
plot(sample_number,theta(:,3),"black")
hold on
plot(sample_number,a_3_hat,"green")
legend('a_3 real','a_3 hat')

figure()
plot(sample_number,theta(:,4),"black")
hold on
plot(sample_number,a_4_hat,"green")
legend('a_4 real','a_4 hat')

figure()
plot(sample_number,theta(:,5),"black")
hold on
plot(sample_number,b_1_hat,"green")
legend('b_1 real','b_1 hat')

figure()
plot(sample_number,theta(:,6),"black")
hold on
plot(sample_number,b_2_hat,"green")
legend('b_2 real','b_2 hat')

figure()
plot(sample_number,theta(:,7),"black")
hold on
plot(sample_number,b_3_hat,"green")
legend('b_3 real','b_3 hat')

figure()
plot(sample_number,theta(:,8),"black")
hold on
plot(sample_number,b_4_hat,"green")
legend('b_4 real','b_4 hat')