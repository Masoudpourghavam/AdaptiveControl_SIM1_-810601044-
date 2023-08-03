% Adaptive Control - Simulation 1
% Masoud Pourghavam
% Student Number: 810601044
% Question 2-4

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
a_1 = den(2); a_2 = den(3); a_3 = den(4); a_4 = den(5); ...
    b_1 = num(2); b_2 = num(3); b_3 = num(4); b_4 = num(5);
theta = [a_1 ; a_2 ; a_3 ; a_4 ; b_1 ; b_2 ; b_3 ; b_4];

%% Define a white noise:
white_noise = sqrt(noise_variance)*randn(Samples,1);
white_noise = white_noise - mean(white_noise);
NOISE = white_noise;

u = sqrt(noise_variance)*randn(Samples,1);

y_out = zeros(length(u),1);
y_out(1,1) = 0;
y_out(2,1) = [-y_out(1,1) 0 0 0 u(1,1) 0 0 0]*theta+NOISE(2,1);
y_out(3,1) = [-y_out(2,1) -y_out(1,1) 0 0 u(2,1) u(1,1) 0 0]*theta+NOISE(3,1);
y_out(4,1) = [-y_out(3,1) -y_out(2,1) -y_out(1,1) 0 u(3,1) u(2,1) u(1,1) 0]*theta+NOISE(4,1);
for t = 5:length(y_out)
    y_out(t,1) = [-y_out(t-1,1) -y_out(t-2,1) -y_out(t-3,1) ...
        -y_out(t-4,1) u(t-1,1) u(t-2,1) u(t-3,1) u(t-4,1)]*theta+NOISE(t,1);
end

theta_hat = zeros(q,1); init_value = 2000*eye(q); phi_t = zeros(1,q); PHI_org = zeros(Samples,q); e = 10; LF = 0;

a_1_hat = zeros(Samples,1); a_2_hat = zeros(Samples,1);
a_3_hat = zeros(Samples,1); a_4_hat = zeros(Samples,1);
b_1_hat = zeros(Samples,1); b_2_hat = zeros(Samples,1);
b_3_hat = zeros(Samples,1); b_4_hat = zeros(Samples,1);

sample_number = zeros(Samples,1);

for t = 1:Samples
    sample_number(t,1) = t;
end


for L = 2:Samples
    for i = 1:4
        if L-i<=0
            y_i=0;
            u_i=0;
        else 
            y_i = -y_out(L-i,1);
            u_i = u(L-i,1);
        end
        phi_t(1,i) = y_i;
        phi_t(1,i+4) = u_i;
    end
    init_value = init_value-((init_value*(phi_t')*phi_t*init_value)/(1+phi_t*init_value*(phi_t')));
    gain = init_value*(phi_t');
    theta_hat = theta_hat+gain*(y_out(L,1)-(phi_t*theta_hat));
    PHI_org(L,:) = phi_t;
    
a_1_hat(L,1) = theta_hat(1); a_2_hat(L,1) = theta_hat(2);
a_3_hat(L,1) = theta_hat(3); a_4_hat(L,1) = theta_hat(4);
b_1_hat(L,1) = theta_hat(5); b_2_hat(L,1) = theta_hat(6);
b_3_hat(L,1) = theta_hat(7); b_4_hat(L,1) = theta_hat(8);
    
    LF = LF+((y_out(L,1)-(phi_t*theta_hat))^2); e = LF;
end

y_hat = zeros(Samples,1);
y_hat = PHI_org*theta_hat;

%%  Over parameterization 
q_overparameterized = 9;
phi_t = zeros(1 , q_overparameterized);
init_value = 2000*eye(q_overparameterized);
teta_hat_overparameterized = zeros(q_overparameterized , 1);
PHI = zeros(Samples , q_overparameterized);
LF_overparameterized = 0;
e_overparameterized = 0;

 for L = 2:Samples
    for i = 2:5
        if L-i<=0
            y_i = 0;
            u_i = 0;
        else 
            y_i=-y_out(L-i,1);
            u_i=u(L-i,1);
        end
        phi_t(1,i) = y_i;
        phi_t(1,i+4) = u_i;
    end
    phi_t(1,1) = -y_out(L-1,1);
    init_value = init_value-((init_value*(phi_t')*phi_t*init_value)/(1+phi_t*init_value*(phi_t')));
    gain = init_value*(phi_t');
    teta_hat_overparameterized = teta_hat_overparameterized+gain*(y_out(L,1)-(phi_t*teta_hat_overparameterized));
    PHI(L,:) = phi_t;
    LF_overparameterized = LF_overparameterized+((y_out(L,1)-(phi_t*teta_hat_overparameterized))^2);
    e_overparameterized = LF_overparameterized;
 end

y_hat_overparameterized = zeros(Samples,1);
y_hat_overparameterized = PHI*teta_hat_overparameterized;

%%  Under Parameterization
q_underparameterized = 7;
phi_t = zeros(1,q_underparameterized);
teta_hat_underparameterized = zeros(q_underparameterized,1);
init_value = 1000*eye(q_underparameterized);
PHI_underparameterized = zeros(Samples,q_underparameterized);
LF_underparameterized = 0;
sse_under = 0;

 for L = 2:Samples
    for i = 2:4
        if L-i<=0
            y_i = 0;
            u_i = 0;
        else 
            y_i = -y_out(L-i,1);
            u_i = u(L-i,1);
        end
        phi_t(1,i) = y_i;
        phi_t(1,i+3) = u_i;
    end
    phi_t(1,1) = -y_out(L-1,1);
    init_value = init_value-((init_value*(phi_t')*phi_t*init_value)/(1+phi_t*init_value*(phi_t')));
    gain = init_value*(phi_t');
    teta_hat_underparameterized = teta_hat_underparameterized+gain*(y_out(L,1)-(phi_t*teta_hat_underparameterized));
    PHI_underparameterized(L,:) = phi_t;
    LF_underparameterized = LF_underparameterized+((y_out(L,1)-(phi_t*teta_hat_underparameterized))^2);
    sse_under=LF_underparameterized;
 end
y_hat_underparameterized = zeros(Samples,1);
y_hat_underparameterized = PHI_underparameterized*teta_hat_underparameterized;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure()
plot(sample_number,y_out,"black")
title('Real output')

figure()
plot(sample_number,y_hat_overparameterized,"green")
title('Predicted over parameterized output')

figure
plot(sample_number,y_out,"black")
hold on 
plot(sample_number,y_hat_overparameterized,"green")
xlabel('Samples')
ylabel('Output')
title('Real and predicted over parameterized systems')
legend('Real output','Predicted over parameterized output')

figure()
plot(sample_number,y_out,"black")
title('Real output')

figure()
plot(sample_number,y_hat_underparameterized,"green")
title('Predicted under parameterized output')

figure
plot(sample_number,y_out,"black")
hold on 
plot(sample_number,y_hat_underparameterized,"green")
xlabel('samples')
ylabel('Output')
title('Real and predicted under parameterized systems')
legend('Real output','Predicted under parameterized output')