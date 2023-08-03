% Adaptive Control - Simulation 1
% Masoud Pourghavam
% Student Number: 810601044
% Question 5
%% --------------------------------------------- %%
clear all;
close all;
clc;

%%

k_1 = 0.6;
k_2 = 0.26;
k_3 = 0.4;
b_1 = 0.26;
b_2 = 0.6;
m_1 = 1;
m_2 = 1;
a_12 = 2.2;
a_21 = -4;
a_22 = -0.34;
s = tf('s');
TF = (k_3)/(-(k_3^2)+(m_2*(s^2)+b_2*s+(k_2+k_3))*(m_1*(s^2)+b_1*s+(k_1+k_3)));
N = 1000;
theta = [a_12 ; a_22];
noise_variance = 0.058;
q = 8;
white_noise = sqrt(noise_variance)*randn(N,1);
white_noise = white_noise - mean(white_noise);

u_in = sqrt(noise_variance)*randn(N,1);

y_out = zeros(N,1);
y_out(1,1) = 0;
y_out(2,1) = [0 y_out(1,1)]*theta+white_noise(2,1);
for t=3:N
   y_out(t,1) = [sin(a_12*y_out(t-2,1)+u_in(t-2,1)) y_out(t-1,1)]*theta+white_noise(t,1); 
end

L = 1; e = 10; LF = 0;

sample_number = zeros(N,1);
for t = 1:N
    sample_number(t,1) = t;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = 9; y_hat = zeros(N,M); E = zeros(M,1); AIC = zeros(M,1);

for n = 1:M
P = 2000*eye(2*n);
phi_t = zeros(1,2*n);
theta_hat = zeros(2*n,1);
PHI = zeros(N,2*n);
 for L = 2:N
    for i = 1:n
        if L-i<=0
            y = 0;
            u = 0;
        else 
            y = -y_out(L-i,1);
            u = u_in(L-i,1);
        end
        phi_t(1,i) = y;
        phi_t(1,i+n) = u;
    end
    PHI(L,:) = phi_t;
    P = P-((P*(phi_t')*phi_t*P)/(1+phi_t*P*(phi_t')));
    gain = P*(phi_t');
    theta_hat = theta_hat+gain*(y_out(L,1)-(phi_t*theta_hat));
    a21_hat(L,1) = theta_hat(1);
    a22_hat(L,1) = theta_hat(2);
    LF = LF+0.5*((y_out(L,1)-(phi_t*theta_hat))^2);
    e = LF;
    
    if L>=N
        break
    end
 end
 
y_hat_0 = PHI*theta_hat;
y_hat(:,n) = y_hat_0;
e = y_out-y_hat_0;
V_theta = 0.5*((e')*e);
E(n,1) = V_theta;

AIC(n,1) = N*log(V_theta+2*(2*n));

end

lin_M = linspace(1,M,M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
plot(sample_number,y_out, "black")
xlabel('Samples')
ylabel('Outputs')
title('Output of nonlinear system')

figure
plot(lin_M , AIC , "black");
xlabel('Estimation Number')
ylabel('AIC')

figure
plot(sample_number,y_out, "black")
hold on 
plot(sample_number,y_hat(:,1), "green")
title('P = 1')

figure
plot(sample_number,y_out, "black")
hold on 
plot(sample_number,y_hat(:,2), "green")
title('P = 2')

figure
plot(sample_number,y_out, "black")
hold on 
plot(sample_number,y_hat(:,3), "green")
title('P = 3')

figure
plot(sample_number,y_out, "black")
hold on 
plot(sample_number,y_hat(:,4), "green")
title('P = 4')

figure
plot(sample_number,y_out, "black")
hold on 
plot(sample_number,y_hat(:,5), "green")
title('P = 5')

figure
plot(sample_number,y_out, "black")
hold on 
plot(sample_number,y_hat(:,6), "green")
title('P = 6')

figure
plot(sample_number,y_out, "black")
hold on 
plot(sample_number,y_hat(:,7), "green")
title('P = 7')

figure
plot(sample_number,y_out, "black")
hold on 
plot(sample_number,y_hat(:,8), "green")
title('P = 8')

figure
plot(sample_number,y_out, "black")
hold on 
plot(sample_number,y_hat(:,9), "green")
title('P = 9')