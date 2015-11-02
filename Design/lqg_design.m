%% LQG CONTROLLER
clear all; load('system'); % LOAD SYSTEM

[Ap,Bp,Cp,Dp] = tf2ss(NumP,DenP); %plant
rho = 2e-4; 
RMSd = 4;
RMSn = 1;

Klqg = dlqry(Ap,Bp,Cp,Dp,1,rho); %LQG control law
Lkalman = dlqe(Ap,Bp,Cp,RMSd^2,RMSn^2); %Kalman gain

%% DESIGN REGULATOR (LQ + Kalman)
[Ac,Bc,Cc,Dc] = dreg(Ap,Bp,Cp,Dp,Klqg,Lkalman); 
Kss = ss(Ac,Bc,Cc,Dc);
[NumK, DenK] = ss2tf(Ac,Bc,Cc,Dc);

K = tf(NumK,DenK,-1)
Q = feedback(K,P)

poles_lqg = pole(Q) % Q is AS stable

% check if solution satisfy "Robustness" requirements
T = feedback(K*P*F,1);
S = feedback(1, K*P*F);

figure; title('Sensitivity Functions (LQG)'); 
bode(T), hold on;
bode(S), legend('T','S');
infnorm = hinfnorm(T) 
W_thresh = 1/infnorm % Robust for Dmin < 0.3839

save('lqg_design', 'K','Q','poles_lqg');