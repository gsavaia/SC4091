clear all, close all

%% Symbolic Proof for Stability
syms Nq Dq Nk Dk Nf Df Np Dp %numerator/denominator of Q,K,F and P 
Q = Nq/Dq;
K = Nk/Dk;
P = Np/Dp;
F = Nf/Df;

K = Q/(1-P*Q*F);
L = K*P*F;
S = 1/(1+L);
disp('Sensitivity Function'); pretty(simplify(S)); 
%if Df,Dp,Dq contains stable poles then the closed loop tf will be stable as well

%% LQG controller
clear all; load('systems'); % LOAD SYSTEM

[Ap,Bp,Cp,Dp] = tf2ss(Np,Dp);
rho = 2e-4;
RMSd = 4;
RMSn = 1;

Klqg = dlqry(Ap,Bp,Cp,Dp,1,rho); %LQG control law
Lkalman = dlqe(Ap,Bp,Cp,RMSd,RMSn); %Kalman gain

%design regulator (LQG + Kalman)
[Ac,Bc,Cc,Dc] = dreg(Ap,Bp,Cp,Dp,Klqg,Lkalman); 
Kss = ss(Ac,Bc,Cc,Dc)
[Nk, Dk] = ss2tf(Ac,Bc,Cc,Dc);

K = tf(Nk,Dk,-1)
Q = feedback(K,P)

poles_lqg = pole(Q) % Q is AS stable

%% check if solution satisfy "Robustness" requirements
T = feedback(K*P,F);

figure, bode(T)
infnorm = norm(T,inf)


