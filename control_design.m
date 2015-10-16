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
infnorm = norm(T,inf) % Robust for Dmin < 1/2.138

%% Noise Sensitivity
Nq = Q.num{1};
Dq = Q.den{1};

[Aq,Bq,Cq,Dq] = tf2ss(Nq,Dq);
[Af,Bf,Cf,Df] = tf2ss(Nf,Df);

%Aq = 4x4 %Ap = 2x2 %Af = 1
At = [  Aq,                 Bq*Cf,         zeros(4,4),    zeros(4,1),   zeros(4,2),     -Bq*Df*Cp;
        zeros(1,4),         Af,            zeros(1,4),    zeros(1,1),   zeros(1,2),     -Bf*Cp;
        zeros(4,4),         zeros(4,1),    Aq,            Bq*Cf,        Bq*Df*Cp,       -Bq*Df*Dp*Cp;
        zeros(1,4),         zeros(1,1),    zeros(1,4),    Af,           Bf*Cp,          -Bf*Dp*Cp;
        zeros(2,4),         zeros(2,1),    zeros(2,4),    zeros(2,1),   Ap,             -Bp*Cp;
        zeros(2,4),         zeros(2,1),    zeros(2,4),    zeros(2,1),   zeros(2,2),     Ap             ];
Bt = [  -Bq*Df*Dp,          -Bq*Df;
        -Bf*Dp,             -Bf;
        -Bq*Df*Dp*Dp,       -Bq*Df*Dp;
        -Bf*Dp*Dp,          -Bf*Dp;
        -Bp*Dp,             -Bp;
         Bp,                 zeros(2,1)  ];
    
    
%  At = [ Af,                 Bf*Cq,         zeros(1,2),    zeros(1,1),   zeros(1,4),     -Bf*Dq*Cp;
%         zeros(4,1),         Aq,            zeros(4,2),    zeros(4,1),   zeros(4,4),     -Bq*Cp;
%         zeros(2,1),         zeros(2,4),    Ap,            Bp*Cf,        Bp*Df*Cq,       -Bp*Df*Dp*Cp;
%         zeros(1,1),         zeros(1,4),    zeros(1,2),    Af,           Bf*Cq,          -Bf*Dq*Cp;
%         zeros(4,1),         zeros(4,4),    zeros(4,2),    zeros(4,1),   Aq,             -Bq*Cp;
%         zeros(2,1),         zeros(2,4),    zeros(2,2),    zeros(2,1),   zeros(2,4),     Ap             ];
% Bt = [  -Bf*Dq*Dp,          -Bf*Dq;
%         -Bq*Dp,             -Bq;
%         -Bp*Df*Dq*Dp,       -Bp*Df*Dq;
%         -Bf*Dq*Dp,          -Bf*Dq;
%         -Bq*Dp,             -Bq;
%         Bp,                 zeros(2,1)  ];


W = diag( [RMSd,RMSn] );
X = dlyap(At, Bt*W*Bt');
    
    
    
    
 