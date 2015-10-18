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
% ===> if Df,Dp,Dq contains stable poles then the closed loop tf will be stable as well