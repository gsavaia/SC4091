clear all, close all

load gradient_J
load lqg_design
load system

rho = 2e-4; 
RMSd = 4;
RMSn = 1;
NumQ = Q.num{1};
DenQ = Q.den{1};

[Ap,Bp,Cp,Dp] = tf2ss(NumP,DenP);
[Aq,Bq,Cq,Dq] = tf2ss(NumQ,DenQ);
[Af,Bf,Cf,Df] = tf2ss(NumF,DenF); %filter

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
     
W = diag( [RMSd^2,RMSn^2] ); %RMS NOISE
X = dlyap(At, Bt*W*Bt'); %RMS STATE

tau = eye(5);

%% OBJECTIVE FUNCTION
[J,G] = noise_sensitivity(NumQ,DenQ,Cp,Dp,Cf,Df,X,W,rho,gradient);
for i=1:5
    [Jtau(i)] = noise_sensitivity(NumQ+tau(i,:),DenQ,Cp,Dp,Cf,Df,X,W,rho,gradient);
    Gtau(i) = G'*tau(:,i);
    result(i) = Jtau(i)-J-G'*tau(:,i) >= 0;
end

result,J, Jtau, Gtau
table_report = [Jtau;J*ones(1,5);Gtau]'

%% CONSTRAINT
[J,~,G,~] = robustness_constraint(NumQ,DenQ,P,F,1);
for i=1:5
    [Jtau(i)] = robustness_constraint(NumQ+tau(i,:),DenQ,P,F,1);
    Gtau(i) = G'*tau(:,i);
    result(i) = Jtau(i)-J-G'*tau(:,i) >= 0;
end

result, J, Jtau, Gtau
table_report = [Jtau;J*ones(1,5);Gtau]'

