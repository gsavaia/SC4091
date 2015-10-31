clear all;

load system
load lqg_design

% SPECS
rho = 2e-4; 
RMSd = 4;
RMSn = 1;

%% TRANSFER FUNCTION BETWEEN [U,Y] AND [D,N] (NOISE SENSITIVITY)
[Aq,Bq,Cq,Dq] = tf2ss(Q.num{1},Q.den{1}); % Q
[Af,Bf,Cf,Df] = tf2ss(NumF,DenF); %filter
[Ap,Bp,Cp,Dp] = tf2ss(NumP,DenP); %plant

%Aq = 4x4 %Ap = 2x2 %Af = 1
%Aq, Bq are known, that is they do not contain parameters to be optimized
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

%% ALTERNATIVE CANONICAL FORM
% Ct = [  Cq,                 Dq*Cf,         zeros(1,4),  zeros(1,1),     zeros(1,2),     -Dq*Df*Cp;
%         zeros(1,4),         zeros(1,1),    Cq,          Dq*Cf,          Dq*Df*Cp,       Cp-Dq*Df*Dp*Cp ];
% 
% Dt = [ -Dq*Df*Dp,           -Dq*Df;
%         Dp-Dq*Df*Dp*Dp,     -Dq*Df*Dp ];
    
% THIS REPRESENTATION LEADS TO MATRIXES WHICH CONTAINS THE VARIABLE TO BE
% OPTIMAZED, THEN THE COMPUTATION OF 'X' SHOULD BE INCLUDED IN THE
% OPTIMIZATION FUNCTION.
%
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

%% LYAPUNOV FUNCTION (needed to compute final objective function)
W = diag( [RMSd^2,RMSn^2] ); %RMS NOISE
X = dlyap(At, Bt*W*Bt'); %RMS STATE

%% OPTIMIZATION SET-UP

% Initialize optimization procedure and set options.
opt=optimoptions('fmincon', 'Algorithm','sqp');     % SQP algorithm.
opt=optimoptions(opt,'TolX',1e-12);       % Termination tolerance for parameters.
opt=optimoptions(opt,'TolFun',1e-6);     % Termination tolerance for objective function.
opt=optimoptions(opt,'TolCon',1e-6);     % Termination tolerance for constraints.
opt=optimoptions(opt,'MaxFunEval',5000);
opt=optimoptions(opt,'MaxIter', 5000);

opt=optimoptions(opt,'GradObj','on');    % Use gradient.
opt=optimoptions(opt,'GradConstr','on'); % Use constraint Jacobian.

%opt=optimoptions(opt,'Display','off') ;
%opt=optimoptions(opt,'Display','iter');  % Display results.
%opt=optimoptions(opt,'Diagnostics','on'); % Display diagnostic
%opt=optimoptions(opt,'DerivativeCheck','on','FinDiffType','central');

%% RUN OPT ALGO
NumQ = Q.num{1}; % Nq from LQG will be our starting point
DenQ = Q.den{1};
load gradient_J

Dmin_inv = 0.2:0.05:2.5;

NumQ_minima = zeros(length(Dmin_inv), 5);
J = zeros(1, length(Dmin_inv));
exitflag = zeros(1, length(Dmin_inv));

NumQ_init = ones(1,5);%NumQ;

for i=1:length(Dmin_inv)    
    [NumQ, J(i), exitflag(i), info(i)]=...
          fmincon( @(NumQ) noise_sensitivity(NumQ,DenQ,Cp,Dp,Cf,Df,X,W,rho,gradient),NumQ,... % goal function
                   [],[],[],[],[],[],...                               % linear constr
                   @(NumQ) robustness_constraint(NumQ,DenQ,P,F,Dmin_inv(i)), ...      % non-linear constr
                   opt); % options 
    
    NumQ_minima(i,:) = NumQ;
    
    disp(Dmin_inv(i));
    disp(NumQ);
    %disp(NumQ_init);
    disp(exitflag(i));
    disp('#######################################');
end

% figure, bar(Dmin_inv(exitflag>0),J(exitflag>0));
% 
% figure; title('Sensitivity for different values of Dmin');
% for i=1:length(Dmin_inv(exitflag>0))
%     Q_minima = tf(NumQ_minima(i,:), DenQ, -1);
%     K_minima = feedback(Q_minima, -P*F);
%     bode( feedback(1, K_minima*P*F) ); hold on;
% end
% 
% save('control_design', 'NumQ_minima', 'DenQ', 'Dmin_inv', 'J');