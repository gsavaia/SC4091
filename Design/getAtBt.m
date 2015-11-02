function [At,Bt] = getAtBt(numQ, DenQ, numF, DenF, numP, DenP)    
[Ap,Bp,Cp,Dp] = tf2ss(NumP,DenP);
[Af,Bf,Cf,Df] = tf2ss(NumF,DenF);
[Aq,Bq,Cq,Dq] = tf2ss(NumQ,DenQ);

At = ...
[Aq,         Bq*Cf,      zeros(4,4), zeros(4,1), zeros(4,2),-Bq*Df*Cp;
 zeros(1,4), Af,         zeros(1,4), zeros(1,1), zeros(1,2),-Bf*Cp;
 zeros(4,4), zeros(4,1), Aq,         Bq*Cf,      Bq*Df*Cp,  -Bq*Df*Dp*Cp;
 zeros(1,4), zeros(1,1), zeros(1,4), Af,         Bf*Cp,     -Bf*Dp*Cp;
 zeros(2,4), zeros(2,1), zeros(2,4), zeros(2,1), Ap,        -Bp*Cp;
 zeros(2,4), zeros(2,1), zeros(2,4), zeros(2,1), zeros(2,2), Ap     ];
 
Bt = [  -Bq*Df*Dp,          -Bq*Df;
        -Bf*Dp,             -Bf;
        -Bq*Df*Dp*Dp,       -Bq*Df*Dp;
        -Bf*Dp*Dp,          -Bf*Dp;
        -Bp*Dp,             -Bp;
         Bp,                 zeros(2,1)  ];
     
%%ALTERNATIVE CANONICAL FORM
% Ct = ...
%[Cq,         Dq*Cf,    zeros(1,4), zeros(1,1), zeros(1,2), -Dq*Df*Cp;
% zeros(1,4), zeros(1,1), Cq,       Dq*Cf,     Dq*Df*Cp, Cp-Dq*Df*Dp*Cp];
%
% Dt = [ -Dq*Df*Dp,           -Dq*Df;
%         Dp-Dq*Df*Dp*Dp,     -Dq*Df*Dp ];

% THIS REPRESENTATION LEADS TO MATRIXES WHICH CONTAINS THE VARIABLE TO BE
% OPTIMAZED, THEN THE COMPUTATION OF 'X' SHOULD BE INCLUDED IN THE
% OPTIMIZATION FUNCTION.
%
%  At = ...
%[ Af,         Bf*Cq,      zeros(1,2), zeros(1,1), zeros(1,4) -Bf*Dq*Cp;
%  zeros(4,1), Aq,         zeros(4,2), zeros(4,1), zeros(4,4), -Bq*Cp;
%  zeros(2,1), zeros(2,4), Ap,         Bp*Cf,      Bp*Df*Cq,-Bp*Df*Dp*Cp;
%  zeros(1,1), zeros(1,4), zeros(1,2), Af,         Bf*Cq,     -Bf*Dq*Cp;
%  zeros(4,1), zeros(4,4), zeros(4,2), zeros(4,1), Aq,        -Bq*Cp;
%  zeros(2,1), zeros(2,4), zeros(2,2), zeros(2,1), zeros(2,4), Ap      ];
% Bt = [  -Bf*Dq*Dp,          -Bf*Dq;
%         -Bq*Dp,             -Bq;
%         -Bp*Df*Dq*Dp,       -Bp*Df*Dq;
%         -Bf*Dq*Dp,          -Bf*Dq;
%         -Bq*Dp,             -Bq;
%         Bp,                 zeros(2,1)  ];
     
end