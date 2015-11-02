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
end