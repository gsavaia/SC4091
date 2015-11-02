clear all;

load system
load lqg_design

syms a0 a1 a2 a3 a4 % Optimization Vector

[Ap,Bp,Cp,Dp] = tf2ss(NumP,DenP);
[Af,Bf,Cf,Df] = tf2ss(NumF,DenF);

p = Q.den{1};
p = p(2:size(p,2));

Cq = [a1-a0*p(1), a2-a0*p(2), a3-a0*p(3), a4-a0*p(4)];
Dq = a0;

Ct = [  Cq,                 Dq*Cf,         zeros(1,4),  zeros(1,1),     zeros(1,2),     -Dq*Df*Cp;
        zeros(1,4),         zeros(1,1),    Cq,          Dq*Cf,          Dq*Df*Cp,       Cp-Dq*Df*Dp*Cp ];

Dt = [ -Dq*Df*Dp,           -Dq*Df;
        Dp-Dq*Df*Dp*Dp,     -Dq*Df*Dp ];

vector = [a0, a1, a2, a3, a4];

dCu = double( jacobian(Ct(1,:), vector) )';
dCy = double( jacobian(Ct(2,:), vector) )';
dDu = double( jacobian(Dt(1,:), vector) )';
dDy = double( jacobian(Dt(2,:), vector) )';

gradient.dCy = dCy;
gradient.dCu = dCu;
gradient.dDy = dDy;
gradient.dDu = dDu;

save('gradient_J', 'gradient');