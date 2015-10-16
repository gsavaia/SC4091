function [J,G]=mcd_J(num_q,den_q,Cp,Dp,Xt,Wt,rho)
% MCD_J  Calculates the function to be minimized (J) and its gradient (G)

% State space representations of Q
[Aq,Bq,Cq,Dq]=tf2ss(num_q,den_q);

% Compute noise sensitivity using Lyapunov theory
Ct=[ Cq, zeros(1,6), -Dq*Cp
       zeros(1,4), Cq, Dq*Cp, Cp-Dq*Dp*Cp ];
Dt=[ -Dq*Dp, -Dq
       Dp-Dq*Dp*Dp, -Dq*Dp ];

Zt=Ct*Xt*Ct'+Dt*Wt*Dt';
J=Zt(2,2)+rho*Zt(1,1);

if ( nargout > 1 )
   % Compute gradient of objective function.
   dc1=[[-den_q(2:5); eye(4)], zeros(5,6), [-Cp; zeros(4,2)]];
   dc2=[zeros(5,4), [-den_q(2:5); eye(4)], [Cp; zeros(4,2)], ...
        [-Dp*Cp; zeros(4,2)]];
   dd1=[[-Dp, -1]; zeros(4,2)];
   dd2=[[-Dp*Dp, -Dp]; zeros(4,2)];
   G=2*(dc2*Xt*Ct(2,:)'+dd2*Wt*Dt(2,:)'+...
         rho*(dc1*Xt*Ct(1,:)'+dd1*Wt*Dt(1,:)'));
end;
