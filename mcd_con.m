function [c,ceq,Gc,Gceq]=mcd_con(num_q,den_q,dv_M11)
% MCD_CON  Calculates the constraints (c and ceq) and their Jacobians (GC
%          and GCeq)

ceq=[];   % No nonlinear equality constraints.
Gceq=[];

% Evaluate infinity norm of M11.
wr=[0:0.01:pi]; 
values_M11=freqz(num_q,den_q,wr);
[infn_M11,ind_M11]=max(abs(values_M11));

% Compute constraint for robustness to additive plant perturbation
c=infn_M11-dv_M11;

if ( nargout > 2 )
   % Compute subgradient of the infinity norm of M11
   w11=wr(ind_M11);
   buffer=[1; exp(-j*w11); exp(-j*2*w11); exp(-j*3*w11); exp(-j*4*w11)];
   Gc=real(infn_M11*buffer/(num_q*buffer));
end;
