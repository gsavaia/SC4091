function [J, G] = noise_sensitivity(NumQ, DenQ, Cp, Dp, Cf, Df, X, W, rho, grad)

    [Aq,Bq,Cq,Dq] = tf2ss(NumQ,DenQ);

    Ct = [  Cq,                 Dq*Cf,         zeros(1,4),  zeros(1,1),     zeros(1,2),     -Dq*Df*Cp;
            zeros(1,4),         zeros(1,1),    Cq,          Dq*Cf,          Dq*Df*Cp,       Cp-Dq*Df*Dp*Cp ];

    Dt = [ -Dq*Df*Dp,           -Dq*Df;
            Dp-Dq*Df*Dp*Dp,     -Dq*Df*Dp ];

    ZZ = Ct*X*Ct' + Dt*W*Dt';    

    J = ZZ(2,2) + rho * ZZ(1,1);

    if nargout > 1        
        Cu = Ct(1,:); Du = Dt(1,:);
        Cy = Ct(2,:); Dy = Dt(2,:);
        
        G = 2*( grad.dCy*X*Cy' + grad.dDy*W*Dy' ) + 2*rho*( grad.dCu*X*Cu' + grad.dDu*W*Du' ); 
    end

end