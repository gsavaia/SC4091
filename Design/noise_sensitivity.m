function [J, G] = noise_sensitivity(Nq, Dq, Cp, Dp, Cf, Df, X, W, rho)

    [Aq,Bq,Cq,Dq] = tf2ss(Nq,Dq);

    Ct = [  Cq,                 Dq*Cf,         zeros(1,4),  zeros(1,1),     zeros(1,2),     -Dq*Df*Cp;
            zeros(1,4),         zeros(1,1),    Cq,          Dq*Cf,          Dq*Df*Cp,       Cp-Dq*Df*Dp*Cp ];

    Dt = [ -Dq*Df*Dp,           -Dq*Df;
            Dp-Dq*Df*Dp*Dp,     -Dq*Df*Dp ];

    ZZ = Ct*X*Ct' + Dt*W*Dt';    

    J = ZZ(2,2) + rho * ZZ(1,1);

    if nargout > 1

        Cu = Ct(1,:); Du = Dt(1,:);
        Cy = Ct(2,:); Dy = Dt(2,:);

        dDu = [ 0 0 0 0 0;
                -0.790000000000117 0 0 0 0]';
        dDy = zeros(2,5)';
        
        dCu = [ -1 1 0 0 0;
                -0.372967707790674 0 1 0 0;
                -0.0927460164703470 0 0 1 0;
                -0.0360416683099865 0 0 0 1;
                0.277800000242183 0 0 0 0;
                zeros(7,5);
                -0.790000000000117 0 0 0 0;
                -0.300199999758189 0 0 0 0  ]';
        dCy = [ zeros(5,5);
                -1 1 0 0 0;
                -0.372967707790674 0 1 0 0;
                -0.0927460164703470 0 0 1 0;
                -0.0360416683099865 0 0 0 1;
                0.277800000242183 0 0 0 0;
                0.790000000000117 0 0 0 0;
                0.300199999758189 0 0 0 0;
                zeros(2,5)                  ]';
        
        G = 2*( dCy*X*Cy' + dDy*W*Dy' ) + 2*rho*( dCu*X*Cu' + dDu*W*Du' ); 
    end

end