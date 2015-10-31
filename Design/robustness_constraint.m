function [c, ceq, Gc, Gceq] = robustness_constraint(NumQ,DenQ,P,F,Dmin_inv)

    Q = tf(NumQ,DenQ,-1);
    [PQFnorm,w0] = hinfnorm(P*Q*F,1e-6);
    
    c = PQFnorm - Dmin_inv; % |PQ| < 1/Dmin
    ceq = [];

    if( nargout > 2 )
        zeta = [1; exp(-1j*w0); exp(-2j*w0); exp(-3j*w0); exp(-4j*w0)];
        
        Gc = real( ( PQFnorm * zeta ) / (NumQ * zeta) );
        Gceq = [];
    end
end