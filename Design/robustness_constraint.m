function [c, ceq, Gc, Gceq] = robustness_constraint(NumQ,DenQ,P,Dmin_inv)

    Q = tf(NumQ,DenQ,-1);

    [PQnorm,w0] = hinfnorm(P*Q);
    
    c = PQnorm - Dmin_inv; % |PQ| < 1/Dmin
    ceq = [];

    if( nargout > 2 )
        delay = [1; exp(-1j*w0); exp(-2j*w0); exp(-3j*w0); exp(-4j*w0)];
        
        Gc = real( ( PQnorm * delay ) / (Nq * delay) );
        Gceq = [];
    end
end