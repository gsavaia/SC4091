function [c, ceq, Gc, Gceq] = robustness_constraint(NumQ,DenQ,P,F,Dmin_inv)

    Q = tf(NumQ,DenQ,-1);
    
   % wr = 0:0.1:pi;
   % NumPQF = -conv( conv(P.num{1},NumQ), F.num{1} );
   % DenPQF = conv( conv(P.den{1},DenQ), F.den{1} );
    
   % values = freqz( NumPQF, DenPQF, wr);
   % [PQFnorm,w0] = max( abs( values ) );
   % w0 = wr(w0);
   
   [PQFnorm,w0] = hinfnorm(Q,1e-6);
    
    c = PQFnorm - Dmin_inv; % |PQ| < 1/Dmin
    ceq = [];

    if( nargout > 2 )
        zeta = [1; exp(-1j*w0); exp(-2j*w0); exp(-3j*w0); exp(-4j*w0)];
        
        Gc = real( ( PQFnorm * zeta ) / (NumQ * zeta) );
        Gceq = [];
    end
end