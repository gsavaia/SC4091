function [c, ceq, Gc, Gceq] = robustness_constraint(NumQ,DenQ,P,F,Dmin_inv)

    Q = tf(NumQ,DenQ,-1);

    [PQFnorm,w0] = hinfnorm(P*Q*F);
    
%     wr = 0:0.1:pi;
%     values = freqz(conv(NumQ,P.num{1}),conv(DenQ,P.den{1}),wr);
%     [PQnorm,w0] = max( abs( values ) );
    
    c = PQFnorm - Dmin_inv; % |PQ| < 1/Dmin
    ceq = [];

    if( nargout > 2 )
        delay = [1; exp(-1j*w0); exp(-2j*w0); exp(-3j*w0); exp(-4j*w0)];
        
        Gc = real( ( PQFnorm * delay ) / (NumQ * delay) );
        Gceq = [];
    end
end