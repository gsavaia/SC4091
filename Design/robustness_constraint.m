function [c, ceq] = robustness_constraint(Nq,Dq,P,Dmin)

Q = tf(Nq,Dq,-1);

c = norm(P*Q, inf) - 1/Dmin; % |PQ| < 1/Dmin
ceq = [];

end