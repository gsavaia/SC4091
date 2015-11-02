syms a b1 b2 e0 e1 f q

F_inv =  (1+f*q)/(e0+e1*q);
P     = -(q*(1+q*a))/(1+b1*q+b2*q*q);

Jsym = jacobian(F_inv+P, [a b1 b2 e0 e1 f]);
pretty( Jsym )