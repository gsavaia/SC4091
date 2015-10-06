syms a b1 b2 e0 e1 f q

Y=(1+f*q)/(e0+e1*q);
 
U=-(q*(1+q*a))/(1+b1*q+b2*q*q);
 
Ye1=diff(Y,e1)
Ye0=diff(Y,e0)
Yb1=diff(Y,b1)
Yb2=diff(Y,b2)
Ya=diff(Y,a)
Yf=diff(Y,f)
 
Ue1=diff(U,e1)
Ue0=diff(U,e0)
Ub1=diff(U,b1)
Ub2=diff(U,b2)
Ua=diff(U,a)
Uf=diff(U,f)