%% 
clear all, close all
load('data_group_35');

%% Plot experimental data
figure, plot(u, 'blue'); hold on; plot(ys, 'red');

%% Find Minima
a = 0.1;
b1 = 0.1;
b2 = 0.1;
e0 = 0.2;
e1 = 0.1;
f = 0.1;

x0 = [a;b1;b2;e0;e1;f];
ub = [inf,inf,inf,inf,inf,inf];
lb = [-inf,-inf,-inf,-inf,-inf,-inf];
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'MaxFunEvals', 2000, 'DerivativeCheck', 'on', 'Jacobian', 'on');


[x,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqnonlin( @(x)cost_function(x,u,ys), x0, lb, ub, options);
disp(x);

%% Check Gradient = 0 for x = x_minima
grad=JACOBIAN'*RESIDUAL;

GRAD=sum(grad.^2)
RES=sum(RESIDUAL.^2)

EXITFLAG

%% System
NumP = [0 1 x(1)];
DenP = [1 x(2) x(3)];
P = tf( NumP, DenP, -1)

NumF = [x(4) x(5)];
DenF = [1 x(6)];
F = tf( NumF, DenF, 1)

save('system', 'P', 'NumP', 'DenP', 'F', 'NumF', 'DenF');

L = F*P

estimation_error = sum( (ys - lsim(L,u)).^2 )