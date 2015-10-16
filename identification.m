%% 
clear all, close all
load('data_group_35');

%% PLOT DATA
figure, plot(u, 'blue'); hold on; plot(ys, 'red');

a = 0.1;
b1 = 0.1;
b2 = 0.1;
e0 = 0.2;
e1 = 0.1;
f = 0.1;

x0 = [a;b1;b2;e0;e1;f];
ub = [inf,inf,inf,inf,inf,inf];
lb = [-inf,-inf,-inf,-inf,-inf,-inf];
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'MaxFunEvals', 2000, 'DerivativeCheck', 'off', 'Jacobian', 'off');

%% Find Minima
x = lsqnonlin( @(x)cost_function(x,u,ys), x0, lb, ub, options);
disp(x);

%% Check Gradient = 0 for x = x_minima
grad = gradient_cost_function(x, u, ys);

%% System
Np = [0 1 x(1)];
Dp = [1 x(2) x(3)];
P = tf( Np, Dp, -1)

Nf = [x(4) x(5)];
Df = [1 x(6)];
F = tf( Nf, Df, 1)

L = F*P

estimation_error = sum( (ys - lsim(L,u)).^2 );

disp(estimation_error);
%%
% a = 0.1;
% b1 = 0.1;
% b2 = 0.1;
% e0 = 0.2;
% e1 = 0.1;
% f = 0.1;
%     0.3800
%    -0.9500
%     0.6500
%     0.7900
%     0.4200
%     0.1800
%%
% a = 0.1;
% b1 = -0.1;
% b2 = 0.1;
% e0 = 0.2;
% e1 = -0.1;
% f = 0.1;
    
%     0.5316
%    -0.9500
%     0.6500
%     0.7900
%     0.3002
%     0.1800