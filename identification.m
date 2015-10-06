%% 

clear all, close all
load('data_group_35');

a = 0.1;
b1 = -0.1;
b2 = 0.1;
e0 = 0.2;
e1 = -0.1;
f = 0.1;

x0 = [a;b1;b2;e0;e1;f];
ub = [inf,inf,inf,inf,inf,inf];
lb = [-inf,-inf,-inf,-inf,-inf,-inf];
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'MaxFunEvals', 2000, 'DerivativeCheck', 'off', 'Jacobian', 'off');

x = lsqnonlin( @(x)cost_function(x,u,ys), x0, lb, ub, options)
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