clear all, close all, clc

%% 
syms a b1 b2 e0 e1 f z;

A_s = expand( (1+f*z)*(1+b1*z+b2*z^2) )
B_s = expand( (e0+e1*z)*(1+a*z) )
C_s = expand( (e0+e1*z)*(1+b1*z+b2*z^2) )

%% load data
load('data_group_35.mat', 'u')
load('data_group_35.mat', 'ys')

%% identification
data = iddata(ys,u, 1);

system = armax(data, [3 3 3 1])
