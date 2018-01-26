clc
clear all;
close all;

a = 15; %m
b = 24.5; %m
sigma_a = 0.03; %m
sigma_b = 0.04; %m
r_ab = 0.3;
sigma_ab = sigma_a*sigma_b*r_ab
sigma_ba = sigma_ab

%% Functional model
%Calculation of the a%rea of the rectangle
A = a*b;

% Design Matrix
syms a0 b0
J = jacobian([a0*b0], [a0, b0])
F = subs(J, [b0,a0], [b,a]);

% Stochastic Model
sigma_ll = [[sigma_a.^2 sigma_ab]; [sigma_ba sigma_b.^2]]
% v = [10^(-18) 10^(-18)];
% SIGMA_ll2 = diag(v)

%VCMatrix
sigma_xx =F*sigma_ll*F';
sigma_A = sqrt(sigma_xx)