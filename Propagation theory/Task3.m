%Several unknowns with uncorrelated obs
clear all
close all
clc

v0=15; %m/s
sigma_v0=0.1; %m/s
a=45*pi()/180; %rad
sigma_a=0.08*pi()/180; %rad
g=9.81; %m/s2

%% Calculation of the Height H of the center of the ring of fire   
H = v0^2*(sin(a))^2/(2*g);

% Calculation of the Lenght L of the center of the safety net 
L = v0^2*sin(2*a)/g;

%% Calculation of the standard deviation of the Height H and the Lenght L
syms v00 a0
J = jacobian([v00^2*(sin(a0))^2/(2*g), v00^2*sin(2*a0)/g], [v00, a0])
F = subs(J, [v00, a0], [v0, a])

%weights
v = [(sigma_v0)^2  (sigma_a)^2];
SIGMA_ll = diag(v);
SIGMA_xx = F*SIGMA_ll*F';

sigma_H= sqrt(SIGMA_xx(1,1))
sigm_L= sqrt(SIGMA_xx(2,2))