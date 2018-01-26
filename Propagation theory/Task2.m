%Task 2
clear all
close all
clc

c = 299792458; %m
t0 = 0.8147236863; %s
t1 = 0.8147240201; %s
sigma_t = 10^(-9); %s


%% Calculation of the distance S between EDM and the reflector 
S = (t1-t0)*c/2;

%% Calculation of the standard deviation of the distance S
% Functional model : S=(t1-t0)*c/2 

%F = zeros(1, 2)
syms t00 t11
J = jacobian([(t11-t00)*c/2], [t00, t11])
F = subs(J, [t00, t11], [t0,t1]);

%weights
SIGMA_ll = zeros(2,2);
for i=1:2
   SIGMA_ll(i,i)=10^(-18) %s^2
end    
% v = [10^(-18) 10^(-18)];
% SIGMA_ll2 = diag(v)

%VCM
Sigma_xx = F*SIGMA_ll*F'

%Standard deviation
sigma_S = sqrt(Sigma_xx)

%% 
% for i=1:2
%     syms F SIGMA_ll
%     sigma_l=sqrt(symsum(F(i,:)^2*(SIGMA_ll(i,i)^2)))
% end
sigma_l=0.001;
sigma_t1 = sqrt(2)*sigma_l/c



