%Task 4
clear all
close all
clc

%Vector of observations 
a1 = 35.1550*pi()/200; %rad
%a1 = subs(a1)
a2 = 55.1200*pi()/200; %rad
%a2 = subs(a2)
s1 = 20.005; %m
s2 = 30.001; %m
t1 = 9.7; %s
t2 = 23.1; %s

x1 = s1*sin(a1)
x2 = s2*sin(a2)
y1 = s1*cos(a1)
y2 = s2*cos(a2)
v=sqrt((x2-x1)^2+(y2-y1)^2)/(t2-t1)

%omali kinisi..provlepsi tritis thesis
% azim23 = atan((x2-x1)/(y2-y1)) !!!orizontia pros katakorufi diastasi
% s23 = v*(t3-t2)
% x3=x2+s23*sin(azim23)
% y3=y2+s23*cos(azim23)

sigma_a = 0.001*pi/200
sigma_s = 0.001
sigma_t = 0.1
     
syms a11 a22 s11 s22 t11 t22
J1 = jacobian([s11*sin(a11) s22*sin(a22) s11*cos(a11) s22*cos(a22) t11 t22], [a11, a22, s11, s22, t11, t22])
F1 = subs(J1, [a11, a22, s11, s22, t11, t22], [a1, a2, s1, s2, t1, t2]);
     
syms x11 x22 y11 y22 t11 t22
J2 = jacobian([sqrt((x22-x11)^2+(y22-y11)^2)/(t22-t11)], [x11, x22, y11, y22, t11, t22])
F2= subs(J2, [x11, x22, y11, y22, t11, t22], [x1, x2, y1, y2, t1, t2]);

F = F2*F1

%VC Matrix of the observations
Sigma_LL = diag([sigma_a sigma_a sigma_s sigma_s sigma_t sigma_t].^2)

Sigma_xx = F*Sigma_LL*F'
sigma_v = sqrt(Sigma_xx)












