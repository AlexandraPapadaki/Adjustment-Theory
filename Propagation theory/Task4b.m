%Task 4b
clear all;
close all;
clc;
% Calculation for the third position
a1 = 35.1550*pi()/200; %rad
a2 = 55.1200*pi()/200; %rad
s1 = 20.005; %m
s2 = 30.001; %m
t1 = 9.7; %s
t2 = 23.1; %s
t3 = 30; %s
x1 = s1*sin(a1);
x2 = s2*sin(a2);
y1 = s1*cos(a1);
y2 = s2*cos(a2);
v=sqrt((x2-x1)^2+(y2-y1)^2)/(t2-t1);

azim23 = direction(x2,y2,x1,y1,0)
s23 = v*(t3-t2)
x3=x2+s23*sin(azim23)
y3=y2+s23*cos(azim23)



%% standar deviations of the observations
sigma_a = 0.001*pi/200;
sigma_s = 0.001;
sigma_t = 0.1;

%% F1
syms a11 a22 s11 s22 t11 t22 t33
J1 = jacobian([s11*sin(a11) s22*sin(a22) s11*cos(a11) s22*cos(a22) t11 t22 t33], [a11, a22, s11, s22, t11, t22, t33])
F1 = subs(J1, [a11, a22, s11, s22, t11, t22, t33], [a1, a2, s1, s2, t1, t2, t3]);

%% F2
syms x11 x22 y11 y22 t11 t22 t33
J2 = jacobian([atan((x22-x11)/(y22-y11)) x22 y22 sqrt((x22-x11)^2+(y22-y11)^2)/(t22-t11) t22 t33], [x11, x22, y11, y22, t11, t22, t33])
F2= subs(J2, [x11, x22, y11, y22, t11, t22, t33], [x1, x2, y1, y2, t1, t2, t3]);

%% F3
syms azim2323 x22 y22 vv t22 t33
J3 = jacobian([azim2323 vv*(t33-t22) x22 y22], [azim2323, x22, y22, vv, t22, t33])
F3= subs(J3, [azim2323, x22, y22, vv, t22, t33], [azim23, x2, y2, v, t2, t3]);
% F3= [[ 1, 0, 0,       0,  0, 0]
% [ 0, 0, 0, t3 - t2, -v, v]
% [ 0, 1, 0,       0,  0, 0]
% [ 0, 0, 1,       0,  0, 0]];

%% F4
syms azim2323 s2323 x22 y22 
J4 = jacobian([x22+s2323*sin(azim2323) y22+s2323*cos(azim2323)], [azim2323, s2323, x22, y22])
F4= subs(J4, [azim2323, s2323, x22, y22], [azim23, s23, x2, y2]);
% F4 = [[  s23*cos(azim23), sin(azim23), 1, 0]
% [ -s23*sin(azim23), cos(azim23), 0, 1]];

F = F4*F3*F2*F1;

% %VC Matrix of the observations
Sigma_LL = diag([sigma_a sigma_a sigma_s sigma_s sigma_t sigma_t sigma_t].^2);

Sigma_xx = F*Sigma_LL*F';

sigma_x3 = sqrt(Sigma_xx(1,1))
sigma_y3 = sqrt(Sigma_xx(2,2))
