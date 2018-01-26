%--------------------------------------------------------------------------
%   
%   VCE
% 
%   Author         : Georgios Malissiovas, Alexandra Papadaki
%   Version        : June 29, 2016
%   Last changes   : June 29, 2016
%
%--------------------------------------------------------------------------

clc;
clear all;
close all;
format long g;

%% Load the necessary data
directions=load('Directions.txt');
distances=load('Distances.txt');
points=load('Points.txt');

%Initial values for VCE
alpha_dist=1;
alpha_dir=1;

s_0_dist = 10^10;
s_0_dir = 10^10;
delta = 10^-5;
while abs(1-s_0_dist)>delta || abs(1-s_0_dir)>delta
    
% Stochastic model
% ------------------------------------------------------------------------------
% VC Matrix of the observations
S_LL = diag([alpha_dist*0.1^2*ones(5,1);
            alpha_dir*(0.001*pi/200)^2*ones(9,1)]);

% Theoretical standard deviation
sigma_0 = 1;

% Cofactor and Weight matrix
Q_LL = 1/sigma_0^2*S_LL;
P = inv(Q_LL);

% ------------------------------------------------------------------------------



% Points 6 and 9 are fixed/control points
% ------------------------------------------------------------------------------
results=network_adjustment_2D(directions,distances,points, sigma_0, Q_LL, P);
% ------------------------------------------------------------------------------



% Variance Components estimation
% ------------------------------------------------------------------------------
% Redundancy Matrix
QvvP = results.Q_vv*P;

% Distances
v_dist = results.v(1:5);
r_dist = sum(diag(QvvP(1:5,1:5)));
P_dist = P(1:5,1:5);

s_0_dist = sqrt((v_dist'*P_dist*v_dist)/r_dist)
alpha_dist = s_0_dist^2*alpha_dist;

% Directions
v_dir = results.v(6:end).*pi/200;
r_dir = sum(diag(QvvP(6:end,6:end)));
P_dir = P(6:end,6:end);

s_0_dir = sqrt((v_dir'*P_dir*v_dir)/r_dir)
alpha_dir = s_0_dir^2*alpha_dir;
% ------------------------------------------------------------------------------

end
% Stoxos einai : s_0_dist+s_0_dir = 1
%Run with variance component and without to decide which is best

