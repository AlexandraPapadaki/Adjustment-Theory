%--------------------------------------------------------------------------
%   
%   Template for non-linear functional models
% 
%   Author         : Sven Weisbrich, Alexandra Papadaki
%   Version        : April 13, 2012
%   Last changes   : 2015
%
%--------------------------------------------------------------------------

clc;
clear all;
close all;

%--------------------------------------------------------------------------
%   Observations and initial values for unknowns
%--------------------------------------------------------------------------
data=load('Homework3-2.txt');
x=data(:,1);
y=data(:,2);

figure
plot(x,y,'.')
hold on

%Vector of Observations
% L =>complete
L = y;

%Initial values for unknowns
% X_0 =>complete
a = 2.2;
% a = 10^(log(L(1)/2)/x(1));
b = 2.3;
c = 0;
X_0 = [a; b; c];

%Number of observations
% no_n =>complete
no_n = length(L);

%Number of unknowns
% no_u =>complete
no_u = length(X_0);

%Redundancy
% r =>complete
r = no_n-no_u;
%--------------------------------------------------------------------------
%  stochastic model
%--------------------------------------------------------------------------
%VC Matrix of the observations
% S_LL =>complete

%Theoretical standard deviation

%Cofactor matrix of the observations

%Weight matrix
% P =>complete
P =  eye(no_n);

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
%break-off conditions
% epsilon =>complete
% delta =>complete
%max_x_hat =>complete 
%max_Phi =>complete 
epsilon = 10^-5;
delta = 10^-12;
max_x_hat = Inf;
max_Phi = Inf;

%Number of iterations
iteration=0;

% y=a*b^x+c
%Iteration
while max_x_hat>epsilon || max_Phi>delta   % if sth is wrong and the iteration is endless...look at the design matrix for mistake
    %Plots
%     plot(x,L,'.')
%     hold on
%     plot(x,a*b.^x+c,'r')
%     hold off
    
    %Vector of reduced observations
%     L_0 =>complete
%     l =>complete
L_0 = a*b.^x+c;
l = L-L_0; 

    %Designmatrix
%     A =>complete
syms a0 b0 c0
J = jacobian(a0*b0.^x+c0, [a0, b0, c0]);
A = subs(J,[a0,b0,c0],[a,b,c]);

    %Normal matrix
%     N =>complete
N = A'*P*A;

    %Vector of absolute values
%     n =>complete
n = A'*P*l;

    %Inversion of normal matrix / Cofactor matrix of the unknowns
%     Q_xx =>complete
Q_xx = inv(N);

    %Solution of normal equation
%     x_hat =>complete
x_hat = Q_xx*n;

    %Adjusted unknowns
%     X_hat =>complete
X_hat = X_0 + x_hat;  

    %Update
a=X_hat(1);
b=X_hat(2);
c=X_hat(3);

%     X_0 =>complete
X_0 = X_hat;  

    %Check 1
%     max_x_hat =>complete
max_x_hat = max(abs(x_hat));

    %Vector of residuals
    % v =>complete
v = A*x_hat-l;

    %Vector of adjusted observations
    % L_hat =>complete
L_hat = L+v;

    %Check 2
    %=>complete
Phi =a*b.^x+c;
max_Phi = max(abs(L_hat-Phi));

    %Update number of iterations
    iteration=iteration+1;

end
figure
plot(A);
bar(v);

figure
plot(x,L,'.')
hold on
plot(x,a*b.^x+c,'r')
hold off
%Empirical reference standard deviation
%s_0 =>complete

%VC matrix of adjusted unknowns
% S_XX_hat =>complete

%Standard deviation of the adjusted unknowns
% s_X =>complete

%Cofactor matrix of adjusted observations
% Q_LL_hat =>complete

%VC matrix of adjusted observations
% S_LL_hat =>complete

%Standard deviation of the adjusted observations
% s_L_hat =>complete

%Cofactor matrix of the residuals
% Q_vv =>complete

%VC matrix of residuals
% S_vv =>complete

%Standard deviation of the residuals
% s_v =>complete




    
