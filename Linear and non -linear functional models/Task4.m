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
data1=load('Task4_Fixed.txt');
data2=load('Task4_obs.txt');
x(:,1)=data1;
y(:,1)=data2;

%Vector of Observations
% L =>complete
L = y;

% Initial values for unknowns
% X_0 =>complete
a=3.1749
b=2.2
X_0 = [a; b];

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
% sigma_0 =>complete

%Cofactor matrix of the observations
% Q_LL =>complete

%Weight matrix
% P =>complete
P = eye(no_n);

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
%break-off conditions
% epsilon =>complete
% delta =>complete
%max_x_hat =>complete 
%max_Phi =>complete 
epsilon = 10^-7;
delta = 10^-12;
max_x_hat = Inf;
max_Phi = Inf;

%Number of iterations
iteration=0;

%Iteration
while max_x_hat>epsilon || max_Phi>delta   % if sth is wrong and the iteration is endless...look at the design matrix for mistake
    plot(x,L,'.')
    hold on
    plot(x,b*sqrt(1-x.^2/a^2),'r')
%     hold on
%     plot(x,-b*sqrt(1-x.^2/a^2),'r')
    hold off
    
    %Vector of reduced observations
%     L_0 =>complete
%     l =>complete

L_0(1,1)= b*sqrt(1-x(1,1)^2/a^2);
L_0(2,1)= -b*sqrt(1-x(2,1)^2/a^2);
L_0(3,1)= -b*sqrt(1-x(3,1)^2/a^2);
L_0(4,1)= b*sqrt(1-x(4,1)^2/a^2);
L_0(5,1)=- b*sqrt(1-x(5,1)^2/a^2);

l = L-L_0;

    %Designmatrix
%     A =>complete
syms a0 b0
J = jacobian([b0*sqrt(1-x(1,1)^2/a0^2),-b0*sqrt(1-x(2,1)^2/a0^2),-b0*sqrt(1-x(3,1)^2/a0^2),b0*sqrt(1-x(4,1)^2/a0^2),-b0*sqrt(1-x(5,1)^2/a0^2)],[a0, b0]);
A = subs(J,[a0,b0],[a,b]);

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
a=X_hat(1,1);
b=X_hat(2,1);

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
    
Phi(1,1) = b*sqrt(1-x(1,1)^2/a^2);
Phi(2,1) = -b*sqrt(1-x(2,1)^2/a^2);
Phi(3,1) = -b*sqrt(1-x(3,1)^2/a^2);
Phi(4,1) = b*sqrt(1-x(4,1)^2/a^2);
Phi(5,1) = -b*sqrt(1-x(5,1)^2/a^2);

max_Phi = max(abs(L_hat-Phi));

    %Update number of iterations
    iteration=iteration+1;

end

figure
bar(v)

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




    
