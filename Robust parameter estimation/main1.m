%--------------------------------------------------------------------------
%   
%   Robust Parameter Estimation
% 
%   Author         : Sven Weisbrich, Alexandra Papadaki
%
%--------------------------------------------------------------------------

clc;
clear all;
close all;
format long g

% Break-off condition and initial settings
epsilon = 1e-9;
c = 1e-9;
max_delta_v = 10e10;

%Number of iterations
iteration   = 0;

%--------------------------------------------------------------------------
%   Observations
%--------------------------------------------------------------------------

L=load('testseries.txt');

%Number of observations
no_n=length(L);

%Initial values for residuals
v=ones(no_n,1);

%Number of unknowns
no_u=1;

%Redundancy
r=no_n-no_u;

%--------------------------------------------------------------------------
%  Initial stochastic model
%--------------------------------------------------------------------------
%Weight matrix
P= eye(no_n, no_n);

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------

% Designmatrix
A= ones(no_n, 1);

while max_delta_v>epsilon
    
    v_approx = v;
    
    %Normal matrix
    N=A'*P*A;

    %Vector of the right hand side
    n=A'*P*L;

    %Inversion of normal matrix
    Q_xx=inv(N);
    
    %Solution of normal equation
    X_hat=Q_xx*n;
      
    %Residuals
    v=A*X_hat-L;
    
    %Update of the stochastic model
    P= diag(1./(abs(v)+c));
    
    %Check
    max_delta_v= max(abs(v-v_approx));
    
    %Update number of iterations
    iteration=iteration+1;

end

%Vector of adjusted observations
L_hat=L+v;

%Empirical reference standard deviation
s_0=sqrt(v'*P*v/r);

%VC matrix of adjusted unknowns
S_XX_hat=s_0^2*Q_xx;

%Standard deviation of the adjusted unknows
s_X=sqrt(diag(S_XX_hat));



