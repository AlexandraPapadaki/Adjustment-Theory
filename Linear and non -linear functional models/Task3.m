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

%Vector of Observations
L = [1.16; 15.15];
ro = 8.93;

%Initial values for unknowns
V = L(2)/ro;
X_0=V;

%Number of observations
no_n =length(L);

%Number of unknowns
no_u =length(X_0);

%Redundancy
r =no_n-no_u;

%--------------------------------------------------------------------------
%  stochastic model
%--------------------------------------------------------------------------
%VC Matrix of the observations
S_LL = diag([0.005 0.05].^2);

%Theoretical standard deviation
sigma_0 =1;

%Cofactor matrix of the observations
Q_LL =1/sigma_0^2*S_LL;

%Weight matrix
P =inv(Q_LL);

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
%break-off conditions
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!! chose small values !!!!!!!!!!!!
epsilon = 10^-5;
delta = 10^-12;
max_x_hat = Inf;
max_Phi = Inf;

%Number of iterations
iteration=0;

%Iteration
while max_x_hat>epsilon || max_Phi>delta
    % Plot results of every iteration for check
    plot(ro,L,'.')
    hold on
    plot(ro,L(2),'.')
    hold off    
  
    %Vector of reduced observations
    %V = abs(V);
    L_0 =[V^(1/3); ro*V];
    l = L-L_0;

    %Designmatrix
    syms Vol
    a = Vol^(1/3);
    m = ro*Vol;
%     a = nthroot(Vol,3);
    J = jacobian([a m], Vol);          
    A = subs(J,Vol,V);
    
    %Normal matrix
    N =A'*P*A;
    
    %Vector of absolute values
    n =A'*P*l;

    %Inversion of normal matrix / Cofactor matrix of the unknowns
    Q_xx =inv(N);
    
    %Solution of normal equation
    x_hat =Q_xx*n;
    
    %Adjusted unknowns
    X_hat =X_0+x_hat;
    
    %Update
V=X_hat;
    X_0 =X_hat;
    
    %Check 1
    max_x_hat = max(abs(x_hat));

    %Vector of residuals
    v = A*x_hat-l;   %%%%%small L ...not 1

    %Vector of adjusted observations
    L_hat = L+v;

    %Check 2
    Phi(1,1) = V^(1/3);
    Phi(2,1) = ro*V;
    max_Phi=max(abs(L_hat-Phi));
    
    %Update number of iterations
    iteration=iteration+1
    
end
figure
plot(A);
bar(v)


% %Empirical reference standard deviation
% s_0 =
% %VC matrix of adjusted unknowns
% S_XX_hat =
% 
% %Standard deviation of the adjusted unknowns
% s_X =
% 
% %Cofactor matrix of adjusted observations
% Q_LL_hat =
% 
% %VC matrix of adjusted observations
% S_LL_hat =
% 
% %Standard deviation of the adjusted observations
% s_L_hat =
% 
% %Cofactor matrix of the residuals
% Q_vv =>complete
% 
% %VC matrix of residuals
% S_vv =>complete
% 
% %Standard deviation of the residuals
% s_v =>complete




    
