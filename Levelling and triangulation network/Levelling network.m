%--------------------------------------------------------------------------
%   
%   Template for linear functional models
% 
%   Author         : Sven Weisbrich, Alexandra Papadaki
%   created        : April 13, 2012
%   Last changes   : December 15, 2015
%
%--------------------------------------------------------------------------

clc;
clear all;
close all;

%--------------------------------------------------------------------------
%   Observations etc. 
%--------------------------------------------------------------------------
%Vector of Observations
L = [10.509; 5.360; -8.523; -7.348; -3.167; 15.881];

%Fixed
H_A=374.027;

%Benchmarks
benchmarks = [-H_A; 0; 0; H_A; 0; -H_A];

%% difference from other afdj problems...diferentiation of the obs
    L_dash = L - benchmarks;

    %Number of observations
    no_o = length(L);

%Number of unknowns
no_u = 3;

    %Redundancy
    r = no_o-no_u;

%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------
%Uncorrelated
%VC Matrix of the observations
S_LL = diag([0.006 0.004 0.005 0.003 0.004 0.012].^2);

%Theoretical standard deviation
sigma_0 = 1;

    %Cofactor matrix of the observations
    Q_LL = 1/sigma_0^2*S_LL;

    %Weight matrix
    P = inv(Q_LL);
    
% %Uncorrelated and equally weighted/unknown/same std
% P = eye(no_o);

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
%Designmatrix
A = [1 0 0; -1 1 0; 0 -1 1; 0 0 -1; -1 0 1; 0 1 0];

%
    %Normal matrix
    N = A'*P*A;

    %Vector of absolute values
    n= A'*P*L_dash;

    %Inversion of normal matrix / Cofactor matrix of the unknowns
    Q_xx = inv(N);

    %Adjusted unknowns
    X_hat = Q_xx*n;

    %Vector of residuals
    v= A*X_hat-L_dash;

    %Vector of adjusted observations
    L_dash_hat = L_dash+v;
    L_hat=L+v;

    %Final Check
    if (L_dash_hat-A*X_hat)<10^-15
        disp('everything fine!')
    else
        disp('We have a problem!')
    end

    %Empirical reference standard deviation
    s_0 = sqrt(v'*P*v/r);

    %VC matrix of adjusted unknowns
    S_XX_hat = s_0^2*Q_xx;

    %Standard deviation of the adjusted unknowns
    s_X_hat = sqrt(diag(S_XX_hat));

    %Cofactor matrix of adjusted observations
    Q_LL_hat = A*Q_xx*A';

    %VC matrix of adjusted observations
    S_LL_hat = s_0^2*Q_LL_hat;

    %Standard deviation of the adjusted observations
    s_L_hat = sqrt(diag(S_LL_hat));

% If Q_ll available 
    %Cofactor matrix of the residuals
    Q_vv = Q_LL - Q_LL_hat;

    %VC matrix of residuals
    S_vv = s_0^2*Q_vv;

    %Standard deviation of the residuals
    s_v  = sqrt(diag(S_vv));

%Compare the residuals
figure
bar(v)

%Compare teh adjusted unknowns
figure
bar(X_hat)

%Compare the observations
figure
bar(L_hat)   
