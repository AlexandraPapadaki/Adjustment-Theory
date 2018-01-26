%--------------------------------------------------------------------------
%   
%   Template for non-linear functional models
% 
%   Author         : Sven Weisbrich, Alexandra Papadaki
%
%--------------------------------------------------------------------------

clc;
clear all;
close all;

%--------------------------------------------------------------------------
%   Observations and initial values for unknowns
%--------------------------------------------------------------------------
%Vector of Observations
L = [6049; 4736.83; 5446.49];

%Fixed Coordinates
x1 = 4527.15;
y1 = 865.4;
x2 = 2047.25;
y2 = 2432.55;
x3 = 27.15;
y3 = 2865.22;

%Initial values for unknowns
x0 = 3000;
y0 = 4400;
X_0 = [x0; y0];

    %Number of observations
    no_o = length(L);

    %Number of unknowns
    no_u = length(X_0);

    %Redundancy
    r = no_o-no_u;

%--------------------------------------------------------------------------
%  stochastic model
%--------------------------------------------------------------------------
%Uncorrelated
%VC Matrix of the observations
%accuracy m+ppm
m=1 %mm
ppm=2
    s_LL = m*10^(-3)+ppm*L*10^(-6);
    S_LL = diag(s_LL.^2); %teleia

    %Theoretical standard deviation
    sigma_0 = 1;

    %Cofactor matrix of the observations
    Q_LL = 1/sigma_0^2*S_LL;

    %Weight matrix
    P = inv(Q_LL);

% %Uncorrelated and equally weighted/unknown std
% P = eye(no_o);

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
    %break-off conditions
    epsilon=10^-5;
    delta=10^-12;
    max_x_hat=Inf;
    max_Phi=Inf;

    %Number of iterations
    iteration=0;

%Iteration
while max_x_hat>epsilon ||  max_Phi>delta
   
    %Vector of reduced observations
L_0(1)= sqrt((x1-x0)^2 + (y1-y0)^2);
L_0(2)= sqrt((x2-x0)^2 + (y2-y0)^2);
L_0(3)= sqrt((x3-x0)^2 + (y3-y0)^2);

    l = L-L_0';
    
    %Design matrix equations from L_0
syms x00 y00
J = jacobian([sqrt((x1-x00)^2 + (y1-y00)^2) sqrt((x2-x00)^2 + (y2-y00)^2) sqrt((x3-x00)^2 + (y3-y00)^2)], [x00 y00])
A = subs(J, [x00 y00], [x0 y0]);   

% A(1,1) = -(x1-x0)/L(1);
% A(1,2) = -(y1-y0)/L(1);
% 
% A(2,1) = -(x2-x0)/L(2);
% A(2,2) = -(y2-y0)/L(2);
% 
% A(3,1) = -(x3-x0)/L(3);
% A(3,2) = -(y3-y0)/L(3);

    %Normal matrix
    N = A'*P*A;

    %Vector of absolute values
    n = A'*P*l;

    %Inversion of normal matrix / Cofactor matrix of the unknowns
    Q_xx = inv(N);
    
    %Solution of normal equation
    x_hat = Q_xx*n;
    
    %Adjusted unknowns
    X_hat = X_0 + x_hat;  
    
    %Update
x0=X_hat(1);
y0=X_hat(2);

    %X_0 =>complete
    X_0 = X_hat;  

    %Check 1
    max_x_hat=max(abs(x_hat));
    
    %Update number of iterations
    iteration=iteration+1;

    %Vector of residuals
    v = A*x_hat-l;

    %Vector of adjusted observations
    L_hat = L+v;

    %Check 2
Phi(1) = sqrt((x1-x0)^2 + (y1-y0)^2);
Phi(2) = sqrt((x2-x0)^2 + (y2-y0)^2);
Phi(3) = sqrt((x3-x0)^2 + (y3-y0)^2);

    max_Phi = max(abs(L_hat-Phi'));
end
figure
plot(A(1:length(A)));
bar(v)


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

