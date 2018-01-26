%--------------------------------------------------------------------------
%   
%   Fitting a straight line in 2D
% 
%   Author         : Georgios Malissiovas, Alexandra Papadaki
%   Version        : July 06, 2016
%   Last changes   : July 06, 2016
%
%--------------------------------------------------------------------------

clc;
clear all;
close all;

%--------------------------------------------------------------------------
%   Observations and initial values for unknowns
%--------------------------------------------------------------------------

% Measured 2D coordinates
y=[0.1; 1.1; 1.8; 2.4];
x=[1.0; 2.0; 3.0; 4.0];
L = [y;x];

%Initial values for the unknowns
a=1;
b=2;

X_0=[a;b];

%Number of condition equations
no_n=length(x);

%Number of unknowns
no_u=length(X_0);

%Redundancy
r=no_n-no_u;

%Initial values for the residuals
v_y=zeros(no_n, 1);
v_x=zeros(no_n, 1);

v=[v_y;v_x];

%--------------------------------------------------------------------------
%  stochastic model
%--------------------------------------------------------------------------
%Weight matrix
Sigma_ll= diag([0.02^2;0.01^2;0.04^2;0.02^2;0.02^2;0.01^2;0.04^2;0.02^2]);
sigma_0 = 1;
Q_ll = 1/sigma_0^2 * Sigma_ll;
P = Q_ll^-1;

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
%break-off condition
delta=10^-15;
epsilon =10^-15 ;

max_X_hat=10^100;
max_Psi = 10^100;

%Number of iterations
iteration=0;

%Iteration
while max_X_hat>delta || max_Psi>epsilon

    %Condition equations Psi_i
    Psi = y+v_y-(a*(x+v_x))-b;

    %Designmatrix A
    A = [-(x+v_x) -ones(no_n,1)];
    
    %Designmatrix B
    B1 = eye(no_n,no_n); %ws pros vy
    B2 = -a*eye(no_n,no_n); %ws pros vx
    B=[B1 B2];
    
    %Vector of misclosures
    w = -B*v+Psi;
    
    %Normal matrix
    N= [B*Q_ll*B' A;A' zeros(no_u,no_u)];

    %Vector of the right hand side
    n= [-w;zeros(size(A',1),1)];
    
    %Solution of normal equation
    x_hat=N\n;
    
    %Adjusted unknowns
    X_hat=x_hat(end-1:end);
    
    %Lagrange multipliers
    k = x_hat(1:end-2);
    
    %Update of the unknowns
    X_0=X_hat+X_0;
    
    a=X_0(1);
    b=X_0(2);
    
    %Residuals
    v=Q_ll*B'*k;
    
    %Update of the residuals
    v_y= v(1:no_n,1);
    v_x= v(no_n+1:end,1);
    
    %Check 1
    max_X_hat = max(abs(X_hat));
    %Check 2
    max_Psi = max(abs(Psi));
    
    %Update number of iterations
    iteration=iteration+1;

end

%VC matrix of adjusted unknowns
Q22 = -(A'*(B*Q_ll*B')^-1*A)^-1;
Q12 = -(B*Q_ll*B')\(A*Q22);
Q21 = Q12';
Q11 = (B*Q_ll*B')\(eye(no_n)-A*Q21);

Q_XX = -Q22;
 
%Vector of adjusted observations
L_hat= L+v;

%Empirical reference standard deviation 
s_0= sqrt(k'*B*Q_ll*B'*k/r); %same as s_02 but written in a way appropriate for GHM. Better because in GHM we dont have to calc design matrix

s_02= sqrt(v'*P*v/r);s_0^2

%VC matrix of adjusted unknowns
S_XX_hat=s_0^2*Q_XX;

%Standard deviation of the adjusted unknowns
s_X=sqrt(diag(S_XX_hat));

%Cofactor matrix of the residuals
Q_vv= Q_ll*B'*Q11*B*Q_ll;

%VC matrix of residuals
S_vv=s_0^2*Q_vv;

%Standard deviation of the residuals
s_v=sqrt(diag(S_vv));

%Cofactor matrix of adjusted observations
Q_LL_hat= Q_ll-Q_vv;

%VC matrix of adjusted observations
S_LL_hat= s_0^2*Q_LL_hat;

%Standard deviation of the adjusted observations
s_L_hat= sqrt(diag(S_LL_hat));

% Plot the adjusted line
figure()
hold on 
grid on
plot(x,y,'r*')
plot(L_hat(5:end),L_hat(1:4))
xlabel('x - direction (m)')
ylabel('y - direction (m)')
axis equal
hold off

