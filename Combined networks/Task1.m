%--------------------------------------------------------------------------
%   
%   Template for non-linear functional models
% 
%   Author : Sven Weisbrich, Alexandra Papapdaki
%
%--------------------------------------------------------------------------

clc;
clear all;
close all;

%--------------------------------------------------------------------------
%   Observations and initial values for unknowns
%--------------------------------------------------------------------------
%Vector of Observations
L =[206.9094; 
    46.5027; 
    84.6449; 
    115.5251; 
    155.5891];

%gon into rad
L = L*pi/200;

%Fixed Coordinates
y1= 682.415;
x1= 321.052;
y2= 203.526;
x2= 310.527;
y4= 251.992;
x4= 506.222;
y5= 420.028;
x5= 522.646;
y6= 594.553;
x6= 501.494;
y=[y1 y2 y4 y5 y6]';
x=[x1 x2 x4 x5 x6]';

%Fixed coordinates for the constraint
y7= 485.959;
x7= 219.089;

%Initial values for unknowns
x3=300;
y3=450;
w3=0;

X_0=[x3; y3; w3];

%Number of observations
no_n = length(L);

%Number of unknowns
no_u = length(X_0);

%Number of constraints
no_b = 1;

%Redundancy
r= no_n - no_u + no_b;

%--------------------------------------------------------------------------
%  stochastic model
%--------------------------------------------------------------------------
%VC Matrix of the observations
sigma_r = 0.001*pi/200;
S_LL=sigma_r^2*eye(no_n);

%Theoretical standard deviation
sigma_0=sigma_r;

%Cofactor matrix of the observations
Q_LL = 1/sigma_0^2*S_LL;

%Weight matrix
P = inv(Q_LL);
% P = eye(no_n);

%--------------------------------------------------------------------------
%  Constraints
%--------------------------------------------------------------------------
C= [2*x3-2*x7, 2*y3-2*y7];

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
%break-off conditions
epsilon=10^-5;
delta=10^-12;
max_x_hat=Inf;
check2=Inf;

%Number of iterations
iteration=0;

%Iteration
while max_x_hat>epsilon && check2>delta

    %Vector of reduced observations
    L_0=zeros(5,1);
        
    for i=1:5 
       L_0(i,1)= atan2((y(i,1)-y3),(x(i,1)-x3))-w3;
%          L_0(i,1) = direction (y(i,1), x(i,1), y3, x3, w3);
    end
  
    l = L-L_0;
    
    %Designmatrix
    A=zeros(5,3);
    for i=1:5
        A(i,1)=dr_dx_origin(x3,y3,x(i,1),y(i,1));
        A(i,2)=dr_dy_origin(x3,y3,x(i,1),y(i,1));  
        A(i,3)=-1;
    end
           
    %Normal matrix
    N=A'*P*A;
    N_ext =[N C'; C 0];

    %Vector of misclosure
    w= (x3-x7)^2+(y3-y7)^2-625;

    %Vector of absolute values
    n = A'*P*l;
    n_ext = [n; -w];

    %Inversion of normal matrix / Cofactor matrix of the unknowns
    Q_xx_ext =inv(N_ext);
    
    %Solution of normal equation
    x_hat =Q_xx_ext*n_ext;
    
    %Adjusted unknowns
    X_hat =X_0+x_hat(1:no_u);
    
    %Update
    x3=X_hat(1);
    y3=X_hat(2);
    w3=X_hat(3);
    
    X_0=X_hat
    
    %Check 1
    max_x_hat = max(abs(x_hat(1:no_u)));
    
    %Update number of iterations
    iteration=iteration+1;
    
    %Vector of residuals
    v = A*x_hat(1:no_u)-l

    %Vector of adjusted observations
    L_hat = L+v;
    
    %Check 2
    Psi=zeros(5,1);
       
    for i=1:5
        Psi(i,1)= atan2((y(i,1)-y3),(x(i,1)-x3))-w3;
%        Psi(i,1)= direction (y(i,1), x(i,1), y3, x3, w3);
    end
    
    check2=max(abs(L_hat-Psi));

end

figure
bar(v)

A'*P*v;
v'*P*v;
-l'*P*v;

v
L_hat
X_hat
iteration

%Cofactormatrix of the unknowns
Q_xx = Q_xx_ext(1:no_u,1:no_u);

%Empirical reference standard deviation
s_0 = sqrt(v'*P*v/r);

%VC matrix of adjusted unknowns
S_XX_hat =s_0^2*Q_xx;

%Standard deviation of the adjusted unknows
s_X =diag(sqrt(S_XX_hat));

%Cofactor matrix of adjusted observations
Q_LL_hat = A*Q_xx*A';

%VC matrix of adjusted observations
S_LL_hat = s_0^2*Q_LL_hat ;

%Standard deviation of the adjusted observations
s_L_hat = diag(sqrt(S_LL_hat));

%Cofactor matrix of the residuals
Q_vv = Q_LL - Q_LL_hat;

%VC matrix of residuals
S_vv = s_0^2*Q_vv;

%Standard deviation of the residuals
s_v = diag(sqrt(S_vv));



    