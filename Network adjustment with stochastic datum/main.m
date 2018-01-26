%--------------------------------------------------------------------------
%   
%   ADJUSTMENT CALCULATION II - Exercise 6
% 
%   Author         : Sven Weisbrich
%   Version        : June 11, 2013
%   Last changes   : June 12, 2014
%
%--------------------------------------------------------------------------

clc;
clear all;
close all;

% Break-off conditions
epsilon=1e-5;
max_x_hat=10^10;
% Check for the linearization error
delta=10^-12;
check2=Inf;
%delta!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! graps to (=0)

%Loading the observations
distances=load('Distances.txt');
directions=load('Directions.txt');

%Observationvector
L_N=[distances(:,3); directions(:,3)*pi/200];

%Loading the Coordinates
fixed=load('ControlPoints.txt');
new=load('NewPoints.txt');

%Initial values for the unknowns
x6=fixed(1,3);
y6=fixed(1,2);

x9=fixed(2,3);
y9=fixed(2,2);

x1=new(1,3);
y1=new(1,2);

x15=new(2,3);
y15=new(2,2);

w1=0;
w6=0;
w9=0;
w15=0;

X_0=[x1 y1 x15 y15 w1 w6 w9 w15 x6 y6 x9 y9]';

%Extension of the observations => observed unknowns
L_A = [x6 y6 x9 y9]';

%Observations
L = [L_N; L_A];

%stochastic model for the observations
S_LL_N = [0.1^2*ones(5,1); (0.001*pi/200)^2*ones(9,1)];

%stochastic model for the observed unknowns
S_LL_A = 0.01^2*ones(4,1);

%Stochastic model 
S_LL = diag([S_LL_N; S_LL_A]);

sigma_0=0.001;
Q_LL=1/sigma_0^2*S_LL;
P=inv(Q_LL);

%Redundancy
r = 18-12;

while max_x_hat>epsilon || check2>delta
	
	%Designmatrices
	A_N = compute_A_N(X_0);
    A_A = compute_A_A(X_0);
    
	%Design matrix A
	A = [A_N A_A; zeros(4,8) eye(4)];
	
	%reduced vector of observations
    L_0=[];
    L_0(1) = sqrt((x1-x6)^2+(y1-y6)^2);
    L_0(2) = sqrt((x1-x9)^2+(y1-y9)^2);
    L_0(3) = sqrt((x9-x6)^2+(y9-y6)^2);
    L_0(4) = sqrt((x1-x15)^2+(y1-y15)^2);
    L_0(5) = sqrt((x9-x15)^2+(y9-y15)^2);
    
    L_0(6) = atan2(y6-y1,x6-x1)-w1;
    L_0(7) = atan2(y15-y1,x15-x1)-w1;
    
    L_0(8) = atan2(y1-y6,x1-x6)-w6;
    L_0(9) = atan2(y9-y6,x9-x6)-w6;
    
    L_0(10) = atan2(y15-y9,x15-x9)-w9;
    L_0(11) = atan2(y1-y9,x1-x9)-w9;
    L_0(12) = atan2(y6-y9,x6-x9)-w9;
    
    L_0(13) = atan2(y1-y15,x1-x15)-w15;
    L_0(14) = atan2(y9-y15,x9-x15)-w15;
    
    %Check for negative angles
    for i=1:9
        if L_0(i+5)<0
            L_0(i+5)=L_0(i+5)+2*pi;
        end 
    end
    
	%Extension of the observations => observed unknowns
	L_0 = [L_0 x6 y6 x9 y9];
	
    l=L-L_0';

	%Normal matrix
    N=(A'*P*A);
    Q_xx=inv(N);
	
	%Solution of the normal equation system
    x_hat=Q_xx*(A'*P*l)

    %Update
    X_hat=X_0+x_hat;
    X_0=X_hat;
    
    
    x1=X_0(1);
    y1=X_0(2);

    x15=X_0(3);
    y15=X_0(4);

    w1=X_0(5);
    w6=X_0(6);
    w9=X_0(7);
    w15=X_0(8);
	
	% Update for the observed unknowns
	x6 = X_0(9)
    y6 = X_0(10)
    x9 = X_0(11)
    y9 = X_0(12)
    

	%Check
    max_x_hat=max(abs(x_hat));


end


% Vector of residuals
v=A*x_hat-l;

%Vector of adjusted observations
L_hat=L+v;

%Empirical reference standard deviation
s_0=sqrt(v'*P*v/r);

%VC matrix of adjusted unknowns
S_XX_hat=s_0^2*Q_xx;

%Standard deviation of the adjusted unknows
s_X=sqrt(diag(S_XX_hat));

%Cofactor matrix of adjusted observations
Q_LL_hat=A*Q_xx*A';

%VC matrix of adjusted observations
S_LL_hat=s_0^2*Q_LL_hat;

%Standard deviation of the adjusted observations
s_L_hat=sqrt(diag(S_LL_hat));

%Cofactor matrix of the residuals
Q_vv=Q_LL-Q_LL_hat;

%VC matrix of residuals
S_vv=s_0^2*Q_vv;

%Standard deviation of the residuals
s_v=sqrt(diag(S_vv));