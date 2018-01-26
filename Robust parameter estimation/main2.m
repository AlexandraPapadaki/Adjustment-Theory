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
c = 1e-5;
max_delta_v = 10e10;

%Number of iterations
iteration = 0;

%--------------------------------------------------------------------------
%   Observations
%--------------------------------------------------------------------------
data=load('StraightLine.txt');

x=data(:,1);
y=data(:,2);

%Number of observations
no_n=length(data);

%Initial values for residuals
v=ones(no_n,1);

%Number of unknowns
no_u=2;

%Redundancy
r=no_n-no_u;

%--------------------------------------------------------------------------
%  Initial stochastic model
%--------------------------------------------------------------------------
%Weight matrix
P= eye(no_n);

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------

% Designmatrix
A= [x ones(no_n,1)];

%Vector of observations
L= y;

while max_delta_v>epsilon
    
    v_approx=v;
    
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

% Plot of results
figure
hold on

x=0:0.1:10;

y_L1=X_hat(1)*x+X_hat(2);

% Plot of the robust estimator
plot(x,y_L1,'g')

%Result of the L2 estimator including blunders
X_hat_L2=(A'*A)\(A'*L);

y_L2=X_hat_L2(1)*x+X_hat_L2(2);

% Plot of the L2 estimator
plot(x,y_L2,'--r')

% legend('Robust Estimator','L2-Estimator')

%Plot of the points
plot(data(:,1),y,'o')
