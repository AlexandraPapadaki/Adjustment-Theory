function [X_hat, s_X, v, s_v, L_hat, s_L_hat]=freenet_2D_pseudoinv(directions,distances,points)

%--------------------------------------------------------------------------
%   Observations and initial values for unknowns
%--------------------------------------------------------------------------
L=[distances(:,3); directions(:,3)*pi/200];

% Initial values for the unknowns 
X_0 = reshape([points(:,3) points(:,2)]',8,1); 
X_0=[X_0; zeros(4,1)];

%Number of observations
no_n=length(L);

%Number of unknowns
no_u=length(X_0);

%Number of constraints
no_b = 3;

%Redundancy
r=no_n-no_u+no_b;

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
% Stochastic model
% VC Matrix of the observations
S_LL = diag([0.1^2*ones(5,1); (0.001*pi/200)^2*ones(9,1)]);

%Theoretical standard deviation
sigma_0=1;

%Cofactor matrix of the observations
Q_LL=1/sigma_0^2*S_LL;

%Weight matrix
P=inv(Q_LL);

%break-off condition
epsilon=10^-5;
delta=10^-10;
max_x_hat=Inf;

%Number of iterations
iteration=0;

%Iteration
while max_x_hat>epsilon

    % Approx. Unknowns
    x1=X_0(1);
    y1=X_0(2);
    x6=X_0(3);
    y6=X_0(4);
    x9=X_0(5);
    y9=X_0(6);
    x15=X_0(7);
    y15=X_0(8);
    w1=X_0(9);
    w6=X_0(10);
    w9=X_0(11);
    w15=X_0(12);
    %Vector of reduced observations
    L_0_dist =[sqrt((x1-x6)^2+(y1-y6)^2);
            sqrt((x1-x9)^2+(y1-y9)^2);
            sqrt((x6-x9)^2+(y6-y9)^2);
            sqrt((x1-x15)^2+(y1-y15)^2);
            sqrt((x9-x15)^2+(y9-y15)^2)];
            
    L_0_dir = [atan((y6-y1)/(x6-x1))-w1;
            atan((y15-y1)/(x15-x1))-w1+pi;
            atan((y1-y6)/(x1-x6))-w6+pi;
            atan((y9-y6)/(x9-x6))-w6;
            atan((y15-y9)/(x15-x9))-w9+pi;
            atan((y1-y9)/(x1-x9))-w9+pi;
            atan((y6-y9)/(x6-x9))-w9+pi;
            atan((y1-y15)/(x1-x15))-w15+2*pi;
            atan((y9-y15)/(x9-x15))-w15];
    L_0 = [L_0_dist; L_0_dir];

    l=L-L_0;

    %Designmatrix
    A = compute_Jacobian(X_0);
    
    %Normal matrix
    N=A'*P*A;

    %Vector of absolute values
    n=A'*P*l;

    % Total trace minimization - Pseudoinverse
    Q_xx = pinv(N);
    
    %Solution of normal equation
    x_hat=Q_xx*n;
    
    %Adjusted unknowns
    X_hat=X_0+x_hat;
    
    %Update
    X_0=X_hat;
    
    %Check 1
    max_x_hat=max(abs(x_hat(1:end-3)));
    
    %Update number of iterations
    iteration=iteration+1;

end

%Vector of residuals
v=A*x_hat-l;

% Sum of squared residuals
vTPv=v'*P*v;

%Vector of adjusted observations
L_hat=L+v;

%Empirical reference standard deviation
s_0=sqrt(vTPv/r);

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

%% Transforming angles from rad to gon and correct negative angles

% Transform angles from rad to gon
X_hat(end-3:end)=X_hat(end-3:end)*200/pi;
s_X(end-3:end)=s_X(end-3:end)*200/pi;

v(6:end) = v(6:end)*200/pi; 
s_v(6:end) = s_v(6:end)*200/pi; 
L_hat(6:end) = L_hat(6:end)*200/pi; 
s_L_hat(6:end) = s_L_hat(6:end)*200/pi; 

% Correct negative angles
X_hat(end-3)=X_hat(end-3)+400;
X_hat(end-2)=X_hat(end-2)+400;
X_hat(end)=X_hat(end)+400;