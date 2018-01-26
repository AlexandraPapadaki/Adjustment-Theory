%--------------------------------------------------------------------------
%   
%   Parametric adjustment
% 
%   Author         :  Papadaki Alexandra, Bourou Stavroula
%
%--------------------------------------------------------------------------
clc;
clear all;
close all;
format long g;

%% Load the necessary data
directions=load('Directions.txt');
distances=load('Distances.txt');
points=load('Points.txt');

%--------------------------------------------------------------------------
% Observations and initial values for unknowns
%--------------------------------------------------------------------------
%Vector of Observations
L=[distances(:,3); directions(:,3)*pi/200];

% Initial values for the unknowns 
points_unknown=points(1:2,:);
X_0 = reshape([points_unknown(:,3) points_unknown(:,2)]',4,1); 
X_0 = [X_0; zeros(4,1)];

%Number of observations
no_n=length(L);

%Number of unknowns
no_u=length(X_0);

%Redundancy
r=no_n-no_u;

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
% Stochastic model
% VC Matrix of the observations
S_LL = diag([0.001^2*ones(5,1); (0.003*pi/200)^2*ones(9,1)]);

%Theoretical standard deviation
sigma_0=1;

%Cofactor matrix of the observations
Q_LL=1/sigma_0^2*S_LL;

%Weight matrix
P=inv(Q_LL);

%break-off condition
epsilon=10^-5;
delta=10^-10;
max_x_hat=10^100;

%Number of iterations
iteration=0;

%fixed points
x9=points(3,3);
y9=points(3,2);
x15=points(4,3);
y15=points(4,2);

%Iteration
while max_x_hat>epsilon 
% Approx. Unknowns
    x1=X_0(1);
    y1=X_0(2);
    x6=X_0(3);
    y6=X_0(4);
    
    w1=X_0(5);
    w6=X_0(6);
    w9=X_0(7);
    w15=X_0(8);
    
    %Vector of reduced observations
    L_0_dist =[sqrt((x1-x6)^2+(y1-y6)^2);
               sqrt((x1-x9)^2+(y1-y9)^2);
               sqrt((x6-x9)^2+(y6-y9)^2);
               sqrt((x1-x15)^2+(y1-y15)^2);
               sqrt((x9-x15)^2+(y9-y15)^2)];
           
    L_0_dir =[atan2((y6-y1),(x6-x1))-w1;
              atan2((y15-y1),(x15-x1))-w1;
              atan2((y1-y6),(x1-x6))-w6;
              atan2((y9-y6),(x9-x6))-w6; 
              atan2((y15-y9),(x15-x9))-w9; 
              atan2((y1-y9),(x1-x9))-w9; 
              atan2((y6-y9),(x6-x9))-w9;
              atan2((y1-y15),(x1-x15))-w15;
              atan2((y9-y15),(x9-x15))-w15];
          
    %Check for negative angles
    for i=1:9
        if L_0_dir(i)<0
            L_0_dir(i)=L_0_dir(i)+2*pi;
        end 
    end
    
     L_0 = [L_0_dist; L_0_dir];    
    
    l=L-L_0;
    
    %Designmatrix
    syms x1 y1 x6 y6 w1 w6 w9 w15
    A = Jacobian_computation_parametric(x1, y1, x6, y6, w1, w6, w9, w15, X_0,points);
    
    %Normal matrix
    N=A'*P*A;
    
    %Vector of absolute values
    n=A'*P*l;
    
    %Inversion of normal matrix / Cofactor matrix of the unknowns
    Q_xx=inv(N);
	
	%Solution of the normal equation system
    x_hat=Q_xx*n;
    
    %Update
    X_hat=X_0+x_hat;
    X_0=X_hat;
   
    %Check
    max_x_hat=max(abs(x_hat));
    
    %Update number of iterations
    iteration=iteration+1; 
end
x1=X_0(1);
y1=X_0(2);
x6=X_0(3);
y6=X_0(4);
    
w1=X_0(5);
w6=X_0(6);
w9=X_0(7);
w15=X_0(8);

% Vector of residuals
v=A*x_hat-l;
vTPv=v'*P*v;
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

X_hat(5:end)=X_hat(5:end)*200/pi;
X_hat(7,1)=X_hat(7,1)+400;
s_X(5:end)=s_X(5:end)*200/pi;
%% Internal reliability
EV=100*diag(Q_vv*P);  %...redundancy number ROUND THESE NUMBERS INTO INTEGER IN THE REPORT
NV=abs(v)./(sigma_0*sqrt(diag(Q_vv))); %standardised residuals
GF=-v./(diag(Q_vv*P));
GRZW=sigma_0*4.13./(sqrt(diag(Q_vv*P).*diag(P)));

%External reliability
EGK=(1-diag(Q_vv*P)).*GRZW;
EP=(1-diag(Q_vv*P)).*GF;

T_chi2=((s_0^2)*r)/sigma_0^2;
trshld=chi2inv(0.95,r);
if T_chi2<trshld
        disp('all good')
else disp('Problem')
end

L_cor =[L(1:5);L(6:end)*200/pi] ; 
v_cor =[v(1:5);v(6:end)*200/pi] ;
L_hat_cor = [L_hat(1:5);L_hat(6:end)*200/pi] ;
s_v_cor = [s_v(1:5);s_v(6:end)*200/pi] ;
s_L_hat_cor = [s_L_hat(1:5);s_L_hat(6:end)*200/pi] ; 

GF_cor = [GF(1:5);GF(6:end)*200/pi] ; 
GRZW_cor = [GRZW(1:5);GRZW(6:end)*200/pi] ; 
EGK_cor = [EGK(1:5);EGK(6:end)*200/pi] ; 
EP_cor = [EP(1:5);EP(6:end)*200/pi] ; 

% mat = [L_cor v_cor L_hat_cor s_v_cor s_L_hat_cor EV NV GF_cor GRZW_cor EGK_cor EP_cor]
mat = [X_hat s_X];
nv = latex(mat, '%0.5f');

figure;
bar(v),title('Residuals')
xlabel('n')
ylabel('residuals(m)')