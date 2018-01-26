%--------------------------------------------------------------------------
%   
%   Observation and adjustment of a combined horizontal network
%   Free net adjustment
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

%Datum [1 6 9 15]
datum=[1 1 1 1]';%datum defining points

%--------------------------------------------------------------------------
%   Observations and initial values for unknowns
%--------------------------------------------------------------------------
L=[distances(:,3); directions(:,3)*pi/200];

% Initial values for the unknowns 
X_0 = reshape([points(:,3) points(:,2)]',8,1); 
X_0 = [X_0; zeros(4,1)];

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
S_LL = diag([0.001^2*ones(5,1); (0.003*pi/200)^2*ones(9,1)]);

%Theoretical standard deviation
sigma_0=1;

%Cofactor matrix of the observations
Q_LL=1/sigma_0^2*S_LL;

%Weight matrix
P=inv(Q_LL);

%Centroid
yc = sum(points(:,2).*datum)/sum(datum);
xc = sum(points(:,3).*datum)/sum(datum);

%Reduced coordinates
y_dash=points(:,2)-yc;
x_dash=points(:,3)-xc;

%Constraint
B=repmat(eye(2),1,4);
B=[B; reshape([y_dash -x_dash]',1,8)];

%Delete entries which are not contributing to the datum
D = (ones(3,1)*reshape((datum*ones(1,2))',1,8));
B=B.*D;

%including the omegas
B=[B zeros(3,4)];

%break-off condition
epsilon=10^-5;
delta=10^-10;
max_x_hat=Inf;

%Number of iterations
iteration=0;

%Iteration
while max_x_hat>epsilon
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
    syms x1 y1 x6 y6 x9 y9 x15 y15 w1 w6 w9 w15
    A = Jacobian_computation(x1, y1, x6, y6, x9, y9, x15, y15, w1, w6, w9, w15, X_0);
    
     %Normal matrix
    N=A'*P*A;
    
    %Vector of absolute values
    n=A'*P*l;
    
    %Extension
    Next=[N B'; B zeros(3)];
    next=[n;zeros(3,1)];
    
    %Inversion of normal matrix / Cofactor matrix of the unknowns
    Q = inv(Next);
    Q_xx=Q(1:end-3,1:end-3);
    
    %Solution of normal equation
    x_hat_ext = Q*next;
    x_hat = x_hat_ext(1:end-3);
    k = x_hat_ext(end-2:end);
    
    %Adjusted unknowns
    X_hat=X_0+x_hat;
    
    %Update
    X_0=X_hat;
    
    %Check 1
    max_x_hat=max(abs(x_hat));
    
    %Update number of iterations
    iteration=iteration+1;
end

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

T_chi2=((s_0^2)*r)/sigma_0^2;
trshld=chi2inv(0.95,r);
if T_chi2<trshld
        disp('all good')
else disp('eliminate Blunder')
end

%% Transforming angles from rad to gon and correct negative angles
% Transform angles from rad to gon
X_hat(9:end)=X_hat(9:end)*200/pi;
s_X(9:end)=s_X(9:end)*200/pi;

%% Internal reliability
EV=100*diag(Q_vv*P);  %...redundancy number ROUND THESE NUMBERS INTO INTEGER IN THE REPORT
NV=abs(v)./(sigma_0*sqrt(diag(Q_vv))); %standardised residuals
GF=-v./(diag(Q_vv*P));
GRZW=sigma_0*4.13./(sqrt(diag(Q_vv*P).*diag(P)));

%External reliability
EGK=(1-diag(Q_vv*P)).*GRZW;
EP=(1-diag(Q_vv*P)).*GF;

L_cor =[L(1:5);L(6:end)*200/pi] ; 
v_cor =[v(1:5);v(6:end)*200/pi] ;
L_hat_cor = [L_hat(1:5);L_hat(6:end)*200/pi] ;
s_v_cor = [s_v(1:5);s_v(6:end)*200/pi] ;
s_L_hat_cor = [s_L_hat(1:5);s_L_hat(6:end)*200/pi] ; 

GF_cor = [GF(1:5);GF(6:end)*200/pi] ; 
GRZW_cor = [GRZW(1:5);GRZW(6:end)*200/pi] ; 
EGK_cor = [EGK(1:5);EGK(6:end)*200/pi] ; 
EP_cor = [EP(1:5);EP(6:end)*200/pi] ; 

mat = [L_cor v_cor L_hat_cor s_v_cor s_L_hat_cor EV NV GF_cor GRZW_cor EGK_cor EP_cor];
mat = [X_hat s_X];
nv = latex(mat, '%0.5f');

%plot the adjusted network
figure; hold all;
plot([y1; y6],[x1; x6], 'color','k');
plot([y1; y15],[x1; x15],  'color','k');
plot([y6; y9], [x6; x9],  'color','k');
plot([y9; y15],[x9; x15],  'color','k');
plot([y1; y9], [x1; x9], 'color','k');

% calculate confidence ellipse and error ellipse
for i=0:3
    j= 1+i*2; 
    qxx = Q_xx(j,j);
    qyy = Q_xx(j+1,j+1);
    qxy = Q_xx(j,i+j);
    qyx = Q_xx(j+1,j);

    sx_i(i+1) = s_0*sqrt(qxx);      
    sy_i(i+1) = s_0*sqrt(qyy);          
    s_pn(i+1) = sqrt(sx_i(i+1)^2+sy_i(i+1)^2);
    w(i+1)  = sqrt((qxx-qyy)^2 + 4*(qxy^2));
    Af(i+1) = sqrt(0.5*(s_0^2)*(qxx+qyy+w(i+1)));
    Bf(i+1) = sqrt(0.5*(s_0^2)*(qxx+qyy-w(i+1)));
    phi(i+1) = 0.5*atan2( (2*qxy),(qxx - qyy) );
    
    Ak(i+1) = sqrt(2*finv(0.95,2,r)*Af(i+1)^2);
    Bk(i+1) = sqrt(2*finv(0.95,2,r)*Bf(i+1)^2);
    
    err_ell=ellipse(Af((i+1))*8000, Bf((i+1))*8000, X_0(j+1), X_0(j),phi((i+1)));
    con_ell=ellipse(Ak((i+1))*8000, Bk((i+1))*8000, X_0(j+1), X_0(j),phi((i+1)));
    
    % plot the ellipses
    plot(err_ell(:,1),err_ell(:,2),'-', 'color','r');
    plot(con_ell(:,1),con_ell(:,2),'-', 'color','c');
end

title('the adjusted coordinates of the network and their ellipse');
xlabel('Easting[m]');
ylabel('Northing[m]');
hold off

figure;
bar(v),title('Residuals')
xlabel('n')
ylabel('residuals(m)')