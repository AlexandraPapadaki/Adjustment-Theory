%--------------------------------------------------------------------------
%   
%   Plane fitting
% 
%   Author         : Sven Weisbrich, Alexandra Papadaki
%   Version        : July 05, 2013
%   Last changes   : July 12, 2016
%
%--------------------------------------------------------------------------

clc;
clear all;
close all;

%--------------------------------------------------------------------------
%   Observations and initial values for unknowns
%--------------------------------------------------------------------------

data=load('plane.txt');

x=data(:,2);
y=data(:,3);
z=data(:,4);

%values for unknowns
 %=> complete
nx = 1;
ny = 1;
nz = 1;
d = 1;

n_v= [nx;ny;nz];

X_0=[n_v;d];

%Initial values for the residuals
vx= zeros(length(x),1);
vy= zeros(length(y),1);
vz= zeros(length(z),1);

v=[vx;vy;vz];

%Number of condition equations
no_n=length(data);

%Number of unknowns
no_u=length(X_0);

%Redundancy
r= no_n - no_u +1; %plus 1 restriction(nx^2+ny^2+nz^2=1)

%--------------------------------------------------------------------------
%  stochastic model
%--------------------------------------------------------------------------
%Cofactor matrix
Q_ll= eye(3*no_n, 3*no_n);
%Weight matrix
P = Q_ll^-1;

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
%break-off condition
epsilon=10^-5;
max_psi=10^100;

%Number of iterations
iteration=0;

%Iteration
while max_psi>epsilon

    %Condition equations Psi_i
    Psi= nx.*(x+vx)+ny.*(y+vy)+nz.*(z+vz)-d;

    %Designmatrix A
    A=[x+vx, y+vy, z+vz, -ones((no_n),1)]; 
   
    %Designmatrix B
    B1 = nx*eye(no_n, no_n);
    B2 = ny*eye(no_n, no_n);
    B3 = nz*eye(no_n, no_n);
    B = [B1, B2, B3];
    
    %Designmatrix C
    C= [2*nx, 2*ny, 2*nz, 0];
    
    %Vector of misclosures
    c1 = zeros(no_n,1);
    c2 = 1;
    w1= -B*v+Psi-c1;
    w2= nx^2+ny^2+nz^2-c2;
    
    %Normal matrix
    N= [-A'*(B*Q_ll*B')^-1*A C';C 0];

    %Vector of the right hand side
    n= [A'*(B*Q_ll*B')^-1*w1;  -w2];

    %Inversion of normal matrix
    N_inv=inv(N);
    
    %Solution of normal equation
    x_hat=N_inv*n;
    
    %Adjusted unknowns
    X_hat= X_0 + x_hat(1:end-1); %end-1 because end=k2
    
    %Update of the unknowns
    X_0=X_hat;
    
	nx = X_0(1);
    ny = X_0(2);
    nz = X_0(3);
    d = X_0(4);
    
    %Lagrangian Multipliers
    k1= (B*Q_ll*B')^-1*(-A*x_hat(1:end-1)-w1);
    k2= x_hat(end);
    
    %Residuals
    v= Q_ll*B'*k1;
    
    %Update of the residuals
    vx= v(1:no_n);
    vy= v(no_n+1:2*no_n);
    vz= v(2*no_n+1:3*no_n);
    
    %Check 1
    max_psi=max(abs(Psi));
    
    %Update number of iterations
    iteration=iteration+1;

end

Q_xx= -N_inv(1:no_u,1:no_u);
 
%Vector of adjusted observations
L_hat= [x;y;z]+v;

%Empirical reference standard deviation
s_0=sqrt(v'*inv(Q_ll)*v/r);

%VC matrix of adjusted unknowns
S_XX_hat=s_0^2*Q_xx;

%Standard deviation of the adjusted unknows
s_X=sqrt(diag(S_XX_hat));

%Cofactor matrix of the residuals
Q_vv= Q_ll*B'*(B*Q_ll*B')^-1*B*Q_ll;

%VC matrix of residuals
S_vv=s_0^2*Q_vv;

%Standard deviation of the residuals
s_v=sqrt(diag(S_vv));

%Cofactor matrix of adjusted observations
Q_LL_hat=Q_ll-Q_vv;

%VC matrix of adjusted observations
S_LL_hat=s_0^2*Q_LL_hat;

%Standard deviation of the adjusted observations
s_L_hat=sqrt(diag(S_LL_hat));

%Plot the adjusted plane 
figure()
hold on
grid on

%plot 3D measured points
ply = [y;y(1,1)];
plx = [x;x(1,1)];
plz = [z;z(1,1)];
plot3(plx,ply,plz, 'r*')

points = [x, y, z];

%Create 3D coords of the plane
[x_plane, y_plane] = ndgrid(9.5:0.5:19.5,50:0.5:60);
z_plane = (-nx*x_plane-ny*y_plane+d)/nz;

%plot surface
surf(x_plane, y_plane, z_plane)

title('Plane fitting')
xlabel('X [m]')
ylabel('y [m]')
zlabel('z [m]')
hold off

