%--------------------------------------------------------------------------
%   
%   Singular matrices
% 
%   Author         : Alexandra Papadaki
%   Version        : May 25, 2016
%
%--------------------------------------------------------------------------
clc;
clear all;
close all;
%--------------------------------------------------------------------------
% Calculate the rank of the design matrix A
%--------------------------------------------------------------------------

A=[-9 -10 -7 -1 -8;
   -3 5 -5 8 3; 
   8 0 8 -8 0;
   -16 0 -10 16 6;
   3 8 0 5 5;
   -13 2 -7 15 8;
   4 2 10 -2 8;
   15 8 4 -7 -3;
   2 6 0 4 4;
   7 2 -1 -5 -6];

r =rank(A);
rank_def_A=size(A,2)-r;

%--------------------------------------------------------------------------
% Calculate the normal matrix N
%--------------------------------------------------------------------------
N=A'*A;
deter_N=det(N)
% rank of N
r_N=rank(N);
rank_def_N=size(N,2)-r_N;
%--------------------------------------------------------------------------
% Calculate the 3 different g-inverses of the matrix N
%--------------------------------------------------------------------------
% First g-inverse
B1=randi([-500 500], size(N,2),rank_def_N);
N_ext_1=[N B1;B1' zeros(rank_def_N)];
N_ext_1_inv=N_ext_1^(-1);
Q_11=N_ext_1_inv(1:size(N,2),1:size(N,2));
N1=N*Q_11*N;
t1=trace(Q_11)
 
% Second g-inverse
B2=randi([-500 500], size(N,2),rank_def_N);
N_ext_2=[N B2;B2' zeros(rank_def_N)];
N_ext_2_inv=N_ext_2^(-1);
Q_11=N_ext_2_inv(1:size(N,2),1:size(N,2));
N2=N*Q_11*N;
t2=trace(Q_11)

% Third g-inverse
B3=randi([-500 500], size(N,2),rank_def_N);
N_ext_3=[N B3;B3' zeros(rank_def_N)];
 N_ext_3_inv=N_ext_3^(-1);
 Q_11=N_ext_3_inv(1:size(N,2),1:size(N,2));
 N3=N*Q_11*N;
 t3=trace(Q_11)
 
%--------------------------------------------------------------------------
% Calculate the eigenvalues and eigenvectors of N
%--------------------------------------------------------------------------
[V U] = eig(N);

G = [V(:,1) V(:,2)];
N_ext_g=[N G; G' zeros(rank_def_N)];
N_extr_g_inv=N_ext_g^(-1);

Q_11=N_extr_g_inv(1:size(N,2),1:size(N,2));

tg=trace(Q_11) % it is the minimum

%--------------------------------------------------------------------------
% Calculate the pseudo inverse of N
%--------------------------------------------------------------------------
Np=pinv(N); % pinv does the above, adding eigenvectors and inverting them
tp=trace(Np) % it is the same with tg



