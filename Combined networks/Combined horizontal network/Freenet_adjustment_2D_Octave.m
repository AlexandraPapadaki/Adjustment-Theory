%--------------------------------------------------------------------------
%   
%   SELECTED SECTIONS OF ADJUSTMENT CALCULATION 
%   Template for non-linear functional models
% 
%   Author         : Georgios Malissiovas
%   Version        : June 7, 2016
%   Last changes   : June 7, 2016
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




%% Adjustment results for Task 1.1
%--------------------------------------------------------------------------
fprintf('\n\nTask 1.1 \n\n')
% Datum definition - Points 6 and 9 contribute to the datum
datum=[0 1 1 0]';
[X_hat, s_X, v, s_v, L_hat, s_L_hat, k]=freenet_2D(directions,distances,points,datum);

% Adjustment results and st devs
Est_par_task11 = [X_hat ];
Adj_res_task11 = [L_hat];
Lagrange_mult_task11 = k;

% Display the adjustment results
fprintf('Adjustment results of Task 1.1 \n')
num2str(Est_par_task11, '%.4f')
num2str(Adj_res_task11, '%.4f') 

% Estimated coordinates from Task 1.1
Y_task11 = X_hat(1:2:8);
X_task11 = X_hat(2:2:8);


%% Adjustment results for Task 1.2
%--------------------------------------------------------------------------
fprintf('\n\n\n\n\nTask 1.2 \n\n')
% Datum definition - Points 1 and 6 contribute to the datum
datum=[1 1 0 0]';
[X_hat, s_X, v, s_v, L_hat, s_L_hat, k]=freenet_2D(directions,distances,points,datum);

% Adjustment results and st devs
Est_par_task12 = [X_hat ];
Adj_res_task12 = [L_hat];
Lagrange_mult_task12 = k;

% Display the adjustment results
fprintf('Adjustment results of Task 1.2 \n')
num2str(Est_par_task12, '%.4f')
num2str(Adj_res_task12, '%.4f') 

% Estimated coordinates from Task 1.2
Y_task12 = X_hat(1:2:8);
X_task12 = X_hat(2:2:8);


%% Adjustment results for Task 1.3
%--------------------------------------------------------------------------
fprintf('\n\n\n\n\nTask 1.3 \n\n')
% Part a
% Datum definition - Total trace minimization
[X_hat, s_X, v, s_v, L_hat, s_L_hat]=freenet_2D_pseudoinv(directions,distances,points);

% Adjustment results and st devs
Est_par_task13a = [X_hat ];
Adj_res_task13a = [L_hat ];

% Display the adjustment results
fprintf('Adjustment results of Task 1.3a \n')
num2str(Est_par_task13a, '%.4f')
num2str(Adj_res_task13a, '%.4f') 



% Part b
% Datum definition - All points contribute to the datum 
datum=[1 1 1 1]';
[X_hat, s_X, v, s_v, L_hat, s_L_hat, k]=freenet_2D(directions,distances,points,datum);

% Adjustment results and st devs
Est_par_task13b = [X_hat ];
Adj_res_task13b = [L_hat];
Lagrange_mult_task13b = k;

% Display the adjustment results
fprintf('Adjustment results of Task 1.3b \n')
num2str(Est_par_task13b, '%.4f')
num2str(Adj_res_task13b, '%.4f') 

% Estimated coordinates from Task 1.3
Y_task13 = X_hat(1:2:8);
X_task13 = X_hat(2:2:8);
Lagrange_mult_task13 = k;




%% Plot for Task 1.4
%--------------------------------------------------------------------------
figure('Name',sprintf('Task 1.4'));
hold all;
plot([X_task11;X_task11(1);X_task11(3)],[Y_task11;Y_task11(1);Y_task11(3)],'ro-')
plot([X_task12;X_task12(1);X_task12(3)],[Y_task12;Y_task12(1);Y_task12(3)],'b*--')
plot([X_task13;X_task13(1);X_task13(3)],[Y_task13;Y_task13(1);Y_task13(3)],'m+-')
legend(('Task1.1'),('Task1.2'),('Task1.3'));
grid on
xlabel('Easting'), ylabel('Northing')
title('\it{Comparison of different datum definitions}','FontSize',12)
hold off;




%% Task 2 - Calculate and interpret the Lagrange multipliers
%--------------------------------------------------------------------------
% Display the Lagrange multipliers
fprintf('\n\n\n\n\nTask 2 \n\n')
fprintf('\nLagrance multipliers from Task 1.1\n');disp(Lagrange_mult_task11);
fprintf('\nLagrance multipliers from Task 1.2\n');disp(Lagrange_mult_task12);    
fprintf('\nLagrance multipliers from Task1.3_beforeb\n');disp(Lagrange_mult_task13);




%% Task 3 - S transformation
%--------------------------------------------------------------------------
fprintf('\n\n\n\n\nTask 3 \n\n')

X_hat = Est_par_task13b(:,1);
X_hat(9:end) = X_hat(9:end)*pi/200;

new_datum = [1 1 0 0]';

[X_hat_transformed] = S_Transformation(X_hat, new_datum, points);

% Check for the S-transformation (Difference between trasnformed and estimated from task 1.2)
X_hat_transformed
Estimated_param = Est_par_task12(:,1)
Check_diff = X_hat_transformed-Est_par_task12(:,1)





