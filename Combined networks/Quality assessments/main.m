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
results=freenet_2D(directions,distances,points,datum);










