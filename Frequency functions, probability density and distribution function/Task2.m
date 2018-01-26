% Task 2
clc;
clear all;
close all;

level=load('levelling.txt')

%Vector of differences
deltaH=level(:,1)-level(:,2)
n=length(level)

%Standard deviation of a single observation
s_l=sqrt(deltaH'*deltaH/(2*n))

%Standard deviation for the arithmetic mean
s_l_mean=s_l/sqrt(2)

fclose all;