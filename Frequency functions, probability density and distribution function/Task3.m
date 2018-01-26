% Task 3
clc;
clear all;
close all;
format long;

sigma1=sqrt(0.01)
sigma2=sqrt(0.09)
sigma12=0.0135;
r12=sigma12/(sigma1*sigma2)
r21=r12


sigma1=0.01;
sigma2=0.025;
r12=0.85;
r21=r12;
VCM1=[[sigma1^2;r21*sigma2*sigma1] [r12*sigma2*sigma1;sigma2^2]]

sigma1=0.001*pi/200
sigma2=0.025;
r12=-0.5;
r21=r12;
VCM2=[[sigma1^2;r21*sigma2*sigma1] [r12*sigma2*sigma1;sigma2^2]]

fclose all;