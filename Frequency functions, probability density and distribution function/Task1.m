% Task 1
clc
clear all
close all

%load distance measurements
l=load ('distances.txt');
%Number of measurements
n=length(l);

%Absolute frequency
bins = round (sqrt(n));
figure (1)
hist(l,bins);

%Name axes
xlabel ('Distance')
ylabel ('Number of observations')

hold on
[y,x]=hist(l,bins);
plot(x,y,'r');

%Relative Frequency
figure (2)
bar(x,y/n*100);

xlabel ('Distance')
ylabel ('Number of observations %')

hold on
plot(x,y/n*100,'r');

%Cumulative Frequency
figure (3)
y_cum=cumsum(y/n*100);
bar(x,y_cum)
xlabel('Distance')
ylabel('Number of observations %')

hold on
plot(x,y_cum,'r');

%Mean Value 
l_mean = sum(l)/n

%Variance for single observation
e=ones(n,1)
v=l_mean*e-l;
sl_2=v'*v/(n-1)

%Standard deviation for single observation
%sl2=sqrt(sl_2)
sl=std(l,0)

%Variance for arithmetic mean
sl_mean_2=sl_2/n

%Standar deviation for arithmetic mean
sl_mean=sl/sqrt(n)

%How often do you have to measure this distance 
%with the previously determined standard deviation 
%in order to obtain an accuracy for the arithmetic mean 
%of  s_l_mean?0.1 mm?
nmin=(sl_2/(0.0001)^2)
nmin2=ceil(sl_2/(0.0001)^2)


fclose all;








