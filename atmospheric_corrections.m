clear all
close all
clc

%% test for temperature determination
a = 1/273.16;
h = 0.4;
p = 1009;
t = load('air_temperature.txt');

n_t = size(t',1);

deltaD_t = zeros(n_t,1);
for i = 1:n_t
    x = 7.5*t(1,i)/(237.3+t(1,i))+0.7857;
    deltaD_t(i,1) = 283.04-( 0.29195*p/(1+a*t(1,i)) -4.126*10^(-4)*h/(1+a*t(1,i))*10^x );
end

figure;
plot(t , deltaD_t)
xlabel('Air Temperature [^oC]')
ylabel('Atmospheric correction[ppm]')

%% test for air pressure determination
a = 1/273.16;
h = 0.4;
p = load('air_pressure.txt');
t = 25;

x = 7.5*t/(237.3+t)+0.7857;

n_p = size(p',1);

deltaD_p = zeros(n_p,1);
for i = 1:n_p
    deltaD_p(i,1) = 283.04-( 0.29195*p(1,i)/(1+a*t) -4.126*10^(-4)*h/(1+a*t)*10^x );
end

figure;
plot(p , deltaD_p)
xlabel('Air Pressure [mbar]')
ylabel('Atmospheric correction[ppm]')

%% test for relative humidity determination
a = 1/273.16;
h = load('relative_humidity.txt');
p = 1009;
t = 25;

x = 7.5*t/(237.3+t)+0.7857;

h = h/100;

n_h = size(h',1);

deltaD_h = zeros(n_h,1);
for i = 1:n_h
    deltaD_h(i,1) = 283.04-( 0.29195*p/(1+a*t) -4.126*10^(-4)*h(1,i)/(1+a*t)*10^x );
end

figure;
plot(h*100 , deltaD_h)
xlabel('Relative Humidity [%]')
ylabel('Atmospheric correction[ppm]')