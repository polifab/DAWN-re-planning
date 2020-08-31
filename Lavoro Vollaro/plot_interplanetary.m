clc
clear all
close all

%Sun gravitational parameter
global mu
mu = 1.327*(10^11); %[km^3/s^2]

%Earth elements and sv in the eliocentric system
% with respect to the sun
[Earth_coe0, Earth_r0, Earth_v0, Earth_julianday] = planet_elements_and_sv(3,2007,9,27,0,0,0);

%Mars elements and sv in the eliocentric system
% with respect to the sun
[Mars_coe0, Mars_r0, Mars_v0, Mars_julianday] = planet_elements_and_sv(4,2007,9,27,0,0,0);

[Earth_coef, Earth_rf, Earth_vf, Earth_juliandayf] = planet_elements_and_sv(3,2009,2,17,0,0,0);
[Mars_coef, Mars_rf, Mars_vf, Mars_juliandayf] = planet_elements_and_sv(4,2009,2,17,0,0,0);

Earth_park(Earth_r0)
orbit = orbitAttempt;

alpha = 0:pi/200:2*pi;
phi = -pi/2:pi/100:pi/2;

plot_orbit(Earth_coe0(3),Earth_coe0(4),Earth_coe0(7),Earth_coe0(7)*sqrt(1-Earth_coe0(2)^2),Rotz(Earth_coe0(3))*Rotx(Earth_coe0(4))*[0;sqrt(Earth_coe0(7)^2-(Earth_coe0(7)*sqrt(1-Earth_coe0(2)^2))^2);0])
xlabel('x')
ylabel('y')
zlabel('z')

plot_orbit(Mars_coe0(3),Mars_coe0(4),Mars_coe0(7),Mars_coe0(7)*sqrt(1-Mars_coe0(2)^2),Rotz(Earth_coe0(3))*Rotx(Earth_coe0(4))*[0;sqrt(Earth_coe0(7)^2-(Earth_coe0(7)*sqrt(1-Earth_coe0(2)^2))^2);0])

plot3(Earth_r0(1),Earth_r0(2),Earth_r0(3),'bo')
plot3(Earth_rf(1),Earth_rf(2),Earth_rf(3),'bx')
plot3(Mars_r0(1),Mars_r0(2),Mars_r0(3),'ro')
plot3(Mars_rf(1),Mars_rf(2),Mars_rf(3),'rx')
