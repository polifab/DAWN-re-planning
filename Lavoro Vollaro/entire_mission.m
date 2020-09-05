clc
clear all
close all
%% Data

%Sun gravitational parameter
global mu
mu = 1.327*(10^11); %[km^3/s^2]

%SOI
Earth_SOI = 9.24*10^5; %[km]
Mars_SOI = 5.74*10^5; %[km]

%% Initial and final configurations
%Sun position for reference
[Sun_coe0, Sun_r0, Sun_v0, Sun_julianday] = planet_elements_and_sv(3,2007,9,27,0,0,0);
%Earth elements and sv in the eliocentric system
% with respect to the sun
[Earth_coe0, Earth_r0, Earth_v0, Earth_julianday] = planet_elements_and_sv(3,2007,9,27,0,0,0);

%Mars elements and sv in the eliocentric system
% with respect to the sun
[Mars_coe0, Mars_r0, Mars_v0, Mars_julianday] = planet_elements_and_sv(4,2007,9,27,0,0,0);

[Earth_coef, Earth_rf, Earth_vf, Earth_juliandayf] = planet_elements_and_sv(3,2009,2,17,0,0,0);
[Mars_coef, Mars_rf, Mars_vf, Mars_juliandayf] = planet_elements_and_sv(4,2009,2,17,0,0,0);

%% Entire mission plot
[xx,yy,zz] = sphere(10);

figure2()

%Parking orbit around Earth
% park_orbit(3,Earth_r0,200)

hold on

%Interplanetary orbit
[body_pos0, sp_v0, body_posf, spacecraft_vf,tof] = ...
                gen_orbit(3,4,[2007 9 27 0 0 0],[2009 2 17 0 0 0]);
orbit = intpl_orbit(12,tof,Earth_r0,sp_v0);

%Planet orbits
plot_orbit(3)
plot_orbit(4)

% %SOIs
% surface(Earth_r0(1)+Earth_SOI*xx, Earth_r0(2)+Earth_SOI*yy,...
%         Earth_r0(3)+Earth_SOI*zz,'FaceColor','none','EdgeColor','b')
% surface(Mars_rf(1)+Mars_SOI*xx, Mars_rf(2)+Mars_SOI*yy,...
%         Mars_rf(3)+Mars_SOI*zz,'FaceColor','none','EdgeColor','b')

%Planet positions
plot3(Earth_r0(1),Earth_r0(2),Earth_r0(3),'bo')
plot3(Earth_rf(1),Earth_rf(2),Earth_rf(3),'bx')
plot3(Mars_r0(1),Mars_r0(2),Mars_r0(3),'ro')
plot3(Mars_rf(1),Mars_rf(2),Mars_rf(3),'rx')

xlabel('x')
ylabel('y')
zlabel('z')
view(-10,45)
grid

%% Earth close-up
figure2()

park_orbit(3,Earth_r0,200)
hold on
track = [[orbit(1,1);orbit(2,1)/1000] ...
    [orbit(1,2);orbit(2,2)/1000] ...
    [orbit(1,3);orbit(2,3)]];
plot3(track(:,1),track(:,2),track(:,3),'k')
hyperbola_attempt;
xlabel('x')
ylabel('y')
zlabel('z')
xlim([ 1.49696*10^8, 1.49714*10^8])
ylim([8.916*10^6, 8.932*10^6])
zlim([-1*10^4, 10^4])
view(-10,45)
grid

%% Earth SOI close-up
figure2()

park_orbit(3,Earth_r0,200)
hold on
surface(Earth_r0(1)+Earth_SOI*xx, Earth_r0(2)+Earth_SOI*yy,...
        Earth_r0(3)+Earth_SOI*zz,'FaceColor','none','EdgeColor','b')
plot3(track(:,1),track(:,2),track(:,3),'k')



xlabel('x')
ylabel('y')
zlabel('z')
xlim([ 1.486*10^8, 1.508*10^8])
ylim([8*10^6, 10*10^6])
zlim([-10*10^5, 10*10^5])
view(-10,45)
grid