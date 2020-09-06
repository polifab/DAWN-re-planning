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

colors = ["g"          %green
          "m"          %magenta
          "b"          %blue
          "r"          %red
          "#A2142F"    %darker red
          "#7E2F8E"    %purple
          "#4DBEEE"    %darker cyan
          "c"          %(bright) cyan
          "#D95319"    %orange
          "#77AC30"    %darker green
          "#EDB120"    %ochre
          "#D95319"];  %orange, not visible due to Sun orbit dimensions

%% Initial and final configurations
%Sun position for reference
% [Sun_coe0, Sun_r0, Sun_v0, Sun_julianday] = planet_elements_and_sv(3,2007,9,27,0,0,0);
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
[body_pos0, sp_v0, body_posf, spacecraft_vf,tof, orb_elem] = ...
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

park_orbit(3,Earth_r0,200,0)
hold on
track = [[orbit(1,1);orbit(2,1)] ...
    [orbit(1,2);orbit(2,2)] ...
    [orbit(1,3);orbit(2,3)]];
plot3(track(:,1),track(:,2),track(:,3),'k')

% hyperbola_attempt;
escape_hyp(3,4,orbit(1:2,1:3),[2007 9 27 0 0 0],200, orb_elem);

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

% park_orbit(3,Earth_r0,200,0)
hold on
surface(Earth_r0(1)+Earth_SOI*xx, Earth_r0(2)+Earth_SOI*yy,...
        Earth_r0(3)+Earth_SOI*zz,'FaceColor','none','EdgeColor','b')
plot3(track(:,1),track(:,2),track(:,3),'k')

hyperbola = escape_hyp(3,4,orbit(1:2,1:3),[2007 9 27 0 0 0],...
                        200, orb_elem);

xlabel('x')
ylabel('y')
zlabel('z')
xlim([ 1.486*10^8, 1.508*10^8])
ylim([8*10^6, 10*10^6])
zlim([-10*10^5, 10*10^5])
view(-10,45)
grid

%% Mars - Vesta travel
[Vesta_coe0, Vesta_r0, Vesta_v0, Vesta_julianday0] = planet_elements_and_sv(10,2009,2,17,0,0,0);
[Vesta_coef, Vesta_rf, Vesta_vf, Vesta_juliandayf] = planet_elements_and_sv(10,2011,7,16,0,0,0);
[Mars_coe2, Mars_r2, Mars_v2, Mars_julianday2] = planet_elements_and_sv(4,2011,7,16,0,0,0);
figure(1)

%Interplanetary orbit
[body_pos2, sp_v2, body_pos2, spacecraft_v2,tof2, orb_elem2] = ...
                gen_orbit(4,10,[2009 2 17 0 0 0],[2011 7 16 0 0 0]);
orbit = intpl_orbit(12,tof2,Mars_rf,sp_v2);

plot_orbit(10)

plot3(Mars_r2(1),Mars_r2(2),Mars_r2(3),'rd')
plot3(Vesta_r0(1),Vesta_r0(2),Vesta_r0(3),'o','Color',colors(10))
plot3(Vesta_rf(1),Vesta_rf(2),Vesta_rf(3),'x','Color',colors(10))

%% Vesta - Ceres travel
[Ceres_coe0, Ceres_r0, Ceres_v0, Ceres_julianday0] = planet_elements_and_sv(11,2012,9,5,0,0,0);
[Ceres_coef, Ceres_rf, Ceres_vf, Ceres_juliandayf] = planet_elements_and_sv(11,2015,3,6,0,0,0);
[Vesta_coe1, Vesta_r1, Vesta_v1, Vesta_julianday1] = planet_elements_and_sv(10,2012,9,5,0,0,0);
[Vesta_coe2, Vesta_r2, Vesta_v2, Vesta_julianday2] = planet_elements_and_sv(10,2015,3,6,0,0,0);
figure(1)

%Interplanetary orbit
[body_pos3, sp_v3, body_pos3, spacecraft_v3,tof3, orb_elem3] = ...
                gen_orbit(10,11,[2012 9 5 0 0 0],[2015 3 6 0 0 0]);
orbit = intpl_orbit(12,tof3,Vesta_r1,sp_v3);

plot_orbit(11)

plot3(Vesta_r1(1),Vesta_r1(2),Vesta_r1(3),'d','Color',colors(10))
plot3(Vesta_r2(1),Vesta_r2(2),Vesta_r2(3),'+','Color',colors(10))
plot3(Ceres_r0(1),Ceres_r0(2),Ceres_r0(3),'o','Color',colors(11))
plot3(Ceres_rf(1),Ceres_rf(2),Ceres_rf(3),'x','Color',colors(11))

%% Gathering positions for a TBD animation
positions = [hyperbola;
             orbit(2:end,1:3)];