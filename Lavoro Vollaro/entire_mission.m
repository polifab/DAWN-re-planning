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

Epark_radius = 200; %[km]
Epark_inclination = 0; %[rad]
Vpark_radius = 200; %[km]
Vpark_inclination = 0; %[rad]
Cpark_radius = 200; %[km]
Cpark_inclination = 0; %[rad]

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
      
[xx,yy,zz] = sphere(10);

%% Initial and final configurations
%Sun position for reference
% [Sun_coe0, Sun_r0, Sun_v0, Sun_julianday] = planet_elements_and_sv(3,2007,9,27,0,0,0);

%Earth elements and sv in the eliocentric system
% with respect to the sun
[Earth_coe0, Earth_r0, Earth_v0, Earth_julianday] = planet_elements_and_sv(3,2007,9,27,0,0,0);
[Earth_coe1, Earth_r1, Earth_v1, Earth_juliandayf] = planet_elements_and_sv(3,2009,2,17,0,0,0);

%Mars elements and sv in the eliocentric system
% with respect to the sun
[Mars_coe0, Mars_r0, Mars_v0, Mars_julianday] = planet_elements_and_sv(4,2007,9,27,0,0,0);
[Mars_coe1, Mars_r1, Mars_v1, Mars_juliandayf] = planet_elements_and_sv(4,2009,2,17,0,0,0);
[Mars_coe2, Mars_r2, Mars_v2, Mars_julianday2] = planet_elements_and_sv(4,2011,7,16,0,0,0);

%Vesta elements and sv in the eliocentric system
% with respect to the sun
[Vesta_coe1, Vesta_r1, Vesta_v1, Vesta_julianday0] = planet_elements_and_sv(10,2009,2,17,0,0,0);
[Vesta_coe2, Vesta_r2, Vesta_v2, Vesta_juliandayf] = planet_elements_and_sv(10,2011,7,16,0,0,0);
[Vesta_coe3, Vesta_r3, Vesta_v3, Vesta_julianday3] = planet_elements_and_sv(10,2012,9,5,0,0,0);
[Vesta_coe4, Vesta_r4, Vesta_v4, Vesta_julianday4] = planet_elements_and_sv(10,2015,3,6,0,0,0);

%Vesta elements and sv in the eliocentric system
% with respect to the sun
[Ceres_coe3, Ceres_r3, Ceres_v3, Ceres_julianday3] = planet_elements_and_sv(11,2012,9,5,0,0,0);
[Ceres_coe4, Ceres_r4, Ceres_v4, Ceres_julianday4] = planet_elements_and_sv(11,2015,3,6,0,0,0);

%% Entire mission plot
%Earth - Mars travel
figure()

hold on

%Interplanetary orbit
fprintf('\n\n EARTH TO MARS \n\n')
[body_pos1, sp_v1, body_posf1, spacecraft_vf1,tof1, orb_elem1] = ...
                gen_orbit(3,4,[2007 9 27 0 0 0],[2009 2 17 0 0 0]);
EM_orbit = intpl_orbit(tof1,Earth_r0,sp_v1);

%Planet orbits
plot_orbit(3,2007)
plot_orbit(4,2009)

%Planet positions
plot3(Earth_r0(1),Earth_r0(2),Earth_r0(3),'o','Color',colors(3))
plot3(Earth_r1(1),Earth_r1(2),Earth_r1(3),'x','Color',colors(3))
plot3(Mars_r0(1),Mars_r0(2),Mars_r0(3),'o','Color',colors(4))
plot3(Mars_r1(1),Mars_r1(2),Mars_r1(3),'x','Color',colors(4))

xlabel('x')
ylabel('y')
zlabel('z')
view(-10,45)
grid

%Mars - Vesta travel

%Interplanetary orbit
fprintf('\n\n MARS TO VESTA \n\n')
[body_pos2, sp_v2, body_posf2, spacecraft_vf2, tof2, orb_elem2] = ...
                gen_orbit(4,10,[2009 2 17 0 0 0],[2011 7 16 0 0 0]);
MV_orbit = intpl_orbit(tof2,Mars_r1,sp_v2);

plot_orbit(10,2011)

plot3(Mars_r2(1),Mars_r2(2),Mars_r2(3),'rd')
plot3(Vesta_r1(1),Vesta_r1(2),Vesta_r1(3),'x','Color',colors(10))
plot3(Vesta_r2(1),Vesta_r2(2),Vesta_r2(3),'*','Color',colors(10))

%Vesta - Ceres travel

%Interplanetary orbit
fprintf('\n\n VESTA TO CERES \n\n')
[body_pos3, sp_v3, body_posf3, spacecraft_vf3,tof3, orb_elem3] = ...
                gen_orbit(10,11,[2012 9 5 0 0 0],[2015 3 6 0 0 0]);
VC_orbit = intpl_orbit(tof3,Vesta_r3,sp_v3);

plot_orbit(11,2015)

plot3(Vesta_r3(1),Vesta_r3(2),Vesta_r3(3),'d','Color',colors(10))
plot3(Vesta_r4(1),Vesta_r4(2),Vesta_r4(3),'+','Color',colors(10))
plot3(Ceres_r3(1),Ceres_r3(2),Ceres_r4(3),'d','Color',colors(11))
plot3(Ceres_r4(1),Ceres_r4(2),Ceres_r4(3),'+','Color',colors(11))

%% Earth close-up
figure()

park_orbit(3,Earth_r0,Epark_radius,Epark_inclination)
hold on
track = [[EM_orbit(1,1);EM_orbit(2,1)] ...
    [EM_orbit(1,2);EM_orbit(2,2)] ...
    [EM_orbit(1,3);EM_orbit(2,3)]];
plot3(track(:,1),track(:,2),track(:,3),'k')

% hyperbola_attempt;
escape_hyp(3,4,EM_orbit(1:2,1:3),[2007 9 27 0 0 0],Epark_radius,0, orb_elem1);

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

% park_orbit(3,Earth_r0,Epark_radius,0)
hold on
surface(Earth_r0(1)+Earth_SOI*xx, Earth_r0(2)+Earth_SOI*yy,...
        Earth_r0(3)+Earth_SOI*zz,'FaceColor','none','EdgeColor','b')
plot3(track(:,1),track(:,2),track(:,3),'k')

hyperbola = escape_hyp(3,4,EM_orbit(1:2,1:3),[2007 9 27 0 0 0],...
                        Epark_radius,0, orb_elem1);

xlabel('x')
ylabel('y')
zlabel('z')
xlim([ 1.486*10^8, 1.508*10^8])
ylim([8*10^6, 10*10^6])
zlim([-10*10^5, 10*10^5])
view(-10,45)
grid

%% Mars close-up
% figure2()
% 
% hold on
% track = [[orbit(1,1);orbit(2,1)] ...
%     [orbit(1,2);orbit(2,2)] ...
%     [orbit(1,3);orbit(2,3)]];
% plot3(track(:,1),track(:,2),track(:,3),'k')
% 
% % hyperbola_attempt;
% escape_hyp(3,4,orbit(1:2,1:3),[2007 9 27 0 0 0],park_radius,0, orb_elem);
% 
% xlabel('x')
% ylabel('y')
% zlabel('z')
% xlim([ 1.49696*10^8, 1.49714*10^8])
% ylim([8.916*10^6, 8.932*10^6])
% zlim([-1*10^4, 10^4])
% view(-10,45)
% grid

%% Mars SOI close-up
% figure2()
% 
% % park_orbit(3,Earth_r0,park_radius,0)
% hold on
% surface(Earth_r0(1)+Earth_SOI*xx, Earth_r0(2)+Earth_SOI*yy,...
%         Earth_r0(3)+Earth_SOI*zz,'FaceColor','none','EdgeColor','b')
% plot3(track(:,1),track(:,2),track(:,3),'k')
% 
% hyperbola = escape_hyp(3,4,EM_orbit(1:2,1:3),[2007 9 27 0 0 0],...
%                         park_radius,0, orb_elem);
% 
% xlabel('x')
% ylabel('y')
% zlabel('z')
% xlim([ 1.486*10^8, 1.508*10^8])
% ylim([8*10^6, 10*10^6])
% zlim([-10*10^5, 10*10^5])
% view(-10,45)
% grid

%% Vesta close-up
figure2()

park_orbit(10,Vesta_r3,Vpark_radius,Vpark_inclination)
hold on
Vtrack = [[VC_orbit(1,1);VC_orbit(2,1)] ...
    [VC_orbit(1,2);VC_orbit(2,2)] ...
    [VC_orbit(1,3);VC_orbit(2,3)]];
plot3(Vtrack(:,1),Vtrack(:,2),Vtrack(:,3),'k')

escape_hyp(10,11,VC_orbit(1:2,1:3),[2012 9 5 0 0 0],Vpark_radius,...
                Vpark_inclination, orb_elem3);

xlabel('x')
ylabel('y')
zlabel('z')
xlim([-2.046542*10^8, -2.04653*10^8])
ylim([-2.49873*10^8, -2.498718*10^8])
zlim([3.22814*10^7, 3.22822*10^7])
view(-10,45)
grid

%% Ceres close-up
figure2()

park_orbit(11,Ceres_r4,Cpark_radius,Cpark_inclination)
hold on
Ctrack = [[VC_orbit(end-1,1);VC_orbit(end,1)] ...
    [VC_orbit(end,2);VC_orbit(end,2)] ...
    [VC_orbit(end-1,3);VC_orbit(end,3)]];
plot3(Ctrack(:,1),Ctrack(:,2),Ctrack(:,3),'k')

xlabel('x')
ylabel('y')
zlabel('z')
xlim([-3.579936*10^8, -3.579918*10^8])
ylim([1.15647*10^8, 1.15649*10^8])
zlim([6.93861*10^7, 6.93874*10^7])
view(-10,45)
grid

%% Gathering positions for a TBD animation
positions = [hyperbola;
             EM_orbit(2:end,1:3)];