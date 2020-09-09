%% Preparations
clc
clear all
close all
%% Constants

%Sun gravitational parameter
global mu
mu = 1.327565122000000e+11; %[km^3/s^2]

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

%% Data

%Distances from the Sun
Ceres_to_Sun = 446145795; %[km]
Vesta_to_Sun = 353*10^6; %[km]

%Masses
Vesta_mass = 2.589*10^20; %[kg]
Ceres_mass = 8.958*10^20; %[kg]
Sun_mass = 1.989*10^30; %[kg]

%SOI: (m_planet/m_Sun)^(2/5) * distance_from_Sun
Earth_SOI = 9.24*10^5; %[km]
Mars_SOI = 5.74*10^5; %[km]
Ceres_SOI = (Ceres_mass/Sun_mass)^(2/5)*Ceres_to_Sun; %[km]
Vesta_SOI = (Vesta_mass/Sun_mass)^(2/5)*Vesta_to_Sun; %[km]

%Parking orbits
Epark_radius = 200; %[km]
Epark_inclination = 0; %[rad]
Vpark_radius = 210; %[km]
Vpark_inclination = 0; %[rad]
Cpark_radius = 365; %[km]
Cpark_inclination = 0; %[rad]

Vesta_hamo = 670; %[km]
Ceres_hamo = 1320; %[km]

%% Planets configurations
%{
    The mission is subdivided in three parts:
    - Earth to Mars;
    - Mars to Vesta;
    - Vesta to Ceres.
    Variables for each part have the same name but different suffix.
    
    Planet configurations have the following numerical indications:
    - 0: at Earth departure (to Mars) -> 27/9/07
    - 1: at Mars arrival/departure (from Earth, to Vesta) -> 17/2/09
    - 2: at Vesta arrival (from Mars) -> 16/7/11
    - 3: at Vesta departure (to Ceres) -> 5/9/12
    - 4: at Ceres arrival (from Vesta) -> 6/3/15
%}

%Earth
[Earth_coe0, Earth_r0, Earth_v0, ~] =...
                        planet_elements_and_sv(3,2007,9,27,0,0,0);
[Earth_coe1, Earth_r1, Earth_v1, ~] =...
                        planet_elements_and_sv(3,2009,2,17,0,0,0);

%Mars
[Mars_coe0, Mars_r0, Mars_v0, ~] =...
                        planet_elements_and_sv(4,2007,9,27,0,0,0);
[Mars_coe1, Mars_r1, Mars_v1, ~] =...
                        planet_elements_and_sv(4,2009,2,17,0,0,0);
[Mars_coe2, Mars_r2, Mars_v2, ~] =...
                        planet_elements_and_sv(4,2011,7,16,0,0,0);

%Vesta
[Vesta_coe1, Vesta_r1, Vesta_v1, ~] =...
                        planet_elements_and_sv(10,2009,2,17,0,0,0);
[Vesta_coe2, Vesta_r2, Vesta_v2, ~] =...
                        planet_elements_and_sv(10,2011,7,16,0,0,0);
[Vesta_coe3, Vesta_r3, Vesta_v3, ~] =...
                        planet_elements_and_sv(10,2012,9,5,0,0,0);
[Vesta_coe4, Vesta_r4, Vesta_v4, ~] =...
                        planet_elements_and_sv(10,2015,3,6,0,0,0);

%Ceres
[Ceres_coe3, Ceres_r3, Ceres_v3, ~] =...
                        planet_elements_and_sv(11,2012,9,5,0,0,0);
[Ceres_coe4, Ceres_r4, Ceres_v4, ~] =...
                        planet_elements_and_sv(11,2015,3,6,0,0,0);

%% Interplanetary plot
%% Earth - Mars travel
if exist('figure2') == 0  %#ok<*EXIST>
    figure()
else
    figure2()
end

hold on

%Interplanetary orbit
fprintf('\n\n EARTH TO MARS \n\n')
[body_pos1, sp_v1, body_posf1, sp_vf1,tof1, orb_elem1] = ...
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
%% Mars - Vesta travel

%Interplanetary orbit
fprintf('\n\n MARS TO VESTA \n\n')
[body_pos2, sp_v2, body_posf2, sp_vf2, tof2, orb_elem2] = ...
                gen_orbit(4,10,[2009 2 17 0 0 0],[2011 7 16 0 0 0]);
MV_orbit = intpl_orbit(tof2,Mars_r1,sp_v2);

%Planet orbits
plot_orbit(10,2011)

%Planet positions
plot3(Mars_r2(1),Mars_r2(2),Mars_r2(3),'rd')
plot3(Vesta_r1(1),Vesta_r1(2),Vesta_r1(3),'x','Color',colors(10))
plot3(Vesta_r2(1),Vesta_r2(2),Vesta_r2(3),'*','Color',colors(10))
%% Vesta - Ceres travel

%Interplanetary orbit
fprintf('\n\n VESTA TO CERES \n\n')
[body_pos3, sp_v3, body_posf3, sp_vf3,tof3, orb_elem3] = ...
                gen_orbit(10,11,[2012 9 5 0 0 0],[2015 3 6 0 0 0]);
VC_orbit = intpl_orbit(tof3,Vesta_r3,sp_v3);

%Planet orbits
plot_orbit(11,2015)

%Planet positions
plot3(Vesta_r3(1),Vesta_r3(2),Vesta_r3(3),'d','Color',colors(10))
plot3(Vesta_r4(1),Vesta_r4(2),Vesta_r4(3),'+','Color',colors(10))
plot3(Ceres_r3(1),Ceres_r3(2),Ceres_r4(3),'d','Color',colors(11))
plot3(Ceres_r4(1),Ceres_r4(2),Ceres_r4(3),'+','Color',colors(11))

title("Interplanetary orbits")

%% Earth close-up
if exist('figure2') == 0
    figure()
else
    figure2()
end

park_orbit(3,Earth_r0,Epark_radius,orb_elem1(4));
park_orbit(3,Earth_r0,Epark_radius,Epark_inclination);
hold on

%           DEBUG
% track1 = [[EM_orbit(1,1);EM_orbit(2,1)] ...
%     [EM_orbit(1,2);EM_orbit(2,2)] ...
%     [EM_orbit(1,3);EM_orbit(2,3)];...
%     [EM_orbit(end-1,1);EM_orbit(end,1)] ...
%     [EM_orbit(end-1,2);EM_orbit(end,2)] ...
%     [EM_orbit(end-1,3);EM_orbit(end,3)]];
% plot3(track1(1:2,1),track1(1:2,2),track1(1:2,3),'k')

escape_hyp(3,EM_orbit(1:2,1:3),[2007 9 27 0 0 0],Epark_radius,...
                                0, orb_elem1, norm(sp_vf1));

xlabel('x')
ylabel('y')
zlabel('z')
xlim([ 1.49696*10^8, 1.49714*10^8])
ylim([8.916*10^6, 8.932*10^6])
zlim([-1*10^4, 10^4])
view(-10,45)
grid
title("Earth close-up")

%% Earth SOI close-up
if exist('figure2') == 0
    figure()
else
    figure2()
end

hold on
surface(Earth_r0(1)+Earth_SOI*xx, Earth_r0(2)+Earth_SOI*yy,...
        Earth_r0(3)+Earth_SOI*zz,'FaceColor','none','EdgeColor',colors(3))

%           DEBUG
% plot3(track1(1:2,1),track1(1:2,2),track1(1:2,3),'k')

hyperbola = escape_hyp(3,EM_orbit(1:2,1:3),[2007 9 27 0 0 0],...
                        Epark_radius,0, orb_elem1, norm(sp_vf1));

xlabel('x')
ylabel('y')
zlabel('z')
xlim([ 1.486*10^8, 1.508*10^8])
ylim([8*10^6, 10*10^6])
zlim([-10*10^5, 10*10^5])
view(-10,45)
grid
title("Earth SOI close-up")

%% Mars close-up
if exist('figure2') == 0
    figure()
else
    figure2()
end

body_sphere(4,Mars_r1);
hold on

%           DEBUG
% track2 = [[MV_orbit(1,1);MV_orbit(2,1)] ...
%     [MV_orbit(1,2);MV_orbit(2,2)] ...
%     [MV_orbit(1,3);MV_orbit(2,3)];...
%     [MV_orbit(end-1,1);MV_orbit(end,1)] ...
%     [MV_orbit(end-1,2);MV_orbit(end,2)] ...
%     [MV_orbit(end-1,3);MV_orbit(end,3)]];
% plot3(track2(1:2,1),track2(1:2,2),track2(1:2,3),'k')

xlabel('x')
ylabel('y')
zlabel('z')
xlim([ 9.3478*10^7, 9.3488*10^7])
ylim([-1.8886*10^8, -1.88851*10^8])
zlim([-6.257*10^6, -6.248*10^6])
view(-10,45)
grid
title("Mars close-up")

%% Mars SOI close-up
if exist('figure2') == 0
    figure()
else
    figure2()
end

body_sphere(4,Mars_r1)
hold on
surface(Mars_r1(1)+Mars_SOI*xx, Mars_r1(2)+Mars_SOI*yy,...
        Mars_r1(3)+Mars_SOI*zz,'FaceColor','none','EdgeColor',colors(4))

%           DEBUG
% plot3(track1(3:4,1),track1(3:4,2),track1(3:4,3),'k')
% plot3(track2(:,1),track2(:,2),track2(:,3),'k')

xlabel('x')
ylabel('y')
zlabel('z')
xlim([ 9.28*10^7, 9.42*10^7])
ylim([-1.896*10^8, -1.88*10^8])
zlim([-7*10^6, -5.5*10^6])
view(-10,45)
grid
title("Mars SOI close-up")

%% Vesta close-up demo (arrival)
if exist('figure2') == 0
    figure()
else
    figure2()
end

park_orbit(10,Vesta_r2,Vesta_hamo,orb_elem2(4));
hold on

%           DEBUG
% track3 = [[VC_orbit(1,1);VC_orbit(2,1)] ...
%     [VC_orbit(1,2);VC_orbit(2,2)] ...
%     [VC_orbit(1,3);VC_orbit(2,3)];...
%     [VC_orbit(end-1,1);VC_orbit(end,1)] ...
%     [VC_orbit(end-1,2);VC_orbit(end,2)] ...
%     [VC_orbit(end-1,3);VC_orbit(end,3)]];
% plot3(track3(1:2,1),track3(1:2,2),track3(1:2,3),'k')

capture_hyp(10,MV_orbit(end-1:end,1:3),[2011 7 16 0 0 0],...
    Vesta_hamo,orb_elem2(4),orb_elem2,sp_vf2);

xlabel('x')
ylabel('y')
zlabel('z')
xlim([-2.045565*10^8, -2.04554*10^8])
ylim([-2.4995*10^8, -2.499475*10^8])
zlim([3.22723*10^7, 3.2273*10^7])
view(-10,45)
grid
title("Vesta close-up (arrival)")

%% Vesta SOI close-up (arrival)
if exist('figure2') == 0
    figure()
else
    figure2()
end

park_orbit(10,Vesta_r2,Vesta_hamo,orb_elem2(4));
hold on
surface(Vesta_r2(1)+Vesta_SOI*xx, Vesta_r2(2)+Vesta_SOI*yy,...
     Vesta_r2(3)+Vesta_SOI*zz,'FaceColor','none','EdgeColor',colors(10))

%           DEBUG
% plot3(track2(3:4,1),track2(3:4,2),track2(3:4,3),'k') 
% plot3(track3(1:2,1),track3(1:2,2),track3(1:2,3),'k')

capture_hyp(10,MV_orbit(end-1:end,1:3),[2011 7 16 0 0 0],...
                        Vesta_hamo,orb_elem2(4),orb_elem2,sp_vf2);

xlabel('x')
ylabel('y')
zlabel('z')
% xlim([-2.047*10^8, -2.0461*10^8])
% ylim([-2.4992*10^8, -2.4982*10^8])
% zlim([3.22*10^7, 3.235*10^7])
view(-10,45)
grid
title("Vesta SOI close-up (arrival)")

%% Vesta close-up demo (departure)
if exist('figure2') == 0
    figure()
else
    figure2()
end

park_orbit(10,Vesta_r3,Vpark_radius,orb_elem2(4));
park_orbit(10,Vesta_r3,Vpark_radius,orb_elem3(4));
hold on

%           DEBUG
% track3 = [[VC_orbit(1,1);VC_orbit(2,1)] ...
%     [VC_orbit(1,2);VC_orbit(2,2)] ...
%     [VC_orbit(1,3);VC_orbit(2,3)];...
%     [VC_orbit(end-1,1);VC_orbit(end,1)] ...
%     [VC_orbit(end-1,2);VC_orbit(end,2)] ...
%     [VC_orbit(end-1,3);VC_orbit(end,3)]];
% plot3(track3(1:2,1),track3(1:2,2),track3(1:2,3),'k')

escape_hyp(10,VC_orbit(1:2,1:3),[2012 9 5 0 0 0],Vpark_radius,...
                              orb_elem3(4), orb_elem3,sp_vf3);

xlabel('x')
ylabel('y')
zlabel('z')
xlim([-2.046542*10^8, -2.04653*10^8])
ylim([-2.49873*10^8, -2.498718*10^8])
zlim([3.22814*10^7, 3.22822*10^7])
view(-10,45)
grid
title("Vesta close-up (departure)")

%% Vesta SOI close-up demo (departure)
if exist('figure2') == 0
    figure()
else
    figure2()
end

body_sphere(10,Vesta_r3)
hold on
surface(Vesta_r3(1)+Vesta_SOI*xx, Vesta_r3(2)+Vesta_SOI*yy,...
      Vesta_r3(3)+Vesta_SOI*zz,'FaceColor','none','EdgeColor',colors(10))
  
%           DEBUG
% plot3(track2(3:4,1),track2(3:4,2),track2(3:4,3),'k') 
% plot3(track3(1:2,1),track3(1:2,2),track3(1:2,3),'k')

escape_hyp(10,VC_orbit(1:2,1:3),[2012 9 5 0 0 0],...
            Vpark_radius, orb_elem3(4), orb_elem3,sp_vf3);

xlabel('x')
ylabel('y')
zlabel('z')
xlim([-2.047*10^8, -2.0461*10^8])
ylim([-2.4992*10^8, -2.4982*10^8])
zlim([3.22*10^7, 3.235*10^7])
view(-10,45)
grid
title("Vesta SOI close-up (departure)")

%% Ceres close-up demo
if exist('figure2') == 0
    figure()
else
    figure2()
end

park_orbit(11,Ceres_r4,Ceres_hamo,orb_elem3(4));
hold on

%           DEBUG
% Ctrack = [[VC_orbit(end-1,1);VC_orbit(end,1)] ...
%     [VC_orbit(end-1,2);VC_orbit(end,2)] ...
%     [VC_orbit(end-1,3);VC_orbit(end,3)]];
% plot3(Ctrack(:,1),Ctrack(:,2),Ctrack(:,3),'k')

capture_hyp(11,VC_orbit(end-1:end,1:3),[2015 3 6 0 0 0],...
                    Ceres_hamo,orb_elem3(4),orb_elem3,sp_vf3);

xlabel('x')
ylabel('y')
zlabel('z')
xlim([-3.57995*10^8, -3.5799*10^8])
ylim([1.156455*10^8, 1.156505*10^8])
zlim([6.93861*10^7, 6.93874*10^7])
view(-10,45)
grid
title("Ceres close-up")

%% Ceres SOI close-up
if exist('figure2') == 0
    figure()
else
    figure2()
end

park_orbit(11,Ceres_r4,Ceres_hamo,orb_elem3(4));
hold on
surface(Ceres_r4(1)+Ceres_SOI*xx, Ceres_r4(2)+Ceres_SOI*yy,...
      Ceres_r4(3)+Ceres_SOI*zz,'FaceColor','none','EdgeColor',colors(11))

%           DEBUG
% plot3(track3(3:4,1),track3(3:4,2),track3(3:4,3),'k')

capture_hyp(11,VC_orbit(end-1:end,1:3),[2015 3 6 0 0 0],...
                    Ceres_hamo,orb_elem3(4),orb_elem3,sp_vf3);

xlabel('x')
ylabel('y')
zlabel('z')
xlim([ -3.5808*10^8, -3.579*10^8])
ylim([1.1556*10^8, 1.1574*10^8])
zlim([6.928*10^7, 6.948*10^7])
view(-10,45)
grid
title("Ceres SOI close-up")

%% Gathering positions for a TBD animation
positions = [hyperbola;
             EM_orbit(2:end,1:3)];