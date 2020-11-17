%% Preparations
clc
clear all
close all
tic

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

 masses = 10^24 * [0.330104
                      4.86732
                      5.97219
                      0.641693
                      1898.13
                      568.319
                      86.8103
                      102.410
                      0.01309
                      0.000259
                      0.0009393
                      1989100]; %[kg]

	radii = [2439.7
             6051.8 
             6371
             3389.5
             69911
             58232
             25362
             24622
             1151
             262.7
             476.2
             695508]; %[km] 

	distances = [57909227
                 108209475
                 149598262
                 227943824
                 778340821
                 1426666422
                 2870658186
                 4498396441
                 5906440628
                 353649000000
                 413690250];%[km]
    
    G = 6.6742e-20; %[km^3/kg/s^2]

%% Data

%Distances from the Sun
Vesta_to_Sun = distances(10);%353*10^6; %[km]
Ceres_to_Sun = distances(11);%446145795; %[km]

%Masses
Vesta_mass = masses(10);%2.589*10^20; %[kg]
Ceres_mass = masses(11);%8.958*10^20; %[kg]
Sun_mass = masses(12);%1.989*10^30; %[kg]

%SOI: (m_planet/m_Sun)^(2/5) * distance_from_Sun
Earth_SOI = (masses(3)/Sun_mass)^(2/5)*distances(3);%9.24*10^5; %[km]
Mars_SOI = (masses(4)/Sun_mass)^(2/5)*distances(4);%5.74*10^5; %[km]
Ceres_SOI = (Ceres_mass/Sun_mass)^(2/5)*Ceres_to_Sun; %[km]
Vesta_SOI = (Vesta_mass/Sun_mass)^(2/5)*Vesta_to_Sun; %[km]

%Parking orbits (not considering body radius)
Epark_radius = 200; %[km]
Epark_inclination = 0; %[rad]

Vesta_hamo = 670; %[km]
Vesta_lamo = 180; %[km]
Ceres_hamo = 1320; %[km]
Ceres_lamo = 700; %[km]

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
    - 4: at Ceres arrival (from Vesta) -> 5/3/15
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
                        planet_elements_and_sv(10,2015,3,5,0,0,0);

%Ceres
[Ceres_coe3, Ceres_r3, Ceres_v3, ~] =...
                        planet_elements_and_sv(11,2012,9,5,0,0,0);
[Ceres_coe4, Ceres_r4, Ceres_v4, ~] =...
                        planet_elements_and_sv(11,2015,3,5,0,0,0);

%% Interplanetary plot
%% Earth - Mars travel
if exist('figure2') == 0  %#ok<*EXIST>
    figure()
else
    figure2()
end

title("Interplanetary orbits")
hold on

%Interplanetary orbit
fprintf('\n\n EARTH TO MARS \n\n')
[body_pos1, sp_v1, body_posf1, sp_vf1,tof1, orb_elem1] = ...
                gen_orbit(3,4,[2007 9 27 0 0 0],[2009 2 17 0 0 0],0);
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
% [body_pos2, sp_v2, body_posf2, sp_vf2, tof2, orb_elem2] = ...
%                 gen_orbit(4,10,[2009 2 17 0 0 0],[2011 7 16 0 0 0],0);
           
%[dep_r, dep_v, arr_r, arr_v, flight, orb_oe]
[body_pos21, sp_todawn, body_posf21, sp_vf21, tof21, orb_elem21] = ...
                gen_orbit(4,10,[2009 2 17 0 0 0],[2011 7 16 0 0 0],1);
[body_pos22, sp_fromdawn, body_posf22, sp_vf22, tof22, orb_elem22] = ...
                gen_orbit(4,10,[2009 2 17 0 0 0],[2011 7 16 0 0 0],2);
MD_orbit = intpl_orbit(tof21,Mars_r1,sp_todawn);
DV_orbit = intpl_orbit(tof22,body_posf21,sp_fromdawn);
MV_orbit = [MD_orbit;DV_orbit];
            
% MV_orbit = intpl_orbit(tof2,Mars_r1,sp_v2);

%Planet orbits
plot_orbit(10,2011)

%Planet positions
plot3(Mars_r2(1),Mars_r2(2),Mars_r2(3),'rd')
plot3(Vesta_r1(1),Vesta_r1(2),Vesta_r1(3),'x','Color',colors(10))
plot3(Vesta_r2(1),Vesta_r2(2),Vesta_r2(3),'*','Color',colors(10))
%% Vesta - Ceres travel

%Interplanetary orbit
fprintf('\n\n VESTA TO CERES \n\n')
[body_pos3, sp_v3, body_posf3, sp_vf3, tof3, orb_elem3] = ...
                gen_orbit(10,11,[2012 9 5 0 0 0],[2015 3 5 0 0 0],0);
VC_orbit = intpl_orbit(tof3,Vesta_r3,sp_v3);

%Planet orbits
plot_orbit(11,2015)

%Planet positions
plot3(Vesta_r3(1),Vesta_r3(2),Vesta_r3(3),'d','Color',colors(10))
plot3(Vesta_r4(1),Vesta_r4(2),Vesta_r4(3),'+','Color',colors(10))
plot3(Ceres_r3(1),Ceres_r3(2),Ceres_r4(3),'d','Color',colors(11))
plot3(Ceres_r4(1),Ceres_r4(2),Ceres_r4(3),'+','Color',colors(11))

%% Earth close-up
if exist('figure2') == 0
    figure()
else
    figure2()
end

xlabel('x')
ylabel('y')
zlabel('z')
xlim([1.49698*10^8, 1.49715*10^8])
ylim([8.916*10^6, 8.932*10^6])
zlim([-1*10^4, 10^4])
view(-10,45)
grid
title("Earth close-up")
hold on

Earth_esc = escape_hyp(3,EM_orbit(1:2,1:3),[2007 9 27 0 0 0],...
                                         Epark_radius, orb_elem1, sp_v1);

[park_E0,t_E0] = park_out(3, Earth_r0, Epark_radius, orb_elem1,...
                  Earth_esc(1,1:3), [2007,9,1,0,0,0], [2007,9,27,0,0,0]);

% park_orbit(3,Earth_r0,Epark_radius,Epark_inclination,orb_elem1(3));
% park_orbit(3,Earth_r0,Epark_radius,orb_elem1(4),orb_elem1(3));

%           DEBUG
% Etrack = [[EM_orbit(1,1);EM_orbit(2,1)] ...
%     [EM_orbit(1,2);EM_orbit(2,2)] ...
%     [EM_orbit(1,3);EM_orbit(2,3)];...
%     [EM_orbit(end-1,1);EM_orbit(end,1)] ...
%     [EM_orbit(end-1,2);EM_orbit(end,2)] ...
%     [EM_orbit(end-1,3);EM_orbit(end,3)]];
% plot3(Etrack(1:2,1),Etrack(1:2,2),Etrack(1:2,3),'k')

%% Earth SOI close-up
if exist('figure2') == 0
    figure()
else
    figure2()
end

xlabel('x')
ylabel('y')
zlabel('z')
xlim([ 1.486*10^8, 1.508*10^8])
ylim([8*10^6, 10*10^6])
zlim([-10*10^5, 10*10^5])
view(-10,45)
grid
title("Earth SOI close-up")
hold on

surface(Earth_r0(1)+Earth_SOI*xx, Earth_r0(2)+Earth_SOI*yy,...
        Earth_r0(3)+Earth_SOI*zz,'FaceColor','none','EdgeColor',colors(3))

%           DEBUG
% Etrack = [[EM_orbit(1,1);EM_orbit(2,1)] ...
%     [EM_orbit(1,2);EM_orbit(2,2)] ...
%     [EM_orbit(1,3);EM_orbit(2,3)];...
%     [EM_orbit(end-1,1);EM_orbit(end,1)] ...
%     [EM_orbit(end-1,2);EM_orbit(end,2)] ...
%     [EM_orbit(end-1,3);EM_orbit(end,3)]];
% plot3(Etrack(1:2,1),Etrack(1:2,2),Etrack(1:2,3),'k')

hyperbola = escape_hyp(3,EM_orbit(1:2,1:3),[2007 9 27 0 0 0],...
                        Epark_radius, orb_elem1, norm(sp_vf1));
                    
park_out(3, Earth_r0, Epark_radius, orb_elem1,...
                  Earth_esc(1,1:3), [2007,9,1,0,0,0], [2007,9,27,0,0,0]);

%% Mars close-up
if exist('figure2') == 0
    figure()
else
    figure2()
end

xlabel('x')
ylabel('y')
zlabel('z')
xlim([ 9.3478*10^7, 9.3488*10^7])
ylim([-1.8886*10^8, -1.88851*10^8])
zlim([-6.257*10^6, -6.248*10^6])
view(-10,45)
grid
title("Mars close-up")
hold on

body_sphere(4,Mars_r1);

%           DEBUG
Etrack = [[EM_orbit(1,1);EM_orbit(2,1)] ...
    [EM_orbit(1,2);EM_orbit(2,2)] ...
    [EM_orbit(1,3);EM_orbit(2,3)];...
    [EM_orbit(end-1,1);EM_orbit(end,1)] ...
    [EM_orbit(end-1,2);EM_orbit(end,2)] ...
    [EM_orbit(end-1,3);EM_orbit(end,3)]];
Mtrack = [[MD_orbit(1,1);MD_orbit(2,1)] ...
    [MD_orbit(1,2);MD_orbit(2,2)] ...
    [MD_orbit(1,3);MD_orbit(2,3)];...
    [MD_orbit(end-1,1);MD_orbit(end,1)] ...
    [MD_orbit(end-1,2);MD_orbit(end,2)] ...
    [MD_orbit(end-1,3);MD_orbit(end,3)]];
plot3(Mtrack(1:2,1),Mtrack(1:2,2),Mtrack(1:2,3),'k')
plot3(Etrack(3:4,1),Etrack(3:4,2),Etrack(3:4,3),'k')

% flyby(planet_id,theta_inf,altitude,flag,year,month,day,hour,minute,second)
[par,fly] = flyby(4,91,512,1,2009,02,17,0,0,0);
plot3(fly(:,1), fly(:,2),fly(:,3),'g-')

%% Mars SOI close-up
if exist('figure2') == 0
    figure()
else
    figure2()
end

xlabel('x')
ylabel('y')
zlabel('z')
xlim([ 9.28*10^7, 9.42*10^7])
ylim([-1.896*10^8, -1.88*10^8])
zlim([-7*10^6, -5.5*10^6])
view(-10,45)
grid
title("Mars SOI close-up")
hold on

body_sphere(4,Mars_r1)

surface(Mars_r1(1)+Mars_SOI*xx, Mars_r1(2)+Mars_SOI*yy,...
        Mars_r1(3)+Mars_SOI*zz,'FaceColor','none','EdgeColor',colors(4))

%           DEBUG
% Mtrack = [[MD_orbit(1,1);MD_orbit(2,1)] ...
%     [MD_orbit(1,2);MD_orbit(2,2)] ...
%     [MD_orbit(1,3);MD_orbit(2,3)];...
%     [MD_orbit(end-1,1);MD_orbit(end,1)] ...
%     [MD_orbit(end-1,2);MD_orbit(end,2)] ...
%     [MD_orbit(end-1,3);MD_orbit(end,3)]];
% plot3(Mtrack(1:2,1),Mtrack(1:2,2),Mtrack(1:2,3),'k')
% plot3(Etrack(3:4,1),Etrack(3:4,2),Etrack(3:4,3),'k')

%% Vesta close-up (arrival)
if exist('figure2') == 0
    figure()
else
    figure2()
end

xlabel('x')
ylabel('y')
zlabel('z')
xlim([1.964214*10^8, 1.964233*10^8])
ylim([-2.676874*10^8, -2.676854*10^8])
zlim([-1.58848*10^7, -1.5884*10^7])
view(-10,45)
grid
title("Vesta close-up (arrival)")
hold on

Vesta_cap = capture_hyp(10,MV_orbit(end-1:end,1:3),[2011 7 16 0 0 0],...
                        Vesta_hamo,orb_elem22,sp_fromdawn);

[park_V2,t_V2] = park_in(10, Vesta_r2, Vesta_hamo, orb_elem22,...
                Vesta_cap(end,1:3), [2011,7,16,0,0,0], [2011,8,1,0,0,0]);


r1 = (1.0e+08 * [ 1.964218321663772  -2.676871375732791  -0.158842884222893]);
r2 = (1.0e+08 * [ 1.964225414871708  -2.676859686053302  -0.158844026563676]);
tf = 12000;
grade = 'pro';
[orb_change_park, t_change_park, deltav_park] = ...
                        park_orbit_change(10, Vesta_r2, r1, r2, tf, grade);

%           DEBUG
% Vtrack = [[MV_orbit(1,1);MV_orbit(2,1)] ...
%     [MV_orbit(1,2);MV_orbit(2,2)] ...
%     [MV_orbit(1,3);MV_orbit(2,3)];...
%     [MV_orbit(end-1,1);MV_orbit(end,1)] ...
%     [MV_orbit(end-1,2);MV_orbit(end,2)] ...
%     [MV_orbit(end-1,3);MV_orbit(end,3)]];
% plot3(Vtrack(end-1:end,1),Vtrack(end-1:end,2),Vtrack(end-1:end,3),'k')

%% Vesta SOI close-up (arrival)
if exist('figure2') == 0
    figure()
else
    figure2()
end
xlabel('x')
ylabel('y')
zlabel('z')
xlim([1.55*10^8,2.4*10^8])
ylim([-3.1*10^8,-2.28*10^8])
zlim([-5.6*10^7,2.4*10^7])
view(-10,45)
grid
title("Vesta SOI close-up (arrival)")
hold on

surface(Vesta_r2(1)+Vesta_SOI*xx, Vesta_r2(2)+Vesta_SOI*yy,...
     Vesta_r2(3)+Vesta_SOI*zz,'FaceColor','none','EdgeColor',colors(10))

%           DEBUG
% Vtrack = [[MV_orbit(1,1);MV_orbit(2,1)] ...
%     [MV_orbit(1,2);MV_orbit(2,2)] ...
%     [MV_orbit(1,3);MV_orbit(2,3)];...
%     [MV_orbit(end-1,1);MV_orbit(end,1)] ...
%     [MV_orbit(end-1,2);MV_orbit(end,2)] ...
%     [MV_orbit(end-1,3);MV_orbit(end,3)]];
% plot3(Vtrack(3:4,1),Vtrack(3:4,2),Vtrack(3:4,3),'k') 

capture_hyp(10,MV_orbit(end-1:end,1:3),[2011 7 16 0 0 0],...
                        Vesta_hamo,orb_elem22,sp_fromdawn);

park_in(10, Vesta_r2, Vesta_hamo, orb_elem22,...
                Vesta_cap(end,1:3), [2011,7,16,0,0,0], [2011,7,16,2,0,0]);

%% Vesta close-up (departure)
if exist('figure2') == 0
    figure()
else
    figure2()
end

xlabel('x')
ylabel('y')
zlabel('z')
xlim([2.144949*10^8, 2.144969*10^8])
ylim([3.143305*10^8, 3.143325*10^8])
zlim([-3.55521*10^7, -3.55515*10^7])
view(-10,45)
grid
title("Vesta close-up (departure)")
hold on

% park_orbit(10,Vesta_r3,Vesta_hamo,orb_elem22(4),orb_elem22(3));
% park_orbit(10,Vesta_r3,Vesta_lamo,orb_elem3(4),orb_elem3(3));

%           DEBUG
% Vtrack2 = [[VC_orbit(1,1);VC_orbit(2,1)] ...
%     [VC_orbit(1,2);VC_orbit(2,2)] ...
%     [VC_orbit(1,3);VC_orbit(2,3)];...
%     [VC_orbit(end-1,1);VC_orbit(end,1)] ...
%     [VC_orbit(end-1,2);VC_orbit(end,2)] ...
%     [VC_orbit(end-1,3);VC_orbit(end,3)]];
% plot3(Vtrack2(1:2,1),Vtrack2(1:2,2),Vtrack2(1:2,3),'k')

Vesta_esc = escape_hyp(10,VC_orbit(1:2,1:3),[2012 9 5 0 0 0],...
                                        Vesta_lamo, orb_elem3,sp_vf3);

[park_V3,t_V3] = park_out(10, Vesta_r3, Vesta_lamo, orb_elem3,...
                  Vesta_esc(1,1:3), [2012,9,4,22,0,0], [2012,10,5,0,0,0]);

%% Vesta SOI close-up (departure)
if exist('figure2') == 0
    figure()
else
    figure2()
end

xlabel('x')
ylabel('y')
zlabel('z')
xlim([1.72*10^8, 2.55*10^8])
ylim([2.75*10^8, 3.55*10^8])
zlim([-7.6*10^7, 0.7*10^7])
view(-10,45)
grid
title("Vesta SOI close-up (departure)")
hold on

surface(Vesta_r3(1)+Vesta_SOI*xx, Vesta_r3(2)+Vesta_SOI*yy,...
      Vesta_r3(3)+Vesta_SOI*zz,'FaceColor','none','EdgeColor',colors(10))
  
%           DEBUG
% Vtrack2 = [[VC_orbit(1,1);VC_orbit(2,1)] ...
%     [VC_orbit(1,2);VC_orbit(2,2)] ...
%     [VC_orbit(1,3);VC_orbit(2,3)];...
%     [VC_orbit(end-1,1);VC_orbit(end,1)] ...
%     [VC_orbit(end-1,2);VC_orbit(end,2)] ...
%     [VC_orbit(end-1,3);VC_orbit(end,3)]];
% plot3(Vtrack2(1:2,1),Vtrack2(1:2,2),Vtrack2(1:2,3),'k')

escape_hyp(10,VC_orbit(1:2,1:3),[2012 9 5 0 0 0],...
            Vesta_lamo, orb_elem3, sp_v3);

park_out(10, Vesta_r3, Vesta_lamo, orb_elem3,...
                  Vesta_esc(1,1:3), [2012,9,4,22,0,0], [2012,10,5,0,0,0]);

%% Ceres close-up 
if exist('figure2') == 0
    figure()
else
    figure2()
end

xlabel('x')
ylabel('y')
zlabel('z')
xlim([3.4862*10^7, 3.48656*10^7])
ylim([-4.271015*10^8, -4.270978*10^8])
zlim([-1.97804*10^7, -1.97794*10^7])
view(-10,45)
grid
title("Ceres close-up")
hold on

%           DEBUG
% Ctrack = [[VC_orbit(end-1,1);VC_orbit(end,1)] ...
%     [VC_orbit(end-1,2);VC_orbit(end,2)] ...
%     [VC_orbit(end-1,3);VC_orbit(end,3)]];
% plot3(Ctrack(:,1),Ctrack(:,2),Ctrack(:,3),'k')

Ceres_cap = capture_hyp(11,VC_orbit(end-1:end,1:3),[2015 3 5 0 0 0],...
                    Ceres_hamo,orb_elem3,sp_vf3);
                
[park_C4,t_C4] = park_in(11, Ceres_r4, Ceres_hamo, orb_elem3,...
                Ceres_cap(end,1:3), [2015,3,5,0,0,0], [2015,4,1,0,0,0]);

%% Ceres SOI close-up 
if exist('figure2') == 0
    figure()
else
    figure2()
end

xlabel('x')
ylabel('y')
zlabel('z')
xlim([ 3.478*10^7, 3.495*10^7])
ylim([-4.2718*10^8,-4.2702*10^8])
zlim([-1.9865*10^7,-1.969*10^7])
view(-10,45)
grid
title("Ceres SOI close-up")
hold on

surface(Ceres_r4(1)+Ceres_SOI*xx, Ceres_r4(2)+Ceres_SOI*yy,...
      Ceres_r4(3)+Ceres_SOI*zz,'FaceColor','none','EdgeColor',colors(11))

park_in(11, Ceres_r4, Ceres_hamo, orb_elem3,...
                Ceres_cap(end,1:3), [2015,3,5,0,0,0], [2015,4,1,0,0,0]);

%           DEBUG
% Ctrack = [[VC_orbit(end-1,1);VC_orbit(end,1)] ...
%     [VC_orbit(end-1,2);VC_orbit(end,2)] ...
%     [VC_orbit(end-1,3);VC_orbit(end,3)]];
% plot3(Ctrack(:,1),Ctrack(:,2),Ctrack(:,3),'k')

capture_hyp(11,VC_orbit(end-1:end,1:3),[2015 3 5 0 0 0],...
                    Ceres_hamo,orb_elem3,sp_vf3);
                
%% Gathering positions for a TBD animation
positions = [hyperbola;
             EM_orbit(2:end,1:3)];
         
toc