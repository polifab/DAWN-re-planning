%% General parameters - To edit
% animation rates
time_pause = 0;		% time [s] after each drawing. Set to zero to avoid pausing
fr_skip = 0 +1;		% frame skip between each drawing

%% Constants
global mu
mu = 1.327565122000000e+11; %[km^3/s^2]

%Distances from the Sun
Ceres_to_Sun = 446145795;	%[km]
Vesta_to_Sun = 353*10^6;	%[km]

%Masses
Vesta_mass = 2.589*10^20;	%[kg]
Ceres_mass = 8.958*10^20;	%[kg]
Sun_mass = 1.989*10^30;		%[kg]

%SOI: (m_planet/m_Sun)^(2/5) * distance_from_Sun
Earth_SOI = 9.24*10^5;	%[km]
Mars_SOI = 5.74*10^5;	%[km]
Ceres_SOI = (Ceres_mass/Sun_mass)^(2/5)*Ceres_to_Sun; %[km]
Vesta_SOI = (Vesta_mass/Sun_mass)^(2/5)*Vesta_to_Sun; %[km]

%Parking orbits (not considering body radius)
Epark_radius = 200;		%[km]
Epark_inclination = 0;	%[rad]

Vesta_hamo = 670;	%[km]
Vesta_lamo = 180;	%[km]
Ceres_hamo = 1320;	%[km]
Ceres_lamo = 700;	%[km]

% time_vector, a row for each day of the mission
time_vector = ymd_gen([2007, 9, 27],[2015, 3, 6]);
n_days = size(time_vector,1);

% radii of planets
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
             695508];

%% Mission info
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

%% Earth Close by

[body_pos1, sp_v1, body_posf1, sp_vf1,tof1, orb_elem1] = ...
                gen_orbit(3,4,[2007 9 27 0 0 0], [2009 2 17 0 0 0],0);
[EM_orbit, t_EM] = intpl_orbit(tof1, Earth_r0, sp_v1);		

orbit_E = park_orbit(3, Earth_r0, Epark_radius, orb_elem1(4), orb_elem1(3));

% removing some points
orbit_E(1:floor(size(orbit_E,1)*4/5),:) = [];

hyperbola_E = escape_hyp(3,EM_orbit(1:2,1:3),[2007 9 27 0 0 0],...
                        Epark_radius, orb_elem1, norm(sp_vf1));

spcr_E = [];
spcr_E = cat(1, spcr_E, orbit_E);
spcr_E = cat(1, spcr_E, hyperbola_E);

%% Mars fly by

%spcr_M = ;

%% Vesta close by

%intro
[body_pos3, sp_v3, body_posf3, sp_vf3, tof3, orb_elem3] = ...
                gen_orbit(10,11,[2012 9 5 0 0 0],[2015 3 5 0 0 0],0);
[body_pos21, sp_todawn, body_posf21, sp_vf21, tof21, orb_elem21] = ...
                gen_orbit(4,10,[2009 2 17 0 0 0],[2011 7 16 0 0 0],1);
[body_pos22, sp_fromdawn, body_posf22, sp_vf22, tof22, orb_elem22] = ...
                gen_orbit(4,10,[2009 2 17 0 0 0],[2011 7 16 0 0 0],2);
MD_orbit = intpl_orbit(tof21,Mars_r1,sp_todawn);
DV_orbit = intpl_orbit(tof22,body_posf21,sp_fromdawn);
MV_orbit = [MD_orbit;DV_orbit];
% to ceres
VC_orbit = intpl_orbit(tof3,Vesta_r3,sp_v3);
			
			
% arrival
%park_orbit(10,Vesta_r2,Vesta_hamo,orb_elem22(4),orb_elem22(3));
orbit_V_arr			= park_orbit(10,Vesta_r2,Vesta_lamo,orb_elem3(4),orb_elem3(3));
hyperbola_V_arr		= capture_hyp(10,MV_orbit(end-1:end,1:3),[2011 7 16 0 0 0],...
									 Vesta_hamo,orb_elem22,sp_fromdawn);



% change park orbits?						 



% departure
%park_orbit(10,Vesta_r3,Vesta_hamo,orb_elem22(4),orb_elem22(3));
orbit_V_dep			= park_orbit(10,Vesta_r3,Vesta_lamo,orb_elem3(4),orb_elem3(3));
hyperbola_V_dep		= escape_hyp(10,VC_orbit(1:2,1:3),[2012 9 5 0 0 0],Vesta_lamo,...
                                             orb_elem3,sp_vf3);
							
% reduce dimensions								 
hyperbola_V_arr = flip(hyperbola_V_arr,1);
hyperbola_V_arr(1:floor(size(hyperbola_V_arr,1)*28.5/29),:) = [];		
hyperbola_V_dep(floor(size(hyperbola_V_dep,1)*1/39):end,:) = [];

% stacking
spcr_V = [];
spcr_V = cat(1, spcr_V, hyperbola_V_arr);
spcr_V = cat(1, spcr_V, orbit_V_arr);
spcr_V = cat(1, spcr_V, orbit_V_dep);
spcr_V = cat(1, spcr_V, hyperbola_V_dep);

%% Ceres close by

% arrival
hyperbola_C = capture_hyp(11,VC_orbit(end-1:end,1:3),[2015 3 5 0 0 0],...
							Ceres_hamo,orb_elem3,sp_vf3);
orbit_C		= park_orbit(11,Ceres_r4,Ceres_hamo,orb_elem3(4),orb_elem3(3));

% stacking
spcr_C = [];
spcr_C = cat(1, spcr_C, hyperbola_C);
spcr_C = cat(1, spcr_C, orbit_C);


%% Animations Parameters
% view angles
View = [120 45];				% initial view
spinlon = 50;
spinlat = -30;

% init of status msg
status_msg = ['Nothing'];

% axis lim in plots
lim_gap = 200;
scale_lims = 1;

% colours of objects, rgb
% nice link to get rbg colour:
% https://www.rapidtables.com/web/color/RGB_Color.html
col = [	"g"		     %green
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

col_bkgnd	= [0.3,	0.3,	0.3];		% Background
col_grid	= [255,	255, 255]	/255;	% Grid
col_sun		= [255, 204, 0]		/255;	% Sun
col_earth	= col(3);	%[0,	102, 204]	/255;	% Earth
col_mars	= col(4);	%[255,	102, 0]		/255;	% Mars
col_vesta	= col(10);	%[0,	102, 204]	/255;	% Vesta
col_ceres	= col(11);	%[102,	153, 153]	/255;	% Ceres
col_spcr	= [0,	196, 255]	/255;	% spacecraft

% width of lines
planet_linewidth = 0.1;
width_spcr  = 1;

% dawn model
% original dimensions: 1.64 × 19.7 × 1.77 m 
dim_dawn = [1.64 19.7 1.77];
scale_stl = dim_dawn(2)/120;
scale_dawn = 100;
col_dawn_face = [0,0,0]./255;
col_dawn_edge = [100,100,100]./255;

%% Movie parameters

% init of frame counter
k = 1;
movie_fps = 20;


