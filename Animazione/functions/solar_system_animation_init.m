%% General parameters - To edit
% animation rates
time_pause = 0;		% time [s] after each drawing. Set to zero to avoid pausing
fr_skip = 0 +1;	% frame skip between each drawing

% view angles
View = [30 10];				% initial view
spinlon = 120;				% how much view angle change long [grad]
spinlat = 20;				% how much view angle change lat[grad]

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

%% Plots Parameters

% init of status msg
status_msg = ['Nothing'];

% axis lim in plots
xy_lim = 5e8;	%lim in xy plane
z_lim = 1e8;		%lim in z coord

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

% dimension of objects
dim_sun		= 30;
dim_earth	= 8;
dim_mars	= 6;
dim_vesta	= 5;
dim_ceres	= 5;
dim_spcr	= 3;

% width of lines
planet_linewidth = 0.1;
width_spcr  = 1;

%% Movie parameters

% init of frame counter
k = 1;
movie_fps = 20;