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
% Epark_radius = 200;		%[km]
% Epark_inclination = 0;	%[rad]
Epark_radius = (165+320)/2;			%[km] extremes of Curiosity park orbit
Epark_inclination = 28.28/180*pi;	%[rad]

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
year = [88 
            225 
            365 
            687 
            4330 
            10748 
            30666 
            60148 
            90560 
            1340 
            1690 
            25];

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

%% intl orbits (only needed ones)

%mars-vesta intl orbit
[body_pos21, sp_todawn, body_posf21, sp_vf21, tof21, orb_elem21] = ...
                gen_orbit(4,10,[2009 2 17 0 0 0],[2011 7 16 0 0 0],1);
[body_pos22, sp_fromdawn, body_posf22, sp_vf22, tof22, orb_elem22] = ...
                gen_orbit(4,10,[2009 2 17 0 0 0],[2011 7 16 0 0 0],2);
MD_orbit = intpl_orbit(tof21,Mars_r1,sp_todawn);
DV_orbit = intpl_orbit(tof22,body_posf21,sp_fromdawn);
MV_orbit = [MD_orbit;DV_orbit];

%vesta ceres inlt orbit
[body_pos3, sp_v3, body_posf3, sp_vf3, tof3, orb_elem3] = ...
                gen_orbit(10,11,[2012 9 5 0 0 0],[2015 3 5 0 0 0],0);
VC_orbit = intpl_orbit(tof3,Vesta_r3,sp_v3);

%% Earth Close by


[body_pos1, sp_v1, body_posf1, sp_vf1,tof1, orb_elem1] = ...
                gen_orbit(3,4,[2007 9 27 0 0 0], [2009 2 17 0 0 0],0);
			
[EM_orbit, t_EM] = intpl_orbit(tof1, Earth_r0, sp_v1);		

Earth_esc = escape_hyp(3,EM_orbit(1:2,1:3),[2007 9 27 0 0 0],...
                                         Epark_radius, orb_elem1, sp_v1);

%change of inclination of parking orbit
orb_elem1_park = orb_elem1;
orb_elem1_park(4) =  Epark_inclination;

[park_E1,t_E1] = park_out(3, Earth_r0, Epark_radius, orb_elem1,...
                  Earth_esc(1,1:3), [2007,9,18,0,0,0], [2007,9,27,0,0,0]);

[park_E0,t_E0] = park_out(3, Earth_r0, Epark_radius, orb_elem1_park,...
                  park_E1(1,1:3), [2007,9,1,0,0,0], [2007,9,20,0,0,0]);

%--------------------------------------------------------------------------
%find cross point  
n = 60;
pos1 = park_E1(1:n,:);
pos0 = park_E0(end-n+1:end,:);

mindiff_true = 1e10;
index_cross = 0;

% scan
for i = 1:n
	temp = pos0(i,:);
	norm_diff =  row_norm2(pos1 - temp);
	
	mindiff = mink(norm_diff,1);
	if mindiff < mindiff_true
		mindiff_true = mindiff;
		index_cross = i;
	end
end

index_crosspoint_0 = find(row_norm2(park_E0) == row_norm2(pos0(index_cross,:)),1,'last');

index_crosspoint_1 = find(row_norm2(park_E1 - pos0(index_cross,:)) < 100 ,1,'first');

park_E0(index_crosspoint_0:end,:) = [];
t_E0(index_crosspoint_0:end,:) = [];
park_E1(1:index_crosspoint_1-1,:) = [];

% stacking		  
spcr_E = [];
spcr_E = cat(1, spcr_E, park_E0);
spcr_E = cat(1, spcr_E, park_E1);
spcr_E = cat(1, spcr_E, Earth_esc);

%--------------------------------------------------------------------------
% time vector and stacking

% days_Earth = ymd_gen([2007,9,1], [2007,9,27]);
days_Earth = ymd_gen([2007,9,1], [2008,9,27]);
% increased number of days around the earth to include the time for escape hyp. 

% creating more time samples
t_E_sofar = [t_E0; t_E0(end) + t_E1];
% best fit (least square error) poly order 1
poly_tE = polyfit(1:size(t_E_sofar,1), t_E_sofar, 1);
t_esc_E = polyval(poly_tE, (size(t_E_sofar,1)+1):size(spcr_E,1) )';


% stacking
time_seconds_E = [];
time_seconds_E = cat(1, time_seconds_E, t_E0);
time_seconds_E = cat(1, time_seconds_E, t_E0(end) + t_E1);
time_seconds_E = cat(1, time_seconds_E, t_esc_E);

close all;

%% Mars fly by

% flyby parameters
delta_ip = 38;		% [deg]
altitude = 512; 
flag = 1;	% Leading-side planetary flyby, planet position on hyperbola 
			% primary focus; 
date_M = [2009 2 17 0 0 0];

ang_M = [0,0,0]; % flyby angles, alpha_x, beta_y, gamma_z


[~, flyby_M]= flyby(	4,	delta_ip, altitude, 1, ...
						date_M(1), date_M(2), date_M(3), date_M(4), date_M(5), date_M(6),...
						ang_M(1), ang_M(1), ang_M(1));

spcr_M = flip(flyby_M - Mars_r1,1);

%--------------------------------------------------------------------------
% time

t_tot = 10*24*60*60;
time_seconds_M = linspace(1, t_tot, size(spcr_M,1) );
% flyby perigee at [2009 2 17]
% adding more days after for safety bounds
days_Mars = ymd_gen([2009,2,13], [2009,2,28]);

%{ 
%Old time management reling on 2*r_soi/v_rel = t_tot
Sp_vf1 = [17.4416249980260  17.7099375279402  0.536070520792872];

%Mars_SOI = (masses(4)/Sun_mass)^(2/5)*distances(4);	%5.77*10^5; %[km]
% ref_length = Mars_SOI;
ref_length = norm(spcr_M(1,:),2) * 2;

t_tot = ref_length/ ( norm(Sp_vf1 - Mars_v1,2)) ;
time_seconds_M = linspace(1, t_tot, size(spcr_M,1) );

% flyby happens at [2009 2 17]
% ndays = floor(t_tot/60/60); % only 45~ hours
% adding more days after for safety bounds
days_Mars = ymd_gen([2009,2,16], [2009,2,19]);

% debug
figure(1)
clf
% plot3(flyby_M_red(:,1),flyby_M_red(:,2),flyby_M_red(:,3));
% plot3(flyby_M(:,1),flyby_M(:,2),flyby_M(:,3));
plot3(spcr_M(:,1),spcr_M(:,2),spcr_M(:,3));
hold on
for i = 1:1e2:size(spcr_M)
	text(spcr_M(i,1),spcr_M(i,2),spcr_M(i,3),num2str(i));
end
grid on
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')

close all;


%}

% redefining mu
global mu
mu = 1.327565122000000e+11; %[km^3/s^2]
%% Vesta close by

%{
% BY Simone
% arrival
Vesta_cap = capture_hyp(10,MV_orbit(end-1:end,1:3),[2011 7 16 0 0 0],...
                        Vesta_hamo,orb_elem22,sp_fromdawn);

[park_V2,t_V2] = park_in(10, Vesta_r2, Vesta_hamo, orb_elem22,...
                Vesta_cap(end,1:3), [2011,7,16,0,0,0], [2011,8,1,0,0,0]);
% departure
Vesta_esc = escape_hyp(10,VC_orbit(1:2,1:3),[2012 9 5 0 0 0],...
                                        Vesta_lamo, orb_elem3,sp_vf3);
[park_V3,t_V3] = park_out(10, Vesta_r3, Vesta_lamo, orb_elem3,...
                  Vesta_esc(1,1:3), [2012,9,4,22,0,0], [2012,10,5,0,0,0]);

%}
%%%BY Francesco
% arrival
Vesta_cap = capture_hyp(10,MV_orbit(end-1:end,1:3),[2011 7 16 0 0 0],...
                        Vesta_hamo,orb_elem22,sp_fromdawn);

[park_V2,t_V2] = park_in(10, Vesta_r2, Vesta_hamo, orb_elem22,...
                Vesta_cap(end,1:3), [2011,7,16,0,0,0], [2012,2,4,22,0,0]);
% departure
Vesta_esc = escape_hyp(10,VC_orbit(1:2,1:3),[2012 9 5 0 0 0],...
                                        Vesta_lamo, orb_elem3,sp_vf3);
[park_V3,t_V3] = park_out(10, Vesta_r3, Vesta_lamo, orb_elem3,...
                  Vesta_esc(1,1:3), [2012,2,4,22,0,0], [2012,10,5,0,0,0]);
%%%

%--------------------------------------------------------------------------
% change parking orbit
r1 = (1.0e+08 * [ 1.964218321663772  -2.676871375732791  -0.158842884222893]);
r2 = (1.0e+08 * [ 1.964225414871708  -2.676859686053302  -0.158844026563676]);
tf = 12000;

grade = 'pro';
[orb_change_park, t_change_park, deltav_park] = ...
                        park_orbit_change(10, Vesta_r2, r1, r2, tf, grade);

% change park orbit adjustment
temporary = zeros(size(park_V2,1),1);
for i = 1:size(park_V2,1)
	temporary(i) = norm(park_V2(i,:) - orb_change_park(1,:) ,2);
end
index_last = find(temporary-mink(temporary,1) < 1e-5, 1, 'last');
%removing the samples after the start of changing orbit
park_V2(index_last+1:end,:) = [];

temporary = zeros(size(park_V3,1),1);
for i = 1:size(park_V3,1)
	temporary(i) = norm(park_V3(i,:)- Vesta_r3 + Vesta_r2 - orb_change_park(end,:), 2);
end
index_first = find(temporary-mink(temporary,1) < 1e-5, 1, 'first');

%removing the samples before the end point of changing orbit
park_V3(1:index_first,:) = [];


% reduce hyperbola capture and exit dimensions
des_length = 500;	% desired length of hyperbolae
Vesta_cap(1:(end-des_length),:) = [];
Vesta_esc(des_length:end,:) = [];

% stacking
spcr_V = [];
spcr_V = cat(1, spcr_V, Vesta_cap			-Vesta_r2);
spcr_V = cat(1, spcr_V, park_V2				-Vesta_r2);
spcr_V = cat(1, spcr_V, orb_change_park		-Vesta_r2);
spcr_V = cat(1, spcr_V, park_V3				-Vesta_r3);
spcr_V = cat(1, spcr_V, Vesta_esc			-Vesta_r3);

%--------------------------------------------------------------------------
% time vector and stacking

% best fit (least square error) poly order 1
poly_tV_cap = polyfit(1:size(t_V2,1), t_V2, 1);
t_V_cap_noshift = polyval(poly_tV_cap, - size(Vesta_cap,1):0 )';
t_V_cap = t_V_cap_noshift + abs(t_V_cap_noshift(1));

poly_tV_esc = polyfit(1:size(t_V3,1), t_V3, 1);
t_V_esc = polyval(poly_tV_esc, (size(t_V3,1)+1):(size(t_V3,1) + size(Vesta_esc,1)) )';

% final stacking
time_seconds_V = [];
time_seconds_V = cat(1, time_seconds_V, t_V_cap);
time_seconds_V = cat(1, time_seconds_V, t_V_cap(end)	+ t_V2);
time_seconds_V = cat(1, time_seconds_V, t_V_cap(end)	+ t_V2(end)		+ t_change_park);
time_seconds_V = cat(1, time_seconds_V, t_V_cap(end)	+ t_V2(end)		+ t_change_park(end)	+ t_V3);
time_seconds_V = cat(1, time_seconds_V, t_V_cap(end)	+ t_V2(end)		+ t_change_park(end)	+ t_V_esc);

% ndays_cap_hyp  =  abs(mink(floor(t_V_cap_noshift/24/3600),1));
% so days_vector should start 4 days before ceres arrival ([2011 7 16).
% length(ymd_gen([2011 7 13], [2011 7 16]))
days_Vesta = ymd_gen([2011 7 13], [2013 7 16]);
% increased number of days after the departure ([2012,10,5])to generate more than enough
% days.

%% Ceres close by

Ceres_cap = capture_hyp(11,VC_orbit(end-1:end,1:3),[2015 3 5 0 0 0],...
                    Ceres_hamo,orb_elem3,sp_vf3);
                
[park_C4,t_C4] = park_in(11, Ceres_r4, Ceres_hamo, orb_elem3,...
                Ceres_cap(end,1:3), [2015,3,5,0,0,0], [2015,4,1,0,0,0]);

%stacking
spcr_C = [];
spcr_C = cat(1, spcr_C, Ceres_cap - Ceres_r4);
spcr_C = cat(1, spcr_C, park_C4 - Ceres_r4);

%--------------------------------------------------------------------------
% time vector and stacking

% best fit (least square error) poly order 1
poly_tC = polyfit(1:size(t_C4,1), t_C4, 1);
t_cap_C_noshift = polyval(poly_tC, -size(Ceres_cap,1):0 )';
% adding abs(t_cap_C_noshift(1)) to shift up correctly
t_cap_C = t_cap_C_noshift + abs(t_cap_C_noshift(1));

% final stacking
time_seconds_C = [];
time_seconds_C = cat(1, time_seconds_C, t_cap_C);
time_seconds_C = cat(1, time_seconds_C, t_C4 + abs(t_cap_C(1)) );

% ndays_cap_hyp  =  abs(mink(floor(t_cap_C_noshift/24/3600),1));
% so days_vector should start 19 days before ceres arrival ([2015,3,5]).
% length(ymd_gen([2015,2,15], [2015,3,5]))
days_Ceres = ymd_gen([2015,2,15], [2016,3,15]);
% increased number of days after the arrival to generate more than enough
% days.

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

col_bkgnd	= [0,	0,	0];		% Background
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
col_dawn_face = [50,50,50]./255;
col_dawn_edge = [130,130,130]./255;

%% Movie parameters

% init of frame counter
k = 1;
movie_fps = 20;


