% This script calculates a vector containing each row position of
% Dawn_spacecraft with respect to solar system

% init
pos_spcr = zeros(n_days, 3);


%% Earth to Mars
% from day 1 (27/9/07) to Mars arrival/departure (17/2/09)
em_days = datenum([2009 2 17]) - datenum([2007 9 27]);

[body_pos1, sp_v1, body_posf1, sp_vf1, tof1, orb_elem1] = ...
                gen_orbit2(3, 4, [2007 9 27 0 0 0], [2009 2 17 0 0 0]);
EM_orbit = intpl_orbit2(tof1, Earth_r0, sp_v1);
em_pos = EM_orbit(:,1:3);
q = floor(size(em_pos,1)/em_days); % rate of samples for each day
% UHM UHM UHM UHM CHECKCHECK
for i = 1:em_days
	pos_spcr(i,:) = em_pos(q * i,:);
end

%% Mars to Vesta
% from (17/2/09) to (16/7/11)
mv_days = datenum([2011 7 16])- datenum([2009 2 17]);
[body_pos2, sp_v2, body_posf2, sp_vf2, tof2, orb_elem2] = ...
                gen_orbit2(4,10,[2009 2 17 0 0 0],[2011 7 16 0 0 0]);
MV_orbit = intpl_orbit2(tof2,Mars_r1,sp_v2);
mv_pos = MV_orbit(:,1:3);
q = floor(size(mv_pos,1)/mv_days); % rate of samples for each day
for i = 1:mv_days
	pos_spcr(i+em_days,:) = mv_pos(q * i,:);
end

%% Vesta Park orbit
% from (16/7/11) to (5/9/12)
v_days = datenum([2012 9 5])- datenum([2011 7 16]);

for i = (em_days + mv_days + 1):(em_days + mv_days + v_days)
	
	year_now	= time_vector(i,1);
	month_now	= time_vector(i,2);
	day_now		= time_vector(i,3);
	
	[~, tmp_Vesta_now, ~, ~] = planet_elements_and_sv2(10, ...
								year_now, month_now, day_now, 0, 0, 0);
	
	pos_spcr(i,:) = tmp_Vesta_now;
end

%% Vesta to Ceres
% from (5/9/12) to (6/3/15)
vc_days = datenum([2015 3 6])- datenum([2012 9 5]);

[body_pos3, sp_v3, body_posf3, sp_vf3, tof3, orb_elem3] = ...
                gen_orbit2(10,11,[2012 9 5 0 0 0],[2015 3 5 0 0 0]);
VC_orbit = intpl_orbit2(tof3,Vesta_r3,sp_v3);
vc_pos = VC_orbit(:,1:3);
q = floor(size(vc_pos,1)/vc_days); % rate of samples for each day
for i = 1:vc_days
	pos_spcr(i+em_days+mv_days+v_days,:) = vc_pos(q * i,:);
end

%tmp
pos_spcr(end,:) = pos_spcr(end-1,:);

%% tests
%sanity check
%(em_days + mv_days + v_days + vc_days) - (n_days-1);


%% days

day_left_earth	= 1;
day_mars		= 1 + em_days;
day_vesta		= 1 + em_days + mv_days;
day_left_vesta	= 1 + em_days + mv_days + v_days;
day_ceres		= 1 + em_days + mv_days + v_days + vc_days;
