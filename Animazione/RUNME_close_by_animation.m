%% Script description
% This script executes the animations of each close by of the Dawn Mission
% The animations works by plotting the position of the spacecraft and the 
% "close" planet and displaying the current day.

% Attention, this script requires:
% 1. "M-Files for Orbital Mechanics for  Engineering Students, 3e" folder in Matlab-path
% 2. "lavoro vollaro" folder in Matlab-path

%% TO DO list
%{

Close by Mars

%}

%% intro & editable parameters
clc; clear;

% addpath
addpath('functions');	
addpath('models');

% init
animation_close_by_init;		% init all par, constants, plot parameters 
								% and all needed vectors
clc;
close all;

% Editable
movie_mode		= -1;	%  1	for movie writing
						%  2	for writing HD movie
						%  0	for only matlab animation
						% -1	for fullscreen animation
						
initial_spin	= 1;	% 1 for initial inspection of spacecraft, 0 otherwise
only_closeby	= 1;	% case option to animate only a specified close by.
						%		0		all
						%		1		Earth
						%		2		Mars
						%		3		Vesta
						%		4		Ceres
						
%% Earth Close by
if only_closeby == 1 || only_closeby == 0

clc
disp('Earth Close by animation started!')

figh = figure(1);
clf
if movie_mode == 2 || movie_mode == -1
	if movie_mode == 2
		warning('Note that a 1080p resolution is needed for movie_mode = 2')
	end
	figh.WindowState = 'maximize';
end
% figure settings
lighting phong;
set(gcf, 'Renderer', 'zbuffer');
set(gca, 'color', col_bkgnd)
ax = gca;
ax.GridColor = col_grid; 
axis equal
axis tight
hold on
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
grid on

% Edit inputs to iterations
y = spcr_E-Earth_r0;				% positions of spacecraft
id_planet = 3;
scale_dawn = 50;					% to edit!

% init positions
pos_planet = [0;0;0];
pos_dawn_init = y(1,:)';			% column vector of the initial dawn position
body_sphere(id_planet, pos_planet);

% init dawn model
stl_offset = [0; 0; 0];
p_dawn = patch(stlread('Dawn_19.stl'));							% EDITARE
T_dawn = hgmat(eye(3)*scale_stl* scale_dawn, pos_dawn_init);	% EDITARE
t_dawn = hgtransform('Parent',gca);
set(t_dawn,'Matrix',T_dawn);
set(p_dawn,'Parent',t_dawn);
set(p_dawn, 'facec', col_dawn_face);            % Set the face colour
set(p_dawn, 'EdgeColor', col_dawn_edge);		% Set the edge colour

% iterations
n = size(y,1)-200;
k = 1;
hold on
for i = 1:fr_skip:n
% for i = 1:100:n
	
	% initial inspection of dawn spacecraft
	if i == 1 && initial_spin
		
		spinlon_insp = 120;
		spinlat_insp = 30;
		
		n_ii = 300;		% half of number of frames in which we do a close inspection
		for ii = 1:n_ii
			%view angle
			view1 = View(1, 1) + spinlon_insp * ii/n_ii;
			view2 = View(1, 2) + spinlat_insp* ii/n_ii;
			view(view1, view2);
			
			%limits
			
			minx1 = pos_dawn_init(1) - scale_dawn*dim_dawn(2)/2;
			maxx1 = pos_dawn_init(1) + scale_dawn*dim_dawn(2)/2;
			xlim([minx1,maxx1])
			miny1 = pos_dawn_init(2) - scale_dawn*dim_dawn(1)/2;
			maxy1 = pos_dawn_init(2) + scale_dawn*dim_dawn(1)/2;
			ylim([miny1,maxy1])
			minz1 = pos_dawn_init(3) - scale_dawn*dim_dawn(3)/2;
			maxz1 = pos_dawn_init(3) + scale_dawn*dim_dawn(3)/2;
			zlim([minz1,maxz1])
			
			% title
			title('Spacecraft Inspection')
			
			% pause
			if time_pause ~= 0
				pause(time_pause)
			end
			drawnow
			%--------------------------------------------------------------
			% write video
			if movie_mode == 1 || movie_mode == 2
				%movieVector(k) = getframe(figh);
				%open to fullscreen figure(1) before using the next line
				movieVector(k) = getframe(figh, [10,10,1910,960]);
				k = k+1;
			end
		end
		
		
		for ii = n_ii:(2*n_ii)
			view1 = View(1, 1) + spinlon_insp - spinlon_insp * (ii-n_ii)/n_ii;
			view2 = View(1, 2) + spinlat_insp - spinlat_insp* (ii-n_ii)/n_ii;
			view(view1, view2);
			
			% lims
			minx = minx1 * (2*n_ii-ii)/n_ii + ...
					min(min(y(1:i,1))-lim_gap , -radii(id_planet)-lim_gap) * (ii-n_ii)/n_ii;
			maxx = maxx1 * (2*n_ii-ii)/n_ii + ...
					max(max(y(1:i,1))+lim_gap , +radii(id_planet)+lim_gap) * (ii-n_ii)/n_ii;
				
			miny = miny1 * (2*n_ii-ii)/n_ii + ...
					min(min(y(1:i,2))-lim_gap , -radii(id_planet)-lim_gap) * (ii-n_ii)/n_ii;
			maxy = maxy1 * (2*n_ii-ii)/n_ii + ...
					max(max(y(1:i,2))+lim_gap , +radii(id_planet)+lim_gap) * (ii-n_ii)/n_ii;
			xlim(scale_lims*[minx, maxx]);
			ylim(scale_lims*[miny, maxy]);
			
			minz = minx1 * (2*n_ii-ii)/n_ii + ...
					(-radii(id_planet)-lim_gap) * (ii-n_ii)/n_ii;
			maxz = maxx1 * (2*n_ii-ii)/n_ii + ...
					(radii(id_planet)+lim_gap) * (ii-n_ii)/n_ii;
				
			zlim([minz, maxz]);
			
			% title
			title('Spacecraft Inspection')
			
			% pause
			if time_pause ~= 0
				pause(time_pause)
			end
			drawnow
			%--------------------------------------------------------------
			% write video
			if movie_mode == 1 || movie_mode == 2
				%movieVector(k) = getframe(figh);
				%open to fullscreen figure(1) before using the next line
				movieVector(k) = getframe(figh, [10,10,1910,960]);
				k = k+1;
			end
			
		end
	end
		
	%----------------------------------------------------------------------
	% dawn update
	pos_now = y(i,:);
	T_dawn_now = hgmat(eye(3)*scale_stl*scale_dawn, pos_now');
	set(t_dawn,'Matrix',T_dawn_now);
	
	
	% trajectory "so far" of spacecraft
	traj =	plot3(  y(1:i,1), y(1:i,2),y(1:i,3),...
            'color', col_spcr, 'LineWidth', width_spcr);
	
		
	%----------------------------------------------------------------------
	% end stuff	
	% Title = ['Earth Close Up (frame ', num2str(i), ' of ', num2str(n), ')'];
	
	nday = floor(time_seconds_E(i)/24/60/60)+1;
	
	date = datetime(days_Earth(nday,:));
	Title = ['Earth Close Up. ' datestr(date) ', (day ' num2str(nday), ...
			' of ', num2str(floor(time_seconds_E(end)/24/60/60)), ').'];
	title(Title)
	
	% axis lims
	minx = min(min(y(1:i,1))-lim_gap , -radii(id_planet)-lim_gap);
	maxx = max(max(y(1:i,1))+lim_gap , +radii(id_planet)+lim_gap);
	xlim(scale_lims*[minx, maxx]);
	
	miny = min(min(y(1:i,2))-lim_gap , -radii(id_planet)-lim_gap);
	maxy = max(max(y(1:i,2))+lim_gap , +radii(id_planet)+lim_gap);
	ylim(scale_lims*[miny, maxy]);
	
	minz = min(min(y(1:i,3))-lim_gap, -radii(id_planet)-lim_gap);
	maxz = max(max(y(1:i,3))+lim_gap, +radii(id_planet)+lim_gap);
	zlim([minz, maxz]);
	
	
	if spinlon == 0 && spinlat == 0
		if i == 1
			view(View(1, 1), View(1, 2));
		end
	else
		view(View(1, 1) + spinlon * i/n, View(1, 2) + spinlat* i/n);
	end
	
	
	if time_pause ~= 0
		pause(time_pause)
	end
	
	drawnow
	%----------------------------------------------------------------------
	% write video
	if movie_mode == 1 || movie_mode == 2
		if movie_mode == 1
			movieVector(k) = getframe(figh);
		elseif movie_mode == 2
			movieVector(k) = getframe(figh, [10,10,1910,960]);
		end		
		k = k+1;
	end
	
end

%--------------------------------------------------------------------------
% Video stuff
if movie_mode == 1 || movie_mode == 2
	movie = VideoWriter('movie_Earthcloseby', 'MPEG-4');
	movie.FrameRate = movie_fps;

	open(movie);
	writeVideo(movie, movieVector);
	close(movie);
	disp('Earth Close by movie saved!')

end

end

%% Mars FlyBy
if only_closeby == 2 || only_closeby == 0

	% to do

end

%% Vesta Close By
if only_closeby == 3 || only_closeby == 0
clc
disp('Vesta Close by animation started!')

figh = figure(1);
clf
if movie_mode == 2 || movie_mode == -1
	if movie_mode == 2
		warning('Note that a 1080p resolution is needed for movie_mode = 2')
	end
	figh.WindowState = 'maximize';
end
% figure settings
lighting phong;
set(gcf, 'Renderer', 'zbuffer');
set(gca, 'color', col_bkgnd)
ax = gca;
ax.GridColor = col_grid; 
axis equal
axis tight
hold on
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
grid on

% Edit inputs to iterations
y = spcr_V;				% positions of spacecraft relative to planet 
						% already done in init.m
id_planet = 10;
scale_dawn = 35;		

% init positions
pos_planet = [0;0;0];
pos_dawn_init = y(1,:)';			% column vector of the initial dawn position
body_sphere(id_planet, pos_planet);

% init dawn model
stl_offset = [0; 0; 0];
p_dawn = patch(stlread('Dawn_19.stl'));							% EDITARE
T_dawn = hgmat(eye(3)*scale_stl* scale_dawn, pos_dawn_init);	% EDITARE
t_dawn = hgtransform('Parent',gca);
set(t_dawn,'Matrix',T_dawn);
set(p_dawn,'Parent',t_dawn);
set(p_dawn, 'facec', col_dawn_face);            % Set the face colour
set(p_dawn, 'EdgeColor', col_dawn_edge);		% Set the edge colour

% iterations
n = size(y,1);
k = 1;
hold on
for i = 1:fr_skip:n
% for i = 27450:fr_skip:n							% to get change of orbit
% for i = (n - size(Vesta_esc,1)-100):fr_skip:n		% to get the exit
	if i ~= 1
		delete(traj)
	end
			
	%----------------------------------------------------------------------
	% dawn update
	pos_now = y(i,:);
	T_dawn_now = hgmat(eye(3)*scale_stl*scale_dawn, pos_now');
	set(t_dawn,'Matrix',T_dawn_now);
	
	
	% trajectory "so far" of spacecraft
	traj =	plot3(  y(1:i,1), y(1:i,2),y(1:i,3),...
            'color', col_spcr, 'LineWidth', width_spcr);
	
		
	%----------------------------------------------------------------------
	% end stuff	
	% Title = ['Vesta Close Up (frame ', num2str(i), ' of ', num2str(n), ')'];
	
	nday = floor(time_seconds_V(i)/24/60/60)+1;
	
	date = datetime(days_Vesta(nday,:));
	Title = ['Vesta Close Up. ' datestr(date) ', (day ' num2str(nday), ...
			' of ', num2str(floor(time_seconds_V(end)/24/60/60)), ').'];
	title(Title)
	
	
	%----------------------------------------------------------------------
	% axis lims
	phase_time_advance = 200;
	if i <= size(Vesta_cap,1) - phase_time_advance
		% arriving
		phase_i = 1;
		lim_gap = 1000;
			
		minx = min(min(y(phase_i:i,1))-lim_gap, -radii(id_planet)-lim_gap);
		maxx = max(max(y(phase_i:i,1))+lim_gap, +radii(id_planet)+lim_gap);
		xlim(scale_lims*[minx, maxx]);

		miny = min(min(y(phase_i:i,2))-lim_gap, -radii(id_planet)-lim_gap);
		maxy = max(max(y(phase_i:i,2))+lim_gap, +radii(id_planet)+lim_gap);
		ylim(scale_lims*[miny, maxy]);

		minz = min(min(y(phase_i:i,3))-lim_gap, -radii(id_planet)-lim_gap);
		maxz = max(max(y(phase_i:i,3))+lim_gap, +radii(id_planet)+lim_gap);
		zlim(scale_lims*[minz, maxz]);
			
	elseif i <= size(Vesta_cap,1)
		% arriving zoom
		i_0 = size(Vesta_cap,1) - phase_time_advance;
		i_1 = size(Vesta_cap,1);
		i_tot = i_1 - i_0;
		ii = (i - i_0) / i_tot;
		
		% calculate initial lims
		lim_gap = 1000;
		minx_0 = min(min(y(1:i_0,1))-lim_gap, -radii(id_planet)-lim_gap);
		maxx_0 = max(max(y(1:i_0,1))+lim_gap, +radii(id_planet)+lim_gap);

		miny_0 = min(min(y(1:i_0,2))-lim_gap, -radii(id_planet)-lim_gap);
		maxy_0 = max(max(y(1:i_0,2))+lim_gap, +radii(id_planet)+lim_gap);

		minz_0 = min(min(y(1:i_0,3))-lim_gap, -radii(id_planet)-lim_gap);
		maxz_0 = max(max(y(1:i_0,3))+lim_gap, +radii(id_planet)+lim_gap);

		% calculate final lims
		lim_gap = 2000;
		orb_samples = 200;	% additional samples to get good orbiting lims
		i_2 = i_1 + orb_samples;
		minx_1 = min(min(y(i_1:i_2,1))-lim_gap, -radii(id_planet)-lim_gap);
		maxx_1 = max(max(y(i_1:i_2,1))+lim_gap, +radii(id_planet)+lim_gap);

		miny_1 = min(min(y(i_1:i_2,2))-lim_gap, -radii(id_planet)-lim_gap);
		maxy_1 = max(max(y(i_1:i_2,2))+lim_gap, +radii(id_planet)+lim_gap);

		minz_1 = min(min(y(i_1:i_2,3))-lim_gap, -radii(id_planet)-lim_gap);
		maxz_1 = max(max(y(i_1:i_2,3))+lim_gap, +radii(id_planet)+lim_gap);

		% sets lims
		p = 1; %power of the index advancement, 1 for linear interpolation
		xlim(scale_lims*[minx_1 * ii^p + minx_0 * (1-ii)^p, maxx_1 * ii^p + maxx_0 * (1-ii)^p])
		ylim(scale_lims*[miny_1 * ii^p + miny_0 * (1-ii)^p, maxy_1 * ii^p + maxy_0 * (1-ii)^p])
		zlim(scale_lims*[minz_1 * ii^p + minz_0 * (1-ii)^p, maxz_1 * ii^p + maxz_0 * (1-ii)^p])
		
	elseif i <= n - size(Vesta_esc,1)
		% orbiting
		phase_i = size(Vesta_cap,1); 
		
		% lim_gap zoom
		zoom_dur = 400;		% zoom duration
		ii = (i-phase_i)/zoom_dur;
		p = 1; %power of the index advancement, 1 for linear interpolation
		if ii >= 1
			% zoom already ended
			lim_gap = 200;
		else
			% zoom 
			lim_gap = 200 * ii^p + 2000 * (1-ii)^p;
		end
		
		% sets lims
		minx = min(min(y(phase_i:i,1))-lim_gap, -radii(id_planet)-lim_gap);
		maxx = max(max(y(phase_i:i,1))+lim_gap, +radii(id_planet)+lim_gap);
		xlim(scale_lims*[minx, maxx]);

		miny = min(min(y(phase_i:i,2))-lim_gap, -radii(id_planet)-lim_gap);
		maxy = max(max(y(phase_i:i,2))+lim_gap, +radii(id_planet)+lim_gap);
		ylim(scale_lims*[miny, maxy]);

		minz = min(min(y(phase_i:i,3))-lim_gap, -radii(id_planet)-lim_gap);
		maxz = max(max(y(phase_i:i,3))+lim_gap, +radii(id_planet)+lim_gap);
		zlim([minz, maxz]);
		
	else
		% departing
		phase_i = n - size(Vesta_esc,1);

		minx = min(min(y(phase_i:i,1))-lim_gap, -radii(id_planet)-lim_gap);
		maxx = max(max(y(phase_i:i,1))+lim_gap, +radii(id_planet)+lim_gap);
		xlim(scale_lims*[minx, maxx]);

		miny = min(min(y(phase_i:i,2))-lim_gap, -radii(id_planet)-lim_gap);
		maxy = max(max(y(phase_i:i,2))+lim_gap, +radii(id_planet)+lim_gap);
		ylim(scale_lims*[miny, maxy]);

		minz = min(min(y(phase_i:i,3))-lim_gap, -radii(id_planet)-lim_gap);
		maxz = max(max(y(phase_i:i,3))+lim_gap, +radii(id_planet)+lim_gap);
		zlim([minz, maxz]);
		
	end
		
	if spinlon == 0 && spinlat == 0
		if i == 1
			view(View(1, 1), View(1, 2));
		end
	else
		view(View(1, 1) + spinlon * i/n, View(1, 2) + spinlat* i/n);
	end
	
	if time_pause ~= 0
		pause(time_pause)
	end
	
	drawnow
	%----------------------------------------------------------------------
	% write video
	if movie_mode == 1 || movie_mode == 2
		if movie_mode == 1
			movieVector(k) = getframe(figh);
		elseif movie_mode == 2
			movieVector(k) = getframe(figh, [10,10,1910,960]);
		end		
		k = k+1;
	end
	
end

%--------------------------------------------------------------------------
% Video stuff
if movie_mode == 1 || movie_mode == 2
	movie = VideoWriter('movie_Vestacloseby', 'MPEG-4');
	movie.FrameRate = movie_fps;

	open(movie);
	writeVideo(movie, movieVector);
	close(movie);
	disp('Vesta Close by movie saved!')

end

end

%% Ceres Close By
if only_closeby == 4 || only_closeby == 0
clc
disp('Ceres Close by animation started!')

figh = figure(1);
clf
if movie_mode == 2 || movie_mode == -1
	if movie_mode == 2
		warning('Note that a 1080p resolution is needed for movie_mode = 2')
	end
	figh.WindowState = 'maximize';
end
% figure settings
lighting phong;
set(gcf, 'Renderer', 'zbuffer');
set(gca, 'color', col_bkgnd)
ax = gca;
ax.GridColor = col_grid; 
axis equal
axis tight
hold on
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
grid on

% Edit inputs to iterations
y = spcr_C;				% positions of spacecraft relative to planet 
						% already done in init.m
id_planet = 11;
scale_dawn = 50;					% to edit!

% init positions
pos_planet = [0;0;0];
pos_dawn_init = y(1,:)';			% column vector of the initial dawn position
body_sphere(id_planet, pos_planet);

% init dawn model
stl_offset = [0; 0; 0];
p_dawn = patch(stlread('Dawn_19.stl'));							% EDITARE
T_dawn = hgmat(eye(3)*scale_stl* scale_dawn, pos_dawn_init);	% EDITARE
t_dawn = hgtransform('Parent',gca);
set(t_dawn,'Matrix',T_dawn);
set(p_dawn,'Parent',t_dawn);
set(p_dawn, 'facec', col_dawn_face);            % Set the face colour
set(p_dawn, 'EdgeColor', col_dawn_edge);		% Set the edge colour


% iterations
n = size(y,1);
k = 1;
hold on
% for i = 1:fr_skip:n
for i = 1:3:n
	if i ~= 1
		delete(traj)
	end
			
	%----------------------------------------------------------------------
	% dawn update
	pos_now = y(i,:);
	T_dawn_now = hgmat(rotz(pi/6)*scale_stl*scale_dawn, pos_now');
	set(t_dawn,'Matrix',T_dawn_now);
	
	
	% trajectory "so far" of spacecraft
	traj =	plot3(  y(1:i,1), y(1:i,2),y(1:i,3),...
            'color', col_spcr, 'LineWidth', width_spcr);
	
	%----------------------------------------------------------------------
	% end stuff	
	% Title = ['Ceres Close Up (frame ', num2str(i), ' of ', num2str(n), ')'];
	
	nday = floor(time_seconds_C(i)/24/60/60)+1;
	
	date = datetime(days_Ceres(nday,:));
	Title = ['Ceres Close Up. ' datestr(date) ', (day ' num2str(nday), ...
			' of ', num2str(floor(time_seconds_C(end)/24/60/60)), ').'];
	title(Title)
	
	%----------------------------------------------------------------------
	% axis lims
	phase_time_advance = 200;
	if i <= size(Ceres_cap,1) - phase_time_advance
		% arriving
		phase_i = 1;
		lim_gap = 1000;
			
		minx = min(min(y(phase_i:i,1))-lim_gap, -radii(id_planet)-lim_gap);
		maxx = max(max(y(phase_i:i,1))+lim_gap, +radii(id_planet)+lim_gap);
		xlim(scale_lims*[minx, maxx]);

		miny = min(min(y(phase_i:i,2))-lim_gap, -radii(id_planet)-lim_gap);
		maxy = max(max(y(phase_i:i,2))+lim_gap, +radii(id_planet)+lim_gap);
		ylim(scale_lims*[miny, maxy]);

		minz = min(min(y(phase_i:i,3))-lim_gap, -radii(id_planet)-lim_gap);
		maxz = max(max(y(phase_i:i,3))+lim_gap, +radii(id_planet)+lim_gap);
		zlim(scale_lims*[minz, maxz]);
			
	elseif i <= size(Ceres_cap,1)
		% arriving zoom
		i_0 = size(Ceres_cap,1) - phase_time_advance;
		i_1 = size(Ceres_cap,1);
		i_tot = i_1 - i_0;
		ii = (i - i_0) / i_tot;
		
		% calculate initial lims
		lim_gap = 1000;
		minx_0 = min(min(y(1:i_0,1))-lim_gap, -radii(id_planet)-lim_gap);
		maxx_0 = max(max(y(1:i_0,1))+lim_gap, +radii(id_planet)+lim_gap);

		miny_0 = min(min(y(1:i_0,2))-lim_gap, -radii(id_planet)-lim_gap);
		maxy_0 = max(max(y(1:i_0,2))+lim_gap, +radii(id_planet)+lim_gap);

		minz_0 = min(min(y(1:i_0,3))-lim_gap, -radii(id_planet)-lim_gap);
		maxz_0 = max(max(y(1:i_0,3))+lim_gap, +radii(id_planet)+lim_gap);

		% calculate final lims
		lim_gap = 2000;
		orb_samples = 200;	% additional samples to get good orbiting lims
		i_2 = i_1 + orb_samples;
		minx_1 = min(min(y(i_1:i_2,1))-lim_gap, -radii(id_planet)-lim_gap);
		maxx_1 = max(max(y(i_1:i_2,1))+lim_gap, +radii(id_planet)+lim_gap);

		miny_1 = min(min(y(i_1:i_2,2))-lim_gap, -radii(id_planet)-lim_gap);
		maxy_1 = max(max(y(i_1:i_2,2))+lim_gap, +radii(id_planet)+lim_gap);

		minz_1 = min(min(y(i_1:i_2,3))-lim_gap, -radii(id_planet)-lim_gap);
		maxz_1 = max(max(y(i_1:i_2,3))+lim_gap, +radii(id_planet)+lim_gap);

		% sets lims
		p = 1; %power of the index advancement, 1 for linear interpolation
		xlim(scale_lims*[minx_1 * ii^p + minx_0 * (1-ii)^p, maxx_1 * ii^p + maxx_0 * (1-ii)^p])
		ylim(scale_lims*[miny_1 * ii^p + miny_0 * (1-ii)^p, maxy_1 * ii^p + maxy_0 * (1-ii)^p])
		zlim(scale_lims*[minz_1 * ii^p + minz_0 * (1-ii)^p, maxz_1 * ii^p + maxz_0 * (1-ii)^p])
		
	else
		% orbiting
		phase_i = size(Ceres_cap,1); 
		
		% lim_gap zoom
		zoom_dur = 400;		% zoom duration
		ii = (i-phase_i)/zoom_dur;
		p = 1; %power of the index advancement, 1 for linear interpolation
		if ii >= 1
			% zoom already ended
			lim_gap = 200;
		else
			% zoom 
			lim_gap = 200 * ii^p + 2000 * (1-ii)^p;
		end
		
		% sets lims
		minx = min(min(y(phase_i:i,1))-lim_gap, -radii(id_planet)-lim_gap);
		maxx = max(max(y(phase_i:i,1))+lim_gap, +radii(id_planet)+lim_gap);
		xlim(scale_lims*[minx, maxx]);

		miny = min(min(y(phase_i:i,2))-lim_gap, -radii(id_planet)-lim_gap);
		maxy = max(max(y(phase_i:i,2))+lim_gap, +radii(id_planet)+lim_gap);
		ylim(scale_lims*[miny, maxy]);

		minz = min(min(y(phase_i:i,3))-lim_gap, -radii(id_planet)-lim_gap);
		maxz = max(max(y(phase_i:i,3))+lim_gap, +radii(id_planet)+lim_gap);
		zlim([minz, maxz]);		
	end
	
	%----------------------------------------------------------------------
	% end stuff
	if spinlon == 0 && spinlat == 0
		if i == 1
			view(View(1, 1), View(1, 2));
		end
	else
		view(View(1, 1) + spinlon * i/n, View(1, 2) + spinlat* i/n);
	end
	
	if time_pause ~= 0
		pause(time_pause)
	end
	
	drawnow
	%----------------------------------------------------------------------
	% write video
	if movie_mode == 1 || movie_mode == 2
		if movie_mode == 1
			movieVector(k) = getframe(figh);
		elseif movie_mode == 2
			movieVector(k) = getframe(figh, [10,10,1910,960]);
		end		
		k = k+1;
	end
	
end

%--------------------------------------------------------------------------
% Video stuff
if movie_mode == 1 || movie_mode == 2
	movie = VideoWriter('movie_Cerescloseby', 'MPEG-4');
	movie.FrameRate = movie_fps;

	open(movie);
	writeVideo(movie, movieVector);
	close(movie);
	disp('Ceres Close by movie saved!')

end

end
