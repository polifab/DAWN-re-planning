%% Script description
% This script executes the animations of each close by of the Dawn Mission
% The animations work by plotting each day the position of the spacecraft
% and the "close" planet. Each animation starts from SOI of such planet.


% Attention, this script requires:
% 1. "M-Files for Orbital Mechanics for  Engineering Students, 3e" folder in Matlab-path
% 2. "lavoro vollaro" folder in Matlab-path


%% TO DO list

% impostare giorni --> guarda vettore t in uscita

% fix spinlat e spinlon

% cosa fare tra le park orbits di vesta e ceres?

% caricare 3D model vesta e ceres

% Close by:
% mars
% vesta		% fare limiti intelligenti!!!!
% ceres		% fare limiti intelligenti!!!!

%% intro 
clc; clear;

animation_close_by_init;		% init all par, constants, plot parameters 
								% and all needed vectors
clc;

% EDITABLE
movie_mode		= 0;	% 1 for movie writing, 2 for HD movie, 0 for only matlab animation
initial_spin	= 1;	% 1 for initial inspection of spacecraft, 0 otherwise
only_closeby	= 3;	% case option to animate only a specified close by.
						%		0		all
						%		1		Earth
						%		2		Mars
						%		3		Vesta
						%		4		Ceres


% addpath
addpath('functions');	
addpath('models');

% prompt figure(1) HD
figh = figure(1);
clf
if movie_mode == 2
	warning('Note that a 1080p resolution is needed for movie_mode = 2')
	figh.WindowState = 'maximize';
end

%% Earth Close by
if only_closeby == 1 || only_closeby == 0

clc
disp('Earth Close by animation started!')

figh = figure(1);
clf 
% figure settings
lighting phong;
set(gcf, 'Renderer', 'zbuffer');
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
		
			% pause
			if time_pause ~= 0
				pause(time_pause)
			end
			drawnow
			%--------------------------------------------------------------
			% write video
			if movie_mode
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
			
			
			% pause
			if time_pause ~= 0
				pause(time_pause)
			end
			drawnow
			%--------------------------------------------------------------
			% write video
			if movie_mode
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
	%Title = ['Earth Close Up (frame ', num2str(i), ' of ', num2str(n), ')'];
	Title = ['Earth Close Up'];
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
	if movie_mode ~= 0
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
if movie_mode
	movie = VideoWriter('movie_earthcloseby', 'MPEG-4');
	movie.FrameRate = movie_fps;

	open(movie);
	writeVideo(movie, movieVector);
	close(movie);
	disp('Earth Close by movie saved!')

end

end

%% Mars FlyBy
if only_closeby == 2 || only_closeby == 0

%TO DO
end

%% Vesta Close By
if only_closeby == 3 || only_closeby == 0
clc
disp('Vesta Close by animation started!')

figh = figure(1);
clf 
% figure settings
lighting phong;
set(gcf, 'Renderer', 'zbuffer');
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
for i = 1:2:n
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
	Title = ['Vesta Close Up'];
	title(Title)
	
	phase_time_advance = 200;
	phase_time_wait = 200;
	
	% axis lims
	if  i <= size(Vesta_cap,1) - phase_time_advance
		% arriving
		close_index = max(i-10,1);
		close_gap = 2e3;
		minx = min(y(close_index:i,1))-close_gap;
		maxx = max(y(close_index:i,1))+close_gap;
		xlim(scale_lims*[minx, maxx]);

		miny = min(y(close_index:i,2))-close_gap;
		maxy = max(y(close_index:i,2))+close_gap;
		ylim(scale_lims*[miny, maxy]);

		minz = min(y(close_index:i,3))-close_gap;
		maxz = max(y(close_index:i,3))+close_gap;
		zlim([minz, maxz]);
		
	elseif i < (size(y,1)-size(Vesta_esc,1)) - phase_time_advance
		%orbiting
		phase_i = size(Vesta_cap,1); % "true" initial i of this phase
		
		% temporary?
		if i < phase_i
			phase_i = phase_i - phase_time_advance;
			lim_gap = 1000;
			
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
			
			
		elseif i < phase_i + phase_time_wait
			
			
			i_0 = phase_i - phase_time_advance;
			i_1 = phase_i + phase_time_wait;
			ii	= (phase_i + phase_time_wait - i)/phase_time_wait; 
			
			
			% calculate initial lims
			lim_gap = 400;
			minx_0 = min(min(y(i_0,1))-lim_gap, -radii(id_planet)-lim_gap);
			maxx_0 = max(max(y(i_0,1))+lim_gap, +radii(id_planet)+lim_gap);
			
			miny_0 = min(min(y(i_0,2))-lim_gap, -radii(id_planet)-lim_gap);
			maxy_0 = max(max(y(i_0,2))+lim_gap, +radii(id_planet)+lim_gap);

			minz_0 = min(min(y(i_0,3))-lim_gap, -radii(id_planet)-lim_gap);
			maxz_0 = max(max(y(i_0,3))+lim_gap, +radii(id_planet)+lim_gap);
			
			% calculate final lims
			lim_gap = 200;
			minx_1 = min(min(y(phase_i:i_1,1))-lim_gap, -radii(id_planet)-lim_gap);
			maxx_1 = max(max(y(phase_i:i_1,1))+lim_gap, +radii(id_planet)+lim_gap);
			
			miny_1 = min(min(y(phase_i:i_1,2))-lim_gap, -radii(id_planet)-lim_gap);
			maxy_1 = max(max(y(phase_i:i_1,2))+lim_gap, +radii(id_planet)+lim_gap);
			
			minz_1 = min(min(y(phase_i:i_1,3))-lim_gap, -radii(id_planet)-lim_gap);
			maxz_1 = max(max(y(phase_i:i_1,3))+lim_gap, +radii(id_planet)+lim_gap);
			
			% sets lims
			xlim(scale_lims*[minx_0 * ii + minx_1 * (1-ii), maxx_0 * ii + maxx_1 * (1-ii)])
			ylim(scale_lims*[miny_0 * ii + miny_1 * (1-ii), maxy_0 * ii + maxy_1 * (1-ii)])
			zlim([minz_0 * ii + minz_1 * (1-ii), maxz_0 * ii + maxz_1 * (1-ii)])
		else
			lim_gap = 400;
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
		
		
	else
		
		phase_i = (size(y,1)-size(Vesta_esc,1)); % initial i of this phase
		%departing
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
	if movie_mode ~= 0
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
if movie_mode
	movie = VideoWriter('movie_vestacloseby', 'MPEG-4');
	movie.FrameRate = movie_fps;

	open(movie);
	writeVideo(movie, movieVector);
	close(movie);
	disp('Earth Close by movie saved!')

end

end

%% Ceres Close by
if only_closeby == 4 || only_closeby == 0
	
	
	
end
