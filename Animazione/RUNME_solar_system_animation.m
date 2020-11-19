%% Script description
% This script executes the animation of the Dawn Mission
% The animation works by plotting each day the position of planets, Sun,
% Ceres, Vesta and the spacecraft. Also, relevant points are plotted
% aswell.

% Attention: this script requires "M-Files for Orbital Mechanics for 
% Engineering Students, 3e" folder in Matlab-path

%% intro
clc; clear;

movie_mode		= 0;	% 1 for movie writing, 2 for HD movie, 0 for only matlab animation
addpath('functions');

solar_system_animation_init;	% init all par, constants and plot parameters
pos_spcr_solar;					% computes the position of the spacecraft 
								% for each day of the mission

%% Solar System Plot and Animation
figh = figure(1);
clf 
set(gcf, 'Renderer', 'zbuffer');
set(gca, 'color', col_bkgnd)
ax = gca;
ax.GridColor = col_grid; 
hold on

% Orbits (NOT ACCURATE, SHOULD BE UPDATED EACH YEAR) 
% TEMPORARY
plot_orbit2(3,	time_vector(1,1), planet_linewidth)		%Earth
plot_orbit2(4,	time_vector(1,1), planet_linewidth)		%Mars
plot_orbit2(10, time_vector(end,1), planet_linewidth)	%Vesta
plot_orbit2(11, time_vector(end,1), planet_linewidth)	%Ceres

for d = 1:3:n_days
% for d = 1:fr_skip:n_days
% for d = (n_days-10):fr_skip:n_days
	%----------------------------------------------------------------------
	%reset figure
	if d ~= 1
		delete(p_Sun)
		delete(p_Earth)
		delete(p_Mars)
		delete(p_Vesta)
		delete(p_Ceres)
		delete(p_spcr)
	end
	
	% now
	year_now	= time_vector(d,1);
	month_now	= time_vector(d,2);
	day_now		= time_vector(d,3);
	
	%Sun (id = 12)
	[~, Sun_now, ~, ~] = planet_elements_and_sv(12, ...
								year_now, month_now, day_now, 0, 0, 0);
							
	%Earth (id = 3)
	[~, Earth_now, ~, ~] = planet_elements_and_sv(3, ...
								year_now, month_now, day_now, 0, 0, 0);
	%Mars (id = 4)
	[~, Mars_now, ~, ~] = planet_elements_and_sv(4, ...
								year_now, month_now, day_now, 0, 0, 0);
	%Vesta (id = 10)
	[~, Vesta_now, ~, ~] = planet_elements_and_sv(10, ...
								year_now, month_now, day_now, 0, 0, 0);
	
	% dirty trick for reaching ceres in the animation.
	early = floor(33 * d/n_days);
	year_now	= time_vector(d-early,1);
	month_now	= time_vector(d-early,2);
	day_now		= time_vector(d-early,3);
	
	%Ceres (id = 11)
	[~, Ceres_now, ~, ~] = planet_elements_and_sv(11, ...
								year_now, month_now, day_now, 0, 0, 0);
	
	% SpaceCraft
	spcr_now = pos_spcr(d,:);
	
	%----------------------------------------------------------------------						
	%%plot
	p_Sun		= plot3(Sun_now(1),Sun_now(2),Sun_now(3),...
						'o','Color',col_sun, 'MarkerSize', dim_sun,...
						'MarkerFaceColor',col_sun);
	
	p_Earth		= plot3(Earth_now(1),Earth_now(2),Earth_now(3),...
						'o','Color',col_earth, 'MarkerSize', dim_earth,...
						'MarkerFaceColor',col_earth);
					
	p_Mars		= plot3(Mars_now(1),Mars_now(2),Mars_now(3),...
						'o','Color',col_mars, 'MarkerSize', dim_mars,...
						'MarkerFaceColor',col_mars);
	p_Vesta		= plot3(Vesta_now(1),Vesta_now(2),Vesta_now(3),...
						'o','Color',col_vesta, 'MarkerSize', dim_vesta,...
						'MarkerFaceColor',col_vesta);				
	p_Ceres		= plot3(Ceres_now(1),Ceres_now(2),Ceres_now(3),...
						'o','Color',col_ceres, 'MarkerSize', dim_ceres,...
						'MarkerFaceColor',col_ceres);
	p_spcr		= plot3(spcr_now(1),spcr_now(2),spcr_now(3),...
						'o','Color',col_spcr, 'MarkerSize', dim_spcr,...
						'MarkerFaceColor',col_spcr);
	%sofar
	pos_tillnow = pos_spcr(1:d,:);
	plot3(pos_tillnow(:,1),pos_tillnow(:,2),pos_tillnow(:,3),...
			'-','Color', col_spcr,'LineWidth', width_spcr);
	
	
	
	
	%----------------------------------------------------------------------
	% end stuff				
	axis equal;
    grid on;
	hold on
	if spinlon == 0 && spinlat == 0
		if d == 1
			view(View(1, 1), View(1, 2));
		end
	else
		view(View(1, 1) + spinlon * d/n_days, View(1, 2) + spinlat* d/n_days);
	end
	xlabel('X [km]');
	ylabel('Y [km]');
	zlabel('Z [km]');
	date = datetime(time_vector(d,:));
	
	
	if d <= day_mars
		status_msg = ['Left Earth! To Mars...'];
	elseif d <= day_vesta 
		status_msg = ['Mars Flyby done! To Vesta...'];
	elseif d <= day_left_vesta
		status_msg = ['Park-orbit Vesta.'];
	elseif d < day_ceres
		status_msg = ['Left Vesta! To Ceres...'];
	else 
		status_msg = ['Ceres Reached.'];
	end
	
	Title = [status_msg '   '  datestr(date), ' (day ', num2str(d), ' of ', num2str(n_days), ').'];
	title(Title)
	
	% end stuff	
	axis equal;
	xlim([-xy_lim, xy_lim]);
	ylim([-xy_lim, xy_lim]);
	zlim([-z_lim, z_lim]);
	drawnow;
	
	if time_pause ~= 0
		pause(time_pause)
	end
	
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

%% Video stuff
if movie_mode
	movie = VideoWriter('movie_heliocentric', 'MPEG-4');
	movie.FrameRate = movie_fps;

	open(movie);
	writeVideo(movie, movieVector);
	close(movie);
end