year_now = 2010;
month_now = 10;
day_now = 1;

%Vesta (id = 10)
	[~, Vesta_now, ~, ~] = planet_elements_and_sv(10, ...
								year_now, month_now, day_now, 0, 0, 0);
	%Ceres (id = 11)
	[~, Ceres_now, ~, ~] = planet_elements_and_sv(11, ...
								year_now, month_now, day_now, 0, 0, 0);
					
Vesta_now/1e6
%Ceres_now/1e6

%now2
year_now2 = 2011;
month_now2 = 11;
day_now2 = 1;

%Vesta (id = 10)
	[~, Vesta_now2, ~, ~] = planet_elements_and_sv(10, ...
								year_now2, month_now2, day_now2, 0, 0, 0);
	%Ceres (id = 11)
	[~, Ceres_now2, ~, ~] = planet_elements_and_sv(11, ...
								year_now2, month_now2, day_now2, 0, 0, 0);
					
Vesta_now2/1e6