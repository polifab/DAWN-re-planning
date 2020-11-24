function [flyby_parameters, G1, tof_sec] = flyby(planet_id,delta_ip,altitude,flag,year,month,day,hour,minute,second,alpha_x,beta_y,gamma_z)
global mu
%example id_planet = 4; delta_ip = 38; altitude = 512; flag = 1; year = 2009
%        month = 02, day = 17, hour = 12, minute = 00, second = 00 alpha_x
%        = 0, beta_y =0, gamma_z = 0
%flyby(4,129,512,1,2009,02,0,17,0,0,0,0,0);

%% Parameters required from the function
%Id_planet	-> Code that indicates the planet where the flyby will be done
%delta_ip	-> turning angle, indicates how many degrees the speed vector of 
%            the probe rotates during the flyby
%Altitude	-> Minimum distance from the planet
%flag		->	 1 = Leading-side planetary flyby, planet position on hyperbola primary focus; 
%				-1 = Trailing-side planetary flyby, planet position on hyperbola secondary focus;
%alpha_x	-> Rotation angle respect X axis 
%beta_y		-> Rotation angle respect Y axis
%gamma_z	-> Rotation angle respect Z axis 

%% Starting parameters
planet = inf_planet(planet_id); %This command gives back a vector that contains planet's data:
                                %planet(1) = the radius 
                                %planet(2) = the radius of the sphere of influence
                                %planet(3) = the name
deg = pi/180;
res = 1000;		% how many points we want to compute for each part of the
				% trajectory

%Angles of hyperbole
    theta_inf = 180-(180-delta_ip)/2;
%     theta = -theta_inf+0.1:0.1:theta_inf-0.1; 
	theta = linspace(-theta_inf+0.1,theta_inf-0.1, res); 
										%Vector that contains angles between 
                                        %-Theta infinity and +Theta
                                        %Infinity
    beta = 180-theta_inf; %Angle between the velocity vector of the planet 
                          %and the asymptote of the hyperbole, this
                          %parameter indicates the entrance angle (-beta)
                          %and exit angle (+beta)
    delta = 180-2*beta; %This parameter indicates how much infinity velocity
                        %vector of spaceship rotates after flyby
e = -1/cos(theta_inf*deg); %Eccentricity
rp = str2double(planet(1)) + altitude; %Minimum distance between the center
                                       %of the planet and spaceship
% a = rp/(e-1); %Semiaxis major of hyperbole
a = rp/(1-e);		%Semiaxis major of hyperbole
b = a*sqrt(e^2-1);	%Semiaxis minor of hyperbole
%% Parameter provided by the function
%flyby_parameters: vettor that contains flyby's parameters in this
%order [planet(3) planet(1) planet(2) altitude theta_inf beta delta e a b]
%% Planet position on the considered date
mu = str2double(planet(3));
[coe, r, v, jd] = planet_elements_and_sv(planet_id, year, month, day, hour, minute, second);

% ---------------------------Da riguardare---------------------------
R_sun_planet = sqrt(r*r');
longitude_angle = asin( r(3) / R_sun_planet );
angle_sun_planet = acos( r(1) / (R_sun_planet*cos(longitude_angle)) );

%% Coordinates vectors 
x = [];     
y = [];     
z = [];

%% Critical points to calculate the trajectory
X_max = -a*cosh(-theta_inf*deg);
Y_max = b*sinh(-theta_inf*deg);

%% Calculate the x,y and z coordinates of Hyperbolic Trajectory
%To calculate the Trajectory the path has been split up in three parts

    %First part of the trajectory: from the enters in the sphere of
    %influence to critical point
    %x_a = (cos(theta_inf*deg)*str2double(planet(2)))-abs(rp)-a:0.1:X_max;
	x_1 = (cos(theta_inf*deg)*str2double(planet(2)))-abs(rp)+a;
    x_a = linspace(x_1, X_max, res);
	y_a = -tan(beta*deg)*(x_a - X_max) - Y_max;
    %Save the coordinates of the second part of trajectory
    x = [x x_a];     
	y = [y y_a];     
    z = [zeros(1,length(x_a))];

    %Second part of the trajectory: when the spaceship trancits close the planet
    i = 1;
    while i <= length(theta)
        q = theta(i)*deg;	%Hyperbolic angle
            X = a*cosh(q); %Calculate the X coordinates
            if q <= 0
                Y =  b*sqrt(X^2/a^2-1); %Calculate the Y negative semi-axis coordinates
            else
                Y = -b*sqrt(X^2/a^2-1); %Calculate the Y positive semi-axis coordinates
            end
            %Save the coordinates of the second part of trajectory
            x = [x X];
            y = [y Y];
            z = [z 0];    
        i = i+1;
    end
    
    %Third part of trajectory: from the critical point to exit from sphere 
    %of influence
% 	x_a = (cos(theta_inf*deg)*str2double(planet(2)))-abs(rp)-a:1:X_max;
	x_1 = (cos(theta_inf*deg)*str2double(planet(2)))-abs(rp)+a;
    x_a = linspace(X_max, x_1, res);
    y_a = tan(beta*deg)*(x_a - X_max) + Y_max;
	
    %Save the coordinates of the third part of trajectory
    x = [x x_a]';    
	y = [y y_a]';    
    z = [z zeros(1,length(x_a))]'; 
    x = x + [abs(e*a)*ones(1,length(x))]';
	
%% Rotation Matrix
% ---------------------------Da riguardare---------------------------
%[R_x R_y R_z] = rotation_matrix(alpha_x,beta_y,gamma_z); 
[R_x R_y R_z] = rotation_matrix(alpha_x,longitude_angle/deg,angle_sun_planet/deg);

%If necessary, rotate the flyby further
%[R_x R_y R_z] = rotation_matrix(alpha_x*deg,beta_y*deg,gamma_z*deg);

%% Plot the Hyperbolic Trajectory 
grid on;
figure(1);
axis equal;
G1 = [flag*x y z];		%Transposes the G1 Matrix 
G1 = G1*R_x*R_y*R_z;	%Executes the rotation 
% Translation to have the Sun in the coordinates (0,0,0) 
    i = 1;
    while i <= length(x)
        j = 1;
        while j <= 3
            G1(i,j) = G1(i,j) + r(j);
            j = j + 1;
        end
        i = i + 1;
    end
hold on;

    figure(1);
    % Plot the trajectory of the spaceship
    plot3(G1(:,1),G1(:,2),G1(:,3));
    % Plot the circle that rapresent the planet
    th = 0:pi/360:2*pi;
    xunit = str2double(planet(1)) * cos(th) + r(1);
    yunit = str2double(planet(1)) * sin(th) + r(2);
    zunit = r(3)*ones(1,length(xunit));
    plot3(xunit, yunit, zunit);
    % Plot the circle that rapresent the radius of influnce of the planet 
    th = 0:pi/360:2*pi;
    xunit = str2double(planet(2)) * cos(th) + r(1);
    yunit = str2double(planet(2)) * sin(th) + r(2);
    zunit = r(3)*ones(1,length(xunit));
    plot3(xunit, yunit, zunit);
    title("Flyby trajectory close " + planet(4));
    % Plot the circle that rapresent the Sun
    sun_posiont = [0 0 0];
    th = 0:pi/360:2*pi;
    xunit = 6.96*10^5 * cos(th) + sun_posiont(1);
    yunit = 6.96*10^5 * sin(th) + sun_posiont(2);
    zunit = zeros(1,length(xunit));
    plot3(xunit, yunit, zunit);
hold off;

%% Polar coordinates of flyby and plot
    rr =[]; %this vector contains the distance from the center of planet and the 
			%point on the spaceship's trajectory 
    i = 1;
    while i <= length(theta)
        vett_r = rp/(1+e*cos(theta(i)*deg)); 
        rr = [rr flag*vett_r]; %save the distance
        i = i + 1;
    end
    % Plot the spaceship's trajectory
    figure('Name','Polar Plot Flyby on ' + planet(4));
    polarplot(theta*deg,rr);
    title("Polar Plot Flyby trajectory close " + planet(4));
    hold on;
    % Plot the planet's influence circle
    t = 0:1:360-1;
    rr = str2double(planet(2))*ones(1,length(t));
    polarplot(t*deg,rr,'-.');
    thetaticks([0:10:360]);
    polarplot(flag*(e*a),'.-k');
    % limit the plot to planet's influence circle
    rlim([0 str2double(planet(2))]);
    hold off;
	
%% Print datas of flyby
v_inf = sqrt(-str2double(planet(3))/a);
    fprintf('\n   Radius of influence of ' + planet(4) + '     = ' + planet(2) + ' (km) %g\n') 
    fprintf('\n   Radius of ' + planet(4) + '                  = ' + planet(1) + ' (km) %g\n')
    fprintf('\n   Altitude                        = ' + string(altitude) + ' (km) %g\n')
    fprintf('\n   True anomaly                    = ' + string(theta_inf) + '%g\n')
    fprintf('\n   Beta angle                      = ' + string(beta) + '%g\n')
    fprintf('\n   Delta angle                     = ' + string(delta) + '%g\n')
    fprintf('\n   Eccentricity                    = ' + string(e) + '%g\n')
    fprintf('\n   V infinity                      = ' + string(v_inf)+ ' (km/s)' + '%g\n') 
    fprintf('\n   Semi major axis                 = ' + string(a) + ' (km) %g\n')
    fprintf('\n   Semi minor axis                 = ' + string(b)+ ' (km)' + '\n')
    
    % Output vector
    flyby_parameters = [planet(3) planet(1) planet(2) altitude theta_inf beta delta e a b];
 
	%% Time of flight
	
	p = a* (1- e^2);
	ni_0 = acos((p/rr(1) -1)/e);
	F_0 = 2*atanh(sqrt((e-1)/(e+1)) * tan(ni_0/2))/deg;
	
% 	F_0 = 2*atanh(sqrt((e-1)/(e+1)) * tan(theta_inf*deg/2));
	F_f = theta(end)*deg;
	
	
	tof = ((flag*a)^3 / mu)^0.5 * ((e * sinh(F_f) - F_f) - (e * sinh(F_0) - F_0) ); % [TU]
	tof_sec = tof * (str2double(planet(1))^3 * 1/mu)^0.5;
	%}
end