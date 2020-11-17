function [flyby_parameters,orb] = ...
 flyby(planet_id,theta_inf,altitude,flag,year,month,day,hour,minute,second)
%   INPUT:
%       Id_planet   - Identifier of the planet where the flyby takes place
%       Theta_inf   - True anomaly, angle between velocity vector of the
%                     planet and escape velocity vector of thespaceship => 
%                     90 < Theta_inf < 180
%       Altitude    - Minimum distance from the planet
%       flag        - 1 = Leading-side planetary flyby
%                    -1 = Trailing-side planetary flyby
%       alpha_x     - Rotation angle rispect X axis 
%       beta_y      - Rotation angle rispect Y axis
%       gamma_z     - Rotation angle rispect Z axis 
%
%   OUTPUT:
%       flyby_parameters - vector containing flyby parameters in this
%                          order:
%                               [planet(3) planet(1) planet(2) 
%                                altitude theta_inf beta delta e a b]
%       orb - point of the hyperbolic trajectory
% 
%example id_planet = 4; theta = 129; altitude = 512; flag = 1; year = 2009
%        month = 02, day = 17, hour = 12, minute = 00, second = 00
%flyby(4,129,512,1,2009,02,17,12,0,0);

    %% Constants
%     global mu
    
    deg = pi/180;
    
    %% Parameters
    planet = inf_planet(planet_id); %vector containing planet data:
                                    %planet(1) = the radius 
                                    %planet(2) = the radius of the sphere of influence
                                    %planet(3) = the name
                                    
    %Angles of hyperbole
    theta = -theta_inf+0.1:0.1:theta_inf-0.1;
    
    beta = 180-theta_inf; %Angle between the velocity vector of the planet 
                          %and the asymptote of the hyperbole; this
                          %parameter indicates the entrance angle (-beta)
                          %and exit angle (+beta)
                          
    delta = 180-2*beta; %Angle of rotation of the velocity vector
    
    e = -1/cos(theta_inf*deg); %Eccentricity
    
    %Minimum distance of the spaceship from the planet
    rp = str2double(planet(1)) + altitude; 
    
    a = rp/(e-1); %Semi-major axis of hyperbole
    b = a*sqrt(e^2-1); %Semi-minor axis of hyperbole
    
    %% Planet position
%     mu = str2double(planet(3));
    [coe, r, v, jd] = ...
                    planet_elements_and_sv(planet_id, year, month, day, ...
                                                    hour, minute, second);
    R_sun_planet = sqrt(r*r');
    longitude_angle = asin( r(3) / R_sun_planet );
    angle_sun_planet = acos( r(1) / (R_sun_planet*cos(longitude_angle)) );
    
    %% Critical points to compute the trajectory
    x = [];     
    y = [];     
    z = [];
    
    X_max = -a*cosh(-theta_inf*deg);
    Y_max = b*sinh(-theta_inf*deg);
    
    %% Compute the coordinates of the hyperbolic trajectory
    %To calculate the Trajectory the path has been split up in three parts

    %First part of the trajectory: from the enterance in the sphere of
    %influence to critical point
    x_a = (cos(theta_inf*deg)*str2double(planet(2)))-abs(rp)-a:0.1:X_max;
    y_a = -tan(beta*deg)*(x_a - X_max) - Y_max;
    
    %Save the coordinates of the second part of trajectory
    x = [x x_a];     
	y = [y y_a];     
    z = [zeros(1,length(x_a))];

    %Second part of the trajectory: when the spaceship transits close to
    %the planet
    i = 1;
    while i <= length(theta)
        q = theta(i)*deg; %Hyperbolic angle
            X = -a*cosh(q); %X coordinates
            if q <= 0
                Y =  b*sqrt(X^2/a^2-1); %negative Y coordinates
            else
                Y = -b*sqrt(X^2/a^2-1); %positive Y coordinates
            end
            %Save the coordinates of the second part of trajectory
            x = [x X];
            y = [y Y];
            z = [z 0];    
        i = i+1;
    end
    
    %Third part of trajectory: from the critical point to the exit from 
    %sphere of influence
    x_a = (cos(theta_inf*deg)*str2double(planet(2)))-abs(rp)-a:1:X_max;
    y_a = tan(beta*deg)*(x_a - X_max) + Y_max;
    
    %Save the coordinates of the third part of trajectory
    x = [x x_a]';    
	y = [y y_a]';    
    z = [z zeros(1,length(x_a))]'; 
    x = x + [abs(e*a)*ones(1,length(x))]';
    
    %% Rotation Matrices
    %R_x = [1 0 0; 0 cos(alpha_x*deg) sin(-alpha_x*deg); 0 sin(alpha_x*deg) cos(alpha_x*deg);];
    R_y = [ cos(longitude_angle)    0   sin(longitude_angle);...
            0                       1   0;...
            -sin(longitude_angle)   0   cos(longitude_angle)];
        
    R_z = [ cos(angle_sun_planet)   -sin(angle_sun_planet)  0;...
            sin(angle_sun_planet)   cos(angle_sun_planet)   0;...
            0                       0                       1];
    
    %% Plotting
    
    %Plot the Hyperbolic Trajectory 
    if exist('figure2') == 0
        figure()
    else
        figure2()
    end
    grid on;
    axis equal;
    hold on;
    title("Flyby trajectory close " + planet(4));
    
    G1 = [flag*x y z]; %Transposes the G1 Matrix 
    G1 = G1*R_y*R_z; %Executes the rotation 
    
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
    
    % Plot the circle that rapresent the Sun
%     sun_posiont = [0 0 0];
%     th = 0:pi/360:2*pi;
%     xunit = 6.96*10^5 * cos(th) + sun_posiont(1);
%     yunit = 6.96*10^5 * sin(th) + sun_posiont(2);
%     zunit = zeros(1,length(xunit));
%     plot3(xunit, yunit, zunit);
    
    hold off;
    
    %% Polar coordinates of flyby and relative plot
    
    rr =[]; %this will contain the distance from the center of the planet
            %to the position of the spaceship along the trajectory 
            
    i = 1;
    while i <= length(theta)
        vett_r = rp/(1+e*cos(theta(i)*deg)); 
        rr = [rr flag*vett_r]; %save the distance
        i = i + 1;
    end
    
    % Plot the trajectory
    if exist('figure2') == 0
        figure()
    else
        figure2()
    end
%     title('Polar Plot Flyby on ' + planet(4));
    polarplot(theta*deg,rr);
    title("Polar Plot: Flyby trajectory close to " + planet(4));
    hold on;
    
    % Plot the planet influence circle
    t = 0:1:360-1;
    rr = str2double(planet(2))*ones(1,length(t));
    polarplot(t*deg,rr,'-.');
    thetaticks([0:10:360]);
    polarplot(0,'.-k');
    
    % limit the plot to the planet influence circle
    rlim([0 str2double(planet(2))]);
    
    hold off;
    %% Print flyby data
    fprintf('\n   Radius of influence of ' + planet(4) + '     = ' + planet(2) + ' (km) %g\n') 
    fprintf('\n   Radius of ' + planet(4) + '                  = ' + planet(1) + ' (km) %g\n')
    fprintf('\n   Altitude                        = ' + string(altitude) + ' (km) %g\n')
    fprintf('\n   True anomaly                    = ' + string(theta_inf) + '%g\n')
    fprintf('\n   Beta angle                      = ' + string(beta) + '%g\n')
    fprintf('\n   Delta angle                     = ' + string(delta) + '%g\n')
    fprintf('\n   Eccentricity                    = ' + string(e) + '%g\n')
    fprintf('\n   V infinity                      = ' + string(sqrt(str2double(planet(3))/a))+ ' (km/s)' + '%g\n') 
    fprintf('\n   Semi major axis                 = ' + string(a) + ' (km) %g\n')
    fprintf('\n   Semi minor axis                 = ' + string(b)+ ' (km)' + '\n')
    
    %% Setting output parameters
    flyby_parameters = [planet(3) planet(1) planet(2) altitude theta_inf beta delta e a b];
    orb = G1;
    
end