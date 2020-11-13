function [flyby_parameters] = ...
                flyby(id_planet,theta_inf,altitude,alpha_x,beta_y,gamma_z)
% FLYBY help section
%
% INPUT:
%   id_planet - Code indicating the body where the flyby will be
%               carried on
%   theta_inf - True anomaly at infinity
%   altitude  - Minimum distance of the spacecraft from the planet
%   alpha_x   - Rotation angle respect the X axis 
%   beta_y    - Rotation angle respect the Y axis
%   gamma_z   - Rotation angle respect the Z axis 
%
% OUTPUT:
%   flyby_parameters - vector of flyby parameters:

    %% Starting parameters
    planet = inf_planet(id_planet); %returns planet data:
                                    %planet(1) = radius 
                                    %planet(2) = radius of the sphere of influence
                                    %planet(3) = name
    deg = pi/180;
    
    %Angles of hyperbole
    theta = -theta_inf+0.1:0.1:theta_inf-0.1;
    beta = 180-theta_inf; %Angle between the planet velocity vector
                          %and the asymptote of the hyperbole; this
                          %parameter indicates the entrance angle (-beta)
                          %and exit angle (+beta)
    delta = 180-2*beta; %This parameter indicates how much the infinity 
                        %velocity vector of the spacecraft rotates after 
                        %the flyby manouver
    e = -1/cos(theta_inf*deg); %Eccentricity
    rp = str2double(planet(1)) + altitude; %Minimum distance between the center
                                           %of the planet and spacecraft
    a = rp/(e-1); %Semi-major axis of hyperbole
    b = a*sqrt(e^2-1); %Semi-minor axis of hyperbole

    %% Coordinates vectors 
    x = [];     
    y = [];     
    z = [];

    %% Critical points to calculate the trajectory
    X_max = -a*cosh(-theta_inf*deg);
    Y_max = b*sinh(-theta_inf*deg);

    %% Calculate the x,y and z coordinates of Hyperbolic Trajectory
    %In order to calculate the trajectory the path has been split up 
    %in three parts

        %First part of the trajectory: from the entering point in the 
        %sphere of influence to a critical point
        x_a = (cos(theta_inf*deg)*str2double(planet(2)))-abs(rp)-a:0.1:X_max;
        y_a = -tan(beta*deg)*(x_a - X_max) - Y_max;
        
        %Save the coordinates of the first part of trajectory
        x = [x x_a];     
        y = [y y_a];     
        z = zeros(1,length(x_a));

        %Second part of the trajectory: when the spaceship transits close 
        %the planet
        i = 1;
        while i <= length(theta)
            q = theta(i)*deg; %Hyperbolic angle
                X = -a*cosh(q); %Calculates the X coordinates
                if q <= 0
                    %Calculates the negative Y semi-axis coordinates
                    Y =  b*sqrt(X^2/a^2-1); 
                else
                    %Calculates the positive Y semi-axis coordinates
                    Y = -b*sqrt(X^2/a^2-1); 
                end
                %Save the coordinates of the second part of trajectory
                x = [x X];
                y = [y Y];
                z = [z 0];    
            i = i+1;
        end

        %Third part of trajectory: from the critical point to the exit 
        %from sphere of influence
        x_a = (cos(theta_inf*deg)*str2double(planet(2)))-abs(rp)-a:1:X_max;
        y_a = tan(beta*deg)*(x_a - X_max) + Y_max;
        
        %Save the coordinates of the third part of trajectory
        x = [x x_a];    
        y = [y y_a];    
        z = [z zeros(1,length(x_a))]; 

    %% Rotation Matrices
    R_x = [1                    0                   0; ...
           0                    cos(alpha_x*deg)    sin(-alpha_x*deg); ...
           0                    sin(alpha_x*deg)    cos(alpha_x*deg);];
       
    R_y = [cos(beta_y*deg)      0                   sin(beta_y*deg); ...
           0                    1                   0; ...
           -sin(beta_y*deg)     0                   cos(beta_y*deg)];
       
    R_z = [cos(gamma_z*deg)     -sin(gamma_z*deg)   0; ...
           sin(gamma_z*deg)     cos(gamma_z*deg)    0; ...
           0                    0                   1];

    %% Plot the Hyperbolic Trajectory 
%     figure2();
%     hold on;
%     grid on;
%     view(-10,45)
    
    G1 = [x' y' z']; %Transposes the G1 Matrix 
    G1 = G1*R_x*R_y*R_z; %Executes the rotation 
    
    % Plots the trajectory of the spacecraft
    plot3(G1(:,1),G1(:,2),G1(:,3));
    
    % Plots the circle that represents the planet
    x = -abs(rp)-abs(a);
    y = 0;
    th = 0:pi/360:2*pi;
    xunit = str2double(planet(1)) * cos(th) + x;
    yunit = str2double(planet(1)) * sin(th) + y;
    plot(xunit, yunit);
    
    % Plot the circle that represents the radius of influence of the planet 
    x = -abs(rp)-abs(a);
    y = 0;
    th = 0:pi/360:2*pi;
    xunit = str2double(planet(2)) * cos(th) + x;
    yunit = str2double(planet(2)) * sin(th) + y;
    plot(xunit, yunit);
    
%     title("Flyby trajectory around " + planet(4));
%     hold off;
%     axis equal;

    %% Polar coordinates of flyby and plot
    %Vector containing the distance from the center of the planet to the 
    %point on the spacecraft trajectory 
    r =[]; 
    
    i = 1;
    while i <= length(theta)
        vett_r = rp/(1+e*cos(theta(i)*deg)); 
        r = [r vett_r]; %save the distance
        i = i + 1;
    end
    
    % Plots the spacecraft trajectory
    figure('Name','Polar Plot Flyby on ' + planet(4));
    polarplot(theta*deg,r);
    title("Polar Plot Flyby trajectory close " + planet(4));
    hold on;
    
    % Plots the planet influence circle
    t = 0:1:360-1;
    rr = str2double(planet(2))*ones(1,length(t));
    polarplot(t*deg,rr,'-.');
    thetaticks([0:10:360]);
    polarplot(-e*a,'.-k')
    
    %Limits the plot to planet influence circle
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

    %% Setting output parameter
    flyby_parameters = G1;
end