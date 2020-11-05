function [flyby_parameters] = flyby(id_planet,theta_inf,altitude,flag,alpha_x,beta_y,gamma_z)

%example id_planet = 4; theta = 129; altitude = 512; flag = 1; alpha_x = 0; beta_y = 0; gamma_z = 0; 
%flyby(4,129,512,1,0,0,0);

%% Parameters required from the function
%Id_planet -> Code that indicates the planet where the flyby will be done
%Theta_inf -> True anomaly, angle between velocity vector of the planet and
%             escape velocity vector of spaceship => 90 < Theta_inf < 180
%Altitude -> Minimum distance from the planet
%flag ->% 1 = Leading-side planetary flyby, planet position on hyperbola primary focus; 
        %-1 = Trailing-side planetary flyby, planet position on hyperbola secondary focus;
%alpha_x -> Rotation angle rispect X axis 
%beta_y  -> Rotation angle rispect Y axis
%gamma_z -> Rotation angle rispect Z axis 
%% Starting parameters
planet = inf_planet(id_planet); %This command gives back a vector that contains planet's data:
                                %planet(1) = the radius 
                                %planet(2) = the radius of the sphere of influence
                                %planet(3) = the name
deg = pi/180;
%Angles of hyperbole
    theta = -theta_inf+0.1:0.1:theta_inf-0.1; %Vetor that contains angles between 
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
a = rp/(e-1); %Semiaxis major of hyperbole
b = a*sqrt(e^2-1); %Semiaxis minor of hyperbole
%% Parameter provided by the function
%flyby_parameters: vettor that contains flyby's parameters in this
%order [planet(3) planet(1) planet(2) altitude theta_inf beta delta e a b]
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
    x_a = (cos(theta_inf*deg)*str2double(planet(2)))-abs(rp)-a:0.1:X_max;
    y_a = -tan(beta*deg)*(x_a - X_max) - Y_max;
    %Save the coordinates of the second part of trajectory
    x = [x x_a];     
	y = [y y_a];     
    z = [zeros(1,length(x_a))];

    %Second part of the trajectory: when the spaceship trancits close the planet
    i = 1;
    while i <= length(theta)
        q = theta(i)*deg; %Hyperbolic angle
            X = -a*cosh(q); %Calculate the X coordinates
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
    x_a = (cos(theta_inf*deg)*str2double(planet(2)))-abs(rp)-a:1:X_max;
    y_a = tan(beta*deg)*(x_a - X_max) + Y_max;
    %Save the coordinates of the third part of trajectory
    x = [x x_a];    
	y = [y y_a];    
    z = [z zeros(1,length(x_a))]; 
%% Rotation Matrix
R_x = [1 0 0; 0 cos(alpha_x*deg) sin(-alpha_x*deg); 0 sin(alpha_x*deg) cos(alpha_x*deg);];
R_y = [cos(beta_y*deg) 0 sin(beta_y*deg); 0 1 0; -sin(beta_y*deg) 0 cos(beta_y*deg)];
R_z = [cos(gamma_z*deg) -sin(gamma_z*deg) 0; sin(gamma_z*deg) cos(gamma_z*deg) 0; 0 0 1];
%% Plot the Hyperbolic Trajectory 
grid on;
figure(1);
axis equal;
G1 = [flag*x' y' z']; %Transposes the G1 Matrix 
G1 = G1*R_x*R_y*R_z; %Executes the rotation 
hold on;
    figure(1);
    % Plot the trajectory of the spaceship
    plot3(G1(:,1),G1(:,2),G1(:,3));
    % Plot the circle that rapresent the planet
    x = -abs(rp)-abs(a);
    y = 0;
    th = 0:pi/360:2*pi;
    xunit = str2double(planet(1)) * cos(th) + flag*x;
    yunit = str2double(planet(1)) * sin(th) + y;
    plot(xunit, yunit);
    % Plot the circle that rapresent the radius of influnce of the planet 
    x = flag*(-abs(rp)-abs(a));
    y = 0;
    th = 0:pi/360:2*pi;
    xunit = str2double(planet(2)) * cos(th) + flag*x;
    yunit = str2double(planet(2)) * sin(th) + y;
    plot(xunit, yunit);
    title("Flyby trajectory close " + planet(4));
hold off;
%% Polar coordinates of flyby and plot
    r =[]; %this vector conteins the distance from the center of planet and the 
           %point on the spaceship's trajectory 
    i = 1;
    while i <= length(theta)
        vett_r = rp/(1+e*cos(theta(i)*deg)); 
        r = [r flag*vett_r]; %save the distance
        i = i + 1;
    end
    % Plot the spaceship's trajectory
    figure('Name','Polar Plot Flyby on ' + planet(4));
    polarplot(theta*deg,r);
    title("Polar Plot Flyby trajectory close " + planet(4));
    hold on;
    % Plot the planet's influence circle
    t = 0:1:360-1;
    rr = str2double(planet(2))*ones(1,length(t));
    polarplot(t*deg,rr,'-.');
    thetaticks([0:10:360]);
    polarplot(-e*a*flag,'.-k')
    % limit the plot to planet's influence circle
    rlim([0 str2double(planet(2))]);
    hold off;
%% Print datas of flyby
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
end