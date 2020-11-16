function [orb, ts] = ...
            park_in(obj_id, body_pos, radius, coe, p_ref, start, finish)
% PARK_IN(obj_id, body_pos, radius, coe, ref, start, finish) computes the 
%   circular parking orbit of a spacecraft around OBJ_ID having the center 
%   in POS and radius RADIUS.
%   The function assumes that the parking orbit follows an arrival
%   hyperbole, takes a reference point P_REF from it and uses the desired
%   orbital elements COE (of which will use the inclination and the RAAN): 
%   it uses P_REF to ensure that's the first point in the parking orbit,
%   "accomodating" in some way the orbit to its coordinates.
%   START and FINISH are the time coordinates for which PARK_IN computes
%   the parking orbit.
% 
%   In case of a parking orbit without the hypothesis of an arrival 
%   hyperbole, it is sufficient to set p_ref = [0 0 0].
%
%   [orb, ts] = PARK_IN(...) returns the points composing the orbit and
%   an array of time istants relative to each point.
%
%   It uses rkf45 to numerically integrate Equation 2.22 in
%   "Orbital Mechanics for Engineering Students" - Howard D. Curtis.
%
%   obj_id   - identifier of the main body in the hypothesis of 2-body
%                 problem:
%                1 = Mercury
%                2 = Venus
%                3 = Earth
%                4 = Mars
%                5 = Jupiter
%                6 = Saturn
%                7 = Uranus
%                8 = Neptune
%                9 = Pluto
%               10 = Vesta
%               11 = Ceres
%               12 = Sun
%
%   pos       - position of the planet in heliocentric components
%
%   radius    - radius of the circular parking orbit
%
%   incl      - inclination of the parking orbit
%
% User M-function required:   rkf45, body_sphere
% User subfunctions implemented: rates, output
    
    %% Constants

    hours = 3600; %[s] to convert from hours to seconds
    G = 6.6742e-20; %[km^3/kg/s^2], gravitational constant

    planets = ["Mercury"
               "Venus  "
               "Earth  "
               "Mars   "
               "Jupiter"
               "Saturn "
               "Uranus "
               "Neptune"
               "Pluto  "
               "Vesta  "
               "Ceres  "
               "Sun    "];
       
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
             695508]; %[km]

    %% Input data
    
    incl = coe(4);
    raan = coe(3);
    
    %Time spent in orbit, from Julian days to seconds
    start_time = J0(start(1),start(2),start(3))*86400 + ...
                 start(4)*3600 + start(5)*60 + start(6); %[s]
    finish_time = J0(finish(1),finish(2),finish(3))*86400 + ...
                 finish(4)*3600 + finish(5)*60 + finish(6); %[s]
	time = finish_time - start_time; %[s]
    
    %Object
    m1 = masses(obj_id); %[kg]
    R  = radii(obj_id); %[km]
    
    %Spacecraft
    m2 = 1000; %[kg]
    
    obj_mu    = G*(m1 + m2); %[km^3/s^2]
    
    %Parking orbit
    park_v0 = sqrt(obj_mu/(R+radius)); %[km/s]
    
    ref_point = p_ref - body_pos;
    ref_raan = deg2rad(atan2d_0_360(ref_point(2),ref_point(1)));
    
    %This or also 2*pi*(R+radius)/sqrt(obj_mu/(R+radius))
    period = 2*pi*(R+radius)/park_v0;
    
    %Computing the point relative to an Argument of Periapsis equal to 0
    TA_0 = wrapTo2Pi(2*pi*start_time/period);
    coe_0 = [(R+radius)*park_v0 0 raan incl 0 TA_0];
    sv_0 = sv_from_coe(coe_0, obj_mu);
    ang_0 = deg2rad(atan2d_0_360(sv_0(2),sv_0(1)));
    
    %Starting point considered for the parking orbit
    TA_start = wrapTo2Pi(2*pi*start_time/period);
    coe_start = [(R+radius)*park_v0 0 raan incl ref_raan-ang_0 TA_start];
    %This does not consider object absolute position
    [r0, v0] = sv_from_coe(coe_start,obj_mu);

%     %Debug
%     pos_start = r0+body_pos;%sv_finish + body_pos;
%     pos_start = sv_start + body_pos;
%     plot3(pos_start(1),pos_start(2),pos_start(3),'bo')
%     plot3(pos_finish(1),pos_finish(2),pos_finish(3),'bo')

    %Time
    t0 = 0; %[s]
    tf = period;%100*hours; %[s]

    %% Numerical integration:
    
    %Initial condition: position, velocity
    y0    = [r0 v0]';
    [t,y] = rkf45(@rates, [t0 tf], y0);

    %% Computations
    div = floor(time/period);
    y_app = [];
    t_app = [];
    ind = 0;
    
    %Here the function checks if the period of time spent in orbit by the
    %spacecraft is bigger or smaller than the time needed for a complete
    %revolution around the body, and it acts accordingly
    if (div > 0) %more than a revolution
        y_app = zeros(size(y,1)*div,3);
        t_app = zeros(size(y,1)*div,1);
        for i = 0:div-1 %for each revolution
            y_app(1+size(y,1)*i:size(y,1)*(i+1),1:3) = y(:,1:3);
            t_app(1+size(y,1)*i:size(y,1)*(i+1)) = t+i*t(end);
        end
        time = time-div*period;
    end
    
    %Checks the length of the last part of the parking orbit
    for i = 1:length(t)
        if (t(i)-t(1)>time)
            ind = i;
            break;
        end
    end
    
    %Selecting the used part of the array for the last part of the orbit
    y_app = cat(1,y_app,y(1:ind,1:3));
    t_app = cat(1,t_app,t(1:ind)+div*t(end));

    %% Setting output parameters
    orb = body_pos + y_app;
    ts = t_app;
    
    %% Output the results
    output
    return

    %% Used functions
    %% ~~~~~~~~~~~~~~~~~~~~~~~~
    function dydt = rates(t,f)
    %{
      This function calculates the acceleration vector using Equation 2.22

      t          - time
      f          - column vector containing the position vector and the
                   velocity vector at time t
      x, y, z    - components of the position vector r
      r          - the magnitude of the the position vector
      vx, vy, vz - components of the velocity vector v
      ax, ay, az - components of the acceleration vector a
      dydt       - column vector containing the velocity and acceleration
                   components
    %}
        
        x    = f(1);
        y    = f(2);
        z    = f(3);
        vx   = f(4);
        vy   = f(5);
        vz   = f(6);

        r    = norm([x y z]);

        ax   = -obj_mu*x/r^3;
        ay   = -obj_mu*y/r^3;
        az   = -obj_mu*z/r^3;

        dydt = [vx vy vz ax ay az]';    
    end %rates
    % ~~~~~~~~~~~~~~~~~~~~~~~

    %% ~~~~~~~~~~~~~~~~~~~~~~~
    function output
    %{
      This function prints computation results to
      the command window and plots the orbit.
    %}

        %% Output to the command window:
        fprintf('\n\n--------------------------------------------------------\n')
        fprintf('\n %s parking orbit\n',planets(obj_id))
        fprintf('\n Radius of the circular orbit: %g (km)\n', R+radius)
        fprintf('\n Inclination of the orbit: %4.2f (°)\n', rad2deg(incl))
        fprintf('\n The initial position is [%g, %g, %g] (km).',...
                                                 r0(1), r0(2), r0(3))
        fprintf('\n   Magnitude = %g km\n', norm(r0))
        fprintf('\n The initial velocity is [%g, %g, %g] (km/s).',...
                                                     v0(1), v0(2), v0(3))
        fprintf('\n   Magnitude = %g km/s\n', norm(v0))
        fprintf('\n Initial time = %g h.\n Final time   = %g h.\n',0,...
                                                            tf/hours) 
        fprintf('\n--------------------------------------------------------\n\n')

        %% Figure plot
        body_sphere(obj_id,body_pos);
        hold on
        if (div>0) %if more than a revolution
            plot3(orb(1:size(y,1)*div,1),...
                  orb(1:size(y,1)*div,2),...
                  orb(1:size(y,1)*div,3),...
                  'r', 'LineWidth', 1)
        end
        plot3(orb(1+size(y,1)*div:end,1),...
              orb(1+size(y,1)*div:end,2),...
              orb(1+size(y,1)*div:end,3),...
              'b', 'LineWidth', 1) %different color to visualize last part
          
%           DEBUG
%         plot3(orb(5,1),orb(5,2),orb(5,3),'go')
%         plot3(orb(10,1),orb(10,2),orb(10,3),'co')
    end %output
    % ~~~~~~~~~~~~~~~~~~~~~~~

end %orbit