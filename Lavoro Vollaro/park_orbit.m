function orb = park_orbit(obj_id,pos,radius,incl)
% PARK_ORBIT(obj_id, pos, dist, incl) computes the circular parking
%   orbit of a spacecraft around OBJ_ID having the center in POS,
%   radius RADIUS and inclination INCL.
%
%   orb = PARK_ORBIT(...) returns the points composing the orbit.
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
    G     = 6.6742e-20; %[km^3/kg/s^2] gravitational constant

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
       
    masses = 10^24 * [0.330
              4.87
              5.97
              0.642
              1898
              568
              86.8
              102
              0.0146
              0.0002589
              0.000947
              1989100]; %[kg]
              
    radii = [2439.5
             6052 
             6378
             3396
             71492
             60268
             25559
             24764
             1185
             262.7
             476.2
             695508]; %[km]

    %% Input data
    %Object
    m1 = masses(obj_id); %[kg]
    R  = radii(obj_id); %[km]
    
    %Spacecraft
    m2 = 1000; %[kg]
    
    obj_mu    = G*(m1 + m2); %[km^3/s^2]
    r0 = (Rotx(incl)*[R+radius; 0; 0])'; %[km,km,km]
    
    %Parking orbit
    Park_v0 = sqrt(obj_mu/r0(1)); %[km/s]
    v0 = (Rotx(incl)*[0; Park_v0; 0])'; %[km/s,km/s,km/s]

    %Time
    t0 = 0; %[s]
    tf = 100*hours; %[s]

    %% Numerical integration:
    
    %Initial condition: position, velocity
    y0    = [r0 v0]';
    [t,y] = rkf45(@rates, [t0 tf], y0);

    %% Output the results:
    output
    
    orb = pos + y(:,1:3);

    return

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
        body_sphere(obj_id,pos);
        hold on
        plot3(  pos(1)+y(:,1),    pos(2)+y(:,2),    pos(3)+y(:,3),...
                                                   'r', 'LineWidth', 1)

    end %output
    % ~~~~~~~~~~~~~~~~~~~~~~~

end %orbit