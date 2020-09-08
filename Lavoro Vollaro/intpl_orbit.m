function orb = intpl_orbit(arr_days,r_init,v_init)
%   INTERPLANETARY_ORBIT(arr_days, r_init, v_init) computes the 
%   interplanetary orbit a spacecraft of 1000kg would travel
%   around the Sun in ARR_DAYS days, starting at position R_INIT
%   with velocity V_INIT.
%   It uses the patched conics method.
%   orb = INTPL_ORBIT(...) returns the points composing the orbit.
%
%   It uses rkf45 to numerically integrate Equation 2.22 in
%   "Orbital Mechanics for Engineering Students" - Howard D. Curtis.
%
%   arr_days  - days needed for the spacecraft to arrive at destination
%   r_init    - initial position of the spacecraft
%   v_init    - initial velocity of the spacecraft
%
% User M-function required:   rkf45
% User subfunctions required: rates, output

    %% Constants
    global mu
    hours = 3600; % [s]
    Sun_radius = 695508; %[km]   

    %% Input data:
    
    %Sun
    R  = Sun_radius; %[km]
    
    r0 = r_init; %[km,km,km]
    v0 = v_init; %[km/s,km/s,km/s]

    t0 = 0; %[s]
    tf = arr_days*24*hours; % [s]

    Sun_mu  = mu; %[km^3/s^2]
    
    %% Numerical integration:
    
    %Initial condition: position, velocity
    y0    = [r0 v0]';
    [t,y] = rkf45(@rates, [t0 tf], y0,1.e-15);

    %% Output the results:
    output
    
    orb = y;

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

        ax   = -Sun_mu*x/r^3;
        ay   = -Sun_mu*y/r^3;
        az   = -Sun_mu*z/r^3;

        dydt = [vx vy vz ax ay az]';    
    end %rates
    % ~~~~~~~~~~~~~~~~~~~~~~~


    %% ~~~~~~~~~~~~~~~~~~~~~~~
    function output
    %{
      This function computes the maximum and minimum radii, the times they
      occur and and the speed at those times. It prints those results to
      the command window and plots the orbit.

      r         - magnitude of the position vector at the times in t
      imax      - the component of r with the largest value
      rmax      - the largest value of r
      imin      - the component of r with the smallest value
      rmin      - the smallest value of r
      v_at_rmax - the speed where r = rmax
      v_at_rmin - the speed where r = rmin

      User subfunction required: light_gray
    %}
        
        for i = 1:length(t)
            r(i) = norm([y(i,1) y(i,2) y(i,3)]);
        end

        [rmax, imax] = max(r);
        [rmin, imin] = min(r);

        v_at_rmax   = norm([y(imax,4) y(imax,5) y(imax,6)]);
        v_at_rmin   = norm([y(imin,4) y(imin,5) y(imin,6)]);

        %% Output to the command window:
        fprintf('\n\n--------------------------------------------------------\n')
        fprintf('\n Spacecraft Orbit\n')
        fprintf(' %s\n', datestr([2007,9,27,0,0,0]))
        fprintf('\n The initial position is [%g, %g, %g] (km).',...
                                                             r0(1), r0(2), r0(3))
        fprintf('\n   Magnitude = %g km\n', norm(r0))
        fprintf('\n The initial velocity is [%g, %g, %g] (km/s).',...
                                                             v0(1), v0(2), v0(3))
        fprintf('\n   Magnitude = %g km/s\n', norm(v0))
        fprintf('\n Initial time = %g h.\n Final time   = %g h.\n',0,tf/hours) 
        fprintf('\n The minimum altitude is %g km at time = %g h.',...
                    rmin-R, t(imin)/hours)
        fprintf('\n The speed at that point is %g km/s.\n', v_at_rmin)
        fprintf('\n The maximum altitude is %g km at time = %g h.',...
                    rmax-R, t(imax)/hours)
        fprintf('\n The speed at that point is %g km/s\n', v_at_rmax)
        fprintf('\n--------------------------------------------------------\n\n')

        %% Plot the results:
        
        hold on
        plot3(  y(:,1),    y(:,2),    y(:,3), 'k', 'LineWidth', 2)

    end %output
    % ~~~~~~~~~~~~~~~~~~~~~~~

end %orbit